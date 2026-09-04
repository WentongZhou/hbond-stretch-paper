# PySCF DM-EDA 详细使用指南

本指南面向 **pyscf-dm-eda** 包的使用者，从安装、单结构 EDA、命令行、能量分项解读，
一直到批量网格计算、断点续算和常见故障排查，按步骤详细讲解。

阅读完本指南后你可以：

- 用 Python API 或命令行对一个 XYZ 结构做 counterpoise 密度矩阵 EDA；
- 正确设置片段、电荷与带符号自旋；
- 读取并解释每个能量分项和数值校验诊断；
- 用 `PySCFGridRunner` 批量扫描一个分子周围的大量探针位置；
- 用分片、锁和检查点安全地并行与续算；
- 根据报错信息快速定位问题。

公式与代码对应关系见 `MATHEMATICAL_PRINCIPLES.md`，本指南偏重操作步骤。

---

## 目录

1. [安装](#1-安装)
2. [核心概念：一次 EDA 内部发生了什么](#2-核心概念一次-eda-内部发生了什么)
3. [五分钟快速开始（命令行）](#3-五分钟快速开始命令行)
4. [Python API 单结构 EDA 分步教程](#4-python-api-单结构-eda-分步教程)
5. [S​​CFConfig 参数表](#5-scfconfig-参数表)
6. [能量分项与输出列含义](#6-能量分项与输出列含义)
7. [结果对象：分量、单位、JSON 与诊断](#7-结果对象分量单位json-与诊断)
8. [常见计算场景](#8-常见计算场景)
9. [批量网格计算：PySCFGridRunner](#9-批量网格计算pyscfgridrunner)
10. [网格检查点、分片与重启](#10-网格检查点分片与重启)
11. [完整实战：球形水二聚体扫描](#11-完整实战球形水二聚体扫描)
12. [故障排查](#12-故障排查)
13. [限制与注意事项](#13-限制与注意事项)

---

## 1. 安装

### 1.1 创建环境

推荐 Python 3.10+。若使用 conda：

```bash
conda create -n edaenv python=3.12 -y
conda activate edaenv
```

也可以使用普通 venv。本仓库开发环境中已有 `edaenv`，其解释器路径为
`~/edaenv/bin/python`；以下命令均假设使用该环境。

### 1.2 安装包

```bash
cd /mnt/d/pyscf-dm-eda
python -m pip install --upgrade pip
python -m pip install -e .
```

需要 D3/D4 经验色散时：

```bash
python -m pip install -e ".[dispersion]"
```

### 1.3 验证安装

```bash
python -c "import pyscf, pyscf_dm_eda; print(pyscf.__version__, pyscf_dm_eda.__version__)"
pyscf-dm-eda --help
pyscf-dm-eda grid --help
```

---

## 2. 核心概念：一次 EDA 内部发生了什么

一次 `PySCFEDA.run()` 会依次完成以下工作：

1. 读取超分子全部原子与基组；
2. 对每个片段分别做 SCF：**片段以外的原子保留为 ghost 原子**，因此所有片段都
   使用与超分子完全相同的 AO 基组（counterpoise 约定），没有 BSSE 修正步骤；
3. 从片段密度构造 promolecule 密度 `P0 = sum_A P_A`；
4. 用片段占据轨道构造 Pauli 正交化密度 `P_Pauli`；
5. 对超分子做一次 SCF，初猜为 promolecule 密度；
6. 在一次批量 J/K 调用中计算所有密度矩阵的库仑/交换势；
7. 按 DM-EDA 定义分解：
   - 静电 = 核-核 + 交叉核-电子 + 交叉库仑；
   - 交换 + Pauli 排斥 = exchange-repulsion；
   - 轨道弛豫/极化 = 超分子 E0 与 Pauli 态 E0 之差；
   - 相关残差 = 泛函（或 HF）残余差；
   - 显式 D3/D4 色散差（若启用）；
8. 做 Mulliken 片段电荷转移，并检查能量闭合、电子数、幂等性等数值诊断。

**必须理解的三个约定：**

| 约定 | 说明 |
|---|---|
| 原子索引 | Python API 从 **0** 开始；命令行从 **1** 开始 |
| `spin` | 带符号 `N_alpha - N_beta`，**不是** `2S+1` |
| 片段电荷/自旋 | 超分子 `charge = Σ 片段 charge`，`spin = Σ 片段 spin`，不满足会直接报错 |

---

## 3. 五分钟快速开始（命令行）

仓库自带 `examples/he2.xyz`。把 He₂ 拆成两个单原子片段：

```bash
cd /mnt/d/pyscf-dm-eda
pyscf-dm-eda examples/he2.xyz \
  --fragment 1 --fragment 2 \
  --method hf --basis sto-3g \
  --unit kcal/mol \
  --output he2_eda.json
```

也可以用模块方式运行：

```bash
python -m pyscf_dm_eda examples/he2.xyz \
  --fragment 1 --fragment 2 \
  --method r2scan --basis def2-svp
```

不写 `--output` 时，结果 JSON 直接打印到 stdout。

**多原子片段**用逗号和闭区间，片段电荷/自旋/标签按相同顺序重复传入：

```bash
pyscf-dm-eda complex.xyz \
  --fragment 1-8 --fragment 9 \
  --fragment-charge 0 --fragment-charge 1 \
  --fragment-spin 0 --fragment-spin 0 \
  --fragment-label molecule --fragment-label probe \
  --charge 1 --spin 0 \
  --method r2scan --basis def2-svp \
  --output complex_eda.json
```

命令行参数说明：

| 参数 | 含义 |
|---|---|
| `--fragment` | 必填，可重复；`1-3,7` 表示原子 1,2,3,7（一基索引） |
| `--fragment-charge/spin/label` | 与 `--fragment` 数量一致，按顺序对应 |
| `--charge / --spin` | 超分子总电荷/带符号自旋；缺省时由片段求和 |
| `--unit` | `hartree`、`kcal/mol`、`kj/mol`、`ev` |
| `--method / --basis` | 方法名 / 基组名，直接传给 PySCF |
| `--output` | 可选 JSON 输出路径 |
| `--strict-validation / --no-strict-validation` | 校验失败抛异常还是降级为 RuntimeWarning |

退出码约定：参数错误 `exit 2`；SCF/数值运行错误打印到 stderr 并 `exit 1`。

---

## 4. Python API 单结构 EDA 分步教程

### 4.1 导入

```python
from pyscf_dm_eda import (
    Atom, FragmentSpec, PySCFEDA, SCFConfig, read_xyz,
)
```

公开 API 还包括：

```python
from pyscf_dm_eda import (
    EDAResult, GridRunSummary, PySCFGridRunner,
    EDAError, EDAValidationError, SCFConvergenceError,
    IncompatibleFragmentError,
)
```

### 4.2 准备原子

三种等价写法：

```python
# 写法 1：Atom 对象
atoms = [
    Atom("O", (0.0, 0.0, 0.0)),
    Atom("H", (0.0, 0.0, 0.9572)),
]

# 写法 2：列表元素为 (symbol, x, y, z)
atoms = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.9572)]

# 写法 3：从 XYZ 文件读取
atoms = read_xyz("water_dimer.xyz")
```

坐标单位由 `SCFConfig(unit=...)` 控制，默认 Angstrom。

### 4.3 定义片段

原子索引从 0 开始，每个原子必须且只能属于一个片段：

```python
fragments = [
    FragmentSpec((0, 1, 2), charge=0, spin=0, label="water_A"),
    FragmentSpec((3, 4, 5), charge=0, spin=0, label="water_B"),
]
```

`label` 可省略，缺省自动命名为 `fragment_1`、`fragment_2`…… 但建议显式命名，
它同时用作 JSON 中 Mulliken 电荷转移的键。

### 4.4 配置计算

```python
config = SCFConfig(
    method="r2scan",      # XC 泛函或 hf/rhf/uhf
    basis="def2-svp",
    grid_level=4,         # DFT 数值积分网格
    conv_tol=1e-9,
    max_cycle=100,
    verbose=0,
)
```

### 4.5 运行并读取能量

```python
eda = PySCFEDA(atoms, fragments, config, charge=0, spin=0)
result = eda.run()

print(result.components("kcal/mol"))
```

`charge` / `spin` 若省略，自动等于所有片段之和；若显式给出但和片段求和不一致，
抛出 `IncompatibleFragmentError`。

从 XYZ 文件一步构造：

```python
eda = PySCFEDA.from_xyz(
    "water_dimer.xyz",
    fragments=fragments,
    config=config,
    charge=0,
    spin=0,
)
result = eda.run()
```

### 4.6 保存结果

```python
result.write_json("water_dimer_eda.json", unit="kcal/mol")
```

如果不需要落盘，也可以直接拿可 JSON 序列化的字典：

```python
payload = result.as_dict("kcal/mol")
print(payload["components"])
```

JSON 中包含：能量分量、Mulliken 电荷转移、全部数值诊断、方法与基组元数据。
它使用原子替换写入，可直接机器读取：

```python
import json
with open("water_dimer_eda.json", encoding="utf-8") as f:
    payload = json.load(f)
print(payload["components"]["Total Interaction energy"])
```

---

## 5. SCFConfig 参数表

```python
SCFConfig(
    method="r2scan",           # XC 泛函；或 "hf", "rhf", "uhf"
    basis="def2-svp",          # 基组；也可传 dict 等 PySCF 接受的形式
    ecp=None,                  # 可选 ECP，如 "def2-ecp"
    dispersion=None,           # None | "d3" | "d3bj" | "d4"
    grid_level=4,              # DFT 网格等级，0..9
    conv_tol=1e-9,             # SCF 收敛阈值
    max_cycle=100,             # 最大 SCF 圈数
    max_memory=4000.0,         # PySCF 内存上限 MB
    density_fit=False,         # True 用 PySCF DF-J/K 近似
    auxbasis=None,             # DF 辅助基组，如 "def2-jkfit"
    unrestricted=None,         # None 自动 / True 强制 / False 强制 restricted
    level_shift=0.0,           # 收敛困难时的 level shift
    damp=0.0,                  # SCF damping
    init_guess="minao",        # 初猜
    newton_fallback=True,      # 不收敛时自动转 Newton SCF
    linear_dep_threshold=1e-9, # AO 线性相关判定
    validation_tol=1e-7,       # 数值校验阈值（a.u.）
    strict_validation=True,    # 校验失败是否抛异常
    verbose=0,
    unit="Angstrom",
    solvent=None,              # None | "cpcm"(="pcm") | "iefpcm" | "ssvpe" | "cosmo" | "smd"
    solvent_eps=None,          # PCM 介电常数；None 时用 solvent_name 或 PySCF 默认水 78.3553
    solvent_name=None,         # PySCF SMD 数据库中的溶剂名，如 "water"、"acetonitrile"
    solvent_options=None,      # 额外设置到 PySCF 溶剂对象的属性，如 {"lebedev_order": 29}
)
```

**隐式溶剂的规则（详见 8.6 节）：**

- `solvent` 选择 PySCF 的 PCM 族或 SMD；ddCOSMO/ddPCM 会被拒绝；
- PCM 模型的介电常数来自 `solvent_eps`，否则来自 `solvent_name` 在 PySCF SMD
  数据库中的条目，两者都不给时用 PySCF 默认值（水）；
- SMD 必须靠 `solvent_name`（默认 `"water"`）取 CDS 描述符，`solvent_eps` 只能
  在给了名字的前提下覆盖介电常数；
- `solvent_options` 的键必须是 PySCF 溶剂对象的可设置属性（例如
  `lebedev_order`、`vdw_scale`、`r_probe`、`radii_table`、`conv_tol`），
  `eps`、`method`、`solvent`、`mol` 等由程序管理的属性会被拒绝。

**方法与色散的规则：**

- `method` 写基础泛函，D3/D4 写进 `dispersion`；
- 不要把 `-d3/-d4` 后缀混进 `method`，会被拒绝；
- `wb97x-d`、`b97-d` 这类参数化本身包含经验色散的泛函会被拒绝，请使用其底层
  XC 泛函并显式指定 `dispersion`；
- 使用 D3/D4 前必须安装 `pyscf-dispersion`，否则在取色散能量时报
  `EDAError`；
- `density_fit=True` 是 PySCF 的 DF-J/K 近似，结果元数据会明确标记
  `"jk_evaluation": "density-fitted"`，不能冒充 exact J/K。

**restricted / unrestricted 选择：**

- `method="hf"` 时按超分子与片段自旋自动选 RHF/UHF；
- `method="uhf"` 强制 UHF；`method="rhf"` 强制 RHF；
- 任意片段开壳层时超分子自动用 UHF/UKS；
- 用 `unrestricted=True` 可对闭壳层体系强制 UHF/UKS；`False` 强制 restricted，
  但不允许非零总自旋。

---

## 6. 能量分项与输出列含义

`result.components()` 的主键：

| 键 | 含义 |
|---|---|
| `Total Interaction energy` | counterpoise 总相互作用能 |
| `Electrostatic Interaction` | 核-核 + 交叉核-电子 + 交叉库仑 |
| `Nuc---Nuc` | 片段间核排斥 |
| `1-electron` | 片段间核-电子吸引 |
| `2-electron` | 片段间库仑排斥 |
| `Exchange Int.` | 片段间 exact exchange（负值） |
| `Repulsion` | Pauli 正交化代价 |
| `Exchange-Repulsion` | `Exchange Int. + Repulsion` |
| `Orbital Relaxation` | 超分子密度相对 Pauli 态的弛豫/极化 |
| `Correlation Interaction` | 不含显式 D3/D4 的泛函残差；HF 为 0 |
| `Dispersion Interaction` | D3/D4 相互作用差（未启用为 0） |
| `Desolvation` | 反应场能量变化 `E_solv(super) - Σ E_solv(A)`（未启用溶剂为 0） |
| `Closure Error` | 总相互作用能与分项之和的差（诊断） |

衍生量：

```text
Corr_Disp = Correlation Interaction + Dispersion Interaction
Steric    = Exchange-Repulsion + Corr_Disp
```

**批量网格 TSV 的列顺序**（固定 kcal/mol）：

```text
Grid_Index Tot Electro Nuc_Nuc 1e 2e Exc_Rep Exc Rep
Orb_Relax Corr Disp Desolv Corr_Disp Steric Mulliken_CT Closure_Error
```

- `Tot` = 总相互作用能；
- `Desolv` = `Desolvation`，气相网格中恒为 0；
- `Mulliken_CT` = 探针片段相对其孤立片段电荷的 Mulliken 电子转移（单位 e）；
  本项目明确不把它伪称为 NBO/NPA；
- `Closure_Error` 应接近 0，默认阈值 `1e-7` a.u.。

---

## 7. 结果对象：分量、单位、JSON 与诊断

### 7.1 单位

```python
for unit in ["hartree", "kcal/mol", "kj/mol", "ev"]:
    c = result.components(unit)
    print(unit, c["Total Interaction energy"])
```

单位名忽略大小写和空格，例如 `"KCAL / MOL"` 也是合法的。

### 7.2 电荷转移与诊断

```python
print(result.fragment_charge_transfer)   # {"water_A": ..., "water_B": ...}
print(result.diagnostics["closure_error_hartree"])
print(result.diagnostics["validation_worst_case"])
print(result.diagnostics["validation_errors"])
print(result.metadata["jk_evaluation"])  # "exact" 或 "density-fitted"
print(result.metadata["pyscf_version"])
```

### 7.3 校验失败时怎么办

默认 `strict_validation=True`，任一残差超过 `validation_tol` 会抛
`EDAValidationError`。异常对象**携带完整结果**，不必重算：

```python
from pyscf_dm_eda import EDAValidationError

try:
    result = eda.run()
except EDAValidationError as exc:
    print(exc)                       # 最差校验项与误差
    result = exc.result              # 完整 EDAResult
    print(result.components("kcal/mol"))
```

适合弥散基组/近线性相关基组时，可先检查分量再决定：

```python
config = SCFConfig(..., strict_validation=False, validation_tol=1e-6)
result = eda.run()   # 改为 RuntimeWarning 并返回结果
```

---

## 8. 常见计算场景

### 8.1 两个水分子

```python
from pyscf_dm_eda import Atom, FragmentSpec, PySCFEDA, SCFConfig

atoms = [
    Atom("O", (0.0000, 0.0000,  0.0000)),
    Atom("H", (0.7569, 0.0000,  0.5859)),
    Atom("H", (-0.7569, 0.0000, 0.5859)),
    Atom("O", (0.0000, 0.0000, -2.9140)),
    Atom("H", (0.0000, 0.0000, -1.9568)),   # O...H = 0.9572 Å，指向受体 O
    Atom("H", (0.9266, 0.0000, -3.1541)),   # 保持 H-O-H = 104.52°
]
fragments = [
    FragmentSpec((0, 1, 2), label="acceptor"),
    FragmentSpec((3, 4, 5), label="donor"),
]
result = PySCFEDA(
    atoms, fragments,
    SCFConfig(method="r2scan", basis="def2-svp", grid_level=4),
).run()
print(result.components("kcal/mol"))
```

### 8.2 开壳层 / broken symmetry

`spin` 是 `N_alpha - N_beta`。两个自由基耦合为 broken-symmetry singlet 时：

```python
fragments = [
    FragmentSpec((0,), spin=1,  label="alpha"),
    FragmentSpec((1,), spin=-1, label="beta"),
]
eda = PySCFEDA(atoms, fragments, config)   # 总 spin 自动为 0
```

> 这是符号约定，不是多重度。给出 `spin=1` 意味着 α 比 β 多一个电子。

### 8.3 带 ECP 的重元素

```python
config = SCFConfig(
    method="pbe0",
    basis="def2-svp",
    ecp="def2-ecp",
)
```

### 8.4 密度拟合（近似 J/K）

```python
config = SCFConfig(
    method="r2scan",
    basis="def2-svp",
    density_fit=True,
    auxbasis="def2-jkfit",
)
```

结果元数据中 `jk_evaluation == "density-fitted"`。

### 8.5 带 D4 色散

```python
config = SCFConfig(
    method="r2scan",
    basis="def2-svp",
    dispersion="d4",
)
```

色散能量在只含真实原子的物理片段几何上计算（ghost 不进入 D3/D4 后端），
输出中 `Dispersion Interaction` 是超分子与片段的色散差。

### 8.6 隐式溶剂（PCM/SMD）

```python
config = SCFConfig(
    method="pbe",
    basis="def2-svp",
    solvent="cpcm",            # 或 "iefpcm" / "ssvpe" / "cosmo" / "smd"
    solvent_name="water",
)
result = PySCFEDA(atoms, fragments, config).run()
values = result.components("kcal/mol")
print(values["Desolvation"])
print(result.diagnostics["solvent"])
```

命令行：

```bash
pyscf-dm-eda water_dimer.xyz --fragment 1-3 --fragment 4-6 \
  --method pbe --basis def2-svp \
  --solvent smd --solvent-name water \
  --solvent-option lebedev_order=29
```

发生了什么：

1. 每个片段 SCF 和超分子 SCF 都附加了 PySCF 的溶剂模型（`mf.PCM()` /
   `mf.SMD()`）。
2. 片段的空腔只由该片段的**真实原子**生成。PySCF 默认会给 ghost 原子一个
   2.0 Å 的球，程序在构造表面、表面核势和 SMD 的 CDS 项时临时换成只含真实
   原子的分子，AO 积分仍用完整 ghost 基组，因此 counterpoise 约定不变。
3. 静电、交换、排斥、极化、相关项照旧用气相算符作用在溶剂化后的密度上；
   反应场能量（SMD 含 CDS）从相关残差中扣除，并作为
   `Desolvation = E_solv(super) - Σ_A E_solv(A)` 单独报告。HF 下
   `Correlation Interaction` 仍为 0，可用来检验溶剂能量是否被完整分离。
4. `diagnostics["solvent"]` 中还给出超分子空腔里 `P0` 与 `P_Pauli` 的反应场
   能量，把 `Desolvation` 拆成 `desolvation_frozen_hartree`
   （空腔合并与相互屏蔽）、`desolvation_pauli_response_hartree` 和
   `desolvation_polarization_response_hartree`（三者之和等于
   `Desolvation`），以及每个空腔的表面点数 `cavity_surface_points`。
5. `metadata["solvent"]` 记录模型、有效介电常数、溶剂名和额外选项。

注意：

- `Desolvation` 通常为正（形成复合物会损失一部分溶剂化能），同时静电项会因
  片段在溶剂中被极化而变得更负，`Orbital Relaxation` 会变小，因为一部分极化
  响应转移到了溶剂响应中；
- SMD 的表面点数约为 C-PCM 的两倍，费用相应更高；可用
  `solvent_options={"lebedev_order": 17}` 之类的设置换取速度；
- ddCOSMO/ddPCM 不支持；没有显式溶剂、非平衡溶剂化或热力学校正。

`examples/water_dimer_solvent_eda.py` 给出气相、C-PCM 和 SMD 的对照表。

---

## 9. 批量网格计算：PySCFGridRunner

`PySCFGridRunner` **不生成网格**，只消费调用方准备好的文件。它把每个组合结构固定
分成“母体 + 探针”两个片段；任意多片段体系请用 `PySCFEDA` 单结构 API。

### 9.1 默认目录约定

```text
work/
├── molecule.xyz                # 母体（前 N 个原子）
├── molecule_filtered.xyz       # 格点文件
└── molecule_probe/
    ├── mol_probe_0.xyz         # 母体原子 + 探针原子的组合 XYZ
    ├── mol_probe_1.xyz
    └── ...
```

### 9.2 三个输入文件的确切格式

**`molecule.xyz`**：标准 XYZ。

**`molecule_filtered.xyz`**：第一行是点数，第二行注释，之后每行至少 4 列：

```text
3
H grid around O: label X Y Z [extra columns...]
H  0.000000  0.000000  1.900000
H  0.850000  0.000000  1.472243
```

- 第 1 列是标签（可以是任意 token，本实现只检查长度）；
- 第 2–4 列是格点坐标；
- 之后的额外列允许携带 `r, theta, phi` 等元数据。

**`mol_probe_<i>.xyz`**：前 N 个原子必须逐元素、逐坐标与 `molecule.xyz`
一致；其余原子属于探针。

### 9.3 单原子探针

单原子探针直接与格点坐标核对，无需额外参数：

```python
runner = PySCFGridRunner(
    molecule_xyz="work/molecule.xyz",
    molecule_charge=0,
    molecule_spin=0,
    probe_charge=1,
    probe_spin=0,
    config=SCFConfig(method="r2scan", basis="def2-svp"),
).run(restart=True)
```

### 9.4 多原子探针：anchor 或 template

多原子探针必须选择一种无歧义校验方式。

**方式 A：`probe_anchor_atom`**——格点就是探针片段的第 `k` 个原子
（零基索引）：

```python
PySCFGridRunner(
    ...,
    probe_anchor_atom=0,   # 探针第 0 个原子位于格点
).run(restart=True)
```

**方式 B：`probe_template_xyz`**——提供以格点为原点的原始探针模板。
程序会检查元素顺序、内部几何、刚体配准 RMSD 和模板原点：

```python
PySCFGridRunner(
    ...,
    probe_template_xyz="work/probe.xyz",
).run(restart=True)
```

> 线形多原子探针绕其分子轴的 Kabsch 旋转不唯一；若模板原点不在分子轴上，
> 请使用 `probe_anchor_atom`，否则可能出现虚假校验失败。

### 9.5 Python 调用

```python
from pyscf_dm_eda import PySCFGridRunner, SCFConfig

summary = PySCFGridRunner(
    molecule_xyz="work/molecule.xyz",
    molecule_charge=0,
    molecule_spin=0,
    probe_charge=0,
    probe_spin=0,
    probe_directory="work/molecule_probe",
    combined_pattern="mol_probe_{index:03d}.xyz",
    grid_xyz="work/molecule_filtered.xyz",
    energy_output="work/eda.tsv",
    xyz_output="work/eda_grid.xyz",
    config=SCFConfig(method="r2scan", basis="def2-svp"),
    probe_anchor_atom=1,
).run(restart=True)

print(summary.completed_points)
print(summary.energy_table)
print(summary.extended_xyz)
print(summary.metadata_file)
```

### 9.6 命令行调用

```bash
pyscf-dm-eda grid work/molecule.xyz \
  --molecule-charge 0 --molecule-spin 0 \
  --probe-charge 0 --probe-spin 0 \
  --probe-anchor-atom 1 \
  --method r2scan --basis def2-svp \
  --progress --restart
```

`--indices` 用于分片，例如 `--indices 0-99` 或 `--indices 0-49,120`。
索引是零基闭区间。

### 9.7 输出文件

- 能量表：默认 `energy_values_pyscf_<stem>.tsv`，固定 kcal/mol；
- 扩展 XYZ：默认 `<stem>_pyscf.xyz`，每行 = 格点行 + 16 个 EDA 数值；
- 元数据：与 TSV 同名的 `.json`；
- 锁文件：`<tsv>.lock`，进程结束后保留，不需要删除。

用 pandas 读取时跳过 `#` 指纹行：

```python
import pandas as pd

table = pd.read_csv("eda.tsv", sep="\t", comment="#")
```

---

## 10. 网格检查点、分片与重启

### 10.1 检查点行为

- 每个格点完成后，TSV 通过临时文件 + `fsync` + 原子替换提交；
- TSV 第一行是 `# pyscf-dm-eda-grid-v1 fingerprint=...`，第二行是列名；
- JSON sidecar 保存代码、PySCF 版本、设置和所有输入文件 SHA-256；
- 任一格点 SCF 失败会抛异常，**不会写伪造的全零结果**；
- `restart=True` 会跳过已完成格点；`restart=False` 初始化并覆盖同名检查点。

### 10.2 fingerprint 的两个注意点

1. fingerprint 包含 `eda.py` 的 SHA-256 和输入/工作目录的**绝对路径**。
   改代码（哪怕只改注释）或移动目录都会使旧检查点失效。这是有意为之的安全取舍。
2. 修改 `molecule.xyz`、格点文件或组合 XYZ 后，旧检查点同样失效。

### 10.3 分片并行

每个分片必须写不同的 `energy_output`（锁按输出文件生效）。可在多进程/多机并行：

```python
shard = PySCFGridRunner(
    ...,
    indices=range(0, 100),
    energy_output="work/eda_000.tsv",
    progress_callback=lambda i, total: print(f"{i}/{total}", flush=True),
).run(restart=True)
```

CLI 对应：

```bash
pyscf-dm-eda grid work/molecule.xyz \
  --molecule-charge 0 --molecule-spin 0 \
  --probe-charge 1 --probe-spin 0 \
  --indices 0-99 \
  --energy-output work/eda_000.tsv \
  --progress --restart
```

多个进程并发时建议设置单线程 BLAS（每个分片一个进程）：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
pyscf-dm-eda grid ... --indices 100-199 --energy-output work/eda_001.tsv --restart
```

### 10.4 性能提示

检查点每完成一点就全量重写已完成行，I/O 按点数平方增长。数千点以内没问题；
数万点请用 `--indices` 分片，避免单输出文件尾部过长。

---

## 11. 完整实战：球形水二聚体扫描

仓库中提供了一个端到端示例，演示“生成球面网格 → PySCFGridRunner 分片计算 →
合并 → 作图 → 写 Markdown”的完整流程：

```text
examples/
├── water_dimer_spherical_eda.py
└── water_dimer_spherical/
    ├── molecule.xyz / molecule_filtered.xyz
    ├── molecule_probe/mol_probe_*.xyz
    ├── energy_values_pyscf_molecule.tsv
    ├── water_dimer_eda_vs_distance.png
    ├── water_dimer_eda_vs_angle.png
    └── water_dimer_spherical_eda.md
```

运行：

```bash
cd /mnt/d/pyscf-dm-eda
WATER_DIMER_EDA_NPROC=8 ~/edaenv/bin/python examples/water_dimer_spherical_eda.py
```

脚本中值得照着改的部分：

- `build_grid_specs()`：构造 `(index, r, theta, phi)` 列表；
- `donor_water()` / `acceptor_water()`：生成探针与母体 XYZ；
- `generate_inputs()`：写出三个输入文件；
- `_grid_worker()`：每个分片调用 `PySCFGridRunner`，`probe_anchor_atom=1`
  表示把探针水的成键 H 当作格点锚原子；
- `run_grid()`：`multiprocessing.Pool` 并行分片，再合并各分片 TSV；
- `load_table()` / `plot_vs_distance()` / `plot_vs_angle()`：读表作图；
- `write_markdown()`：生成含数据表和图的结果 Markdown。

想改成自己的分子，替换 `acceptor_water()` 与 `donor_water()` 即可；若探针不是水，
把 `probe_anchor_atom` 改成实际锚原子索引，或改用 `probe_template_xyz`。

---

## 12. 故障排查

### 12.1 `SCFConvergenceError: SCF did not converge ...`

按顺序尝试：

1. `max_cycle=200` 或更大；
2. `level_shift=0.2`、`damp=0.7`；
3. 换初猜，如 `init_guess="atom"`；
4. 确认 `newton_fallback=True`（默认开启）；
5. 先用小基组/低网格找到波函数，再换目标基组；
6. 对 DFT 提高 `grid_level` 有时也能改善数值噪声。

### 12.2 `EDAValidationError`

读异常携带的结果：

```python
except EDAValidationError as exc:
    print(exc.result.diagnostics["validation_worst_case"])
    print(exc.result.diagnostics["validation_errors"])
```

- `closure` / `frozen identity` 大：通常不是 SCF 问题，应检查是否误用了
  `density_fit` 或非常弥散基组；
- `core partition` 大：片段定义与超分子不一致；
- `electron count` / idempotency 大：SCF 结果质量差；
- 弥散/近线性相关基组可临时 `strict_validation=False` 或放大
  `validation_tol`，但要明确这是自行放宽诊断。

### 12.3 `IncompatibleFragmentError`

常见原因：

- 片段电荷和 ≠ 超分子电荷；
- 片段自旋和 ≠ 超分子自旋；
- Python 索引从 0 开始但写成了 1 基索引；
- 重复/遗漏原子；
- AO 基组近线性相关：调大 `linear_dep_threshold` 或换基组。

### 12.4 `ValueError: Pass the base XC functional ...`

把 `method="r2scan-d4"` 改成 `method="r2scan", dispersion="d4"`。

### 12.5 `wb97x-d / b97-d` 被拒绝

这两个泛函的参数化本身包含经验色散，而本实现会禁用 `mf.disp`。请使用其底层
XC 泛函并显式设置 `dispersion`。

### 12.6 网格重启时报 `Checkpoint fingerprint does not match`

- 移动了工作目录；
- 修改了 `eda.py`；
- 修改了任一输入 XYZ；
- 改了 `combined_pattern`、格点数、`indices` 等设置。

处理：确认要开始新任务时删除旧 TSV/JSON 并用 `restart=False`；否则恢复原始
输入路径与文件。

### 12.7 `Another grid runner holds the output lock`

另一个进程正在写同一个 `energy_output`。分片任务请务必使用不同输出文件。

### 12.8 多原子探针校验失败

- 组合 XYZ 中母体前 N 个原子与 `molecule.xyz` 不完全一致；
- 没有传 `probe_anchor_atom` 或 `probe_template_xyz`；
- 线形探针用了模板方式且模板原点不在分子轴上——改 `probe_anchor_atom`。

---

## 13. 限制与注意事项

- 隐式溶剂仅限 PySCF 的 PCM 族与 SMD，以单一 `Desolvation` 项报告；没有
  ddCOSMO/ddPCM、显式溶剂、非平衡溶剂化或热力学校正；
- 不支持分数占据、smearing、ROHF/ROKS 和非单行列式自旋耦合；
- `Orbital Relaxation` 包含无法唯一分离的 polarization / charge-transfer；
- `density_fit=True` 是 DF-J/K 近似，不是 exact J/K；
- `Correlation Interaction` 是所选泛函的代数残差；其化学解释依赖泛函选择；
- Mulliken 电荷转移对基组敏感，不要当作 NBO/NPA；
- DM-EDA 的定义与其他 EDA scheme 的同名分项不逐项可比；
- PySCF 无 Windows 原生轮子，Windows 请使用 WSL2。

---

## 参考

1. Density-matrix EDA: *J. Chem. Phys.* **160**, 174101 (2024),
   DOI: 10.1063/5.0202787.
2. 公式推导与代码对应：`MATHEMATICAL_PRINCIPLES.md`。
3. 项目概览与批量网格约定：`README.md`。
