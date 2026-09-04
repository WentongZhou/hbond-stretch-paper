# PySCF DM-EDA

一个独立、可安装的 PySCF 密度矩阵能量分解分析（DM-EDA）实现。

项目直接使用 PySCF 的 AO 密度矩阵、重叠矩阵和 J/K 构造器计算能量分项，
不调用或解析其他量化程序，也不依赖任何特定的网格生成、可视化或上层应用。

完整公式、推导和代码对应关系见
[MATHEMATICAL_PRINCIPLES.md](MATHEMATICAL_PRINCIPLES.md)。

## 功能概览

- 2024 density-matrix EDA 的 exact-J/K 定义；
- 完整超分子 AO 空间中的 ghost/counterpoise 片段计算；
- RHF、UHF、RKS、UKS，以及 broken-symmetry 初猜；
- 静电、交换、Pauli 排斥、极化、相关残差和经验色散分项；
- 可选 PySCF 隐式溶剂（C-PCM/IEF-PCM/SS(V)PE/COSMO/SMD）：片段用真实原子
  空腔溶剂化，反应场能量变化作为 `Desolvation` 分项；
- 可选 PySCF density fitting，并在结果中明确标记为近似 J/K；
- JSON 单结构输出；
- 可重启的批量网格计算，带输入哈希、原子写入和进程锁；
- Mulliken 片段电荷转移及严格的能量闭合、电子数和密度矩阵诊断。

## 项目结构

```text
pyscf-dm-eda/
├── pyproject.toml
├── README.md
├── MATHEMATICAL_PRINCIPLES.md
├── LICENSE
├── .github/
│   └── workflows/
│       └── ci.yml
├── examples/
│   └── he2.xyz
├── src/pyscf_dm_eda/
│   ├── __init__.py
│   ├── __main__.py
│   ├── eda.py
│   └── py.typed
└── tests/
    └── test_eda.py
```

发行包、Python 包和命令行入口分别是：

```text
pyscf-dm-eda       # distribution / CLI
pyscf_dm_eda       # Python import
```

## 安装

推荐使用 Python 3.10 或更高版本，并在独立虚拟环境中安装：

```bash
cd pyscf-dm-eda
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

如果需要 D3/D4 色散：

```bash
python -m pip install -e ".[dispersion]"
```

> PySCF 没有原生 Windows 轮子，代码中的 Windows 文件锁分支无法在原生
> Windows 上安装/测试。Windows 用户请使用 WSL2（本项目的开发与 CI 均采用
> Linux）。

## 快速开始

### 命令行

CLI 中原子索引从 **1** 开始。下面把 He₂ 分成两个单原子片段：

```bash
pyscf-dm-eda examples/he2.xyz \
  --fragment 1 --fragment 2 \
  --method hf --basis sto-3g \
  --unit kcal/mol \
  --output he2_eda.json
```

也可以不安装 console script，直接运行模块：

```bash
python -m pyscf_dm_eda examples/he2.xyz \
  --fragment 1 --fragment 2 \
  --method r2scan --basis def2-svp
```

多原子片段可用逗号和闭区间表示，例如 `--fragment 1-3,7`。每个片段的
电荷、自旋和标签按相同顺序重复传入：

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

`--unit` 支持 `hartree`、`kcal/mol`、`kj/mol` 和 `ev`。SCF 选项
`--ecp`、`--max-memory`、`--level-shift`、`--damp`、`--init-guess`、
`--linear-dep-threshold`、`--density-fit`、`--auxbasis`、
`--validation-tol` 与 `--strict-validation/--no-strict-validation` 均已暴露，
隐式溶剂选项见下文[隐式溶剂](#隐式溶剂)，`pyscf-dm-eda --help` 可查看完整列表。参数错误仍由 argparse 以 exit 2 报告；
SCF/数值运行期错误打印到 stderr 并返回 exit 1。

批量网格计算也有命令行入口（`pyscf-dm-eda grid --help`），见下文
[批量网格计算](#批量网格计算)。

### Python API

Python API 中原子索引从 **0** 开始：

```python
from pyscf_dm_eda import FragmentSpec, PySCFEDA, SCFConfig

result = PySCFEDA.from_xyz(
    "complex.xyz",
    fragments=[
        FragmentSpec((0, 1, 2), charge=0, spin=0, label="molecule"),
        FragmentSpec((3,), charge=1, spin=0, label="probe"),
    ],
    config=SCFConfig(
        method="r2scan",
        basis="def2-svp",
        grid_level=4,
        density_fit=False,
        dispersion=None,
    ),
    charge=1,
    spin=0,
).run()

print(result.components("kcal/mol"))
result.write_json("complex_eda.json")
```

## 自旋约定

`spin` 使用 PySCF 的带符号定义：

```text
spin = N_alpha - N_beta
```

它不是多重度 `2S+1`。总体系的 charge 和 spin 必须分别等于所有片段的和。
两个自由基形成 broken-symmetry singlet 时，可给片段相反的自旋投影：

```python
fragments = [
    FragmentSpec((0,), spin=1, label="alpha"),
    FragmentSpec((1,), spin=-1, label="beta"),
]
```

只要任一片段开壳层，自动模式会对超分子使用 UHF/UKS；restricted 与
unrestricted 超分子 SCF 都会以片段 promolecule 密度作为初猜（restricted
情形取 α+β 总密度），以保留片段起点、加快收敛。也可以通过
`SCFConfig(unrestricted=True)` 或 CLI 的 `--unrestricted` 强制使用
unrestricted reference。

## 方法定义

所有片段先在完整超分子 AO 基组中独立计算，其他片段的原子作为 ghost：

```text
P0       = sum_A P_A
P_Pauli  = C0 (C0^T S C0)^(-1) C0^T   # alpha/beta 分别构造
E0[P]    = Tr(P h) + 1/2 Tr(P J[P])
           - 1/2 sum_sigma Tr(P_sigma K[P_sigma])
```

能量分项为：

```text
E_electrostatic = E_nuc-nuc + E_electron-nuclear + E_Coulomb
E_exchange      = fragment-pair exact exchange
E_repulsion     = E0[P_Pauli] - E0[P0]
E_polarization  = E0[P_super] - E0[P_Pauli]
E_correlation   = supermolecule-fragment DFT/HF residual difference
E_dispersion    = E_disp(super) - sum_A E_disp(fragment A)
E_desolvation   = E_solv(super) - sum_A E_solv(fragment A)   # 仅启用隐式溶剂时
```

`Correlation Interaction` 为去除显式 D3/D4（以及反应场能量）后的 residual；
输出中的 `Corr_Disp = Correlation Interaction + Dispersion Interaction` 对应
论文把经验色散包含在内的广义相关项。

## 输出

`EDAResult.components()` 返回以下主分量：

| 键 | 含义 |
|---|---|
| `Total Interaction energy` | counterpoise 总相互作用能 |
| `Electrostatic Interaction` | 核-核、交叉核-电子与交叉 Coulomb 之和 |
| `Exchange Int.` | 片段间 exact exchange |
| `Repulsion` | Pauli 正交化代价 |
| `Exchange-Repulsion` | exchange + repulsion |
| `Orbital Relaxation` | 终态密度相对 Pauli 态的极化/弛豫 |
| `Correlation Interaction` | 不含显式色散的泛函 residual |
| `Dispersion Interaction` | D3/D4 相互作用差值（若启用） |
| `Desolvation` | PCM/SMD 反应场能量变化（若启用隐式溶剂；否则为 0） |
| `Closure Error` | 直接相互作用能与分项和的差 |

JSON 还包含：

- 各片段 `Mulliken` 电荷转移；
- SCF 收敛状态；
- AO overlap 条件数；
- 电子数、幂等性和 Hermitian 误差；
- 逐项 `validation_errors`、最差校验项与阈值；
- 计算方法、基组、PySCF 版本和已知限制。

默认 `strict_validation=True`：任一校验残差超过 `validation_tol` 时抛出
`EDAValidationError`，异常对象携带完整 `EDAResult`，可直接检查所有分项和
诊断信息而无需重算。设置 `strict_validation=False`（CLI 为
`--no-strict-validation`）则降级为 `RuntimeWarning` 并返回该结果；这适合
弥散/近线性相关基组下需要自行判断的情形。

项目不把 Mulliken 电荷标记为 NBO/NPA；PySCF 核心没有 Weinhold NBO 实现。

## D3/D4 色散

色散通过可选的 `pyscf-dispersion` 后端计算：

```python
config = SCFConfig(
    method="r2scan",
    basis="def2-svp",
    dispersion="d4",
)
```

请把基础泛函写在 `method`，把 `d3`、`d3bj` 或 `d4` 单独写在
`dispersion`；不要把 `-d3/-d4` 后缀混入方法名。`wb97x-d`、`b97-d` 这类
**参数化本身包含经验色散**的泛函会被拒绝：libxc 只提供其 XC 部分，而程序
会禁用 `mf.disp`，静默运行会丢掉色散尾项，得到不完整的泛函。片段色散只在
包含真实原子的物理几何上求值，ghost 原子不会传给 D3/D4 后端。

## 隐式溶剂

`SCFConfig(solvent=...)` 把 PySCF 自带的连续介质模型附加到每个片段 SCF 和
超分子 SCF 上：

```python
config = SCFConfig(
    method="pbe",
    basis="def2-svp",
    solvent="cpcm",            # cpcm(=pcm) | iefpcm | ssvpe | cosmo | smd
    solvent_name="water",      # 来自 PySCF SMD 数据库；PCM 只取其介电常数
    # solvent_eps=78.3553,     # 直接给介电常数（PCM）；SMD 需配合 solvent_name
    # solvent_options={"lebedev_order": 29, "vdw_scale": 1.2},
)
```

命令行对应 `--solvent`、`--solvent-eps`、`--solvent-name` 和可重复的
`--solvent-option KEY=VALUE`：

```bash
pyscf-dm-eda water_dimer.xyz --fragment 1-3 --fragment 4-6 \
  --method pbe --basis def2-svp --solvent smd --solvent-name water
```

定义与约定：

- 片段在**只含自身真实原子**的空腔中溶剂化。PySCF 会给 ghost 原子分配
  2.0 Å 的默认半径，程序在构造空腔、表面核势和 SMD 的 CDS 项时改用真实原子
  分子，而 AO 积分仍在完整 ghost 基组中完成；超分子使用全部原子的空腔。
- 静电、交换、排斥、极化、相关项继续用气相算符作用在溶剂化密度上求值；
  反应场能量（SMD 含 CDS）从 residual 中去掉并单独报告为
  `Desolvation = E_solv(super) - Σ_A E_solv(A)`，与 Su 等人 LMO-EDA/GKS-EDA
  溶液版本中 desolvation 项的定义方式一致。HF 下 `Correlation Interaction`
  仍严格为零。
- `diagnostics["solvent"]` 额外给出在超分子空腔中求得的 `E_solv[P0]`、
  `E_solv[P_Pauli]`，把 `Desolvation` 拆成 frozen、Pauli 响应和轨道弛豫响应三
  部分（三者之和等于 `Desolvation`），以及各空腔的表面点数。
- ddCOSMO/ddPCM 无法把空腔限制到真实原子，会被拒绝；不含显式溶剂、非平衡
  溶剂化或热力学校正。

`examples/water_dimer_solvent_eda.py` 对同一水二聚体给出气相、C-PCM 和 SMD
的对照表。

## 批量网格计算

`PySCFGridRunner` 不生成网格或探针位置，只消费调用方准备好的文件。默认约定：

```text
work/
├── molecule.xyz
├── molecule_filtered.xyz
└── molecule_probe/
    ├── mol_probe_0.xyz
    ├── mol_probe_1.xyz
    └── ...
```

`molecule_filtered.xyz` 的每一行至少需要：

```text
label  X  Y  Z
```

每个组合 XYZ 的前 N 个原子必须与 `molecule.xyz` 完全对应，剩余原子属于探针。
网格 runner 固定把结构分成“母体 + 探针”两个片段；任意多片段体系请使用
`PySCFEDA` 单结构 API。
目录、网格文件、组合文件命名和输出路径均可显式配置：

```python
from pyscf_dm_eda import PySCFGridRunner, SCFConfig

summary = PySCFGridRunner(
    molecule_xyz="work/molecule.xyz",
    molecule_charge=0,
    molecule_spin=0,
    probe_charge=1,
    probe_spin=0,
    probe_directory="work/placed_probes",
    grid_xyz="work/points.xyz",
    combined_pattern="complex_{index:05d}.xyz",
    energy_output="work/eda.tsv",
    xyz_output="work/eda_grid.xyz",
    config=SCFConfig(method="r2scan", basis="def2-svp"),
).run(restart=True)

print(summary.completed_points)
```

同一任务也可以用 CLI 子命令运行：

```bash
pyscf-dm-eda grid work/molecule.xyz \
  --molecule-charge 0 --molecule-spin 0 \
  --probe-charge 1 --probe-spin 0 \
  --method r2scan --basis def2-svp \
  --progress --restart
```

单原子探针会直接与网格坐标核对。多原子探针必须选择一种无歧义校验方式：

```python
# 网格点就是探针片段中的第 0 个原子
probe_anchor_atom=0

# 或提供生成组合结构时使用的、以网格锚点为原点的原始模板
probe_template_xyz="probe.xyz"
```

模板方式会检查元素顺序、内部几何、刚体配准 RMSD 和模板原点位置，因此也
支持锚点处没有真实原子的探针。注意：**线形多原子探针**绕分子轴的旋转使
Kabsch 旋转不唯一，如果模板原点不在分子轴上，平移量会随退化方向变化并可能
产生虚假校验失败；此类情形请使用 `probe_anchor_atom`。

### 分片与进度

`indices` 接受零基索引子集，`progress_callback` 在每个格点成功写入检查点后
回调 `(index, total_selected)`。每个分片必须写不同的 `energy_output`
（锁按输出文件生效），可在多进程/多机间并行：

```python
shard = PySCFGridRunner(
    ...,
    indices=range(0, 100),
    energy_output="work/eda_000.tsv",
    progress_callback=lambda i, total: print(f"{i}/{total}", flush=True),
).run(restart=True)
```

CLI 对应 `--indices 0-99` 与 `--progress`。

### 检查点与重启安全

- 每完成一个格点，TSV 都通过临时文件、`fsync` 和原子替换提交；
- TSV 首行携带任务 fingerprint，第二行是列名；
- JSON sidecar 保存代码、PySCF、设置和所有输入文件的 SHA-256；
- restart 同时核对 JSON 与 TSV fingerprint；
- advisory lock 阻止两个进程写入同一个输出；
- SCF 失败会抛出异常，不会写入伪造的全零结果。

两点使用注意事项：

1. fingerprint 包含 `eda.py` 的 SHA-256 和输入/工作目录的绝对路径。修改
   代码（即使只改注释）或移动工作目录会使既有检查点失效，这是有意为之的
   安全性取舍。
2. 当前检查点每完成一点就全量重写已完成行，I/O 按点数平方增长。数千点
   以内没有问题；数万点请用 `indices` 分片，避免单输出文件的检查点提交
   尾部过长。

批量能量固定以 kcal/mol 写入，列顺序为：

```text
Grid_Index Tot Electro Nuc_Nuc 1e 2e Exc_Rep Exc Rep
Orb_Relax Corr Disp Corr_Disp Steric Mulliken_CT Closure_Error
```

默认输出为 `energy_values_pyscf_<stem>.tsv`、同名 `.json` metadata、
`<stem>_pyscf.xyz` 和 `.tsv.lock`。`restart=False` 会初始化并覆盖同名
checkpoint；只有确认要开始新任务时才使用它。

用 pandas 读取 TSV 时，应忽略以 `#` 开头的指纹行：

```python
import pandas as pd

table = pd.read_csv("eda.tsv", sep="\t", comment="#")
```

`.lock` 文件是持久的锁载体；进程结束后锁会自动释放，不需要删除该文件。

## 测试

先以 editable 模式安装，再运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m unittest discover -s tests -v
```

当前测试覆盖：

- HF 与 DFT 能量闭合；
- 静电及 exchange-repulsion 子恒等式；
- HF 零 correlation residual；
- density-fitting 路径与 DF 标记；
- mock 后端下的 D3/D4 差值逻辑；
- opposite-spin broken-symmetry reference（UHF 与 UKS）；
- 单位转换（hartree/kcal/mol/kJ/mol/eV）与非法输入；
- XYZ/JSON 原子写出；
- 校验失败携带完整 `EDAResult`，以及 `strict_validation=False` 降级；
- CLI 单结构端到端、运行期错误退出码、`grid` 子命令端到端；
- fragment/index 解析边界（反向区间、越界、重复）；
- 网格检查点、内容指纹、分片 indices、进度回调、模拟崩溃；
- 并发 writer 锁；
- 多原子模板刚体配准与错配拒绝。

CI 配置见 `.github/workflows/ci.yml`（Linux + Python 3.10/3.11/3.12）。

## 限制

- 隐式溶剂仅限 PySCF 的 PCM 族与 SMD，并以单一 `Desolvation` 项报告；没有
  ddCOSMO/ddPCM、显式溶剂或非平衡溶剂化；
- 不支持分数占据、smearing、ROHF/ROKS 和非单行列式自旋耦合；
- `Orbital Relaxation` 包含不能唯一分离的 polarization/charge-transfer 效应；
- `density_fit=True` 是 PySCF DF-J/K 近似，不是 exact-J/K 结果；
- ECP、真实 D3/D4 后端及大规模多片段体系仍需要更多参考数据回归；
- 这是 DM-EDA 实现，不应假定与其他 EDA scheme 的同名分项逐项相等；
- PySCF 无 Windows 原生轮子，Windows 请使用 WSL2。

## 参考文献

1. Density-matrix EDA: *J. Chem. Phys.* **160**, 174101 (2024),
   [DOI: 10.1063/5.0202787](https://doi.org/10.1063/5.0202787).
2. Generalized Kohn-Sham EDA: *J. Phys. Chem. A* (2014),
   [DOI: 10.1021/jp500405s](https://doi.org/10.1021/jp500405s).

## 许可证

本项目采用 [MIT License](LICENSE)。
