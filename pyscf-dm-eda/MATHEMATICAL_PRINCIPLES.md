# PySCF DM-EDA 的数学原理

本文说明 `pyscf_dm_eda.eda` 实际计算的数学对象、能量分解恒等式以及代码中的
数值检查。所有公式均使用原子单位；实现的主要理论依据是 2024 年提出的
density-matrix energy decomposition analysis（DM-EDA）。

程序内部以 Hartree（$E_h$）保存能量；输出单位换算常数为
$1\ E_h=627.5094740631\ \mathrm{kcal\,mol^{-1}}=2625.4996394799\ \mathrm{kJ\,mol^{-1}}=27.211386245988\ \mathrm{eV}$。

## 1. 适用范围

给定由若干片段组成的超分子，程序计算

$$
\Delta E_{\mathrm{int}}
=E_{\mathrm{super}}-\sum_A E_A
$$

并将其分解为

$$
\boxed{
\Delta E_{\mathrm{int}}
=\Delta E_{\mathrm{ele}}
+\Delta E_{\mathrm{ex}}
+\Delta E_{\mathrm{rep}}
+\Delta E_{\mathrm{pol}}
+\Delta E_{\mathrm{corr}}
+\Delta E_{\mathrm{disp}}
}
$$

其中 `corr` 是不含显式 D3/D4 的泛函 residual。若按 2024 论文把经验色散并入
广义相关项，则应使用

$$
\Delta E_{\mathrm{corr}}^{\mathrm{paper}}
=\Delta E_{\mathrm{corr}}+\Delta E_{\mathrm{disp}},
$$

也就是程序输出中的 `Corr_Disp`。

该分解建立在单行列式 HF/GKS 密度矩阵上。它不是任意其他 EDA 方法中同名
分量的数值定义，也不能仅凭分项名称推断跨程序逐项相等。

## 2. 记号与 AO 度量

使用非正交原子轨道基组 $\{\chi_\mu\}$：

$$
S_{\mu\nu}=\langle\chi_\mu|\chi_\nu\rangle.
$$

自旋密度矩阵记为 $P^\alpha$ 和 $P^\beta$，总密度为

$$
P=P^\alpha+P^\beta.
$$

电子数不是普通矩阵迹，而是 AO 度量下的迹：

$$
N_\sigma=\operatorname{Tr}(P^\sigma S),\qquad
N=N_\alpha+N_\beta.
$$

这里电子数需要 $S$，而一电子能写作 $\operatorname{Tr}(Ph)$、不额外乘 $S$。
两者并不矛盾：$P$ 是在非正交 AO 基中定义的密度矩阵，恒等算符的 AO 矩阵表示
引入重叠度量 $S$，而一般算符的期望值直接由其 AO 积分矩阵与 $P$ 收缩。

本文使用 $\dagger$ 表示共轭转置。当前实际计算通常使用实轨道，但代码仍按
Hermitian 矩阵形式实现。

对限制性参考态，PySCF 给出的总密度被拆成

$$
P^\alpha=P^\beta=\frac{1}{2}P.
$$

对 UHF/UKS，两个自旋通道分别保留。

## 3. Ghost/counterpoise 片段空间

每个片段都在完整超分子 AO 基组中计算。属于片段 $A$ 的原子保留真实核电荷；
其他片段的原子改为 ghost atoms：

- ghost 原子不提供核电荷和核吸引势；
- ghost 原子的基函数仍然存在；
- 所有片段与超分子因此具有相同的 AO 数目、顺序、$S$ 和动能矩阵 $T$。

这使不同片段的密度矩阵可以直接相加，并给出 counterpoise 定义的片段能量。
因此程序的 $\Delta E_{\mathrm{int}}$ 是 CP/ghost 基组约定下的相互作用能；它与
各片段只使用自身基函数的非 CP 结果通常不同。

所有片段均保持在输入复合物中的几何位置，所以这里得到的是**固定几何相互作用能**。
若要计算结合能或解离能，还需另行加入片段几何弛豫、零点能、热修正等；这些
deformation/preparation 项不属于本脚本的分解。

令片段 $A$ 的核吸引算符为 $V_A$，则

$$
h_A=T+V_A,
$$

而超分子的一电子算符为

$$
h_{\mathrm{super}}=T+\sum_A V_A.
$$

代码显式检查

$$
\left\|h_{\mathrm{super}}-T-\sum_A V_A\right\|_{\max}
$$

是否小于数值容差，防止片段 AO/ECP/核势定义不一致。

使用 ECP 时，$V_A$ 包含相应的有效芯势；片段与超分子必须采用完全一致的 ECP
和基组约定。

## 4. HF-like 密度矩阵泛函 $E_0[P]$

DM-EDA 的分阶段能量使用下列 HF-like 泛函：

$$
\boxed{
E_0[P^\alpha,P^\beta;h]
=\operatorname{Tr}(Ph)
+\frac{1}{2}\operatorname{Tr}\!\left(PJ[P]\right)
-\frac{1}{2}\sum_{\sigma\in\{\alpha,\beta\}}
\operatorname{Tr}\!\left(P^\sigma K[P^\sigma]\right)
}
$$

其中 $J[P]$ 是总密度产生的 Coulomb 矩阵，$K[P^\sigma]$ 是同自旋交换矩阵。
核-核排斥能不包含在 $E_0$ 中，需要单独加入。

注意 $E_0$ 的交换部分始终采用完整 HF-like exact exchange。对 DFT 而言，所选
交换-相关泛函与该参考形式之间的差额随后进入 residual，而不是改变这里的系数。

在 AO 指标下，可写成

$$
J[P]_{\mu\nu}
=\sum_{\lambda\kappa}(\mu\nu|\lambda\kappa)P_{\kappa\lambda},
$$

$$
K[P^\sigma]_{\mu\nu}
=\sum_{\lambda\kappa}(\mu\lambda|\nu\kappa)
P^\sigma_{\kappa\lambda}.
$$

默认情况下，程序通过 PySCF `get_jk` 使用 exact J/K。若
`density_fit=True`，相同公式使用 PySCF density-fitted J/K；能量恒等式仍在
同一近似内部闭合，但结果必须理解为 DF-J/K 近似，而不是 exact-J/K。

## 5. Promolecule 密度 $P^0$

每个片段先独立完成 SCF。片段自旋密度的简单和定义未反对称化的 promolecule：

$$
P^{0,\sigma}=\sum_A P_A^\sigma,
\qquad
P^0=P^{0,\alpha}+P^{0,\beta}.
$$

片段内部占据轨道已经正交，但不同片段的占据轨道一般彼此重叠，因此 $P^0$ 通常
不满足 AO 度量下的幂等关系

$$
P^{0,\sigma}SP^{0,\sigma}=P^{0,\sigma}.
$$

$P^0$ 描述尚未施加片段间 Pauli 反对称化的冻结密度叠加。
其非幂等性通常正是片段占据空间相互重叠的预期结果，因此程序不把 $P^0$ 的
幂等性作为验收条件；幂等性检查针对 Pauli 密度和最终单行列式密度。

## 6. Pauli 密度的构造

对每个自旋通道，把所有片段的占据轨道系数拼接：

$$
C_0^\sigma=
\left[C_{1,\mathrm{occ}}^\sigma\;
C_{2,\mathrm{occ}}^\sigma\;\cdots\right].
$$

这些列在 AO 度量中的 Gram 矩阵为

$$
G^\sigma=(C_0^\sigma)^\dagger S C_0^\sigma.
$$

对占据空间作对称正交化：

$$
\widetilde C_0^\sigma
=C_0^\sigma(G^\sigma)^{-1/2}.
$$

于是

$$
(\widetilde C_0^\sigma)^\dagger
S\widetilde C_0^\sigma=I,
$$

Pauli 密度为

$$
\boxed{
P^{\mathrm{Pauli},\sigma}
=C_0^\sigma(G^\sigma)^{-1}(C_0^\sigma)^\dagger
}
$$

并满足

$$
P^{\mathrm{Pauli},\sigma}S
P^{\mathrm{Pauli},\sigma}
=P^{\mathrm{Pauli},\sigma}.
$$

若 $G^\sigma$ 的最小特征值相对于最大特征值过小，说明片段占据空间近线性相关，
程序会停止而不是对病态矩阵求逆。

## 7. 冻结静电与交换作用

### 7.1 核-核项

$$
\Delta E_{\mathrm{nuc}}
=E_{\mathrm{nuc}}^{\mathrm{super}}
-\sum_A E_{\mathrm{nuc}}^A.
$$

它等价于不同片段原子核之间的 Coulomb 排斥和。

### 7.2 交叉核-电子项

$$
\Delta E_{\mathrm{en}}
=\sum_A\sum_{B\ne A}
\operatorname{Tr}(P_A V_B).
$$

每一项表示片段 $A$ 的电子密度与另一片段 $B$ 的核吸引势之间的作用。

### 7.3 片段间 Coulomb 项

$$
\Delta E_J
=\sum_{A<B}\operatorname{Tr}\!\left(P_AJ[P_B]\right).
$$

从 $\frac12\operatorname{Tr}(PJ[P])$ 展开时，$A,B$ 与 $B,A$ 两个交叉项由
$1/2$ 抵消，因此每个无序片段对只计一次。

等价地，若采用论文中的有序求和，可写成

$$
\Delta E_J
=\frac{1}{2}\sum_{A\ne B}\operatorname{Tr}\!\left(P_AJ[P_B]\right).
$$

### 7.4 静电总项

$$
\boxed{
\Delta E_{\mathrm{ele}}
=\Delta E_{\mathrm{nuc}}
+\Delta E_{\mathrm{en}}
+\Delta E_J
}
$$

对应输出中的：

```text
Electrostatic Interaction = Nuc---Nuc + 1-electron + 2-electron
```

### 7.5 片段间 exact exchange

$$
\boxed{
\Delta E_{\mathrm{ex}}
=-\sum_{A<B}\sum_\sigma
\operatorname{Tr}\!\left(P_A^\sigma K[P_B^\sigma]\right)
}
$$

交换只发生在相同自旋通道。与 Coulomb 项相同，代码对每个无序片段对只计算
一次；上式已经包含从 HF-like 泛函展开后得到的正确系数。

对应的有序求和形式为

$$
\Delta E_{\mathrm{ex}}
=-\frac{1}{2}\sum_{A\ne B}\sum_\sigma
\operatorname{Tr}\!\left(P_A^\sigma K[P_B^\sigma]\right).
$$

冻结态恒等式为

$$
\begin{aligned}
&E_0[P^0;h_{\mathrm{super}}]
+E_{\mathrm{nuc}}^{\mathrm{super}}
-\sum_A\left(E_0[P_A;h_A]+E_{\mathrm{nuc}}^A\right)\\
&\qquad=\Delta E_{\mathrm{ele}}+\Delta E_{\mathrm{ex}}.
\end{aligned}
$$

程序把等式左右之差记录为 `frozen_identity_error_hartree`。
由于完整相关项本身按 residual 定义，总闭合可以被这种定义代数保证；因此
`frozen_identity_error_hartree` 对交叉静电/交换项的系数、符号和片段势分割错误
尤其敏感，是比总闭合更直接的冻结项实现诊断。

## 8. Pauli 排斥与轨道极化

施加片段间反对称化但尚未进行最终 SCF 弛豫时，能量变化定义为 Pauli repulsion：

$$
\boxed{
\Delta E_{\mathrm{rep}}
=E_0[P^{\mathrm{Pauli}};h_{\mathrm{super}}]
-E_0[P^0;h_{\mathrm{super}}]
}
$$

“repulsion” 是该 EDA 阶段的约定名称；由于它是两个非自洽密度的 $E_0$ 差，
数学上并没有对任意体系都严格为正的定理保证。

程序另外提供组合量

$$
\Delta E_{\mathrm{ex-rep}}
=\Delta E_{\mathrm{ex}}+\Delta E_{\mathrm{rep}}.
$$

令超分子最终收敛密度为 $P^S$。从 Pauli 态到最终密度的 HF-like 能量降低为

$$
\boxed{
\Delta E_{\mathrm{pol}}
=E_0[P^S;h_{\mathrm{super}}]
-E_0[P^{\mathrm{Pauli}};h_{\mathrm{super}}]
}
$$

输出名称为 `Orbital Relaxation`。这个量包含片段内极化、片段间电荷流动以及它们
的耦合；本密度矩阵定义没有把 polarization 和 charge transfer 唯一分开。

## 9. HF/DFT residual 与色散

PySCF SCF 能量（不含额外经验色散）记为 $E_X^{\mathrm{SCF}}$，其中
$X$ 可以是超分子或某个片段。定义状态 residual：

$$
R_X
=E_X^{\mathrm{SCF}}
-E_{\mathrm{nuc}}^X
-E_0[P_X;h_X].
$$

于是

$$
\boxed{
\Delta E_{\mathrm{corr}}
=R_{\mathrm{super}}-\sum_A R_A
}
$$

对于不使用 density fitting 的 HF，$E^{\mathrm{SCF}}-E_{\mathrm{nuc}}$ 与
$E_0$ 是同一泛函，因此 $R_X$ 和 $\Delta E_{\mathrm{corr}}$ 应在数值误差内为零。

对 KS/GKS 方法，$R_X$ 收集 $E_0$ 未包含的泛函贡献。因为 $E_0$ 采用完整
exact-exchange 形式，residual 会吸收所选纯泛函、杂化泛函、meta-GGA 或
range-separated 泛函相对于该参考形式的差异。分解在代数上闭合，但 residual
的化学解释仍依赖具体泛函。

因此 `Correlation Interaction` 是**广义泛函 residual**，不能无条件解释为
波函数理论中的纯电子相关能。2024 论文对普通纯泛函/全局杂化泛函给出的单一
混合参数解释可直接对应；对 range-separated、meta-GGA 或其他更复杂的
orbital-dependent 形式，当前实现仍能代数闭合，但需要额外理论验证后才能赋予
同样的逐项物理解释。

若启用 D3/D4，程序把它作为后 SCF 的加和修正：

$$
\boxed{
\Delta E_{\mathrm{disp}}
=E_{\mathrm{disp}}^{\mathrm{super}}
-\sum_A E_{\mathrm{disp}}^A
}
$$

电子结构片段仍使用 ghost/counterpoise 基组；但 D3/D4 只在每个片段的真实原子
几何上计算，不把 ghost 原子传给经验色散模型。

此外，`density_fit=True` 仅表示使用 PySCF 的 DF-J/K 数值近似；它不是 2024
论文中 RI/COSX 半数值 `DM-EDA(LR)` 方案的实现。

### 9.1 隐式溶剂（PCM/SMD）与去溶剂化项

`SCFConfig(solvent=...)` 把 PySCF 的连续介质模型（C-PCM、IEF-PCM、SS(V)PE、
COSMO 或 SMD）附加到每一个片段 SCF 和超分子 SCF 上。对状态 $X$，PySCF 的
SCF 能量变为

$$
E_X^{\mathrm{SCF}}
=E_{\mathrm{nuc}}^X+E_{\mathrm{elec}}[P_X;h_X]+E_{\mathrm{solv}}^X[P_X],
$$

其中反应场能量

$$
E_{\mathrm{solv}}^X[P]
=\tfrac12\,\mathbf q_X[P]^{\mathsf T}\,\mathbf v_X[P]
\;\left(+\,E_{\mathrm{CDS}}^X\ \text{仅 SMD}\right).
$$

$\mathbf v_X[P]$ 是溶质（核 + 密度 $P$）在空腔表面离散点上的静电势，
$\mathbf q_X[P]$ 是由 PCM 线性方程 $\mathbf K\mathbf q=\mathbf R\mathbf v$
得到的表观表面电荷（PySCF 采用 SWIG 离散化）。SMD 的 CDS
（cavitation–dispersion–solvent structure）项只依赖几何和溶剂描述符，与
密度无关。

**空腔约定。** PySCF 为 `mol` 中的每个原子生成一个球面，ghost 原子会拿到
半径表中 $Z=0$ 的默认值（2.0 Å）。因此程序在构造片段的空腔、表面核势
$\mathbf v_{\mathrm{nuc}}$ 和 CDS 项时临时改用只含该片段真实原子的分子对象，
而 AO 积分（$\mathbf v[P]$ 与反应场矩阵 $V_{\mathrm{solv}}$）始终在完整 ghost
基组中完成。于是每个片段在**自己的空腔**中溶剂化，超分子在**全部原子的空腔**
中溶剂化；这正是定义溶液中相互作用能
$\Delta E_{\mathrm{int}}^{\mathrm{sol}}=E_{AB}^{\mathrm{sol}}-E_A^{\mathrm{sol}}-E_B^{\mathrm{sol}}$
时的通常约定。

**分解定义。** 第 7、8 节的全部项继续用气相算符 $h$、$J$、$K$ 作用在
**溶剂化的**密度 $P_A$、$P^0$、$P^{\mathrm{Pauli}}$、$P^S$ 上求值；residual
改为

$$
R_X=E_X^{\mathrm{SCF}}-E_{\mathrm{nuc}}^X-E_{\mathrm{solv}}^X[P_X]-E_0[P_X;h_X],
$$

而反应场能量的变化单独构成去溶剂化项

$$
\boxed{
\Delta E_{\mathrm{desolv}}
=E_{\mathrm{solv}}^{\mathrm{super}}[P^S]
-\sum_A E_{\mathrm{solv}}^A[P_A]
}
$$

这与 Su 等人 LMO-EDA/GKS-EDA 溶液版本中 desolvation 项的定义方式一致：静电、
交换、排斥、极化和相关项描述溶剂化片段之间的直接相互作用，
$\Delta E_{\mathrm{desolv}}$ 收集形成复合物时损失（或获得）的溶剂化能。HF 下
$R_X$ 仍严格为零，这是检验 $E_{\mathrm{solv}}$ 被完整分离出来的直接判据。

**诊断性拆分。** 程序还在超分子空腔中求
$E_{\mathrm{solv}}^{\mathrm{super}}[P^0]$ 与
$E_{\mathrm{solv}}^{\mathrm{super}}[P^{\mathrm{Pauli}}]$（SMD 时加同一个 CDS），
把去溶剂化项沿 DM-EDA 的密度序列拆为

$$
\Delta E_{\mathrm{desolv}}
=\underbrace{E_{\mathrm{solv}}^{\mathrm{super}}[P^0]-\sum_A E_{\mathrm{solv}}^A[P_A]}_{\text{frozen：空腔合并与相互屏蔽}}
+\underbrace{E_{\mathrm{solv}}^{\mathrm{super}}[P^{\mathrm{Pauli}}]-E_{\mathrm{solv}}^{\mathrm{super}}[P^0]}_{\text{Pauli 响应}}
+\underbrace{E_{\mathrm{solv}}^{\mathrm{super}}[P^S]-E_{\mathrm{solv}}^{\mathrm{super}}[P^{\mathrm{Pauli}}]}_{\text{轨道弛豫响应}} .
$$

三项之和代数上等于 $\Delta E_{\mathrm{desolv}}$，保存在
`diagnostics["solvent"]` 中，不作为独立分项进入闭合式。

**边界。** ddCOSMO/ddPCM 按 `mol` 的每个原子离散空腔并在原子中心 DFT 网格上
积分密度，无法把空腔限制到真实原子，因此被拒绝。溶剂模型只影响 SCF 密度与
$E_{\mathrm{solv}}$；D3/D4 与 CDS 一样只在真实原子几何上求值。`Desolvation`
是连续介质的电子能变化，不含显式溶剂、热校正或非平衡溶剂化。

## 10. 总能量闭合

程序直接计算的总相互作用能为

$$
\Delta E_{\mathrm{int}}
=\left(E_{\mathrm{super}}^{\mathrm{SCF}}
+E_{\mathrm{disp}}^{\mathrm{super}}\right)
-\sum_A\left(E_A^{\mathrm{SCF}}+E_{\mathrm{disp}}^A\right).
$$

启用隐式溶剂时，$E^{\mathrm{SCF}}$ 已经包含第 9.1 节的 $E_{\mathrm{solv}}$，
因此 $\Delta E_{\mathrm{int}}$ 就是溶液中的（counterpoise）相互作用能。

分项和为

$$
E_{\mathrm{sum}}
=\Delta E_{\mathrm{ele}}
+\Delta E_{\mathrm{ex}}
+\Delta E_{\mathrm{rep}}
+\Delta E_{\mathrm{pol}}
+\Delta E_{\mathrm{corr}}
+\Delta E_{\mathrm{disp}}
+\Delta E_{\mathrm{desolv}},
$$

其中 $\Delta E_{\mathrm{desolv}}$ 在气相计算中恒为零。

闭合误差定义为

$$
\boxed{
\varepsilon_{\mathrm{closure}}
=\Delta E_{\mathrm{int}}-E_{\mathrm{sum}}
}
$$

并以 `closure_error_hartree` 保存。默认情况下，若它或其他核心数值恒等式的误差
超过 `SCFConfig.validation_tol`，计算会抛出 `EDAValidationError`；异常对象携带
完整的 `EDAResult`（全部分项、`validation_errors` 和诊断信息），可直接检查而
无需重算。`SCFConfig(strict_validation=False)` 会把校验失败降级为
`RuntimeWarning` 并返回结果，适合弥散/近线性相关基组下需要自行判断的情形。

同一个 `validation_tol` 被用于若干不同量级和含义的数值残差，这是程序的统一
防护阈值，不是物理误差条、实验不确定度或 EDA 分项精度估计；`validation_errors`
中的逐项值可用于分项诊断。

需要强调：相关项被定义为 residual，因此极小的闭合误差证明实现内部代数一致，
但不能单独证明所选泛函、基组或 EDA scheme 对某个化学问题足够准确。

## 11. 开壳层与 broken-symmetry

PySCF 的 `spin` 定义为

$$
m=N_\alpha-N_\beta,
$$

不是多重度。给定总电子数 $N$：

$$
N_\alpha=\frac{N+m}{2},\qquad
N_\beta=\frac{N-m}{2}.
$$

程序要求

$$
Q_{\mathrm{super}}=\sum_A Q_A,
\qquad
m_{\mathrm{super}}=\sum_A m_A,
$$

并进一步检查 alpha/beta 电子数分别守恒。

对于由两个自由基形成的 broken-symmetry singlet，可以指定

$$
m_1=+1,\qquad m_2=-1,\qquad m_{\mathrm{super}}=0.
$$

即使超分子总 spin 为零，只要任一片段开壳层，自动模式仍使用 UHF/UKS；
restricted 与 unrestricted 超分子 SCF 均以片段 $P^0$ 作为初猜（restricted 情形
取 $P^0_\alpha+P^0_\beta$），从而保留片段起点和 broken-symmetry 起点。

broken-symmetry 行列式通常不是 $\hat S^2$ 的本征函数，结果可能依赖所收敛的
局域极小值并受到自旋污染影响。当前脚本不做自旋投影，也不把
$\langle\hat S^2\rangle$ 作为 EDA 输出；解释开壳层结果时应另行检查 SCF 解。

## 12. Mulliken 片段电荷转移

最终超分子密度的 Mulliken AO population 为

$$
n_\mu=(PS)_{\mu\mu}.
$$

片段 $A$ 的电子布居为

$$
N_A^{\mathrm{Mulliken}}
=\sum_{\mu\in A}(PS)_{\mu\mu}.
$$

其最终片段电荷为

$$
q_A^{\mathrm{final}}
=\sum_{a\in A}Z_a-N_A^{\mathrm{Mulliken}},
$$

程序报告

$$
\boxed{
\Delta q_A=q_A^{\mathrm{final}}-q_A^{\mathrm{formal}}
}
$$

正值表示片段比形式电荷更正，即净失电子；负值表示净得电子。在电荷守恒时
$\sum_A\Delta q_A=0$。

该量只是一种基组依赖的 Mulliken 指标，不是独立能量分项，也不是 NBO/NPA
电荷或严格的 charge-transfer energy。

若使用 ECP，上式中的 $Z_a$ 是 PySCF 采用的有效核电荷而非裸核电荷，所以
Mulliken 电荷也应在同一赝势/价电子约定下解释。

## 13. 数值稳定性与验证

代码除 SCF 收敛外，还验证下列量。

### 13.1 AO 线性相关

对 $S$ 的特征值检查

$$
\lambda_{\min}(S)
>\tau_{\mathrm{lin}}\lambda_{\max}(S).
$$

近线性相关的 ghost/弥散基组会在 SCF 前被拒绝。

### 13.2 电子数

$$
\operatorname{Tr}(P^\sigma S)=N_\sigma
$$

分别对片段和、Pauli 密度及最终超分子密度检查。

### 13.3 Hermiticity

$$
\left\|P^\sigma-(P^\sigma)^\dagger\right\|_{\max}.
$$

### 13.4 度量幂等性

$$
\left\|P^\sigma SP^\sigma-P^\sigma\right\|_{\max}.
$$

该检查适用于整数占据的单行列式密度。分数占据和 smearing 因而不在当前保证
范围内。

### 13.5 两个能量恒等式

- `frozen_identity_error_hartree`：冻结 promolecule 的静电+交换恒等式误差；
- `closure_error_hartree`：完整 EDA 分项和的闭合误差。

这些量只验证数值和代数自洽，不是与论文参考程序、Turbomole、实验值一致性的
替代测试。

## 14. 代码与公式的对应关系

### 14.1 输出量的依赖关系

结果表同时包含基本分项、基本分项的子项和便于观察的派生组合。只有下表中标为
“是”的七个基本分项各计一次，才能组成第 10 节的闭合式；**不能把输出表所有
能量列直接相加**。

| 输出名称 | 定义或角色 | 作为独立项加入闭合式 |
|---|---|:---:|
| `Nuc---Nuc` | $\Delta E_{\mathrm{nuc}}$，静电子项的子项 | 否 |
| `1-electron` | $\Delta E_{\mathrm{en}}$，静电子项的子项 | 否 |
| `2-electron` | $\Delta E_J$，静电子项的子项 | 否 |
| `Electrostatic Interaction` | 上述三个子项之和 $\Delta E_{\mathrm{ele}}$ | 是 |
| `Exchange Int.` | $\Delta E_{\mathrm{ex}}$ | 是 |
| `Repulsion` | $\Delta E_{\mathrm{rep}}$ | 是 |
| `Exchange-Repulsion` | `Exchange Int.` + `Repulsion` | 否 |
| `Orbital Relaxation` | $\Delta E_{\mathrm{pol}}$ | 是 |
| `Correlation Interaction` | $\Delta E_{\mathrm{corr}}$ | 是 |
| `Dispersion Interaction` | $\Delta E_{\mathrm{disp}}$ | 是 |
| `Desolvation` | $\Delta E_{\mathrm{desolv}}$（第 9.1 节；气相为 0） | 是 |
| `Corr_Disp` | `Correlation Interaction` + `Dispersion Interaction` | 否 |
| `Steric` | `Exchange-Repulsion` + `Corr_Disp` | 否 |
| `Total Interaction energy` | 直接 SCF/色散能差；闭合式的目标值 | 结果 |
| `Closure Error` | 目标值减去七个基本分项之和 | 诊断 |
| `Mulliken_CT` | Mulliken 电荷变化，不是能量 | 不适用 |

### 14.2 实现位置

| 代码位置 | 数学作用 |
|---|---|
| `PySCFEDA.run` | 构造 CP 片段、片段 SCF、$P^0$ 初猜和超分子 SCF |
| `PySCFEDA._decompose` | 计算全部交叉项、阶段能量、residual 和闭合误差 |
| `_build_pauli_density` | 构造 $G^\sigma$ 与 $P^{\mathrm{Pauli},\sigma}$ |
| `_e0` | 计算 HF-like 泛函 $E_0[P]$ |
| `_batch_jk` | 在统一超分子 AO 空间批量生成 J/K |
| `_mulliken_fragment_charge_transfer` | 计算 $\Delta q_A$ |
| `_electron_count` | 计算 $\operatorname{Tr}(P^\sigma S)$ |
| `_idempotency_error` | 检查 $P^\sigma SP^\sigma=P^\sigma$ |
| `_solvent_model` | 构造以真实原子为空腔、在 ghost 基组中积分的 PySCF PCM/SMD 对象 |
| `_scf_solvent_energies` | 计算收敛密度的 $E_{\mathrm{solv}}^X[P_X]$（含 SMD 的 CDS） |
| `_solvent_diagnostics` | 在超分子空腔中求 $E_{\mathrm{solv}}[P^0]$、$E_{\mathrm{solv}}[P^{\mathrm{Pauli}}]$ 并拆分 $\Delta E_{\mathrm{desolv}}$ |
| `EDAResult.components` | 单位转换与 `Corr_Disp`、`Steric` 派生量 |

`Steric` 定义为

$$
\mathrm{Steric}
=\Delta E_{\mathrm{ex-rep}}
+\Delta E_{\mathrm{corr}}
+\Delta E_{\mathrm{disp}}.
$$

它与 `Corr_Disp` 一样是便于输出和分析的派生组合，不是新增的独立理论阶段，
也不应视为跨所有 EDA scheme 通用且可独立比较的“位阻能”。

批量网格运行器只是对每个组合结构重复完全相同的单结构 DM-EDA；它不改变任何
数学定义。

## 15. 算法摘要

```text
for each fragment A:
    build A with all other atoms as ghosts
    if a solvent model is requested:
        attach PySCF PCM/SMD with a cavity built from A's real atoms only
    run fragment SCF in the common AO space
    collect P_A^alpha, P_A^beta, occupied orbitals and E_solv[P_A]

P0 = sum_A P_A
choose restricted/unrestricted supermolecule reference
run supermolecule SCF, using P0 as the unrestricted initial guess when needed

construct P_Pauli separately for alpha and beta
batch-build J/K for every fragment, P_Pauli and P_super
evaluate electrostatic and interfragment exchange terms
evaluate repulsion and polarization from E0 differences
evaluate the HF/DFT residual (reaction field removed) and optional real-atom D3/D4 difference
evaluate desolvation = E_solv[P_super] - sum_A E_solv[P_A] when a solvent is attached
verify frozen identity, electron counts, idempotency and total closure
```

## 16. 理论与实现边界

- 所有片段和终态必须可由整数占据的单行列式参考描述；
- 当前没有 ROHF/ROKS、分数占据、smearing 或多参考自旋耦合；
- 当前没有 MP2/CC 相关修正、double-hybrid 或 SCF-MI；
- 隐式溶剂只支持 PySCF 的 PCM 族（C-PCM/IEF-PCM/SS(V)PE/COSMO）与 SMD，并以
  单一 `Desolvation` 项进入闭合式；没有 ddCOSMO/ddPCM、显式溶剂、非平衡溶剂化
  或热力学校正；
- `Orbital Relaxation` 不提供唯一的 polarization/CT 再分解；
- Mulliken 电荷转移不能替代 NBO/NPA；
- DFT residual 的数值和解释依赖泛函；
- density-fitted 与 exact-J/K 是不同数值层级；
- CP/ghost 与非 CP 片段定义会给出不同相互作用能；
- 内部恒等式闭合不等同于对其他 EDA 程序的 golden regression。

本程序直接实现的是 2024 年 AO density-matrix DM-EDA 的 $P^0$、
$P^{\mathrm{Pauli}}$ 和 $P^S$ 路径。2014 年 GKS-EDA 的逐级中间行列式方案仅作为
概念背景：同名分项不代表定义或数值相等。尤其不能据此宣称本程序是 Turbomole
EDA、GAMESS GKS-EDA、LMO-EDA 或 SAPT 的数值克隆。

## 17. 参考文献

1. Density-matrix EDA，*J. Chem. Phys.* **160**, 174101 (2024)，
   [DOI: 10.1063/5.0202787](https://doi.org/10.1063/5.0202787)。
2. Generalized Kohn-Sham EDA 背景，*J. Phys. Chem. A* (2014)，
   [DOI: 10.1021/jp500405s](https://doi.org/10.1021/jp500405s)。
3. 溶液中的 LMO-EDA 与 desolvation 项：P. Su, H. Liu, W. Wu,
   *J. Chem. Phys.* **137**, 034111 (2012)。
4. PCM 综述：J. Tomasi, B. Mennucci, R. Cammi, *Chem. Rev.* **105**, 2999
   (2005)，[DOI: 10.1021/cr9904009](https://doi.org/10.1021/cr9904009)。
5. SMD：A. V. Marenich, C. J. Cramer, D. G. Truhlar, *J. Phys. Chem. B*
   **113**, 6378 (2009)，
   [DOI: 10.1021/jp810292n](https://doi.org/10.1021/jp810292n)。
6. PySCF PCM 所用的 SWIG 表面离散化：A. W. Lange, J. M. Herbert,
   *J. Chem. Phys.* **133**, 244111 (2010)，
   [DOI: 10.1063/1.3511297](https://doi.org/10.1063/1.3511297)。
