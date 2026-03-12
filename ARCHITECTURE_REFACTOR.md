# LatticeQSL 架构重构与当前框架说明

## 本次修改概览

这次修改的目标，是把原本较为分散的 ED 相关实现整理成一个更清晰、可扩展的模块化框架，同时保留现有代码的兼容入口，便于后续逐步迁移算法与模型。

本轮已经完成的核心工作包括：

- 将主包按 `Lattice / Basis / Hamiltonian / Solver / Observable` 五个层次拆分。
- 为有限尺寸晶格、平移对称、反射对称、Sz 守恒建立统一的数据结构。
- 建立基于 `BitBasis.BitStr` 的自旋 `1/2` 基底表示。
- 用稀疏矩阵构造 ED Hamiltonian，避免 Kronecker 积直接展开。
- 接入平移约化基底与 1D 平移 + 反射约化基底。
- 实现约化基底态向父基底的展开接口，便于观测量计算。
- 实现约化密度矩阵与 von Neumann 纠缠熵。
- 补充了一组面向新架构的基础测试，并已通过 `Pkg.test()`。

---

## 当前主模块结构

`src/LatticeQSL.jl` 现在负责统一 include 和 export，整体结构如下：

```text
LatticeQSL
├── Lattice
│   ├── AbstractLattice.jl
│   ├── PredefinedLattices.jl
│   ├── LatticeUtils.jl
│   └── FiniteLattice.jl
├── Basis
│   ├── SpinBasis.jl
│   ├── Symmetry.jl
│   ├── SpinSymmetry.jl
│   ├── TranslationSymmetry.jl
│   ├── PointGroupSymmetry.jl
│   └── SymmetryReducedBasis.jl
├── Hamiltonian
│   ├── SpinOperators.jl
│   ├── AbstractModel.jl
│   ├── Models.jl
│   └── HamiltonianBuilder.jl
├── Solver
│   ├── ExactDiag.jl
│   └── Lanczos.jl
├── Observable
│   ├── GroundState.jl
│   ├── Correlation.jl
│   ├── StructureFactor.jl
│   ├── Entanglement.jl
│   └── WilsonLoop.jl
├── EDdemo.jl
└── Basis.jl
```

其中：

- 新架构代码在上述子目录内。
- `EDdemo.jl` 和旧版 `Basis.jl` 仍作为兼容层保留。
- 后续可以继续把旧逻辑逐步并入新模块，而不必一次性重写全部接口。

---

## 各模块当前职责

### 1. `Lattice/`

负责定义晶格的几何结构与有限簇信息。

#### 已实现内容

- `AbstractLattice{D}`：统一的晶格抽象接口。
- `Lattice`：通用 Bravais 晶格类型。
- 预定义晶格：
  - `HoneycombLattice`
  - `SquareLattice`
  - `TriangularLattice`
  - `ChainLattice`
  - `KagomeLattice`
  - `LiebLattice`
  - `RectangularLattice`
  - `GeneralLattice`
- `FiniteLattice`：
  - 有限尺寸 cluster
  - PBC/OBC
  - site 编号
  - bond 列表
  - 邻接查询
  - 平移群排列 `translation_group`
- `LatticeUtils`：
  - `generate_sites`
  - `AtomList`
  - `MaskedGrid`
  - 坐标平移/缩放/裁剪等工具

#### 当前定位

`Lattice/` 解决的是“几何与索引”的问题，不涉及 Hilbert space 和 Hamiltonian 本身。

---

### 2. `Basis/`

负责 Hilbert space 的表示、量子数扇区、以及对称约化。

#### 已实现内容

- `SpinBasis`：
  - 基于 `BitBasis.BitStr{N, Int}` 的自旋 `1/2` 基底
  - 支持全空间基底
  - 支持固定 `nup` 的 `Sz` 守恒扇区
  - 提供 `state_index`、`flip_spin`、`permute_bits` 等工具
- `Symmetry.jl`：
  - `AbstractSymmetry`
  - `AbstractSpatialSymmetry`
  - `AbstractInternalSymmetry`
- `SpinSymmetry.jl`：
  - `SzConservation`
  - `SpinInversion`
- `TranslationSymmetry.jl`：
  - 平移群对象
  - 动量点列表 `momentum_values`
- `PointGroupSymmetry.jl`：
  - 当前是占位接口，后续用于旋转/镜像/点群不可约表示
- `SymmetryReducedBasis.jl`：
  - 平移约化基底
  - 1D 平移 + 反射约化基底
  - representative 映射
  - phase / normalization 管理
  - `matrix_element_factor`
  - `expand_state`

#### 本轮重点补充

本轮把原来远程分支里 `iso_total2K / iso_K2MSS` 的思路，整理进了 `SymmetryReducedBasis` 统一框架：

- 先按平移轨道找 representative。
- 用 character 选择目标动量扇区。
- 在 1D 情况下进一步用 `BitBasis.breflect` 做反射约化。
- 约化态与父基底之间可通过 `expand_state` 显式互相联系。

#### 当前定位

`Basis/` 解决的是“哪些量子态存在、如何标号、如何按对称性降维”的问题。

---

### 3. `Hamiltonian/`

负责定义模型、局域算符动作，以及稀疏 Hamiltonian 的装配。

#### 已实现内容

- `SpinOperators.jl`：
  - 单点算符：`Sx/Sy/Sz/Sp/Sm`
  - 双点算符：`SzSz`, `SpSm`, `SmSp`
  - 组合项：`apply_heisenberg`, `apply_kitaev_bond`
- `AbstractModel.jl`：
  - 模型统一抽象接口 `AbstractModel`
- `Models.jl`：
  - `HeisenbergModel`
  - `KitaevModel`
  - `IsingModel`
  - `XXZModel`
  - `GeneralSpinModel`
- `HamiltonianBuilder.jl`：
  - 在 `SpinBasis` 上构造稀疏 Hamiltonian
  - 在 `SymmetryReducedBasis` 上构造稀疏 Hamiltonian
  - 提供 `basis_summary`

#### 当前设计思想

- 不直接搭张量积矩阵。
- 对每个基底态逐项施加局域作用。
- 把结果写入 `SparseMatrixCSC`。
- 对约化基底，通过 `matrix_element_factor` 自动补上轨道相位和归一化因子。

#### 当前定位

`Hamiltonian/` 解决的是“模型如何作用在基底上，以及如何生成矩阵”的问题。

---

### 4. `Solver/`

负责谱求解。

#### 已实现内容

- `ExactDiag.jl`
  - `exact_diagonalization`
  - `ground_state_exact`
- `Lanczos.jl`
  - `ground_state`
  - `low_energy_spectrum`

#### 当前定位

- 小体系：完整对角化
- 稀疏大体系：Arpack / Lanczos

这层保持尽量薄，只负责“拿矩阵求谱”。

---

### 5. `Observable/`

负责从本征态或波函数中计算物理量。

#### 已实现内容

- `GroundState.jl`
  - `expectation`
  - `energy_per_site`
- `Correlation.jl`
  - `szsz_correlation`
- `Entanglement.jl`
  - `reduced_density_matrix`
  - `entanglement_entropy`
- `StructureFactor.jl`
  - 当前占位
- `WilsonLoop.jl`
  - 当前占位

#### 本轮重点补充

将原先 `rdm_PXP / ee` 的核心思路，整理成更通用的接口：

- 支持直接在 `SpinBasis` 上计算约化密度矩阵。
- 支持先从 `SymmetryReducedBasis` 展开回父基底，再做 partial trace。
- 统一提供 `entanglement_entropy` 计算 von Neumann entropy。

#### 当前定位

`Observable/` 解决的是“态已知后，如何计算能量、关联、纠缠等量”的问题。

---

## 当前核心数据流

现在一条典型的 ED 工作流可以概括为：

```text
Lattice
  -> FiniteLattice
  -> SpinBasis / SymmetryReducedBasis
  -> Model
  -> build_hamiltonian
  -> exact_diagonalization / ground_state
  -> Observable
```

更具体地说：

1. 先定义晶格几何，例如 `FiniteLattice(HoneycombLattice(), (Lx, Ly))`。
2. 再定义 Hilbert space，例如 `SpinBasis(N, nup)`。
3. 如果需要动量或反射约化，再构造 `SymmetryReducedBasis`。
4. 用 `HeisenbergModel` / `KitaevModel` / `IsingModel` 等模型生成 Hamiltonian。
5. 用 `ExactDiag` 或 `Lanczos` 求基态/低能谱。
6. 用 `Observable` 模块计算关联函数、纠缠熵等。

---

## 当前已经验证的内容

本轮新增/更新测试位于 `test/architecture.jl`，目前已经验证：

- Honeycomb 有限晶格构造正常。
- `SpinBasis` 全空间与固定 `nup` 构造正常。
- 2-site Heisenberg chain 的谱正确。
- 平移约化基底能成功构造，并满足正交性测试。
- 1D 平移 + 反射约化基底能成功构造。
- 受限基底上的约化密度矩阵与纠缠熵结果正确。

测试命令：

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

当前结果：测试已全部通过。

---

## 当前框架的优点

### 1. 模块边界更清晰

几何、基底、模型、求解、观测量已分层，后续扩展不会再把逻辑混在单个文件里。

### 2. 易于扩展新晶格与新模型

- 新晶格：主要补 `PredefinedLattices` + `FiniteLattice` 的 bond 构造。
- 新模型：主要补 `AbstractModel` 的实现。
- 新观测量：直接在 `Observable/` 下扩展。

### 3. 方便做对称约化

目前已经把：

- `Sz`
- translation
- 1D reflection

纳入了统一框架，后面继续接 point group 会更自然。

### 4. 与旧代码兼容过渡

`EDdemo.jl` 和旧 `Basis.jl` 仍保留，因此可以逐步迁移，而不是一次性推倒。

---

## 当前仍是占位或待扩展的部分

下面这些模块已经留出位置，但还没有完全展开：

- `PointGroupSymmetry.jl`
- `Observable/StructureFactor.jl`
- `Observable/WilsonLoop.jl`

另外，当前 `SymmetryReducedBasis` 主要完成了：

- 平移约化
- 1D 平移 + 反射约化

如果后续要做 2D 点群不可约表示分块，还需要继续把具体群作用和 character 表接进去。

---

## 建议的下一步

从当前框架出发，后续最自然的推进方向有：

1. 把指定扇区的 Heisenberg / Kitaev 求谱接口封装成更高层 API。
2. 把 point group symmetry 真正接入 `SymmetryReducedBasis`。
3. 为 `StructureFactor` 和 `WilsonLoop` 填充真实实现。
4. 增加更多 2D cluster 的回归测试。
5. 逐步把旧 `EDdemo.jl` 中仍有用的逻辑迁移进新模块。

---

## 总结

当前 `LatticeQSL` 已经从“若干功能文件并列”的形态，进入了“可持续扩展的 ED 包框架”阶段。

本轮最关键的成果是：

- 建立了清晰的模块边界；
- 把几何、基底、模型、求解、观测量贯通起来；
- 接上了平移/反射约化和纠缠熵；
- 并且通过了基础测试。

这意味着后续继续补 2D 自旋模型、更多对称性、更多 observable 时，已经有比较稳的骨架可沿用。
