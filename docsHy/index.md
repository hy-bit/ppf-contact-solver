# ZOZO Contact Solver 有限元理论文档

> 本文档库整理了本项目（ppf-contact-solver / ZOZO's Contact Solver）中**实际使用**的有限元与接触动力学核心理论。所有公式、模型、算法均严格来源于项目源码实现、配套文章与主论文描述，便于中文读者理解、维护与贡献。

**项目定位**：一个 GPU 加速（单精度 CUDA）的混合维度（壳/体/杆）隐式有限元接触求解器。核心创新包括：
- 仅依赖奇异值导数的闭式正定本征雅可比（支持 3×3 体与 3×2 壳）
- 弹性包含的动态接触刚度（elasticity-inclusive dynamic stiffness）
- 自定义 cubic 势垒 + IPC 风格的穿透-free 接触处理
- 高效的应变限制、塑性、摩擦与多材料混合

**当前状态**：文档基于 2026 年 6 月代码快照（主分支持续演进中）。与 `articles/eigensys.md`（详尽推导）互补；主论文为 ACM TOG 2024 "A Cubic Barrier with Elasticity-Inclusive Dynamic Stiffness"。

**快速导航**
- [主论文与材料](#主论文与技术材料)
- [文档结构](#文档结构)
- [核心架构概览](#核心架构概览-mermaid)
- [源码与文章交叉引用](#源码与文章交叉引用)
- [如何贡献/更新本理论文档](#如何贡献更新本理论文档)

## 主论文与技术材料

- 📘 **A Cubic Barrier with Elasticity-Inclusive Dynamic Stiffness**（ACM TOG Vol.43, No.6）
  - 主论文 PDF 与补充材料见仓库根 `README.md` 的 [Technical Materials 节](../README.md#-%F0%9F%8E%93-technical-materials)
  - 参考分支：`sigasia-2024`（与论文一致，较少维护）
- 🔍 奇异值本征分析详解：[`articles/eigensys.md`](../articles/eigensys.md)（含 hindsight 提及 Poya et al. 2023 先前工作）
- 🧐 事后说明与修正：[`articles/hindsight.md`](../articles/hindsight.md)（quadratic 势垒曲率问题、ACCD 浮点误差、倾斜角度勘误）
- 🐞 Bug 修复与更新：[`articles/bug.md`](../articles/bug.md)（应变限制线搜索改进等）

## 文档结构

1. [01-有限元离散化与变形梯度.md](./01-有限元离散化与变形梯度.md)  
   元素类型（壳/体/杆）、变形梯度 F 计算（`F = dx * inv_rest`）、SVD 预处理、inv_rest 意义、rest 形状。

2. [02-超弹性势能模型.md](./02-超弹性势能模型.md)  
   ARAP / StVK / SNHk / BaraffWitkin / Hook + `detsqr` 体积项 + 弯曲能（dihedral/strand）。每模型精确 Ψ、梯度/二阶表（DiffTable）、适用维度与实现位置。

3. [03-奇异值本征系统与雅可比.md](./03-奇异值本征系统与雅可比.md)  
   项目核心技术。twist/flip/scaling 模态、λ 公式、ε 近似消除奇点、壳 3×2 vs 体 3×3 区别、与 Smith 2019 / Poya 2023 关系、从 DiffTable+SVD 到顶点空间力/滤波 Hessian 的 `convert` 流程。强烈推荐先读 `articles/eigensys.md`。

4. [04-接触势垒与动态刚度.md](./04-接触势垒与动态刚度.md)  
   IPC 风格变体。三种势垒（cubic/quadratic/logarithm）的精确 energy/grad/curv 公式（含 hindsight 修正）。gap/proximity（VF/EE/PE/VV + 碰撞网格）。`compute_stiffness` 的弹性 hess + mass/g² 混合 + 静态侧质量替换。push 单侧约束。

5. [05-摩擦模型与应变限制.md](./05-摩擦模型与应变限制.md)  
   正则化 Coulomb 摩擦（`Friction` 结构体、P 投影、λ 计算、friction_eps 平滑、combine 模式）。应变限制（barrier 复用于 stretch/奇异值、壳与杆特化、stiffness 估计、二分精确线搜索）。

6. [06-碰撞检测、ACCD 与约束.md](./06-碰撞检测、ACCD与约束.md)  
   LBVH 宽阶段 + swept AABB + 窄阶段精确距离。ACCD 迭代 TOI（保守推进）。两阶段 dyn-CSR 接触 Hessian 嵌入。线搜索流程。pins/sphere/floor/collision-mesh 约束、intersection 验证。

7. [07-时间积分、Newton 求解器与管线.md](./07-时间积分、Newton 求解器与管线.md)  
   隐式变分形式（动量惯性项 + 所有势能累加）。target 外推 + gravity + inactive_momentum。Newton 循环（min_newton_steps + final error-reduction + TOI 收缩自适应）。PCG（fixed CSR + dyn CSR + 块对角预条件）。playback、塑性更新、其他能量（air/torque/stitch/inflate 等）。

8. [08-参考文献与出处.md](./08-参考文献与出处.md)  
   主 TOG 论文引用、eigensys.md 中的 Smith/Poya/Stomakhin/Zhu 等、代码 License（Apache 2.0）与作者、hindsight 链接、进一步阅读源码建议。更新文档时的文件检查清单。

## 核心架构概览 (Mermaid)

```mermaid
flowchart TD
    A[输入网格与 rest 形状<br/>tri/tet/rod + inv_rest] --> B[每步预处理<br/>质量/面积/长度/初始角度]
    B --> C{元素类型?}

    C -->|壳 face| D[计算 F 3x2 = dx * inv_rest2x2<br/>svd3x2 / svd3x2_shifted]
    C -->|体 tet| E[计算 F 3x3 = dx * inv_rest3x3<br/>svd3x3_rv]
    C -->|杆 edge| F[边长 l/l0 + 三点角度 theta]

    D & E --> G[超弹性模型<br/>ARAP/StVK/SNHk/BaraffWitkin<br/>+ detsqr 体积]
    F --> H[杆拉伸 Hook + 弯曲 strand]

    G --> I[make_diff_table<br/>deda = ∂Ψ/∂σ, d2ed2a = ∂²Ψ/∂σ²]
    I --> J[eigenanalysis<br/>twist/flip/scale 模式 Q_k<br/>λ_k = ... max(λ,0) 投影]
    J --> K[convert_force/hessian<br/>回顶点空间 dedx / d2edx2]

    H --> L[解析梯度与 Gauss-Newton hess<br/>dihedral_angle]

    K & L --> M[embed_elastic_force_hessian<br/>atomic 写 fixed CSR + force]

    M --> N[动量项 momentum<br/>½m‖y-target‖²/dt² + 空气/风/拉力]
    N --> O[应变限制 strainlimiting<br/>barrier 作用于 stretch/S<br/>复用 eigen 或解析 PSD]
    O --> P[接触势垒 contact<br/>VF/EE/PE/VV + 碰撞网格<br/>cubic/quad/log + compute_stiffness<br/>+ Friction 投影]
    P --> Q[约束 pins/fix/sphere/floor/push/stitch<br/>torque 预计算]

    Q --> R[总梯度 g + Hessian<br/>fixed CSR（弹性）<br/>+ dyn CSR（接触，可变 nnz）<br/>+ 每顶点 diag hess]
    R --> S[PCG 求解<br/>块对角预条件 + Unrolled 3x3]
    S --> T[线搜索<br/>contact::line_search（ACCD）<br/>+ strain/rod 二分 TOI]
    T --> U[接受步长 t<br/>更新位置 + 可选 final error-reduction]
    U --> V[步后塑性更新<br/>inv_rest / rest_angle 缓慢演化]
    V --> W[日志与可视化]

    style J fill:#e6f3ff
    style P fill:#fff4e6
    style S fill:#f0e6ff
```

（上图浓缩了 Newton 一步的主要数据流。实际代码见 `main/main.cu:482` 左右的装配顺序与 `contact/contact.cu:976/1170` 的两阶段嵌入。）

## 源码与文章交叉引用

- 弹性模型与本征：`crates/ppf-cts-solver/src/cpp/energy/model/{arap,stvk,snhk,baraffwitkin,hook,dihedral_angle,detsqr}.hpp` + `eigenanalysis/eigenanalysis.{hpp,cu}` + `utility/utility.{hpp,cu}`
- 接触与势垒：`barrier/{cubic,quadratic,logarithm,barrier}.{hpp,cu}` + `contact/{contact,accd,distance,aabb}.{hpp,cu}` + `lbvh/*`
- 应变限制：`strainlimiting/strainlimiting.{hpp,cu}`
- 求解器管线：`main/main.cu`（advance 循环）、`solver/solver.{hpp,cu}`（PCG）、`energy/energy.{hpp,cu}`（momentum/stitch/embed 入口）
- 详尽理论推导：`articles/eigensys.md`（含 SymPy 验证代码与重推导 Smith 系统）
- 事后勘误与数值洞见：`articles/hindsight.md`（quadratic 曲率、ACCD 浮点、倾斜角度）
- 主论文引用位置：根 `README.md:134-156`

每篇理论文档末尾均有“源码出处”小节，便于对照。

## 如何贡献/更新本理论文档

1. 当修改 `energy/model/*`、`barrier/*`、`contact/*`、`strainlimiting/*`、`utility/*` 或 `main/main.cu` 中的公式/逻辑时，同步检查并更新对应理论文档。
2. 公式必须与至少两个源码位置（.hpp + .cu 或实现调用点）一致；优先引用精确代码片段。
3. 保持简体中文 + MathJax 兼容公式（用 `$$` 或 `$`）。代码块注明文件与行号。
4. 表格优先用于对比（四种弹性模型、三种势垒、壳 vs 体 vs 杆）。
5. 使用 `> [!NOTE]` / `> [!TIP]` 强调数值稳定性要点（eps guard、质量替换、Valanis-Landel 假设、Gauss-Newton for bend 等）。
6. 更新后运行根目录的 lint（`.markdownlint.json` 已配置）或手动检查内部链接。
7. 大型重构后，考虑更新本 `index.md` 的 Mermaid 架构图与“更新检查清单”（见 08-参考文献与出处.md）。

**许可**：本理论文档随项目采用 Apache 2.0 许可。引用时请同时提及主论文与对应源码文件。

---

**开始阅读建议**：先浏览本 index 与 Mermaid 获得全貌；然后阅读 `03-奇异值本征系统与雅可比.md`（项目灵魂）+ `04-接触势垒与动态刚度.md`（核心创新）；最后对照 `articles/eigensys.md` 深入推导与 `main/main.cu` 管线实现。

祝阅读愉快！如有疑问，欢迎在 GitHub Discussions 或 Issues 中讨论（见根 README 相关章节）。