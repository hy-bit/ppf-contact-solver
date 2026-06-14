# 07 时间积分、Newton 求解器与管线

本篇把前几章的所有势能（弹性本征雅可比、接触势垒+动态刚度、摩擦、应变限制、动量等）串成完整的**隐式时间积分管线**。

核心文件：`main/main.cu`（`advance` / `main_helper`）、`energy/energy.{hpp,cu}`（momentum/elastic/stitch 嵌入入口）、`solver/solver.{hpp,cu}`（PCG）、`data.hpp`（ParamSet 与全局状态）。

## 隐式目标函数（所有势能的和）

每 Newton 迭代在当前位置 `eval_x`（记 y）处装配**梯度 g**（“force”缓冲，正方向为下降方向的反）与 **Hessian**，然后解 H d = g，更新 y -= d。

累加的势能（通过 embed_* 贡献 force 与 hess）包括：

- **惯性/动量 (momentum)**：自由顶点 (m/2) ||y - target||² / dt²
- **超弹性 (elastic)**：壳/体/杆的 ARAP/StVK/SNHk/Baraff/Hook + 弯曲 + 体积（detsqr 或内建）
- **缝合 (stitch)**：4 顶点 seam 势
- **应变限制 (strainlimit)**：barrier 作用于 stretch/σ（壳走 eigen，杆解析 PSD）
- **接触势垒 (contact barrier)**：VF/EE/... + 碰撞网格（cubic/quad/log + compute_stiffness）
- **摩擦 (friction)**：正则化 Coulomb，叠加在接触力上
- **空气/风/拉力/扭矩/各向同性空气摩擦/fix_xz/inflate/pressure/push/fix**：其他辅助项
- **约束 (pins/fix/sphere/floor)**：软/硬销定与单侧推斥

**动量项精确定义**（`energy/model/momentum.hpp:14`）：

```cpp
__device__ float energy(float dt, const Vec3f &x, const Vec3f &y) {
    return 0.5f * (x - y).squaredNorm() / (dt * dt);
}
__device__ Vec3f gradient(float dt, const Vec3f &x, const Vec3f &y) {
    return (x - y) / (dt * dt);
}
__device__ Mat3x3f hessian(float dt) { return Mat3x3f::Identity() / (dt * dt); }
```

**target 计算**（线性外推 + 重力，variable-dt 友好）：

```cpp
// main/main.cu 附近 compute_target lambda
if (fix) target = constraint_fix.position;
else {
    float tr = dt / prev_dt, h2 = dt*dt;
    Vec3f y = (x1 - x0)*tr + h2 * gravity;
    target = inactive_momentum ? x1 : x1 + y;
}
```

`inactive_momentum`：跳过动量 + 风，target 直接 = curr（静态求解或姿态匹配有用）。

`playback`：`dt = param.dt * playback`；最终时间推进 `time += prev_dt / playback`（不改变内部物理 dt）。

## Newton 循环结构（伪代码 + 关键点）

`main/main.cu:advance()` 核心（~434 起）：

```cpp
eval_x = curr
compute_target(dt)
toi_advanced = 0; step=1; final_step=false

while (true) {
    if (final_step) { logging "error reduction step"; dt *= toi_advanced; compute_target(dt); }
    else { logging "newton step %u", step; }

    clear(dyn/fixed/diag/force/dx)

    preload fixed dx = eval_x - target   // 固定顶点直接拉到目标
    torque pre-pass (PCA)

    // ===== 矩阵装配 =====
    embed_momentum(...)          // 风/拉力/扭矩/动量/空气/... → diag_hess + force
    embed_elastic(...)           // 杆弯曲 + 杆 Hook + 壳(模型) + 体 + 铰链 → fixed + force
    if (stitch) embed_stitch(...)
    tmp_fixed = fixed
    if (shell) strainlim embed (barrier on sing + eigen 或 杆解析)   // 写 fixed
    if (rod)   rod_strainlim embed
    if (!disable_contact) {
        num_contact += embed_contact_force_hessian(...)   // self VF/EE + 两阶段 dyn CSR
    }
    num_contact += embed_constraint_force_hessian(...)    // pins + sphere/floor + collmesh

    // ===== PCG =====
    success, iter, res = solver::solve(dyn_hess, fixed_hess, diag_hess,
                                       force /*=g*/, cg_tol, cg_max_iter, dx, ...)
    if (!success) return fail

    toi_recale = min(1, prm.max_dx / max|dx|)   // 方向截断，防爆炸
    eval_x -= recale * dx
    (可选 fix_xz 后处理)

    // 保守 AABB 更新（contact 用 line_search_max_t）
    toi = contact::line_search(target, eval_x)          # CCD broad + exact
    if (shell) toi = min(toi, strain::line_search bisection)
    if (rod)   toi = min(toi, rod::line_search bisection)
    if (toi <= eps) return CCD fail

    if (!final) toi_advanced += (1-toi_advanced) * recale * toi
    eval_x = target + toi * (eval_x - target)

    if (final) break
    else if (toi_advanced >= target_toi && step >= min_newton_steps) final=true
    else { step++; compute_target(dt) }   // 下一 Newton 迭代
}

// 成功后
final AABB + check_intersection (显式穿透检测)
prev_dt = dt; time += prev_dt / playback
prev = curr; curr = eval_x
plasticity updates (face/tet/hinge/rod_bend)
return result (带丰富日志)
```

**final_step（误差缩减步）**：
- 积累足够 TOI 后做一次“额外”迭代（dt 收缩到已前进量）。
- 只做误差减小，不再推进时间/TOI。
- 注释明确：“------ error reduction step ------”。

**自适应步长**：主要通过 **TOI 收缩** 实现（contact/strain 强迫小步）。`enable_retry` 仅在日志注释中提及（上层可据 PCG fail 全局缩 dt 重试），核心 advance 本身失败即返回 flag，无内部重试循环。

## PCG 线性求解器

`solver/solver.cu`：

- Hessian 组装后 = **DynCSR**（接触/约束，可变稀疏） + **FixedCSR**（弹性/固定连接 + 转置引用） + **每顶点 diag Mat3x3f**（动量 + 空气 + 拉力等） + 可选标量 D I。
- `apply`：全 3×3 块展开（UnrolledMat3x3f 加速）。
- 预条件：块对角 `inv_diag[i] = invert(A(i,i)+B(i,i)+C[i])`（解析 3×3 det，奇异时退化为 diag(1/m_ii)）。
- CG：标准，初 x=0，r = b - A x（b=force=g）；双精度累加 rz/alpha/beta/err 减 FP 误差；`reresid = err/err0`；NaN 或超 iter 失败。
- `invert` 有 det==0 安全回退。

**力符号约定**：装配的 force 是**总势能梯度 g**；解 H d = g 后 y -= d（能量下降）。

## 其他重要机制

- **playback / gravity / inactive_momentum**：见 target 计算。
- **固定顶点处理**：fix_index >0 时跳过动量，target=pin 位置，装配预设 dx 使更新直接拉到目标；kinematic pin 有 stiffness + tol 保护。
- **Dyn hess 预算**：仅接触贡献；超 1.0 失败（日志 `dyn_consumed` / `max_nnz_row`）。
- **塑性更新**：成功步后（非 Newton 内），见 plasticity/（inv_rest / rest_angle creep）。
- **日志通道**（SimpleLog）：time-per-frame、matrix-assembly、pcg-linsolve、line-search、newton-steps、num-contact、max-sigma、dt、dyn_consumed 等。多数在 Newton 迭代内多次记录（因为内层 Newton 步）。

## 单精度 + GPU 特点

- 几乎全部 float32（仅 time、少数 CG 累加、日志用 double）。
- `DISPATCH_START(N) [...] __device__(unsigned i) ... DISPATCH_END`（utility/dispatcher）。
- 原子写：`atomic_embed_force<N>` / `atomic_embed_hessian`（force 是简单原子加；hess 走 CSR push 优先 fixed，不行进 dyn）。
- 内存：`buffer::MemoryPool` 每迭代临时分配（自动释放）；持久 Vec 用于 AABB/BVH 等。
- 无 cuSPARSE 等，纯自定义内核。

## 源码出处（本篇关键引用）

- 完整循环与装配顺序：`main/main.cu:434`（while Newton）、`482`（momentum+elastic+... 调用）、`688`（line_search 组合）、`822` 附近（最终更新 + plasticity + 日志）
- 动量：`energy/model/momentum.hpp:14-22`
- 目标与 playback：`main/main.cu:399`（compute_target）、`326`（dt = param->dt * playback）
- PCG：`solver/solver.{hpp,cu}`（solve/cg/apply/invert）、`main:510`（调用）
- 其他能量嵌入：`energy/energy.cu`（embed_momentum 153 起、embed_elastic 441 起、embed_stitch、torque）
- 参数：`data.hpp:343`（ParamSet 完整定义：dt, prev_dt, playback, inactive_momentum, cg_tol, cg_max_iter, line_search_max_t, target_toi, min_newton_steps, eiganalysis_eps, friction_eps, ccd_* 等）
- 进一步阅读：`01~06` 各篇（势能来源）、`articles/bug.md`（应变限制 LS 改进、BVH 更新策略、hang 例子步长调整）、`articles/refactor_202510.md`（整体重构历史）

---

**总结**：这是一个为**极端接触 + 高刚度材料 + 混合维度**高度优化的隐式求解器。Newton 外层通过 TOI 积累实现自适应，final error-reduction 步 polish 质量；装配高度模块化（diag 放 momentum 类，fixed 放弹性/缝，dyn 专接触）；大量 FP 防护散布在 SVD、eigen、invert、LS 二分、质量替换、sin eps guard 等处。

文档到此结束主流程。最后一篇 `08-参考文献与出处.md` 将汇总引用、许可与维护建议。