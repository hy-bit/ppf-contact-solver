# 06 碰撞检测、ACCD 与约束

本篇描述项目如何在 GPU 上实现高效、鲁棒的**连续碰撞检测（CCD）** 与**约束力嵌入**，支撑“穿透-free”保证。

核心组件：
- LBVH（宽阶段） + swept AABB
- 窄阶段精确距离（distance.hpp）与 Proximity 构造（VF/EE/PE/VV + 碰撞网格）
- ACCD 迭代保守推进（accd.hpp）
- 两阶段 dyn-CSR 接触 Hessian 嵌入（embed_contact / embed_constraint）
- 独立线搜索 + 事后 intersection 验证
- 各类约束（pins/fix、sphere、floor、collision-mesh push）

## 广阶段：LBVH + Swept AABB

- `lbvh/lbvh.{hpp,cu}`：为 face/edge/vertex 分别构建 LBVH（Morton code + 并行 radix sort + 树构建）。
- 使用 **swept AABB**：AABB 从 x0 扩展到 x0 + extrapolate*(x1-x0)，支持 CCD。
- 碰撞网格（静态）单独 BVH（`build_collision_mesh_bvh`）。
- AABB 带 margin（常用 0.5*ghat + offset）。
- `aabb::query(BVH, aabbs, op, query_aabb)`：小栈（AABB_MAX_QUERY=128）遍历，命中叶时调用 op(leaf)（通常是某个 `*ContactForceHessEmbed` 结构体）。

更新策略：Newton 每步用 `line_search_max_t` 做保守 AABB 更新；背景 CPU 持续重建（不阻塞）。

## 窄阶段与 Proximity

自接触（内部网格）使用统一“barrier edge”形式，通过 Proximity<N> + bary 权重描述多顶点接触：

- **Point-Face (VF)**：`point_triangle_distance_coeff` 得 3-bary（夹取 [0,1]）；Proximity<4>，value = {1, -c0,-c1,-c2}
- **Edge-Edge (EE)**：类似 4-bary，Proximity<4>
- **Point-Edge (PE)**：Proximity<3>
- **Point-Point (VV)**：Proximity<2>

g 代理：后续 `ex = sum w_i * x_i`；`g = ||ex|| - offset` 传给 barrier。

邻接过滤（shared face/edge）决定是否计入摩擦与 count 修正（避免重复计数）。

碰撞网格 / sphere / floor 多使用 **push**（已知 normal 的单侧 cubic）而非通用 edge barrier，stiffness 也简化（normal 方向 + mass/gap²）。

`distance.hpp` 提供精确 bary 系数计算（带 unclassified 夹取）。

## ACCD 连续碰撞（保守迭代）

`accd.hpp` 的 `ccd_helper` + 特化 `point_triangle_ccd`、`edge_edge_ccd` 等：

关键技巧：
- `centerize(x)`：相对运动中心化（大幅改善数值）
- `max_relative_u`：估计最大相对速度
- 迭代：`toi += (d - target) / u_max`
- `eps = ccd_reduction * (current_dist - offset)`
- `target = eps + offset`
- 提前终止条件 + `line_search_max_t` 保护
- `ccd_max_iter` 控制

`contact::line_search`（1684 行起）先做解析约束（球/地板/固定），再对自/碰撞网格做 broad + 每对 CCD，返回 min TOI（归一化）。

## 接触力/ Hessian 嵌入（两阶段）

`embed_contact_force_hessian`（976 行）：

```cpp
for (int stage = 0; stage < 2; ++stage) {
    if (stage == 0) {
        // dry pass: 只检查 fixed 是否存在，否则标记 dyn 需要的 nnz
        dyn_out.start_rebuild_buffer();
    }
    // DISPATCH vert/edge + AABB query(op) 其中 op 是 PointFaceContact... 等结构体
    // 结构体 operator() 里做精确距离检查、构造 prox、调用 embed_contact_force_hess
}
```

`embed_contact_force_hess`（内部）：
- 计算 ghat/offset/friction（per-param combine）
- stiff_k = compute_stiffness（弹性 hess + mass/g²）
- f = stiff_k * barrier_grad_dir (+ optional Friction)
- H = stiff_k * rank1_curv (+ Friction hess)
- extend_by_bary + atomic 到 force / fixed 或 dyn

EE 单独 DISPATCH + count 修正。

结束后把临时 contact_force_vec 原子加回主 force。

**DynCSRMat 预算**：`csrmat_max_nnz`；超预算则 dyn_consumed >1，装配失败。

## 约束力嵌入（pins / sphere / floor / collmesh）

`embed_constraint_force_hessian`（1170 行起）类似流程，但负责：
- 固定 pin / FixPair（fix.hpp，简单弹簧 + stiffness + constraint_tol 保护）
- Sphere / Floor（push 或直接 + 厚度 pass-through）
- Collision-mesh（M2C / C2M / EE，AABB query + push + friction）

Kinematic vs 非 kinematic 有不同 gap / stiffness 处理。

## 线搜索与事后验证

- contact line_search 返回**归一化** TOI（实际接受分数 t / max_t）
- 与 strain/rod TOI 取 min
- `toi <= eps` → CCD 失败
- 成功接受后：最终 AABB + `check_intersection`（记录 IntersectionRecord，原子计数，用于 debug / CI 显式穿透检测）

## 碰撞窗口（可选优化）

`vert/edge/face_collision_active` 标志可剔除不参与碰撞的元素（由上层设置）。

## 与论文 / hindsight 的关系

- 论文核心是 cubic barrier + 动态刚度；CCD/ACCD/LBVH 是支撑工程（常见于 IPC 家族，但本项目有自定义 swept LBVH + 两阶段 dyn CSR + 静态质量替换等细节）。
- hindsight 提到：ACCD 在极紧间隙下单精度可能不准；cubic barrier 通过把距离曲线“推开”危险区来缓解（而非根本解决）。
- 参考分支 `sigasia-2024` 更接近论文时刻的碰撞策略。

## 源码出处（本篇关键引用）

- ACCD：`contact/accd.hpp:39`（ccd_helper + centerize/max_relative_u）
- 距离与 Proximity：`contact/distance.hpp`、`contact/contact.cu` 内大量 `*ContactForceHessEmbed` 结构体（PointFace... 447 起、EdgeEdge 505 起等）
- 两阶段主入口：`contact/contact.cu:976`（embed_contact）、`1170`（embed_constraint）
- 线搜索：`1684`（contact line_search）、`contact.hpp:30`
- AABB/LBVH：`contact/aabb.hpp`（query）、`lbvh/lbvh.*`（build/update/query）
- 约束模型：`energy/model/{fix,push}.hpp`
- 装配与 intersection：`main/main.cu:482`（调用顺序）、`535`（check_intersection）
- 进一步阅读：`04-接触势垒...md`（stiff_k 与 barrier 公式）、`05-...摩擦应变.md`（TOI 如何与 strain 组合）、`07-...管线.md`（Newton 整体 + final_step + dyn 预算监控）

---

本章完成“接触”部分的工程拼图。下一章将把所有势能（弹性 + 接触 + 应变 + 摩擦 + 动量 + ...）放在一起，看 Newton/PCG/自适应 TOI 收缩的完整时间积分管线。