# 术语中英对照表 (Glossary)

> 按英文字母序。随读随补，新术语请注明首次出现的章节。

| English | 中文 | 说明 | 章节 |
|---|---|---|---|
| Apparent derivative | 表观导数 | 只对分量求导，把基当常量 | 2 |
| Articulated-body algorithm (ABA) | 铰接体算法 | $O(n)$ 的正动力学算法 | 7 |
| Articulated-body inertia | 铰接体惯性 $I^A$ | 子树关节**自由**时的等效惯性 | 7 |
| Bias force | 偏置力 | ①单刚体的 $\mathbf{v}\times^*I\mathbf{v}$；②系统的 $C$；③ABA 的 $\mathbf{p}^A$ | 2,3,7 |
| Baumgarte stabilization | Baumgarte 稳定化 | 抑制约束漂移的方法 | 8 |
| Closed-loop system | 闭环系统 | 连通图含回路 | 8 |
| Composite rigid body algorithm (CRBA) | 复合刚体算法 | 构造 $H$ 的 $O(n^2)$ 算法 | 6 |
| Composite rigid body inertia | 复合刚体惯性 $I^c$ | 子树**焊死**时的总惯性 | 6 |
| Connectivity graph | 连通图 | 节点=刚体，边=关节 | 4 |
| Constraint drift | 约束漂移 | 只在加速度层施加约束导致的位置误差累积 | 8 |
| Constraint force | 约束力 | 关节维持约束所需的力，不做功 | 3 |
| Constraint force subspace | 约束力子空间 $T$ | $F^6$ 中与 $\mathrm{range}(S)$ 对偶正交的子空间 | 3 |
| Dual space | 对偶空间 | $F^6=(M^6)^*$ | 2 |
| Fixed base | 固定基座 | body 0 | 4 |
| Floating base | 浮动基座 | 用 6-DoF 自由关节连到惯性系 | 4,9 |
| Force space | 力空间 $F^6$ | 装力、动量的 6D 空间 | 2 |
| Forward dynamics (FD) | 正动力学 | 已知 $\tau$ 求 $\ddot q$ | 6,7 |
| Free vector | 自由向量 | 纯平移的运动向量 / 纯力偶的力向量 | 2 |
| Generalized coordinates | 广义坐标 $q$ | | 3 |
| Generalized force | 广义力 $\tau$ | | 3 |
| Gravity trick | 重力技巧 | 令 $\mathbf{a}_0=-\mathbf{a}_g$ 代替逐体加重力 | 5 |
| Hybrid dynamics | 混合动力学 | 部分关节已知力、部分已知加速度 | 9 |
| Impulsive dynamics | 冲量动力学 | 碰撞瞬间的速度跳变 $H\Delta\dot q=\iota$ | 3,9 |
| Inverse dynamics (ID) | 逆动力学 | 已知 $\ddot q$ 求 $\tau$ | 5 |
| Joint model | 关节模型 | `jcalc` 返回的 $(X_J, S, \mathbf{c}_J)$ | 4 |
| Joint transform | 关节变换 $X_J$ | 依赖 $q_i$ 的那部分变换 | 4 |
| Joint-space inertia matrix (JSIM) | 关节空间惯性矩阵 $H$ | 对称正定，分支诱导稀疏 | 3,6 |
| Kinematic tree | 运动学树 | 无回路的连通图 | 4 |
| Line vector | 线向量 | 纯转动的运动向量 / 沿一条线的纯力 | 2 |
| Loop constraint | 环约束 | 闭环关节给出的约束方程 | 8 |
| Loop Jacobian | 环雅可比 $K$ | $\partial\phi/\partial q$ | 8 |
| Loop joint | 闭环关节 | 生成树之外的边 | 4,8 |
| $LTL$ / $LTDL$ factorization | $LTL$ 分解 | $H=L^{\mathsf T}L$，正则编号下零填充 | 6 |
| Motion space | 运动空间 $M^6$ | 装速度、加速度的 6D 空间 | 2 |
| Motion subspace matrix | 运动子空间矩阵 $S$ | 关节允许的运动方向，$6\times n_i$ | 3,4 |
| Operational space inertia matrix | 操作空间惯性矩阵 $\Lambda$ | $(JH^{-1}J^{\mathsf T})^{-1}$ | 8,9 |
| Outward pass | 外推 | 根 → 叶，`for i = 1 to N` | 4,5 |
| Inward pass | 内推 | 叶 → 根，`for i = N to 1` | 4,5 |
| Parent array | 父节点数组 $\lambda(i)$ | 编码树拓扑 | 4 |
| Plücker coordinates | Plücker 坐标 | 空间向量的标准 6D 坐标 | 2 |
| Plücker basis | Plücker 基 | $M^6$ 与 $F^6$ 的互为对偶的基 | 2 |
| Recursive Newton-Euler algorithm (RNEA) | 递推牛顿-欧拉算法 | $O(n)$ 的逆动力学算法 | 5 |
| Regular numbering | 正则编号 | $\lambda(i)<i$ | 4 |
| Rotor inertia | 转子惯量 | 反射到关节侧为 $r^2I_{\text{rotor}}$ | 9 |
| Schur complement | Schur 补 | $I^A-UD^{-1}U^{\mathsf T}$，ABA 的消元核心 | 7 |
| Screw | 螺旋 | 转动+沿轴平移的组合 | 2 |
| Screw pitch | 螺距 | $h=\omega\cdot v_O/\omega\cdot\omega$ | 2 |
| Spanning tree | 生成树 | 闭环系统的树形骨架 | 4,8 |
| Spatial acceleration | 空间加速度 $\mathbf{a}$ | $\dot{\mathbf v}$，**不是**物质点加速度 | 2 |
| Spatial force | 空间力 $\mathbf{f}$ | $[n_O; f]\in F^6$，即 wrench | 2 |
| Spatial inertia | 空间惯性 $I$ | $M^6\to F^6$ 的映射 | 2 |
| Spatial momentum | 空间动量 $\mathbf{h}$ | $I\mathbf{v}\in F^6$ | 2 |
| Spatial vector | 空间向量 | 6D 的运动量或力量 | 2 |
| Spatial velocity | 空间速度 $\mathbf{v}$ | $[\omega; v_O]\in M^6$，即 twist | 2 |
| Subtree | 子树 $\nu(i)$ | 以 $i$ 为根（含 $i$） | 4 |
| Support | 支撑集 $\kappa(i)$ | 从 $i$ 到根路径上的关节 | 4 |
| Tree transform | 树变换 $X_T$ | 常量，描述连杆几何 | 4 |
| Twist / Wrench | 旋量 / 力旋量 | 螺旋理论中对应 $\mathbf{v}$ / $\mathbf{f}$ | 2 |
| Unilateral constraint | 单边约束 | 接触约束，$\phi\ge0$ 且 $f_n\ge0$ | 9 |

## 容易混淆的术语对

| A | B | 区别 |
|---|---|---|
| 复合刚体惯性 $I^c$ | 铰接体惯性 $I^A$ | 子树**焊死** vs 子树**自由**；$I^A\preceq I^c$ |
| 空间加速度 | 经典加速度 | 差一个 $\omega\times v_O$ |
| $X$ | $X^*$ | 运动 vs 力；$X^*=X^{-\mathsf T}$ |
| $\times$ | $\times^*$ | 作用于运动 vs 作用于力；$\times^*=-(\times)^{\mathsf T}$ |
| $C$（本书，向量） | $C$（多数教材，矩阵） | 本书 $C$ 已含重力，且是向量 |
| 外推 outward | 内推 inward | 根→叶 vs 叶→根 |
| $\lambda(i)$（parent） | $\lambda$（Lagrange 乘子） | 同字母，靠有无下标区分 |

---

## ✍️ 随读补充

| English | 中文 | 说明 | 章节 |
|---|---|---|---|
|  |  |  |  |
