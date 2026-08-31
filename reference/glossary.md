# 术语中英对照表 (Glossary)

> 按英文字母序，标注了原书首次出现的章节。

| English | 中文 | 说明 | 章 |
|---|---|---|---|
| Abstract vector | 抽象向量 | 被坐标向量表示的那个几何/物理对象 | 2.1 |
| Active contacts | 活动接触 | 当前接触的子集，法向构成 range($T$) 的基 | 11.5 |
| Active force subspace | 主动力子空间 $T_a$ | 驱动器等的力方向，$T_a^{\mathsf T}S=\mathbf 1$；**不唯一** | 3.5 |
| Apparent derivative | 表观导数 $\mathring{(\cdot)}$ | 动系观察者看到的变化率 | 2.10 |
| Apparent inverse inertia | 表观逆惯性 $\Phi=S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}$ | 受约束刚体的 | 3.6 |
| Articulated body | 铰接体 | 其动力学效应被计入 $I^A,p^A$ 的（子）系统 | 7.1 |
| Articulated-body algorithm (ABA) | 铰接体算法 | $O(N_B)$ 的正动力学算法 | 7.3 |
| Articulated-body inertia | 铰接体惯性 $I^A$ | 刚体**作为系统一部分**时表现的惯性；**21 个参数**；只映射加速度→力 | 7.1 |
| Axial screw transform | 轴向螺旋变换 | 绕单轴旋转 + 沿同轴平移；**特别便宜** | A.4 |
| Base | 基座 | 连通图中唯一的固定节点，body 0 | 4.1 |
| Base inertia parameters | 基惯性参数 | 可从动力学行为辨识的那个子集 | 9.7 |
| Basic solution | 基本解 | 非零变量个数最少的 LCP/QP 解 | 11.5 |
| Baumgarte stabilization | Baumgarte 稳定化 | 抑制约束漂移；**目的是稳定不是精度** | 8.3 |
| Bias acceleration | 偏置加速度 $b$ | $f=0$ 时的加速度 | 2.15 |
| Bias force | 偏置力 | **产生零加速度所需的力**；$p$、$C$、$p^A$ 都是 | 2.14, 3.1, 7.1 |
| Bias velocity | 偏置速度 $\sigma$ | $\dot q=0$ 时的 $v_J$；仅非定常关节非零 | 3.5 |
| Body coordinates | 体坐标系 | 嵌在刚体中的坐标系；**数值上比绝对坐标系更准** | 4.2, 10.1 |
| Body Jacobian | 体雅可比 $J_i$ | $v_i=J_i\dot q$ | 2.7, 4.1 |
| Branch-induced sparsity | 分支诱导稀疏 | $H$ 中由树拓扑决定的结构性零 | 6.4 |
| Chord | 弦 | 生成树之外的图的边 | 4.1 |
| Classical acceleration | 经典加速度 $\hat a'$ | $[\dot\omega;\ddot r]$；与空间加速度差 $\omega\times v_O$ | 2.11 |
| Closed-loop system | 闭环系统 | 不是运动学树的刚体系统 | 1.1, 4.1 |
| Coefficient of restitution | 恢复系数 $e$ | 碰撞前后分离速度之比 | 11.7 |
| Compliant / soft contact | 柔顺 / 软接触 | 把接触面建成一阶动力系统 | 11.8 |
| Composite-rigid-body algorithm (CRBA) | 复合刚体算法 | 构造 $H$ 的最快算法 | 6.2 |
| Composite-rigid-body inertia | 复合刚体惯性 $I^c$ | 子树**焊死**时的总惯性 | 6.2 |
| Configuration ambiguity | 位形歧义 | 独立变量 $y$ 不唯一确定闭环系统位形 | 8.10 |
| Connectivity graph | 连通图 | 节点=刚体，弧=关节 | 4.1 |
| Constraint force subspace | 约束力子空间 $T=\mathcal S^{\perp}$ | | 3.5 |
| Constraint stabilization | 约束稳定化 | 见 Baumgarte | 8.3 |
| Coordinate vector | 坐标向量 | $\mathbb R^n$ 的元素；需区分时加下划线 | 2.1 |
| Current contacts | 当前接触 | 处于状态 3 且 $\dot\zeta_i=0$ 的接触 | 11.5 |
| Denavit-Hartenberg parameters | DH 参数 | $d_i,\theta_i,a_i,\alpha_i$；$4n+6$ 个 | 4.3 |
| Depth | 深度 $d$ | 树中最深刚体的深度 | 6.4, 10.3 |
| Direct sum | 直和 $V=\mathcal S_1\oplus\mathcal S_2$ | | 3.3 |
| Dual basis / reciprocity condition | 对偶基 / 互易条件 | $d_i\cdot e_j=\delta_{ij}$ | 2.1 |
| Dyad / dyadic | 并矢 / 并矢张量 | $a\,b\cdot$（秩 1）/ 其和 | 2.1 |
| Dynamic equivalence | 动力学等价 | 两个系统有相同的运动方程 | 9.7 |
| Euler parameters | Euler 参数 | 单位四元数；**无奇异但需 4 条额外处理** | 4.5 |
| Floating base | 浮动基 | 经 6-DoF 关节连到固定基座的刚体 | 4.1, 9.3 |
| Forward dynamics (FD) | 正动力学 | 已知 $\tau$ 求 $\ddot q$；主要用于仿真 | 1.1 |
| Free vector | 自由向量 | 纯平移 / 纯力偶；$s=0$；3 个数 | 2.5 |
| GU(n) | — | 效率基准：无分支、全转动、一般参数的链 | 10.3 |
| Handle | 柄 | $I^A,p^A$ 所指的那个刚体 | 7.1 |
| Hybrid dynamics | 混合动力学 | 部分关节已知力、部分已知加速度 | 9.1 |
| Impulse | 冲量 $\iota$ | 力的时间积分；$\iota=I\Delta v$ | 11.7 |
| Inequality / unilateral constraint | 不等式 / 单边约束 | $\phi(q)\ge0$ | 3.4, 11.4 |
| Inverse dynamics (ID) | 逆动力学 | 已知 $\ddot q$ 求 $\tau$ | 1.1 |
| Inverse inertia | 逆惯性 $\Phi$ | 受约束刚体只有它、没有 $I$ | 2.15 |
| Jourdain's principle of virtual power | Jourdain 虚功率原理 | 约束力沿任何相容速度自由度不做功 | 3.2 |
| Joint model | 关节模型 | `jcalc` 返回的 $(X_J,S,c_J)$ | 4.4 |
| Joint polarity | 关节极性 | 前驱→后继的朝向；正向/反向 | 4.1.3 |
| Joint-space inertia matrix (JSIM) | 关节空间惯性矩阵 $H$ | 对称正定、分支诱导稀疏 | 6.1 |
| Kinematic tree | 运动学树 | 不含运动学回路的刚体系统 | 1.1, 4.1 |
| LCP (linear complementarity problem) | 线性互补问题 | $\dot\zeta=M\lambda+d,\ \dot\zeta\ge0,\lambda\ge0,\dot\zeta^{\mathsf T}\lambda=0$ | 11.5 |
| Line vector | 线向量 | 纯转动 / 沿一条线的纯力；$s\cdot s_O=0$；5 个数 | 2.5 |
| Loop closure function | 闭环函数 $\gamma$ | $q=\gamma(y)$；**用它则无闭环误差** | 8.11 |
| Loop Jacobian | 环雅可比 $J_{Ll}=J_{s(k)}-J_{p(k)}$ | | 8.2 |
| Loop joint | 闭环关节 | 弦对应的关节 | 4.1 |
| $LTL$ / $LTDL$ factorization | $LTL$ 分解 | **重排的 Cholesky / $LDL^{\mathsf T}$**；无填充 | 6.4, 6.5 |
| Mobility | 自由度（机构） | $n-r$；恰约束时 $=n_{tot}-6N_L$ | 8.10 |
| Model-based algorithm | 基于模型的算法 | 把 `model` 作为输入 | 1.1 |
| Motion subspace matrix | 运动子空间矩阵 $S$ | $6\times n_f$ | 3.5 |
| Multiple handles | 多柄 | 一个铰接体可有多个柄 | 7.1, 7.5 |
| Operational-space inertia matrix | 操作空间惯性矩阵 $\Lambda$ | $=(JH^{-1}J^{\mathsf T})^{-1}$，即投影法的 $I^A$ | 7.2.1 |
| Outward / inward pass | 外推 / 内推 | 根→叶 / 叶→根 | 5.3 |
| Overconstrained | 过约束 | $r<n_c$；**含平面回路的系统一定是** | 3.2, 8.10 |
| Parent array | 父数组 $\lambda$ | 编码树拓扑 | 4.1 |
| Planar vectors | 平面向量 | $M^3,F^3$；平面刚体惯性 **4 个参数** | 2.16 |
| Plücker coordinates / basis | Plücker 坐标 / 基 | $M^6$ 与 $F^6$ 上互为对偶的基 | 2.2–2.4 |
| Predecessor / successor | 前驱 / 后继 | 关节"从前驱连到后继" | 3.5, 4.1.3 |
| Prescribed motion | 规定运动 | 混合动力学的典型用途 | 9.1 |
| Projection method | 投影法 | 左乘 $G^{\mathsf T}$ 消掉约束力；或 §7.2.1 的 $JH^{-1}J^{\mathsf T}$ | 3.2, 7.2.1 |
| Properly constrained / actuated | 恰约束 / 恰驱动 | $r=n_c$ / $p=\mathrm{rank}(G_u)=n-r$ | 8.10, 8.12 |
| Radius of gyration | 回转半径 | 判断坐标系是否"太远"的尺度 | 10.1 |
| Recurrence relation | 递推关系 | **现代算法高效的根本原因** | 5.2 |
| Recursive Newton-Euler algorithm (RNEA) | 递推牛顿-欧拉算法 | $O(n)$ 逆动力学 | 5.3 |
| Regular numbering | 正则编号 | $\lambda(i)<i$ | 4.1.2 |
| Rheonomic / scleronomic | 非定常 / 定常 | 显含时间 / 不显含 | 3.4 |
| Schur complement | Schur 补 | $I^A-UD^{-1}U^{\mathsf T}$ | 7.2 |
| Screw / screw axis / pitch | 螺旋 / 螺旋轴 / 螺距 | Chasles 定理 | 2.5 |
| Sensitivity | 敏感性 | 力的微小变化引起加速度剧变 | 10.2 |
| Separation velocity / acceleration | 分离速度 / 加速度 $\zeta,\dot\zeta$ | | 11.1 |
| Spanning tree | 生成树 | 含全部节点的树形子图；**应选深度最小的** | 4.1, 10.2 |
| Spatial acceleration | 空间加速度 $\hat a$ | $\dot{\hat v}$；**不是**物质点加速度 | 2.11 |
| Spatial inertia | 空间惯性 $I$ | **10 个参数**；可加 | 2.13 |
| Spatial vector | 空间向量 | 6D 的运动量或力量；**不是欧氏向量** | 2.1 |
| Support / child / subtree set | 支撑集 / 子集合 / 子树集 | $\kappa(i),\mu(i),\nu(i)$ | 4.1.4 |
| Symbolic simplification | 符号化简 | 自动生成定制代码 | 10.4 |
| System model | 系统模型 | 描述**系统本身**（区别于数学模型） | 1.1, 4.6 |
| Torsional friction / restitution | 扭转摩擦 / 恢复 | 因局部变形使点接触变面接触 | 11.7 |
| Tree transform | 树变换 $X_T$ | 常量，描述连杆几何 | 4.2 |
| Truncation error | 截断误差 | 主要来自数值积分 | 10.1 |
| Underactuated / redundantly actuated | 欠驱动 / 冗余驱动 | | 8.12 |
| Variable mobility | 可变自由度 | 闭环系统的 $r$ 随位形变化 | 8.10 |
| Velocity-product acceleration | 速度乘积加速度 $a^{vp}$ | $\ddot q=0$ 时的加速度；算 $C$ 的副产品 | 8.2 |

## 容易混淆的术语对

| A | B | 区别 |
|---|---|---|
| 复合刚体惯性 $I^c$ | 铰接体惯性 $I^A$ | 子树**焊死** vs **自由**；10 vs 21 个参数 |
| $I^A$ | $I^a$ | 用**子**的加速度 vs 用**父**的加速度表达同一个关节力（式 7.25） |
| 空间加速度 $\hat a$ | 经典加速度 $\hat a'$ | 差 $\omega\times v_O$；前者在静系求导，后者在跟随平移的系求导 |
| $X$ / $\times$ | $X^{*}$ / $\times^{*}$ | 运动 vs 力 |
| 本书的 $C$（向量、含重力） | 多数教材的 $C$（矩阵） | 跨教材头号混淆点 |
| 线向量 | 自由向量 | 有作用线 vs 无位置 |
| 欧氏正交补 | **对偶**正交补 | 同一空间内 vs 互为对偶的两个空间之间 |
| 虚功原理 | **Jourdain 虚功率原理** | 虚位移 vs 虚速度；后者更适合非完整约束 |
| 恰约束 | 恰驱动 | $r=n_c$ vs $p=\mathrm{rank}(G_u)=n-r$ |
| $\tau^c$ | $\tau^a$ | 闭环关节的**约束力** vs **主动力** |
| 当前接触 | 活动接触 | 状态 3 的全部 vs 其中构成基的子集 |
| $\lambda(i)$（parent） | $\lambda$（Lagrange 乘子） | 同字母，靠有无下标区分 |

---

## ✍️ 随读补充

| English | 中文 | 说明 | 章 |
|---|---|---|---|
|  |  |  |  |
