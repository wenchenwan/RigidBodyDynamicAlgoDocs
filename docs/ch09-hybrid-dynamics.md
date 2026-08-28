# 第 9 章 混合动力学与其他专题 (Hybrid Dynamics and Other Topics)

> **一句话概括**：ID 和 FD 是两个极端（全部已知力 / 全部已知加速度），
> 混合动力学处理中间情形；外加浮动基、传动、碰撞与接触等工程上绕不开的专题。

## 本章定位

前面几章解决了"标准问题"，这一章处理**实际系统里真正会遇到的复杂情况**。
内容比较杂，各节相对独立，**可以按需选读**。

## 前置依赖

- 第 5 章 RNEA、第 7 章 ABA（尤其是铰接体惯性 $I^A$）
- 第 3 章 §3.7 冲量动力学
- 第 8 章（接触部分与闭环约束同构）

---

## 9.1 混合动力学 (Hybrid Dynamics)

### 问题定义

把关节集合分成两部分：

| 集合 | 已知 | 求 |
|---|---|---|
| $\mathcal{F}$（force-known） | $\tau_i$ | $\ddot q_i$ |
| $\mathcal{A}$（acceleration-known） | $\ddot q_i$ | $\tau_i$ |

- $\mathcal{A}=\varnothing$ → 退化为 **FD**
- $\mathcal{F}=\varnothing$ → 退化为 **ID**

**典型应用**：

- **浮动基机器人**：基座的 6-DoF 关节 $\tau=0$ 已知 → 属于 $\mathcal{F}$；
  被电机严格位置控制的关节 $\ddot q$ 已知 → 属于 $\mathcal{A}$
- **部分位置控制**：机械臂某几个关节被锁定或做位置伺服
- **运动学标定 / 轨迹跟踪**：末端轨迹给定，某些关节的运动被规定
- **人机交互仿真**：人手拖动的那个关节加速度已知

### 矩阵形式

按 $\mathcal{F}/\mathcal{A}$ 分块：

$$
\begin{bmatrix}H_{\mathcal{FF}} & H_{\mathcal{FA}}\\ H_{\mathcal{AF}} & H_{\mathcal{AA}}\end{bmatrix}
\begin{bmatrix}\ddot q_{\mathcal{F}}\\ \ddot q_{\mathcal{A}}\end{bmatrix}
+\begin{bmatrix}C_{\mathcal{F}}\\ C_{\mathcal{A}}\end{bmatrix}
=\begin{bmatrix}\tau_{\mathcal{F}}\\ \tau_{\mathcal{A}}\end{bmatrix}
$$

已知 $\tau_{\mathcal F}$ 和 $\ddot q_{\mathcal A}$，求 $\ddot q_{\mathcal F}$ 和 $\tau_{\mathcal A}$：

$$
\ddot q_{\mathcal F} = H_{\mathcal{FF}}^{-1}\left(\tau_{\mathcal F} - C_{\mathcal F} - H_{\mathcal{FA}}\ddot q_{\mathcal A}\right)
$$

$$
\tau_{\mathcal A} = H_{\mathcal{AF}}\ddot q_{\mathcal F} + H_{\mathcal{AA}}\ddot q_{\mathcal A} + C_{\mathcal A}
$$

代价 $O(n^3)$（或利用稀疏性降低）。

### 铰接体混合动力学算法 —— $O(n)$ 版本

更漂亮的做法：**改造 ABA**。在第 2 趟（内推）里，对每个关节判断类型：

```
第 2 趟内推中，对关节 i：
  if i ∈ F  (力已知，加速度未知):
      照常做 ABA 的消元：
      U_i = I^A_i S_i,  D_i = S_iᵀ U_i,  u_i = τ_i − S_iᵀ p^A_i
      I^a = I^A_i − U_i D_i⁻¹ U_iᵀ          ← 消去该自由度（Schur 补）
      p^a = p^A_i + I^a c_i + U_i D_i⁻¹ u_i

  if i ∈ A  (加速度已知，力未知):
      不做消元，因为 q̈_i 已知，直接把它的效应并进偏置力：
      I^a = I^A_i                            ← 惯性原样传递
      p^a = p^A_i + I^A_i (c_i + S_i q̈_i)    ← 已知加速度贡献到偏置力

第 3 趟外推中：
  if i ∈ F:  q̈_i = D_i⁻¹ (u_i − U_iᵀ a')    ← 求加速度（同 ABA）
  if i ∈ A:  q̈_i 已知
  a_i = ⁱX_λ(i) a_λ(i) + c_i + S_i q̈_i
  if i ∈ A:  τ_i = S_iᵀ (I^A_i a_i + p^A_i)  ← 求关节力（类似 RNEA 的投影）
```

> 💡 **统一视角**：$\mathcal{F}$ 类关节做**消元**（未知量，用 Schur 补消掉），
> $\mathcal{A}$ 类关节做**代入**（已知量，直接并入偏置项）。
> 这正是解线性方程组时"消元 vs 代入"的选择，只不过是在树上做的。
> 从这个角度看，**RNEA 是全代入的极端，ABA 是全消元的极端，混合动力学是两者的插值**。

---

## 9.2 浮动基座系统 (Floating-Base Systems)

**建模**（第 4 章已提过）：在惯性系与基座之间插一个 6-DoF 自由关节。
系统仍是树，所有算法直接可用。

**特点**：

- $n = 6 + n_j$（$n_j$ = 实际关节自由度）
- 前 6 个广义力 $\tau_{1:6} = 0$（没有东西驱动基座）——**这是欠驱动的来源**
- 于是 FD 问题天然是一个混合动力学问题的特例

**分块形式**（人形/四足控制中的标准写法）：

$$
\begin{bmatrix}H_{bb} & H_{bj}\\ H_{jb} & H_{jj}\end{bmatrix}
\begin{bmatrix}\ddot q_{b}\\ \ddot q_{j}\end{bmatrix}
+\begin{bmatrix}C_{b}\\ C_{j}\end{bmatrix}
=\begin{bmatrix}\mathbf{0}\\ \tau_{j}\end{bmatrix}
+\begin{bmatrix}J_{b}^{\mathsf T}\\ J_{j}^{\mathsf T}\end{bmatrix}f_{\text{ext}}
$$

上面 6 行称为**欠驱动动力学方程**，它等价于**动量守恒定律**：
没有外力时系统的空间动量守恒。这是行走控制里 ZMP、
质心动量 (centroidal momentum) 等概念的出处。

**位形表示的坑**（第 4 章 §4.5 提过）：
基座姿态用四元数（4 个数）+ 位置（3 个数）= 7 维 $q_b$，
但速度是 6 维。所以：

$$
\dim q = 7 + n_j \ne \dim\dot q = 6 + n_j
$$

积分时**不能简单地 $q\mathrel{+}=\dot q\,\Delta t$**，
姿态部分必须用四元数积分（或在 $SO(3)$ 上做指数映射更新）。
这是浮动基仿真的头号 bug 来源。

---

## 9.3 传动与齿轮 (Gears and Transmissions)

真实电机有**转子惯量 (rotor inertia)**，经过减速比 $r$ 放大后
在关节侧表现为 $r^2 I_{\text{rotor}}$。

**为什么不能忽略**：谐波减速器 $r$ 可达 100~160，
$r^2$ 就是 $10^4$ 量级。对轻质连杆的机器人，
**反射惯量常常和连杆惯量同量级甚至更大**。忽略它会让模型严重失真。

**最简处理**：在 $H$ 的对角元上加 $r_i^2 I_{\text{rotor},i}$。
这个近似忽略了转子的陀螺效应，对多数应用足够。

**顺带的好处**：加大对角元会**改善 $H$ 的条件数**（第 10 章），
数值上更好解。

**完整处理**：把转子当作独立刚体建模，用齿轮约束连接——
这会引入闭环（第 8 章），代价更高。

---

## 9.4 冲量动力学与碰撞 (Impulsive Dynamics and Collision)

### 单点碰撞

碰撞瞬间的方程（第 3 章 §3.7）：

$$
H(q)\,\Delta\dot q = J^{\mathsf T}\iota
$$

- $q$ 不变，$\dot q$ 跳变
- $C$ 项有限，积分后消失
- $J$：接触点的雅可比，$\iota$：接触冲量

接触点速度的跳变：$\Delta v_c = J\Delta\dot q = \underbrace{JH^{-1}J^{\mathsf T}}_{\Lambda^{-1}}\iota$

**$\Lambda = (JH^{-1}J^{\mathsf T})^{-1}$ 就是操作空间惯性矩阵**——
"从接触点看过去，机器人有多重"。注意它和第 8 章路线 B 里的
$KH^{-1}K^{\mathsf T}$ 是同一个东西。

### 恢复系数

**牛顿恢复模型**：$v_c^{+} = -e\,v_c^{-}$（$e\in[0,1]$）
- $e=0$：完全非弹性（塑性碰撞）
- $e=1$：完全弹性

联立即可解出 $\iota$，进而得到 $\dot q^{+}$。

⚠️ **多点同时碰撞时牛顿模型会产生能量增加**（不物理）。
更严谨的做法用 **Poisson 恢复模型**或基于能量的模型。

### 持续接触

持续接触与闭环约束在数学上**同构**，都是 $J\ddot q = \dots$，
但有两个关键差别：

| | 闭环约束 | 接触约束 |
|---|---|---|
| 类型 | 双边 (bilateral)，$\phi=0$ | **单边** (unilateral)，$\phi\ge0$ |
| 法向力 | 可正可负 | 只能推不能拉：$f_n\ge0$ |
| 互补性 | 无 | $f_n\cdot\phi = 0$（要么接触要么无力） |
| 摩擦 | 无 | 摩擦锥 $\|f_t\|\le\mu f_n$ |

于是求解从**线性方程组**变成 **LCP (线性互补问题)** 或 **QP**：

$$
0 \le f_n \perp (J\ddot q - \dots) \ge 0
$$

这是现代物理引擎（Bullet、MuJoCo、Drake、ODE）的核心。
摩擦锥通常被线性化成多面锥以化为 LP/QP。

〔待核对〕本书对接触的处理深度——Featherstone 主要给框架，
完整的接触求解属于另一个领域，建议配合 Stewart & Trinkle、
Anitescu、Todorov 等人的工作阅读。

---

## 9.5 其他专题〔待补充〕

读的时候按实际内容补：

- [ ] 灵敏度分析 / 动力学导数（现代最优控制、可微仿真需要）
- [ ] 弹性关节 / 柔性连杆
- [ ] 对称性与守恒量
- [ ] 参数辨识 (system identification)：动力学关于惯性参数是线性的

## 易错点与陷阱

1. **浮动基的 $q$ 积分**。姿态必须用四元数/李群积分，
   不能当成普通向量做欧拉积分。别忘了每步归一化四元数。
2. **忘记转子惯量**。轻质高减速比机器人的模型误差会大到无法接受。
3. **碰撞后忘了检查能量**。多点碰撞用牛顿模型可能产生能量，
   仿真会莫名其妙"越弹越高"。
4. **把单边接触当双边处理**。会出现"地面拉住脚"的现象——
   机器人抬腿时地面反而把它拽回去。
5. **混合动力学里 $\mathcal{F}/\mathcal{A}$ 划分错误**。
   同一个关节不能既指定力又指定加速度（过约束），也不能都不指定（欠定）。

## 与其他章的联系

- ← 第 3 章：冲量动力学的基本方程
- ← 第 7 章：铰接体惯性 $I^A$ 是混合动力学 $O(n)$ 算法的基础
- ← 第 8 章：接触约束 ≅ 闭环约束（单边化）
- → 第 10 章：这些扩展的代价分析

---

## ✍️ 我的理解

<!-- 读完后填 -->

## ❓ 疑问与待办

- [ ] 确认本章的实际小节构成（本笔记按主题组织，可能与原书划分不同）
- [ ] 实现浮动基 ABA，验证无外力时空间动量守恒
- [ ] 实现单点碰撞，验证 $e=1$ 时能量守恒、$e=0$ 时接触点法向速度为零
- [ ] 查一下本书对接触/摩擦的处理深度，决定是否需要补充其他资料

## 📌 与原文的出入

<!-- 本章按主题组织，与原书的小节划分可能出入较大，读时重点核对 -->
