# 第 11 章 接触与碰撞 (Contact and Impact)

> **原书 pp. 213–240**，共 9 节。
> ⚠️ **这是全新的一章——我此前的笔记误以为全书只有 10 章，把接触内容塞进了第 9 章。**

> **一句话概括**：接触约束是**单边 (unilateral)** 的——只能推不能拉。
> 这把动力学从**线性方程**变成了**线性互补问题 (LCP)**，
> 并带来一个根本变化：**加速度不再是施加力的线性函数，只是分段线性函数**。

## 原书节次对照

| 节 | 标题 | 页 |
|---|---|---|
| 11.1 | Single Point Contact | 213 |
| 11.2 | Multiple Point Contacts | 216 |
| 11.3 | A Rigid-Body System with Contacts | 219 |
| 11.4 | **Inequality Constraints**（四个状态） | 222 |
| 11.5 | Solving Contact Equations | 224 |
| 11.6 | Contact Geometry | 227 |
| 11.7 | Impulsive Dynamics | 230 |
| 11.8 | Soft Contact | 235 |
| 11.9 | Further Reading | 239 |

## 开篇：接触问题的本质

> 两刚体接触时，受到一个说"两者不许相互穿透"的**接触约束**。
> 若 $\phi$ 度量二者之间的**有符号距离**（重叠时 $\phi<0$），则约束是不等式 $\phi\ge0$。
>
> **施加这个约束的力同样是单边的**：它们能阻止穿透，**但不能阻止分离**——
> **只能排斥，不能吸引**。
>
> 若两物体相遇时的速度与接触约束不相容，就会产生一个**冲量**使速度发生阶跃变化，
> 这个事件称为**碰撞 (impact)**。

---

## 11.1 单点接触

### 设置

刚体 $B$ 与固定表面在单点 $C$ 接触。**假定无摩擦**，且刚体在当前时刻前后
**只受有限力**（无冲量）。

$$
\boxed{\ \textbf{接触法向 }n\in F^6：\text{沿接触法线作用的}\textbf{单位力}\ }
$$

**$n$ 定义为指向表面外侧**，这样 $n$ 的正数倍会把刚体从表面**推开**。

### 三个标量

**接触力**（无摩擦 ⟹ 必沿法向）：

$$
f_c=n\lambda,\qquad \boxed{\lambda\ge0}\tag{11.1}
$$

**分离速度 (separation velocity)** = $C'$（刚体上此刻与 $C$ 重合的点）线速度的法向分量：

$$
\zeta=n\cdot v\tag{11.2}
$$

> 💡 **原书验证这个式子的办法很直观**：把坐标系原点放在 $C$、$z$ 轴沿法向，
> 则 $n$ 的 Plücker 坐标是 $[0\ 0\ 0\ 0\ 0\ 1]^{\mathsf T}$，
> 于是 $n\cdot v=v_{Cz}$——正是 $C'$ 线速度的 $z$ 分量。

### ⭐ 关键论证：当前时刻必有 $\zeta=0$

> **原书的论证**：若刚体只受有限力，则加速度有限。于是：
> - 若当前时刻 $\zeta<0$，刚体**在紧接着的未来会穿透**表面；
> - 若 $\zeta>0$，刚体**在刚刚过去的时刻已经穿透**过表面。
>
> 但**任何时刻都不允许穿透**，所以 $\zeta=0$。

$$
\boxed{\ \textbf{推论：}\zeta\ne0\ \Longrightarrow\ \textbf{存在冲量}\ }
$$

> 🔑 **这条推论把本章分成了两半**：
> $\zeta=0$ 的情形用**接触动力学**（§11.1–11.5），
> $\zeta\ne0$ 的情形用**冲量动力学**（§11.7）。

**分离加速度**：

$$
\dot\zeta=\frac{d}{dt}(n\cdot v)=n\cdot a+\dot n\cdot v\tag{11.3}
$$

$\dot n\cdot v$ 假定已知（它是 $B$ 的位置速度和接触面形状的函数）。
由 $\zeta=0$ 得 $\dot\zeta\ge0$（否则会穿透）。

### ⭐ 互补条件

**接触的精确行为**：当前时刻要么**断开**、要么**持续**。

| 情形 | 分离加速度 | 接触力 |
|---|---|---|
| **断开** | $\dot\zeta>0$（严格为正） | $\lambda=0$（正在失去接触） |
| **持续** | $\dot\zeta=0$ | $\lambda\ge0$ |

这两种行为被下面**一组式子简洁地统一**了：

$$
\boxed{\ \dot\zeta\ge0,\qquad \lambda\ge0,\qquad \dot\zeta\,\lambda=0\ }\tag{11.4}
$$

> 💡 **$\dot\zeta\lambda=0$ 是"互补性"的精髓**：两个非负量的乘积为零
> ⟺ **至少一个为零** ⟺ "要么分开（无力），要么保持接触（无相对加速度）"。

### 解（单点接触可以"看出来"）

由运动方程 $Ia+v\times^{*}Iv=f+n\lambda$（11.5）解出 $a$ 代入 11.3：

$$
\dot\zeta=M\lambda+d\tag{11.7}
$$

$$
M=n\cdot I^{-1}n
\qquad
d=n\cdot I^{-1}(f-v\times^{*}Iv)+\dot n\cdot v\tag{11.8, 11.9}
$$

$$
\boxed{
\begin{aligned}
d\ge0&：\quad \dot\zeta=d,\ \ \lambda=0\quad(\text{接触断开})\\
d<0&：\quad \dot\zeta=0,\ \ \lambda=-d/M\quad(\text{接触持续})
\end{aligned}}\tag{11.10}
$$

（$M$ **严格为正**，因为 $I$ 正定。）

### 🔑 本章最重要的概念结论

> *"Observe that **$d$ is a linear function of $f$, but $\lambda$ is only a
> piecewise-linear function of $d$**. Thus, the acceleration of $B$ is **not a linear
> function of the applied force, but only a piecewise-linear function**.
> In this respect, **contact dynamics differs fundamentally from the dynamics of
> systems with equality constraints**."*

$$
\boxed{\ \text{等式约束：}a\text{ 是 }f\text{ 的}\textbf{线性}\text{函数}\qquad
\text{接触约束：}a\text{ 只是 }f\text{ 的}\textbf{分段线性}\text{函数}\ }
$$

> 💡 **后果**：
> - 不能再用"叠加原理"；
> - 不能直接求偏导（**可微仿真在接触处不可微**，这是现代可微物理引擎的核心难点）；
> - 求解需要**组合搜索**（哪些接触保持、哪些断开），而不只是解线性方程。

---

## 11.2 多点接触

$n_c$ 个接触，各有接触点 $C_i$ 与法向 $n_i$。

$$
f_c=\sum_{i=1}^{n_c}n_i\lambda_i
\qquad
\zeta_i=n_i\cdot v
\qquad
\dot\zeta_i=n_i\cdot a+\dot n_i\cdot v\tag{11.11, 11.12, 11.13}
$$

> 💡 **原书特别提醒**：*"Observe that **each contact has its own separation velocity**,
> even though there is only a single moving body in the system."*

**每个 $\zeta_i$ 都必须为零** ⟹ $n_i\cdot v=0$ 对所有 $i$ 成立。

> ⚠️ **若接触法向张成整个 $F^6$，则唯一解是 $v=0$。
> 但这并不必然意味着刚体不能动**——它可能能加速离开（从而断开）一个或多个接触。

**矩阵形式**：$N=[n_1\ \cdots\ n_{n_c}]$（$6\times n_c$），

$$
f_c=N\lambda
\qquad
\dot\zeta=N^{\mathsf T}a+\dot N^{\mathsf T}v\tag{11.16, 11.17}
$$

$$
\boxed{\ \dot\zeta\ge0,\quad \lambda\ge0,\quad \dot\zeta^{\mathsf T}\lambda=0\ }\tag{11.18}
$$

$$
\dot\zeta=M\lambda+d
\qquad
M=N^{\mathsf T}I^{-1}N\tag{11.20, 11.21}
$$

### ⚠️ 解的存在性与唯一性

$M$ **对称半正定**，**当且仅当接触法向线性无关时正定**。

| 情形 | 解 |
|---|---|
| $d\ge0$ | $\lambda=0$ 是一个可行解 |
| $M$ **正定** | **恰好一个**解 |
| $M$ **奇异** | 可能**一个、没有、或无穷多个**解 |

> **若有多个 $\lambda$ 解，它们给出的 $a$ 是相同的。**
>
> 🔑 **"无解"或"多解"只可能在接触法向线性相关时出现。**

### ⭐ 例 11.1：一个无解的接触系统（原书图 11.3）

**设置**：圆盘在槽中滑动。直线 $L$ 左侧槽壁笔直、宽度恰等于圆盘直径，
圆盘与槽有两个点接触，法向 $n_1$、$n_2$ 满足 $n_2=-n_1$（**线性相关**）。
$L$ 右侧槽壁是向内弯的圆弧，槽变窄，所以圆盘中心 $D$ **不能越过 $L$ 到右边**。

**悖论**：当 $D$ 落在 $L$ 上时，圆盘不能再向右，
**但两个接触法向都竖直朝上下，根本无法影响圆盘的水平运动**。
若施加一个向右推的力——**向右运动必须不发生，但接触力无法阻止它**。

**数学表现**：$n_2=-n_1$ ⟹ $M=\begin{bmatrix}m&-m\\ -m&m\end{bmatrix}$，故
$\dot\zeta_1+\dot\zeta_2=d_1+d_2$。有解要求 $\dot\zeta_1+\dot\zeta_2\ge0$，
所以**无解的条件是 $d_1+d_2<0$**。

> 🔑 **原书给的实践解法**：
> *"the practical solution is to **allow $D$ to move slightly to the right of $L$**,
> so that the two contact normals change direction and become able to resist any
> further motion to the right. ... In both cases, **the simulator must tolerate
> a small amount of penetration** between the disc and the slot."*
>
> **这解释了为什么所有物理引擎都允许少量穿透**——不是实现上的偷懒，
> 而是**刚体接触模型在数学上就可能无解**，必须靠一点穿透让法向转向。

---

## 11.3 含接触的一般刚体系统

### 与单刚体情形的两个差别

1. **运动方程在关节空间中**（$H\ddot q+C=\tau+\tau^c$，11.23）
   ⟹ 必须把接触力和加速度表达式**翻译到关节空间**；
2. **接触可以发生在不同的刚体对之间** ⟹ 必须标明每个接触涉及哪两个刚体。

**术语**（沿用第 3、4 章）：接触 $i$ 的**前驱** $pc(i)$ 与**后继** $sc(i)$；
接触力定义为**从前驱传给后继**；分离加速度定义为**后继相对前驱**的加速度。
**$n_i$ 由前驱指向后继。**

### ⭐ 关节空间的接触法向 $t_i$

由 $\tau=J_i^{\mathsf T}f$：

$$
\boxed{\ t_i=(J_{sc(i)}-J_{pc(i)})^{\mathsf T}n_i\ }\tag{11.25}
$$

> **$t_i$ 可以看作接触 $i$ 的接触法向向量在关节空间中的表达。**

于是一切都变得和 §11.2 平行：

$$
\tau^c=T\lambda
\qquad
\zeta_i=t_i^{\mathsf T}\dot q
\qquad
\dot\zeta=T^{\mathsf T}\ddot q+\dot T^{\mathsf T}\dot q\tag{11.32, 11.28, 11.33}
$$

$$
\dot\zeta=M\lambda+d
\qquad
\boxed{M=T^{\mathsf T}H^{-1}T}
\qquad
d=T^{\mathsf T}H^{-1}(\tau-C)+\dot T^{\mathsf T}\dot q\tag{11.34, 11.35, 11.36}
$$

$\dot T^{\mathsf T}\dot q$ 由**速度乘积加速度** $a^{vp}$（第 8 章的副产品）算出：

$$
\dot t_i^{\mathsf T}\dot q=n_i\cdot(a^{vp}_{sc(i)}-a^{vp}_{pc(i)})
+\dot n_i\cdot(v_{sc(i)}-v_{pc(i)})\tag{11.30}
$$

> 🔑 **原书指出式 11.34 与 11.20 的唯一差别**：
> **式 11.21 的 $M$ 秩至多为 6；式 11.35 的 $M$ 的秩没有上限。**
> 求解过程完全一样。

> 💡 **认出 $M=T^{\mathsf T}H^{-1}T$**：这就是第 7 章 §7.2.1 的
> **操作空间惯性矩阵的逆** $\Lambda^{-1}=JH^{-1}J^{\mathsf T}$，
> 也是第 8 章 §8.5 方法 2 的 $A=KH^{-1}K^{\mathsf T}$。
> **三章里出现的是同一个量**，只是约束从双边变成了单边。

---

## 11.4 ⭐ 不等式约束的四个状态（本章最实用的一张图）

$$
\phi(q)\ge0\tag{11.37}
$$

称为**不等式约束**或**单边约束 (unilateral constraint)**——
"若两面接触，则一个方向的运动被允许、相反方向不被允许"。

**原书图 11.4 的四个状态**：

| 状态 | 条件 | 含义 | 对瞬时动力学的影响 |
|---|---|---|---|
| **1** | $\phi>0$ | 分离，无接触 | **无**（约束**不激活**） |
| **2a** | $\phi=0,\ \dot\phi<0$ | **碰撞的瞬间** | **产生冲量**，速度阶跃；只持续一个瞬间 |
| **2b** | $\phi=0,\ \dot\phi>0$ | 接触但正在飞离 | **无**（同状态 1）；只持续一个瞬间 |
| **3** | $\phi=0,\ \dot\phi=0$ | **持续接触** | $\ddot\phi\ge0$，有接触力 |

**状态转移**：

- **2a 之后**：视碰撞是**弹性**（反弹）还是**塑性**（不反弹）而进入 **2b** 或 **3**
- **2b 只可能出现在**：冲量刚发生之后，或**几何失接触**的瞬间（§11.6）
- **状态 3 的两种去向**：两面加速分开（$\ddot\phi>0$，下一瞬间转入状态 1，接触力为零），
  或保持接触（$\ddot\phi=0$，留在状态 3，接触力阻止 $\ddot\phi<0$）

$$
\boxed{\ \textbf{§11.1–11.3 讨论的全部是状态 3}\ }
$$

### 与前面各节的联系

对 $\phi(q)$ 求导：

$$
\dot\phi=t^{\mathsf T}\dot q
\qquad
\ddot\phi=t^{\mathsf T}\ddot q+\dot t^{\mathsf T}\dot q
\qquad
t^{\mathsf T}=\frac{\partial\phi}{\partial q}\tag{11.38, 11.39, 11.40}
$$

> 🔑 **对照式 11.28、11.29 可见**：
>
> $$
> \zeta_i\leftrightarrow\dot\phi_i,\qquad \dot\zeta_i\leftrightarrow\ddot\phi_i,
> \qquad t_i\leftrightarrow\left(\frac{\partial\phi_i}{\partial q}\right)^{\mathsf T}
> $$
>
> **接触法向在关节空间中就是约束函数的梯度。**

---

## 11.5 求解接触方程

### LCP

$$
\boxed{\ \dot\zeta=M\lambda+d,\quad \dot\zeta\ge0,\quad\lambda\ge0,
\quad\dot\zeta^{\mathsf T}\lambda=0\ }\tag{11.42}
$$

> **这组式子定义了数学上的标准问题：线性互补问题 (linear complementarity problem, LCP)。**
> 文献：Cottle et al. (1992)、Cottle & Dantzig (1968)；
> $M$ 正定的特例算法见 Featherstone (1987)。

### 对称半正定 LCP 的三条性质

$$
\boxed{
\begin{aligned}
&1.\ M\text{ 正定}\Rightarrow\textbf{恰好一个解}；\text{否则可能无解、一解或无穷多解}\\
&2.\ \text{若有解，}\textbf{}T\lambda\text{ 与 }\dot\zeta\text{ 是唯一的}\text{（多解只在 }\lambda\text{ 上差
一个 }\mathrm{null}(T)\text{ 中的量）}\\
&3.\ \text{式 11.42 }\textbf{等价于}\text{一个二次规划}
\end{aligned}}
$$

**等价的 QP**：

$$
\min_{\lambda}\ \tfrac12\lambda^{\mathsf T}M\lambda+\lambda^{\mathsf T}d
\qquad\text{s.t.}\quad \lambda\ge0\tag{11.43}
$$

（式 11.42 正是式 11.43 解的 **Kuhn-Tucker 条件**。）

> 🔑 **原书给的实用理由**：
> *"11.43 has the advantage that **software to solve quadratic programs is
> more widely available than software to solve LCPs**."*
> —— **这就是现代物理引擎多用 QP 求解器的原因。**

**另一种形式**（Lötstedt 1982），**对 $\ddot q$ 而非 $\lambda$ 最小化**：

$$
\min_{\ddot q}\ \tfrac12\big(\ddot q-H^{-1}(\tau-C)\big)^{\mathsf T}H\big(\ddot q-H^{-1}(\tau-C)\big)
\quad\text{s.t.}\quad T^{\mathsf T}\ddot q+\dot T^{\mathsf T}\dot q\ge0\tag{11.44}
$$

> 💡 **原书对它的评价**：*"This formulation will produce a **unique solution**,
> if one exists, and it may have a **computational advantage if $n_c>n$**,
> but **recovering $\lambda$ from $\ddot q$ is not straightforward**.
> Equation 11.44 can be regarded as a **generalization of Gauss' principle of
> least constraint**, making it applicable to rigid-body systems with inequality
> constraints."*
>
> **这正是 MuJoCo 等引擎所用表述的理论来源。**

### ⭐ 表 11.1：接触动力学仿真的通用流程

> **动机**：接触多时，解 LCP/QP 的代价**远大于**用第 8 章的技术施加等价的等式约束。
> **因此要尽量减少求解 LCP/QP 的次数。**

```
1. 向前积分，把「活动接触」当作等式约束，忽略其他所有接触
2. 积分过程中监测两类事件：
     (a) 几何事件（接触的产生与失去）
     (b) 负的接触力
3. 若最近一步中检测到事件，把系统插值回到最早那个事件的时刻
4. 若事件是接触产生（如碰撞）：施加冲量动力学（§11.7），
     并确定新的「当前接触集」与「活动接触集」
5. 若事件是几何失接触：从两个集合中移除失去的接触
6. 若事件是负接触力：建立并求解式 11.42 或 11.43，确定新的两个集合
7. 回到步骤 1
```

### 两个集合的定义

| 集合 | 定义 |
|---|---|
| **当前接触 (current contacts)** | 所有处于**状态 3** 且满足 $\dot\zeta_i=0$ 的接触 |
| **活动接触 (active contacts)** | 当前接触的**子集**，其法向构成 $\mathrm{range}(T)$ 的一组**基**，且包含所有 $\lambda_i>0$ 的接触 |

> 🔑 **这个设计的核心命题**：
> *"If we treat the active contacts as equality constraints, then the following
> statement is true: **for as long as the contact force variables remain non-negative,
> there is no change in the state of contact**."*
>
> **所以只在检测到活动接触上出现负 $\lambda_i$、或发生碰撞时，才需要解 LCP/QP。**
> 其余时间用便宜得多的等式约束方法（第 8 章）。

**基本解 (basic solution)**：LCP/QP 求解器多返回**非零变量个数最少**的解。
基本解的性质是**正接触力对应的法向线性无关**。

> ⚠️ **原书的一个边角情形**：*"On rare occasions, this set will be insufficient to
> form a basis, and will need to be supplemented with contacts satisfying
> $\dot\zeta_i=\lambda_i=0$."*

### 插值积分步（很实用的技巧）

**问题**：步骤 3 要"把系统插值回事件时刻"，怎么做？

**原理**：许多数值积分法的工作方式是对 $\dot y$ 在 $[t_0,t_1]$ 上**拟合一个多项式**
并加上它的积分：

$$
\dot y(t)=p(t),\qquad y_1=y_0+\int_{t_0}^{t_1}p(t)\,dt\tag{11.45}
$$

$p$ 的次数比方法的阶低 1。于是**任意中间时刻的值**可以算出来：

$$
y(t)=y_0+\int_{t_0}^{t}p(t)\,dt\tag{11.46}
$$

令 $\delta=(t-t_0)/h\in[0,1]$：

**Euler 法**：$y_1=y_0+hk$，$k=f(y_0,t_0)$

$$
y(\delta)=y_0+hk\delta\tag{11.49}
$$

**Heun 法**：$k_1=f(y_0,t_0)$，$k_2=f(y_0+\tfrac23hk_1,\ t_0+\tfrac23h)$，
$y_1=y_0+h(\tfrac14k_1+\tfrac34k_2)$

$$
y(\delta)=y_0+h\big((\delta-\tfrac34\delta^{2})k_1+\tfrac34\delta^{2}k_2\big)\tag{11.51}
$$

> 💡 **价值**：这样仿真器就能**在已完成的积分步内部**算出任意中间时刻的状态，
> 从而**把最近一步缩短任意需要的量**——**不用重算这一步**。

---

## 11.6 接触几何

### ⚠️ 多面体不适合做动力学

> **最常用的形状表示是多面体**，围绕它有大量技术积累与软件
> （创建、编辑、显示、高效碰撞检测）。
> **不幸的是，多面体并不适合做动力学。**

**原书图 11.5 的例子**：真球 vs 多面体近似的球，沿斜坡滚下。

| | 行为 |
|---|---|
| 真球 | **平滑滚动** |
| 多面体球 | **每个新顶点撞上斜面时都产生一次碰撞** |

**两难**：

- 当作**塑性**碰撞 ⟹ **每次碰撞都损失能量**
- 当作**弹性**碰撞 ⟹ **多面体开始弹跳**

> 🔑 **这是所有物理引擎都要面对的问题**。实践中的对策通常是：
> 对滚动体用**解析曲面**（球、圆柱、胶囊）而非网格；
> 或引入接触阻尼/软接触（§11.8）来吸收这些伪碰撞。

本节还讨论**几何事件**（接触的产生与失去）以及**如何用点接触表示线接触和面接触**。

---

## 11.7 冲量动力学

### 定义与基本性质

$$
\iota=\int_{t_0}^{t_1}f(t)\,dt=h(t_1)-h(t_0)\tag{11.52, 11.53}
$$

> **刚体在给定时间区间内受到的冲量 = 它在该区间内动量的净变化。**

**刚体动力学关心的是一种特殊冲量**——把硬度推到无穷大时的极限：
接触力发散到无穷、时段收敛到零，但**冲量保持有限**：

$$
\iota=\lim_{\delta t\to0}\int_{t}^{t+\delta t}\frac{f(t)}{\delta t}\,dt\tag{11.54}
$$

**三条基本性质**：

$$
\boxed{
\begin{aligned}
&1.\ \text{它们是}\textbf{力向量}\text{（}F^6\text{ 的元素）}\\
&2.\ \text{它们由刚体间的}\textbf{碰撞（接触产生）}\text{引起}\\
&3.\ \text{它们引起刚体速度的}\textbf{阶跃变化}
\end{aligned}}
$$

### 三个层次的冲量运动方程

$$
\boxed{\ \iota=I\,\Delta v\ }\tag{11.55 单刚体}
$$

$$
\boxed{\ \iota=I^{A}\,\Delta v\ }\tag{11.56 铰接体的柄}
$$

$$
\boxed{\ u=H\,\Delta\dot q\ }\tag{11.57 一般刚体系统}
$$

$$
u=J_i^{\mathsf T}\iota
\qquad
\Delta v_i=J_i\,\Delta\dot q\tag{11.58, 11.59}
$$

> 🔑 **式 11.56 很值得注意**：**铰接体惯性 $I^A$ 也满足冲量方程**。
> 这是因为 $I^A$ 本来就是"加速度→力"的映射，
> 而冲量方程正是它的时间积分形式。第 7 章的 $I^A$ 在这里直接可用。

**推导的两个关键步骤**（原书）：
$I^{-1}$ 能提到积分号外，因为 $I$ 只随位置变化而位置是时间的连续函数；
第二个极限为零，因为被积函数在 $\delta t\to0$ 时保持有限。

### ⚠️ 原书的重要警告

> *"**Caution:** Equations 11.56, 11.57 and 11.58 are **mathematically correct,
> but they do not necessarily offer an accurate prediction of the impulsive behaviour
> of a real multibody system**. This is because real impacts involve physical processes
> that are **not modelled by the rigid-body assumption**, such as **the propagation of
> compression waves through solids and across joint bearings**."*

> 🔑 **这条提醒非常诚实**：多体系统的碰撞，冲量会通过关节轴承以压缩波的形式传播，
> 刚体模型完全没有描述这个过程。
> **所以多体碰撞的仿真结果应当被当作近似，而不是预测。**

### 两体碰撞

$$
\Delta v_1=-I_1^{-1}n\lambda
\qquad
\Delta v_2=I_2^{-1}n\lambda\tag{11.60, 11.61}
$$

$$
\zeta=n\cdot(v_2-v_1)
\qquad
\Delta\zeta=n\cdot(\Delta v_2-\Delta v_1)\tag{11.62, 11.63}
$$

必须满足 $\zeta\le0$ 且 $\zeta+\Delta\zeta\ge0$。**还差一个方程**——由**恢复系数** $e\in[0,1]$ 提供：

$$
e=\frac{\zeta+\Delta\zeta}{-\zeta}
\quad\Longrightarrow\quad
\Delta\zeta=-(1+e)\zeta\tag{11.64}
$$

$$
\boxed{\ \lambda=\frac{-(1+e)\,n\cdot(v_2-v_1)}{n\cdot(I_1^{-1}+I_2^{-1})\,n}\ }\tag{11.65}
$$

### 摩擦（原书列出的五个效应）

允许 Coulomb 摩擦后，一般要考虑：

$$
\boxed{
\begin{aligned}
&1.\ \textbf{切向摩擦}\ (\text{tangential friction})\\
&2.\ \textbf{扭转摩擦}\ (\text{torsional friction})\\
&3.\ \textbf{切向恢复}\ (\text{tangential restitution})\\
&4.\ \textbf{扭转恢复}\ (\text{torsional restitution})\\
&5.\ \text{切向/扭转}\textbf{滑动是否在接触期结束前被制止}
\end{aligned}}
$$

**扭转效应的来源**（原书解释得很好）：

> **扭转效应指反抗绕接触法线方向角运动的摩擦力偶。**
> 它们之所以存在，是因为**真实物体并非真正刚性**——
> 接触面的局部变形使**点接触变成了面接触**。
>
> **因此扭转效应在较软的物体上更明显。** 即便如此，
> **扭转效应通常比切向效应小一个量级，常常可以忽略。**

**切向/扭转恢复**由接触体在切向/扭转方向上**局部弹性变形的回复**引起。

（更多细节见 Brach (1991)。）

---

## 11.8 软接触 (Soft Contact)

### ⭐ 为什么要用软接触（原书给了两类理由）

**理由 1：物理上就需要**——很多系统含软硬混合的物体、
或表面软的硬物体、或整体可视为刚体但接触时有显著局部变形的物体。

**理由 2：⭐ 实现上简单得多**：

$$
\boxed{
\begin{aligned}
&1.\ \textbf{碰撞瞬间没有冲量}\ \Longrightarrow\ \text{不需要冲量动力学计算}\\
&2.\ \textbf{失接触可由位置和速度数据判定}\ \Longrightarrow\ \textbf{不需要解 LCP}\\
&3.\ \text{在柔顺接触上}\textbf{实现 Coulomb 摩擦也更容易}
\end{aligned}}
$$

> 🔑 **原书的结论很实在**：*"For these reasons, it may be **more sensible to use
> compliant contact in place of rigid contact, even in cases where rigid contact is
> a better model** of the physical system, **if implementation effort is an issue**."*

**主要缺点**：

> **引入高频动力学**（柔顺表面里有刚硬的弹簧）。
> 若这些动力学迫使积分例程取更小的步长，**仿真速度就会下降**。

### 建模技巧：把柔顺表面建成一阶动力系统

$$
\boxed{\ \text{柔顺表面 = }\textbf{含弹簧和阻尼器但}\textbf{无质量}\text{的系统}\ }
$$

在这种系统里，**施加力产生速度**，**强加速度遇到反作用力**。

> 🔑 **最大的好处**（原书原话）：
> *"it **separates the rigid-body dynamics from the contact dynamics**:
> the contact forces can be calculated from the position and velocity variables,
> and then supplied to the rigid-body dynamics routine as a collection of
> **known force inputs**. It is therefore possible to use **any standard
> forward-dynamics algorithm** with this contact model, such as those described
> in Chapters 6, 7 and 8."*
>
> **也就是说：软接触把接触问题降级成了"外力输入"，第 6/7/8 章的算法原样可用。**

### 一维模型（原书图 11.10）

刚体 $B$（位置 $p$）+ 无质量表面片 $S$（位置 $z$），$S$ 经弹簧（刚度 $K$）
和阻尼器（阻尼 $D$）连到固定基座。

$S$ 无质量 ⟹ 净力为零：$Kz+D\dot z+f_c=0$，即

$$
\boxed{\ \dot z=-(Kz+f_c)/D\ }\tag{11.73}
$$

> ⚠️ **注意这是一阶微分方程**（含 $z$ 和 $\dot z$，**不含 $\ddot z$**），
> 且**接触力线性关联于 $S$ 的速度，而非加速度**。

**接触力**：

$$
\boxed{\ f_c=\begin{cases}
0&p>z\ (\text{无接触})\\[3pt]
\max(0,\ -Kz-D\dot p)&p=z
\end{cases}\ }\tag{11.74}
$$

> **推导**：若接触要持续，当前时刻必须 $\dot p=\dot z$，
> 所需接触力是 $-Kz-D\dot p$。
> 若这个表达式 $\ge0$，它就是正确的 $f_c$，接触确实持续；
> 若为负，正确的 $f_c$ 是零，且 $\dot p>\dot z$，说明**接触正在断开**。
>
> 💡 **`max(0, ·)` 就是单边性的全部体现**——比 LCP 简单太多。

### 三个扩展

1. **两个运动物体之间的柔顺接触**：把 $p$ 重定义为二者的**相对位置**
2. **两个柔顺表面之间的接触**：每个表面用一对弹簧-阻尼器和**一个状态变量 $z$**；
   **若两表面的 $K/D$ 相同，两个柔度可以合并成一个**
3. ⭐ **免去显式状态变量的技巧**：式 11.73 在 $f_c=0$ 时有简单解析解，故
   $$z=p\ (\text{接触时}),\qquad z=z_0e^{-K(t-t_0)/D}\ (\text{其他时候})$$
   其中 $t_0$、$z_0$ 是最近一次失接触的时刻与当时的 $z$。
   > **这免去了数值积分式 11.73 的需要，可能提高仿真效率。**

也可以引入**非线性**弹簧与阻尼器。

---

## 本章要点回顾

1. **接触约束是单边的**：只能推不能拉，$\phi\ge0$、$\lambda\ge0$。
2. **有限力 ⟹ $\zeta=0$**；**$\zeta\ne0$ 意味着有冲量**——这条把本章分成接触与碰撞两半。
3. **互补条件** $\dot\zeta\ge0,\lambda\ge0,\dot\zeta\lambda=0$ 简洁地统一了"断开"与"持续"。
4. ⭐ **加速度只是施加力的分段线性函数**——这是接触动力学与等式约束动力学的根本差别。
5. **$M$ 正定 ⟺ 接触法向线性无关**；否则可能无解（例 11.1 的槽中圆盘）。
6. ⭐ **无解的存在性解释了为什么所有物理引擎都容忍少量穿透。**
7. **$t_i=(J_{sc}-J_{pc})^{\mathsf T}n_i$ 是关节空间的接触法向**，也就是 $\partial\phi_i/\partial q$。
8. **$M=T^{\mathsf T}H^{-1}T$** 与第 7 章的 $\Lambda^{-1}$、第 8 章的 $KH^{-1}K^{\mathsf T}$ 是同一个量。
9. ⭐ **四个状态**（图 11.4）：1 无接触、2a 碰撞、2b 飞离、3 持续接触；**§11.1–11.3 讲的全是状态 3**。
10. **LCP 等价于 QP**（式 11.43），而 **QP 求解器比 LCP 求解器更常见**。
11. **式 11.44 是 Gauss 最小约束原理向不等式约束的推广。**
12. ⭐ **表 11.1 的七步流程 + 两个集合**：只在出现负 $\lambda$ 或碰撞时才解 LCP。
13. **插值积分步**可以在不重算的情况下回退到事件时刻。
14. **多面体不适合做动力学**（滚动会变成一连串碰撞）。
15. **$\iota=I\Delta v$、$\iota=I^A\Delta v$、$u=H\Delta\dot q$**；⚠️ 但**多体碰撞的刚体预测不可靠**。
16. **软接触把接触降级成外力输入**，第 6/7/8 章算法原样可用；代价是引入高频动力学。

---

## 易错点

1. **以为接触动力学是"带约束的线性问题"**——它是分段线性的，需要组合搜索。
2. **不容忍任何穿透**——数学上可能无解（例 11.1）。
3. **把状态 2b 当成需要处理的接触**——它和状态 1 一样无影响。
4. **每一步都解 LCP**——按表 11.1，只在必要时解。
5. **用多面体近似球做滚动仿真**——会产生一连串伪碰撞。
6. **把多体碰撞的刚体计算结果当成物理预测**——见 §11.7 的警告。
7. **软接触参数取得过硬**——引入高频动力学，逼迫积分器减小步长。
8. **忘了 $n$ 指向表面外侧 / 由前驱指向后继**——符号错会得到"吸引"的接触力。

## 与其他章的联系

- ← 第 3 章 §3.4：不等式约束在约束分类树中的位置
- ← 第 6 章：$H$、$H^{-1}$ 的高效计算（$M=T^{\mathsf T}H^{-1}T$ 要用）
- ← 第 7 章 §7.2.1：$\Lambda^{-1}=JH^{-1}J^{\mathsf T}$ 就是本章的 $M$；式 7.12 的 $I^A$ 用于式 11.56
- ← 第 8 章：接触约束与闭环约束**数学同构**，差别只在单边性；
  $a^{vp}$ 是 §8.6 的副产品；活动接触当作等式约束处理时用第 8 章的方法
- ← 第 9 章：混合动力学；浮动基（接触多发生在浮动基系统上）
- ← 第 10 章：软接触引入的高频动力学 ⟹ 截断误差与步长的权衡

---

## ✍️ 我的理解

<!-- 建议：说清「为什么接触让加速度只是分段线性函数」及其后果 -->

## ❓ 疑问与待办

- [ ] 实现单点接触的式 11.10，验证"接触断开/持续"的分支
- [ ] 复现例 11.1 的槽中圆盘，构造 $d_1+d_2<0$ 的无解情形
- [ ] 用 QP 求解器（式 11.43）实现多点接触，与 LCP 求解器对比
- [ ] 实现表 11.1 的七步流程，统计一次仿真中真正解 LCP 的次数
- [ ] 实现 Euler / Heun 的插值多项式（式 11.49、11.51）
- [ ] 实现 §11.8 的软接触模型，含免状态变量的解析技巧
- [ ] 用两体碰撞公式 11.65 验证 $e=1$ 时能量守恒、$e=0$ 时分离速度为零

## 📌 与原文的出入

<!-- 本笔记按原书 pp.213–240 逐节撰写。
     此前的笔记完全没有本章——误以为全书只有 10 章，
     并把接触/碰撞的零散内容错放在第 9 章。 -->
