# 第 8 章 闭环系统 (Closed Loop Systems)

> **原书 pp. 141–170**，共 13 节——**全书最长的一章**。

> **一句话概括**：含运动学回路的系统。做法是"**生成树 + 闭环约束**"，
> 得到一个 $(n+n_c)$ 维的鞍点方程。回路带来一整套新问题：
> **可变自由度、位形歧义、过约束**，而且**逆动力学的性质发生了根本变化**。

## 原书开篇的定位

> *"The presence of kinematic loops brings **a new level of complexity**:
> new formulations are required; new problems arise; the systems exhibit
> new behaviours; and **the inverse dynamics problem, in particular,
> changes its nature** in the presence of kinematic loops."*

**两大策略**（对应第 3 章 §3.2 的方法 1 和方法 2）：

| 策略 | 做法 | 结果 | 本书 |
|---|---|---|---|
| **1** | 从一组**无约束刚体**出发，同时施加所有关节约束 | **大而稀疏**的矩阵方程 | §8.13（简述） |
| **2** | 从**生成树**出发，施加**闭环约束** | **典型闭环系统的最佳选择** | **本章主体** |
| （3） | §7.5 的装配法 | 本章不考虑 | — |

## 原书节次对照

| 节 | 标题 | 页 |
|---|---|---|
| 8.1 | Equations of Motion | 141 |
| 8.2 | Loop Constraint Equations | 143 |
| 8.3 | Constraint Stabilization | 145 |
| 8.4 | Loop Joint Forces | 148 |
| 8.5 | Solving the Equations of Motion | 149 |
| 8.6 | Algorithm for $C-\tau^a$ | 152 |
| 8.7 | Algorithm for $K$ and $k$ | 154 |
| 8.8 | Algorithm for $G$ and $g$ | 156 |
| 8.9 | Exploiting Sparsity in $K$ and $G$ | 158 |
| 8.10 | **Some Properties of Closed-Loop Systems** | 159 |
| 8.11 | Loop Closure Functions | 161 |
| 8.12 | **Inverse Dynamics** | 164 |
| 8.13 | Sparse Matrix Method | 166 |

## ⚠️ 下标约定（原书专门规定，读本章前必须记住）

$$
\boxed{
\begin{aligned}
i,j&：\text{树关节 / 刚体编号}\quad(1..N_B)\\
l&：\text{回路编号}\quad(1..N_L)\\
k&：\text{闭合回路 }l\text{ 的那个}\textbf{闭环关节}\quad(N_B{+}1..N_J)\\
&\qquad\text{且恒有 }\boxed{k=l+N_B}
\end{aligned}}
$$

**$k$ 与 $l$ 同时出现时，一律假定 $k=l+N_B$。**

---

## 8.1 运动方程

### 四步流程

$$
\boxed{
\begin{aligned}
&1.\ \text{写出}\textbf{生成树}\text{的运动方程}\\
&2.\ \text{加上表示}\textbf{闭环关节对树施力}\text{的项}\\
&3.\ \text{写出闭环关节对树施加的}\textbf{运动约束}\\
&4.\ \text{把两个方程合并}
\end{aligned}}
$$

生成树 $n$ 自由度、闭环关节施加 $n_c$ 个约束 ⟹ **$n+n_c$ 个方程、$n+n_c$ 个未知量**。

$$
n=\sum_{i=1}^{N_B}n_i,\qquad n_c=\sum_{k=N_B+1}^{N_J}n_{ck}\tag{8.1}
$$

### 从树到闭环

生成树：$H\ddot q+C=\tau$（8.2）。闭环系统与生成树的**唯一差别是多了闭环关节**，故

$$
\boxed{\ H\ddot q+C=\tau+\tau^{c}+\tau^{a}\ }\tag{8.3}
$$

| 项 | 含义 | 已知？ |
|---|---|---|
| $\tau^{c}$ | 闭环关节产生的**约束力** | **未知**（用 $n_c$ 个变量表达） |
| $\tau^{a}$ | 闭环关节产生的**主动力** | **已知** |

> ⚠️ **$\tau^a$ 是我此前完全遗漏的一项。** 它来自作用在闭环关节上的
> 弹簧、阻尼器、驱动器等。
>
> **术语**：若某关节上没有这类力，称该关节**被动 (passive)**；
> **若所有闭环关节都被动，则 $\tau^a=0$**，
> 此时可以**把算 $\tau^a$ 的代码整段省掉**。

### 合并成鞍点方程

$$
\tau^{c}=K^{\mathsf T}\lambda\tag{8.5}
\qquad\qquad
K\ddot q=k\tag{8.4}
$$

$$
\boxed{\ \begin{bmatrix}H&K^{\mathsf T}\\ K&0\end{bmatrix}
\begin{bmatrix}\ddot q\\ -\lambda\end{bmatrix}
=\begin{bmatrix}\tau-C+\tau^{a}\\ k\end{bmatrix}\ }\tag{8.6}
$$

**系数矩阵的性质**（原书明确给出）：

- **对称**，但**不正定**（⟹ 不能用 Cholesky）
- **秩 $=n+\mathrm{rank}(K)$**
- 满秩 ⟹ $\ddot q$ 与 $\lambda$ 都唯一
- **不满秩 ⟹ $\lambda$ 的部分元素不定，但 $\ddot q$ 仍然唯一**

---

## 8.2 环约束方程

### 推导

闭环关节 $k$ 的关节速度：

$$
v_{Jk}=v_{s(k)}-v_{p(k)}\tag{8.7}
$$

约束（**为简单起见假定所有闭环关节的偏置速度 $\sigma_k=0$**）：

$$
T_k^{\mathsf T}v_{Jk}=0\tag{8.8}
\quad\Longrightarrow\quad
T_k^{\mathsf T}(v_{s(k)}-v_{p(k)})=0\tag{8.9}
$$

求导得加速度约束：

$$
T_k^{\mathsf T}(a_{s(k)}-a_{p(k)})+\dot T_k^{\mathsf T}(v_{s(k)}-v_{p(k)})=0\tag{8.10}
$$

用体雅可比（第 4 章式 4.6）$v_i=J_i\dot q$ 和

$$
a_i=J_i\ddot q+\dot J_i\dot q=J_i\ddot q+a^{vp}_i\tag{8.13}
$$

> 🔑 **$a^{vp}_i$ = "速度乘积"加速度**：所有树关节加速度变量为零时 body $i$ 的加速度。
> **它是计算 $C$ 时的一个副产品**（§8.6），不需要额外计算。

代入得

$$
T_k^{\mathsf T}(J_{s(k)}-J_{p(k)})\,\ddot q=k_l\tag{8.14}
$$

$$
k_l=-T_k^{\mathsf T}(a^{vp}_{s(k)}-a^{vp}_{p(k)})-\dot T_k^{\mathsf T}(v_{s(k)}-v_{p(k)})\tag{8.15}
$$

### ⭐ 环雅可比与 $\epsilon_{lj}$

$$
\boxed{\ J_{Ll}=J_{s(k)}-J_{p(k)}
=\begin{bmatrix}\epsilon_{l1}S_1&\epsilon_{l2}S_2&\cdots&\epsilon_{lN_B}S_{N_B}\end{bmatrix}\ }\tag{8.18}
$$

**$\epsilon_{lj}$ 怎么取**（原书给了一个非常清楚的几何刻画）：

> **定义回路 $l$ 的"根" = $p(k)$ 与 $s(k)$ 的共同祖先中编号最大的那个刚体。**
> 参与回路 $l$ 的树关节集合 = 两个子集之并：
> 连接**根到 $p(k)$** 的那些，与连接**根到 $s(k)$** 的那些。

$$
\epsilon_{lj}=\begin{cases}
-1&\text{关节 }j\text{ 在"根}\to p(k)\text{"路径上}\\
+1&\text{关节 }j\text{ 在"根}\to s(k)\text{"路径上}\\
0&\text{其他}
\end{cases}
$$

> 💡 **直观理解**：绕回路走一圈，沿一支下行、沿另一支上行，
> 符号相反的两段合起来构成"闭合条件"。$\epsilon$ 就是走向的符号。

于是 $K$ 的分块表达（实际用于计算的形式）：

$$
K_{lj}=\epsilon_{lj}\,T_k^{\mathsf T}S_j\tag{8.20}
$$

> 🔑 **式 8.20 的实用价值**：**$K$ 不需要先构造 $J_{s(k)}$ 和 $J_{p(k)}$ 再相减**——
> 只需沿两条路径走，把 $T_k^{\mathsf T}S_j$ 填到对应位置并带上符号。
> 这就是 §8.9 "利用 $K$ 的稀疏性"的基础：
> **$K$ 的第 $l$ 行只在"参与回路 $l$ 的关节"处非零。**

---

## 8.3 约束稳定化 (Constraint Stabilization)

### 问题

> 式 8.16 理论上正确，但**数值积分时不稳定**。原理上它表现得像
> $$\ddot e=0$$
> （$e$ = 闭环位置误差），但实践中表现得像
> $$\ddot e=\text{noise}$$
> 因此**加速度误差保持很小，却没有任何东西阻止位置和速度误差无界累积**。

### Baumgarte 稳定化

给式 8.4 加一个稳定化项：

$$
K\ddot q=k+k_{stab}\tag{8.21}
$$

使约束方程在数值积分下表现为

$$
\boxed{\ \ddot e+2\alpha\dot e+\beta^{2}e=\text{noise}\qquad(\alpha>0,\ \beta\ne0)\ }
$$

$$
k_{stab}=-2\alpha\begin{bmatrix}T_{N_B+1}^{\mathsf T}v_{JN_B+1}\\ \vdots\\ T_{N_J}^{\mathsf T}v_{JN_J}\end{bmatrix}
-\beta^{2}\begin{bmatrix}\delta_1\\ \vdots\\ \delta_{N_L}\end{bmatrix}\tag{8.22}
$$

$\delta_l$ 度量生成树的位形**违反闭环关节约束的程度**：$\delta_l=T_k^{\mathsf T}d_l$，
$d_l$ 是度量回路 $l$ 位置误差的空间向量。

### ⭐ $T_{stab}$ 怎么选（原书给了非常实用的指导）

$$
\alpha=\beta=1/T_{stab}
$$

> *"A reasonable strategy is to ask **how many snapshots per second of the moving
> system you would need in order not to miss anything important**, and choose a value
> of $T_{stab}$ that is **similar to, or a little shorter than, the period between snapshots**."*

**原书给的具体数字（大型工业机器人）**：

| $T_{stab}$ | 评价 |
|---|---|
| $1$ | **太大** |
| $0.1$ | **合理** |
| $0.01$ | 有点小，但仍可接受 |

> ⚠️ **$T_{stab}$ 选得太小 ⟹ 微分方程变得不必要地刚性 (stiff)。**
> 但原书也说：*"The exact value of $T_{stab}$ is **not critical**—
> changing it by a factor of 2 makes only a small difference."*

### 🔑 一条极重要的观念纠正

> *"There is a temptation to choose $T_{stab}$ as small as possible... in the belief
> that this will maximize the accuracy of the simulation. **This is a bad strategy.**
> **The purpose of constraint stabilization is to achieve stability, not accuracy.**
> If the simulation is not accurate enough, then the best way to improve it is to use
> **a better integration method and/or a shorter integration time step**."*

$$
\boxed{\ \textbf{稳定化的目的是稳定，不是精度。想要精度请换积分器或减小步长。}\ }
$$

### 怎么算 $\delta$

把生成树位形决定的 ${}^{s}X_p$ 分解为

$$
{}^{s}X_p=X_{err}\,X_J\tag{8.25}
$$

$X_J$ **精确满足**关节约束，$X_{err}$ 表示（假定很小的）**约束误差**。
这个分解必须**用关节类型的知识符号地求解**，每种关节类型的表达式不同。

由于 $X_{err}$ 表示小位移，存在运动向量 $d\in M^6$ 使

$$
\mathbf 1-d\times\simeq X_{err}\tag{8.26}
\qquad\qquad
\delta=T^{\mathsf T}d\tag{8.27}
$$

> **$d$ 是"后继坐标系应该在的位置"到"它实际在的位置"的（近似）位移。**

**原书表 8.1** 给出零自由度、转动、移动、圆柱、球关节的 $\delta$ 公式。

> 💡 **表 8.1 后的一条极实用的说明**：
> *"a **zero-DoF joint can serve as a universal loop joint**, since any closed-loop
> system can be modified, **by adding massless bodies and zero-DoF joints**,
> in such a way that every loop-closing joint is a zero-DoF joint."*
>
> **也就是说：只要实现零自由度闭环关节这一种，就够用了**——
> 其余情形都能通过加无质量刚体转化过去。这是很好的实现简化。

---

## 8.4 闭环关节力

若空间力 $f$ 作用在树中 body $i$ 上，其效果等价于关节空间力 $\tau=J_i^{\mathsf T}f$。
故闭环关节 $k$ 对生成树的效果等价于

$$
\tau=(J_{s(k)}^{\mathsf T}-J_{p(k)}^{\mathsf T})f_k=J_{Ll}^{\mathsf T}f_k
$$

分解成主动力与约束力：

$$
\tau^{a}=\sum_{l=1}^{N_L}J_{Ll}^{\mathsf T}f^a_k\tag{8.28}
\qquad
\tau^{c}=\sum_{l=1}^{N_L}J_{Ll}^{\mathsf T}f^c_k\tag{8.29}
$$

用 $f^c_k=T_k\lambda_k$（8.30）代入并与式 8.17 的 $K$ 对比，即得 $\tau^{c}=K^{\mathsf T}\lambda$（8.32）。

> 🔑 **原书对 $\tau^a$ 的处理建议**：
> *"The **best way** to incorporate $\tau^a$ into the equation of motion is to
> **modify the algorithm for calculating $C$ so that it calculates $C-\tau^a$ instead**."*
> —— 见 §8.6。

---

## 8.5 求解运动方程：三种方法

$$
H\ddot q+C=\tau+K^{\mathsf T}\lambda+\tau^{a}\tag{8.34}
\qquad
K\ddot q=k+k_{stab}\tag{8.35}
$$

$$
\begin{bmatrix}H&K^{\mathsf T}\\ K&0\end{bmatrix}\begin{bmatrix}\ddot q\\ -\lambda\end{bmatrix}
=\begin{bmatrix}\tau-C+\tau^{a}\\ k+k_{stab}\end{bmatrix}\tag{8.36}
$$

**系数矩阵**：$(n+n_c)\times(n+n_c)$，对称、**不正定**、秩 $=n+r$（$r=\mathrm{rank}(K)$），
**$r<n_c$ 时奇异**。

**$r<n_c$ 时的两个问题**（原书指出）：

1. $\lambda$ 的部分元素**不定**（但 $\ddot q$ 仍唯一）；
2. 式 8.35 可能**轻微不相容**——$k+k_{stab}$ 可能略微落在 $\mathrm{range}(K)$ 之外。
   起因是 $K,k,k_{stab}$ 的数值误差，以及位置/速度层闭环约束的不精确满足。
   > **此时按最小二乘意义求解是可以接受的。**

### 方法 1：直接解式 8.36

用通用线性方程求解器，或专为**对称不定**矩阵设计的求解器。

> *"This is the **simplest** solution method, but **not the most efficient**.
> It is therefore appropriate whenever **human effort is more important than
> computational efficiency**."*

**剪枝技巧**：若事先知道哪些行线性相关，可以**先剔除 $n_c-r$ 个相关行**再组装，
相当于把 $n_c$ 减小到 $=r$。**剪枝也可以配合方法 2、3 使用。**

### 方法 2：先解 $\lambda$

用 $KH^{-1}$ 乘第一行减去第二行：

$$
A\lambda=b\tag{8.37}
\qquad
\boxed{A=KH^{-1}K^{\mathsf T}}\tag{8.38}
\qquad
b=k+k_{stab}-KH^{-1}(\tau-C+\tau^{a})\tag{8.39}
$$

$A$：$n_c\times n_c$，**对称半正定**，秩 $=r$。

- $r=n_c$ ⟹ $A$ 可逆，$\lambda=A^{-1}b$
- $r<n_c$ ⟹ $A$ 奇异，解有无穷多个：
  $$\lambda=A^{+}b+(\mathbf 1-A^{+}A)z\tag{8.41}$$
  $A^{+}$ 是伪逆，$z$ 任意，$\mathbf 1-A^{+}A$ 把向量投影到 $A$ 的零空间。
  > **$z$ 对 $\ddot q$ 没有影响**，因为 $K^{\mathsf T}(\mathbf 1-A^{+}A)=0$。
  > 所以可以随便取 $z=0$。

**⭐ 原书给的高效计算流程**（用第 6 章的 $LTL$ 分解，因而**利用了 $H$ 的分支稀疏性**）：

```
1. 分解 H = LᵀL
2. τ' = τ − C + τ^a
3. 回代算 Y = L⁻ᵀKᵀ  和  z = L⁻ᵀτ'
4. A = YᵀY,  b = k + k_stab − Yᵀz
5. 解 Aλ = b
6. 用步骤 1 的因子解 H q̈ = τ' + Kᵀλ
```

> 💡 **注意第 3、4 步的技巧**：$A=KH^{-1}K^{\mathsf T}=(L^{-\mathsf T}K^{\mathsf T})^{\mathsf T}(L^{-\mathsf T}K^{\mathsf T})=Y^{\mathsf T}Y$。
> **从不显式求 $H^{-1}$**，而且 $Y^{\mathsf T}Y$ 天然对称半正定。

### 方法 3：先解式 8.35（独立坐标法）

$$
\ddot q=G\ddot y+g\tag{8.42}
\qquad
KG=0\tag{8.43}
\qquad
Kg=k+k_{stab}\tag{8.44}
$$

$\ddot y$ 是 $(n-r)$ 维**独立加速度变量**向量（**通常取 $\ddot q$ 元素的一个子集**）。

左乘 $G^{\mathsf T}$ 消掉 $\lambda$，代入 8.42：

$$
\boxed{\ G^{\mathsf T}HG\,\ddot y=G^{\mathsf T}(\tau-C+\tau^{a}-Hg)\ }\tag{8.45}
$$

> **这可以看作用独立关节加速度变量表达的闭环系统运动方程**（与第 3 章式 3.20 对照）。
> **$G^{\mathsf T}HG$ 对称正定**，因此可以立即解出 $\ddot y$（式 8.46）。

---

## 8.6–8.9 四个配套算法

| 节 | 算法 | 要点 |
|---|---|---|
| **8.6** | $C-\tau^a$ | 改造 RNEA（表 8.2、8.3） |
| **8.7** | $K$ 与 $k$、$k_{stab}$ | 把所有量**变换到基座坐标系**再组合（最简单，未必最高效） |
| **8.8** | $G$ 与 $g$ | 从 $K$ 出发 |
| **8.9** | 利用 $K$、$G$ 的稀疏性 | — |

### §8.6 的算法（表 8.2）：RNEA 的三处改动

| 改动 | 内容 |
|---|---|
| 1 | **去掉 $S_i\ddot q_i$ 项**（算 $C$ 时加速度变量置零） |
| 2 | 变量名 $a_i\to a^{vp}_i$、$\tau_i\to C'_i$（**纯粹是改名**，$C'=C-\tau^a$） |
| 3 | **关节力计算里新增一项**：在减去一般外力的**同一位置**减去闭环关节的主动力 |

$$
f_i=f^B_i-f^x_i-\sum_{k=N_B+1}^{N_J}e_{ik}f^a_k+\sum_{j\in\mu(i)}f_j
$$

$$
e_{ik}=\begin{cases}+1&i=s(k)\\ -1&i=p(k)\\ 0&\text{否则}\end{cases}
$$

> 🔑 **原书的解释**：*"**From the perspective of the spanning tree, the active
> loop-joint forces are indeed external forces**, since they come from something that
> is **not a part of the tree**, and that is exactly how they are treated."*
>
> **这个视角很重要**：闭环关节的主动力对树来说就是外力，因此处理方式完全一样。

**伪代码结构**（表 8.3）：**三个循环**——

1. 外推：算 $v_i$、$a^{vp}_i$，末尾把 $f_i$ 初始化为 $f^B_i-f^x_i$
2. **新增的中间循环**：对每个回路 $l$，减去 $f^a_k$（**注意要做坐标变换**）
3. 内推：算 $C'_i=S_i^{\mathsf T}f_i$ 并累加到父节点

> 💡 **原书的一个假设**：中间循环里用 `jcalc` 算 $T^a_k$ 时，**假定 $T^a_k$ 是常量**，
> 所以 `jcalc` 只需要关节类型就够了。不做这个假定的话，
> 还得算出关节 $k$ 从前驱到后继坐标系的变换传给 `jcalc`。

---

## 8.10 ⭐ 闭环系统的特殊性质（本章最有价值的一节）

> *"Closed-loop systems exhibit a **greater variety of properties and behaviours**
> than kinematic trees, which can **pose new difficulties** for dynamics algorithms."*

### A. 自由度 (Mobility)

$$
\boxed{\ \text{mobility}=n-r\ }\tag{8.52}
$$

（运动学树的自由度就是 $n$。）独立坐标向量 $y$ 的维数就是 $n-r$。

**恰约束 (properly constrained)**：$r=n_c$ 的系统。此时

$$
\text{mobility}=n_{tot}-6N_L\ \ (\text{3D})\tag{8.53}
\qquad
\text{mobility}=n_{tot}-3N_L\ \ (\text{2D})\tag{8.54}
$$

其中 $n_{tot}$ = **所有**关节（不只是树关节）的自由度之和。

> 💡 这就是机构学里的 **Grübler / Kutzbach 公式**。

### B. 可变自由度 (Variable Mobility)

> **运动学树的自由度是固定的；闭环系统的自由度依赖 $r$，而 $r$ 会随位形变化。**

**原书图 8.2(a) 的例子**：只要角度 $\theta\ne0$，机构自由度为 1；
但 $\theta=0$ 时两臂可以独立运动，自由度变成 2；
而在两种运动模式的**边界位形**上，机构的**瞬时自由度是 3**。

用 $n$ 和 $r$ 表达：全程 $n=5$，但 $r$ 分别是 4、3、2。

> ⚠️ **实践含义**：闭环机构的 $K$ 的秩会**在仿真过程中变化**。
> 写的代码若假定 $r$ 固定，在这类位形附近会出问题。

### C. 位形歧义 (Configuration Ambiguity)

> **$q$ 总能唯一确定树的位形，但 $y$ 不一定能唯一确定闭环系统的位形。**

**原书图 8.2(b)**：知道独立变量 $\theta$ 的值**不足以**唯一确定系统位形（有两个解）。

**两种解决办法**：

1. **用 $q$ 而不是 $y$**（$q$ 是生成树的位置向量）—— 前面各节用的就是这个
2. **提供一个把每个 $y$ 映到唯一 $q$ 的函数** —— **§8.11 的闭环函数**

> 💡 对图 8.2(b)，办法 2 只在**事先知道两个位形之一不可达**时适用；
> 否则要么换独立变量，要么退回办法 1。

### D. ⚠️ 过约束 (Overconstraint)

$$
r<n_c\ \Longrightarrow\ \textbf{过约束}
$$

后果：$n_c-r$ 个**冗余约束** ⟹ 运动方程中有 $n_c-r$ 个**不定的约束力**；
且相对恰约束系统有**多余的自由度**：

$$
\text{mobility}=n_{tot}-6N_L+(n_c-r)\tag{8.55}
$$

> 🔑 **原书的一句警告，实践中极其重要**：
>
> *"Overconstrained systems are **actually very common**. For example,
> **any system containing planar kinematic loops will be overconstrained**."*

**原书给的具体例子**：图 8.2(b) 的平面回路含 4 个转动关节，
$n=3$、$n_c=5$，但自由度是 1 ⟹ $r=2$，所以 $n_c-r=3$ 个冗余约束。

**⭐ 实用的转化技巧**（原书给的解法）：

> 过约束系统的动力学比恰约束系统**更难算**，因此**尽可能把前者转化为后者**。
> 做法：**把原来的闭环关节换成施加更少约束的关节**，从而减小 $n_c$。
>
> **例**：把平面回路中闭合用的**转动关节**换成**球-圆柱关节
> (sphere-in-cylinder joint)**，后者只施加 **2 个**约束，
> 于是 $n_c$ 从 5 降到 2，恰好 $n_c=r$。

> 💡 **这是一个非常实用的建模技巧**：平面四杆机构若按空间机构建模，
> 直接用转动关节闭合就会过约束；换成球-圆柱关节就变成恰约束，
> **运动完全一样，但方程好解得多、约束力也变确定了**。

---

## 8.11 闭环函数 (Loop Closure Functions)

设 $y$ 是独立位置变量向量，定义可接受值集合 $C\subseteq\mathbb R^{n-r}$，
假定 $y\in C$ 时 $y$ 唯一确定 $q$。则存在函数 $\gamma$：

$$
q=\gamma(y)\tag{8.56}
\qquad
\dot q=G\dot y,\ G=\frac{\partial\gamma}{\partial y}\tag{8.57, 8.58}
\qquad
\ddot q=G\ddot y+g,\ g=\dot G\dot y\tag{8.59, 8.60}
$$

**四步方法**：定义 $y$ 与 $C$ → 求 $\gamma$ 的表达式 → **符号求导**得 $G$、$g$ → 写代码。

**与 §8.5 方法的对比**：

| | §8.5 的算法 | §8.11 的闭环函数 |
|---|---|---|
| 通用性 | **更通用** | 较低（第 2 步不总是可行） |
| 便利性 | 只需系统模型 | 用户还要额外提供 $G$、$g$ 的表达式或代码 |
| 计算代价 | 较高 | **可以低得多** |
| **闭环误差** | 会有，**需要稳定化** | ✅ **不会发生** |

> 🔑 **最后一条是决定性的优势**（原书解释）：
> 因为**不积分 $\ddot q$ 得到 $\dot q,q$，而是积分 $\ddot y$ 得到 $\dot y,y$**，
> 然后用式 8.56、8.57 从 $y,\dot y$ 算出 $q,\dot q$，
> **这样算出的值自动满足闭环约束**。
> **因此完全不需要约束稳定化项。**

**原书给了一个完整算例**（图 8.3）：4 刚体、5 关节、1 回路的平面机构，
自由度 2；$B_1$ 与 $B_4$ 组成两连杆臂，$B_2$、$B_3$ 是一个线性驱动器的缸体和活塞。

---

## 8.12 ⭐ 逆动力学：性质发生了根本变化

> *"In a kinematic tree, there is a **1:1 relationship** between $\tau$ and $\ddot q$...
> In a closed-loop system, the given value of $\ddot q$ must comply with the
> loop-closure constraints, and **there are infinitely many values of $\tau$
> that produce the same acceleration**."*

从生成树方程出发：

$$
\tau=H\ddot q+C-\tau^{a}-K^{\mathsf T}\lambda\tag{8.61}
$$

两个未知量 $\tau$ 和 $\lambda$。**左乘 $G^{\mathsf T}$ 消掉 $\lambda$**：

$$
\boxed{\ G^{\mathsf T}\tau=G^{\mathsf T}\tau_{ID}\ }\tag{8.62}
\qquad
\tau_{ID}=H\ddot q+C-\tau^{a}=\mathrm{ID}(q,\dot q,\ddot q)
$$

> $G^{\mathsf T}$ 是 $(n-r)\times n$ 矩阵 ⟹ 式 8.62 对 $n$ 维未知向量只施加 $n-r$ 个约束，
> **剩下 $r$ 个自由度**。
>
> $$\boxed{\ \text{有 }\infty^{r}\text{ 个不同的 }\tau\text{ 产生同一个加速度}\ }$$
>
> **想要唯一解，必须再加约束或引入最优性准则。**

### 驱动器力

设 $p$ = 树中被驱动的自由度数，$u$ = $p$ 维驱动力向量。
把 $G$ 按"被驱动/未驱动"分成 $G_u$、$G_0$，则

$$
\boxed{\ G_u^{\mathsf T}u=G^{\mathsf T}\tau_{ID}\ }\tag{8.63}
$$

**三种情形**：

| 情形 | 条件 | 名称 | 含义 |
|---|---|---|---|
| 1 | $\mathrm{rank}(G_u)<n-r$ | **欠驱动 (underactuated)** | 运动自由度多于驱动器能控制的 |
| 2 | $p>\mathrm{rank}(G_u)$ | **冗余驱动 (redundantly actuated)** | 无穷多组 $u$ 产生同一加速度 |
| 3 | $p=\mathrm{rank}(G_u)=n-r$ | **恰驱动 (properly actuated)** | $G_u$ 可逆，**唯一解** |

> 💡 **这三个术语在机器人控制里天天用**，但它们的精确定义就在这里。
> 人形机器人是情形 1（欠驱动，浮动基的 6 个自由度没有驱动器）；
> 多臂协同抓取是情形 2（冗余驱动，内力可以任意分配）。

---

## 8.13 稀疏矩阵法

即 §8.1 提到的**策略 1**：从一组无约束刚体出发，同时施加所有关节约束。
产生**大而稀疏**的矩阵方程（就是第 3 章式 3.74）。

原书对它的评价：对**闭环系统非常有用**，但用在运动学树上**竞争不过**前两类方法。

---

## 本章要点回顾

1. **策略 2（生成树 + 闭环约束）是典型闭环系统的最佳选择。**
2. 运动方程含**三**项：$\tau+\tau^c+\tau^a$——别忘了闭环关节的**主动力** $\tau^a$。
3. 鞍点矩阵**对称但不正定**，秩 $=n+\mathrm{rank}(K)$；$r<n_c$ 时 $\lambda$ 不定但 $\ddot q$ 唯一。
4. $K_{lj}=\epsilon_{lj}T_k^{\mathsf T}S_j$，$\epsilon$ 的符号由"根 → $p(k)$ / $s(k)$"两条路径决定。
5. **Baumgarte 的目的是稳定，不是精度**；$T_{stab}$ 按"每秒需要多少张快照"来选。
6. **零自由度关节可以作为万能闭环关节**（加无质量刚体即可转化）。
7. **三种求解方法**：直接解鞍点 / 先解 $\lambda$（$A=KH^{-1}K^{\mathsf T}$）/ 独立坐标（$G^{\mathsf T}HG$）。
8. **闭环关节的主动力对生成树而言就是外力。**
9. **自由度 $=n-r$**；**可变自由度**、**位形歧义**、**过约束**是闭环特有的三个麻烦。
10. **含平面回路的系统一定过约束**；换成球-圆柱关节可转成恰约束。
11. **闭环函数法完全避免闭环误差**，因为积分的是 $\ddot y$ 而不是 $\ddot q$。
12. **闭环系统的逆动力学有 $\infty^r$ 个解**；欠驱动/冗余驱动/恰驱动的精确定义在 §8.12。

---

## 易错点

1. **忘了 $\tau^a$**（闭环关节的主动力）。
2. **对鞍点矩阵用 Cholesky**——它不正定。
3. **忘记约束稳定化**：短时仿真看不出，长时必然散架。
4. **把 $T_{stab}$ 调到极小以求精度**——这是错误策略，只会让方程变刚性。
5. **没意识到平面回路一定过约束**，然后困惑于 $\lambda$ 为什么算不出唯一值。
6. **假定 $r$ 在仿真中固定**——可变自由度会打破这个假设。
7. **显式求 $H^{-1}$**——用 §8.5 方法 2 的 $LTL$ 流程。
8. **用 $y$ 做位置变量却没处理位形歧义**。

## 与其他章的联系

- ← 第 3 章：§3.2 方法 1/2；式 3.17 = 式 8.6；式 3.20 = 式 8.45
- ← 第 4 章：生成树、弦、$N_L=N_J-N_B$、$X_P$/$X_S$、$p(k)$/$s(k)$
- ← 第 5 章：§8.6 是 RNEA 的改造版
- ← 第 6 章：§8.5 方法 2 用 $LTL$ 分解，**利用了 $H$ 的分支稀疏性**
- ← 第 7 章：§7.4 变体 1 可处理部分闭环系统；§7.5 的装配法是第三条路
- → 第 11 章：接触约束与闭环约束**数学同构**，差别只在**单边性**

---

## ✍️ 我的理解

<!-- 建议：说清「为什么闭环的逆动力学有无穷多解」 -->

## ❓ 疑问与待办

- [ ] 用曲柄滑块机构走一遍 §8.1 的四步流程
- [ ] 手推平面四杆的 $K$（用式 8.20 与 $\epsilon_{lj}$），验证 $r=2<n_c=5$
- [ ] 实现 §8.5 方法 2 的六步流程，与直接解鞍点对拍
- [ ] 对比不加稳定化 / Baumgarte / 闭环函数三种做法的漂移曲线
- [ ] 把平面回路的闭环转动关节换成球-圆柱关节，观察 $n_c$ 从 5 降到 2
- [ ] 验证 §8.12：给定同一个 $\ddot q$，构造两组不同的 $\tau$

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.141–170 逐节核对。
     此前版本遗漏了 τ^a、ε_lj 与回路的根、δ 的计算与表 8.1、
     §8.6-8.9 的四个算法、§8.10 的四个特殊性质、§8.11 闭环函数、
     §8.12 逆动力学的性质变化，现已补上。 -->
