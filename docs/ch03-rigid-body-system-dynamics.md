# 第 3 章 刚体系统动力学 (Dynamics of Rigid Body Systems)

> **原书 pp. 39–64**，共 7 节。与第 2 章并列为"最数学的两章"。

> **一句话概括**：给出刚体系统运动方程的各种标准形式，
> 并系统讲解**如何从更简单的方程把它构造出来**——
> 这些构造方法就是后面各章算法的骨架。

## 本章在全书中的角色

原书开篇讲得很直白：

> *"The purpose of a dynamics algorithm is to evaluate numerically an equation of motion;
> but, before we can evaluate it, we must first know **what that equation is**."*

**关节的本质**（原书的表述值得逐字读）：

> 关节的作用是给它连接的两个刚体施加**运动约束**：某些方向允许相对运动，某些不允许。
> 更精确地说，关节的作用是往系统里**引入一个约束力**。这个力的特殊之处在于：
> **我们不知道它的值，但知道它的效果**——它会取"任何使结果运动满足约束所需的值"。

> 🔑 **本章的一大半内容就是在讲：怎么描述运动约束，怎么把它施加到运动方程上。**

⚠️ **读法提醒**：本章比较抽象，且原书明说"讲得比后续算法所需的最低限度更深"。
第一遍可以读快，**读完第 5–7 章后回头重读会通透很多**。

## 原书节次对照

| 节  | 原书标题                             | 页 | 内容                                         |
| --- | ------------------------------------ | -- | -------------------------------------------- |
| 3.1 | Equations of Motion                  | 40 | 标准形式、广义速度变量、状态空间形式         |
| 3.2 | Constructing Equations of Motion     | 42 | **两种操作、四种构造方法**、隐式/显式约束    |
| 3.3 | Vector Subspaces                     | 46 | 子空间、直和、**正交补（对偶 vs 欧氏）**     |
| 3.4 | Classification of Constraints        | 50 | 等式/不等式、完整/非完整、定常/非定常        |
| 3.5 | Joint Constraints                    | 53 | 关节模型：$S$、$T$、$T_a$、偏置速度 $\sigma$ |
| 3.6 | Dynamics of a Constrained Rigid Body | 57 | **同一问题的三种解法**                       |
| 3.7 | Dynamics of a Multibody System       | 60 | 六步组装出整体方程                           |

所有公式已数值验证：`python3 code/verify_ch03.py`

---

## 3.1 运动方程 (Equations of Motion)

### 标准形式

$$
\boxed{\ H(q)\,\ddot q+C(q,\dot q)=\tau\ }\tag{3.1}
$$

| 符号               | 原书名称                                      | 说明                |
| ------------------ | --------------------------------------------- | ------------------- |
| $q,\dot q,\ddot q$ | 广义位置/速度/加速度变量                      |                     |
| $\tau$             | 广义力                                        |                     |
| $H(q)$             | **广义惯性矩阵 (generalized inertia matrix)** | 只依赖$q$           |
| $C(q,\dot q)$      | **广义偏置力 (generalized bias force)**       | 依赖$q$ 和 $\dot q$ |

**原书对 $C$ 的定义（一句话，非常好记）**：

> *"The bias force is simply the value of $\tau$ that will produce **zero acceleration**."*

它涵盖：科氏力、离心力、**重力**，以及 $\tau$ 之外作用于系统的**任何其他力**。

**术语**：$H$ 与 $C$ 是**系数**，$\tau$ 与 $\ddot q$ 是**变量**。

$$
T=\tfrac12\,\dot q^{\mathsf T}H\dot q\tag{3.2}
$$

### FD / ID

$$
\ddot q=\mathrm{FD}(model,q,\dot q,\tau)\tag{3.3}
$$

$$
\tau=\mathrm{ID}(model,q,\dot q,\ddot q)\tag{3.4}
$$

显然 $\mathrm{FD}=H^{-1}(\tau-C)$、$\mathrm{ID}=H\ddot q+C$。**但原书特别强调**：

> *"the algorithms that implement them **need not necessarily work by calculating
> $C$ or $H$ or $H^{-1}$**."*

> 🔑 **这句话是第 7 章 ABA 的伏笔**：ABA 计算 FD，却**从不构造 $H$**。
> 式 3.3/3.4 的价值在于"清楚显示输入输出，而不暗示任何特定计算方法"。

### 当位置变量不是速度变量的积分时

式 3.1 隐含假设"位置变量是速度变量的积分"。**这并不总成立**，原书给了三种情形：

1. 系统含**非完整约束**（§3.4）⟹ 位置变量比速度变量**多**；
2. 某些运动用**冗余**位置变量描述更好（例如用四元数表示姿态）⟹ 同样更多；
3. 即使个数相同，有时也希望速度变量不是位置变量的导数。

此时把 $\dot q$ 换成另一个速度变量向量 $\alpha$：

$$
H(q)\,\dot\alpha+C(q,\alpha)=\tau\tag{3.5}
$$

并补一个从 $\alpha$ 算 $\dot q$ 的方程（左：矩阵形式；右：函数形式）：

$$
\dot q=Q(q)\,\alpha
\qquad\qquad
\dot q=\mathrm{qdfn}(model,q,\alpha)\tag{3.6, 3.7}
$$

> 💡 式 3.6 揭示**数学结构**（$\dot q$ 线性依赖于 $\alpha$）；
> 式 3.7 只说**输入输出**，不暗示算法。原书反复使用这种"一对写法"。

### 状态空间形式（仿真器要的形式）

$$
x=\begin{bmatrix}q\\ \alpha\end{bmatrix},
\qquad
\dot x=\begin{bmatrix}\dot q\\ \dot\alpha\end{bmatrix}
=\begin{bmatrix}Q\alpha\\ H^{-1}(\tau-C)\end{bmatrix}\tag{3.8}
$$

$$
\dot x=\mathrm{FD}_x(model,x,\tau)
=\begin{bmatrix}\mathrm{qdfn}(model,q,\alpha)\\ \mathrm{FD}(model,q,\alpha,\tau)\end{bmatrix}\tag{3.9}
$$

> 原书说明：**后续大多用式 3.1 这种较简单的系统**，
> 因为大部分结论都能直截了当地推广到式 3.5；
> *"Most dynamics algorithms are indifferent to the distinction between
> Eqs. 3.1 and 3.5, except at the level of programming details."*

**最后一种可能**（原书结尾提到）：有时不便于找出一组显式的独立速度变量。
此时仍用式 3.1/3.5，但 $\dot q,\ddot q$ 要**额外满足一组运动约束**，
方程要配一个约束方程——这就引出了 §3.2。

---

## 3.2 构造运动方程 (Constructing Equations of Motion)

### 两种基本操作

系统的运动方程是一串数学操作的**最终结果**。起点是单个刚体/子系统的运动方程，
操作只有两种：

1. **把方程收集起来**，构成更大子系统的方程；
2. **施加额外的运动约束**。

> 🔑 **"按什么顺序做这两件事" = "什么算法"。** 这是本章最重要的组织思想。

### ⭐ 四种构造方法（对应四类算法）

| 方法  | 做法                                                            | 对应的算法                          | 特点                                                       |
| ----- | --------------------------------------------------------------- | ----------------------------------- | ---------------------------------------------------------- |
| **1** | 收集**所有**刚体，然后施加**所有**约束                          | §3.7 的例子；朴素算法              | 易懂，但产生**大而稀疏**的矩阵，**必须**用稀疏技术否则太慢 |
| **2** | 先得到**生成树**的运动方程，再施加**剩余**约束                  | **闭环系统的标准做法（第 8 章）**   | 生成树有极高效的算法可用                                   |
| **3** | 从单个刚体/子系统开始：加一个刚体 → 施加与之相关的约束 → 重复 | **≈ 铰接体算法 ABA（第 7 章）**    | ABA 用的是铰接体方程而非完整运动方程                       |
| **4** | 子系统**两两合并**，各自内部施加约束，重复                      | **分治算法** (Featherstone 1999a,b) | 适合**并行计算**                                           |

> 💡 **这张表是理解全书结构的钥匙**：第 5–8 章不是四个孤立的技巧，
> 而是"两种操作的四种排列顺序"。
> 特别是**方法 3 = ABA**：这解释了为什么 ABA 要一个刚体一个刚体地往上"装配"。

### 收集方程（式 3.10）

若子系统 $i$ 的方程是 $H_i\ddot q_i+C_i=\tau_i$，则合并后

$$
\begin{bmatrix}H_1&&&\\&H_2&&\\&&\ddots&\\&&&H_N\end{bmatrix}
\begin{bmatrix}\ddot q_1\\ \ddot q_2\\ \vdots\\ \ddot q_N\end{bmatrix}
+\begin{bmatrix}C_1\\ C_2\\ \vdots\\ C_N\end{bmatrix}
=\begin{bmatrix}\tau_1\\ \tau_2\\ \vdots\\ \tau_N\end{bmatrix}\tag{3.10}
$$

**关键洞察（原书特别解释）**：单个刚体也是一个"子系统"，
可以直接用它的空间运动方程 $I_ia_i+p_i=f_i$ 顶替广义坐标形式（用 $I_i$ 代替 $H_i$）。

**为什么可以这样混用**：

> *"$M^6$ and $F^6$ are **special cases** of $M^n$ and $F^n$, and Plücker coordinates
> on $M^6$ and $F^6$ are a special case of generalized coordinates on $M^n$ and $F^n$."*

> 🔑 **空间向量是广义向量的特例**——所以空间速度可以和广义速度自由混用。
> 这一条许可证在第 7 章 ABA 和第 8 章会被反复使用。

### 运动约束的两种描述（式 3.11）

|                     | 位置层        | 速度层           | 加速度层             |
| ------------------- | ------------- | ---------------- | -------------------- |
| **隐式 (implicit)** | $\phi(q)=0$   | $K\dot q=0$      | $K\ddot q=k$         |
| **显式 (explicit)** | $q=\gamma(y)$ | $\dot q=G\dot y$ | $\ddot q=G\ddot y+g$ |

$$
K=\frac{\partial\phi}{\partial q},\quad k=-\dot K\dot q,\quad
G=\frac{\partial\gamma}{\partial y},\quad g=\dot G\dot y
$$

若二者描述同一约束，则

$$
\phi\circ\gamma=0,\qquad KG=0,\qquad Kg=k\tag{3.12}
$$

> 💡 **直观理解**：隐式说"**不许往哪走**"（$K$ 的行 = 被禁止的方向）；
> 显式说"**只许往哪走**"（$G$ 的列 = 被允许的方向）。
> $KG=0$ 就是"允许的方向必须避开被禁止的方向"。

### 约束力与 Jourdain 虚功率原理

带约束的方程：

$$
H\ddot q+C=\tau+\tau_c\tag{3.14}
$$

$\tau_c$ 是**约束力**，未知，但满足一条关键性质——原书原文：

> **The constraint force delivers zero power along every direction of velocity freedom
> that is compatible with the motion constraints.**
> (**Jourdain's principle of virtual power**)

⚠️ **注意是 Jourdain 虚功率原理，不是达朗贝尔虚功原理**。
区别在于：虚功原理谈**虚位移**，Jourdain 原理谈**虚速度/功率**。
处理非完整约束时，Jourdain 形式更合适（非完整约束是速度层的）。

由此立刻得到两个形式（左：隐式，式 3.15；右：显式，式 3.16）：

$$
\tau_c=K^{\mathsf T}\lambda
\qquad\qquad
G^{\mathsf T}\tau_c=0\tag{3.15, 3.16}
$$

**验证（原书给的，一行）**：

- 若 $\tau_c=K^{\mathsf T}\lambda$，则 $\tau_c\cdot\dot q=\lambda^{\mathsf T}K\dot q$，
  对任何满足 $K\dot q=0$ 的 $\dot q$ 都为零 ✓
- 若 $\dot q=G\dot y$，则 $\dot q\cdot\tau_c=\dot y^{\mathsf T}G^{\mathsf T}\tau_c$，
  对所有 $\dot y$ 为零当且仅当 $G^{\mathsf T}\tau_c=0$ ✓

$\lambda$ 的元素就是 **Lagrange 乘子**。

### 施加隐式约束 → KKT 系统

把式 3.15 代入 3.14，再与加速度约束联立：

$$
\boxed{\ \begin{bmatrix}H&K^{\mathsf T}\\ K&0\end{bmatrix}
\begin{bmatrix}\ddot q\\ -\lambda\end{bmatrix}
=\begin{bmatrix}\tau-C\\ k\end{bmatrix}\ }\tag{3.17}
$$

⚠️ **注意原书的符号约定**：矩阵里是 $+K^{\mathsf T}$，未知量里是 $-\lambda$。
（不同文献放法不同，自己内部一致即可。）

**关于秩的重要讨论**（原书这段值得记住）：

- $n=\dim(\ddot q)$：无约束系统的自由度；$n_c=\dim(\lambda)=\dim(\phi)$：施加的约束数
- 系数矩阵的秩 $=n+n_{ic}$，其中 $n_{ic}=\mathrm{rank}(K)$ 是**独立**约束数
- $n_{ic}=n_c$ ⟹ 系数矩阵非奇异，$\ddot q$ 和 $\lambda$ 都能解出
- $n_{ic}<n_c$ ⟹ **仍能解出 $\ddot q$**，但 $\lambda$ 有 $n_c-n_{ic}$ 维不定性

> 🔑 **$n_{ic}<n_c$ 的系统称为"过约束 (overconstrained)"。**
> 原书说得很实在：*"Ideally, we would always have $n_{ic}=n_c$,
> but it is not always possible to guarantee this in practice."*
>
> **物理含义**：过约束系统的运动是确定的，但**内力不确定**——
> 就像一张四条腿的桌子放在平地上，四个支反力的分配是超静定的。
> 做仿真无所谓，做**结构受力分析**就必须处理。

### 施加显式约束 → 投影法

把式 3.14、3.16 和显式加速度约束联立（式 3.18），做两步高斯消元，得到

$$
G^{\mathsf T}HG\,\ddot y=G^{\mathsf T}(\tau-C-Hg)\tag{3.20}
$$

令 $u=G^{\mathsf T}\tau$，则

$$
\boxed{\ H_G\,\ddot y+C_G=u\ }\tag{3.21}
$$

$$
H_G=G^{\mathsf T}HG,\qquad C_G=G^{\mathsf T}(C+Hg)
$$

> 🔑 **结果与式 3.13 具有完全相同的代数形式**——
> 也就是说，**施加显式约束之后，你得到的还是一个"标准"运动方程**，
> 只是维数降低了、变量换成了 $y$。这就是"化归"。

**原书给的四点评注**（都很有价值）：

1. 式 3.19 可以看作式 3.17 的**对偶**（运动变量与力变量互换位置）。
   但**实践中 3.17 比 3.19 有用得多**。
2. 式 3.20 有一条**捷径**：直接把式 3.14 **左乘 $G^{\mathsf T}$**（这就消掉了 $\tau_c$），
   再代入显式加速度约束。这叫**投影法 (projection method)**——
   因为第一步相当于把运动方程投影到 $\mathrm{range}(G^{\mathsf T})$ 上。
3. $H_G$ **继承** $H$ 的对称性和正定性，且 $T=\tfrac12\dot y^{\mathsf T}H_G\dot y$（式 3.22）。
4. $u$ 确实是广义力，因为 $u\cdot\dot y=(G^{\mathsf T}\tau)^{\mathsf T}\dot y=\tau^{\mathsf T}G\dot y=\tau\cdot\dot q$（式 3.23）——**功率不变**。

---

## 3.3 向量子空间 (Vector Subspaces)

### 为什么用子空间语言

原书的理由非常清晰：

> 若约束用隐式描述，则运动被限制在 $\mathcal S=\mathrm{null}(K)$；
> 用显式描述则 $\mathcal S=\mathrm{range}(G)$。
> 若 $K_1,K_2$ 描述同一约束，它们**唯一的共同点**是 $\mathrm{null}(K_1)=\mathrm{null}(K_2)=\mathcal S$。
> **因此捕捉约束本质的是子空间 $\mathcal S$，而不是任何一个特定矩阵。**

> 🔑 **这是本节的中心思想**：$K$、$G$、$S$、$T$ 都不唯一，
> **它们张成/零化的子空间才是物理**。写算法时不要依赖某个特定矩阵的具体取值。

### 基础

- **子空间**：对加法和数乘封闭的子集
- $\mathrm{span}(Y)$：$Y$ 中元素的全部线性组合，必是子空间
- $S\subseteq V$ 用**子集符号**表示"是子空间"（原书说明：没有专门符号，
  要靠上下文区分"子集"和"子空间"）；$S\subset V$ 表示**真**子空间（$\dim S<\dim V$）

**一个有用的计数**（原书给出）：

$$
\boxed{\ n \text{ 维空间中，}m\text{ 维子空间由 }m(n-m)\text{ 个参数唯一确定}\ }
$$

（$m=0$ 或 $m=n$ 时为 0，即唯一；其余情形有无穷多个。）

### 矩阵表示

$\mathcal S$ 的一组基 $\{s_1,\dots,s_m\}$ 拼成 $n\times m$ 矩阵
$S=[s_1\cdots s_m]$，则 $\mathrm{range}(S)=\mathcal S$，任意 $v\in\mathcal S$ 写作 $v=S\alpha$。

**关键区分（原书讲得很好）**：

> $S$ **同时定义了一个子空间和一组基**。在 $S$ 的 $mn$ 个数中，
> $m(n-m)$ 个定义**子空间**，$m^2$ 个定义**基**。

$$
\boxed{\ S'=SA\ (A\text{ 可逆})\ \Longrightarrow\ S'\text{ 与 }S\text{ 描述同一子空间}\ }
$$

> 💡 **实用推论**：任何只依赖子空间的公式，在 $S\to SA$ 下必须**不变**。
> 这是检验公式正确性的好办法——例 3.1 的式 3.28 就用了这个检验。

### 向量的分解与直和

把 $v$ 分解为 $v=v_1+v_2$（$v_1\in\mathcal S_1$，$v_2\in\mathcal S_2$）：

- **可能** ⟺ $v\in\mathrm{span}(\mathcal S_1\cup\mathcal S_2)$
- **唯一** ⟺ $\mathcal S_1\cap\mathcal S_2=\{0\}$

对**任意** $v$ 都既可能又唯一的条件：

$$
\mathcal S_1\cap\mathcal S_2=\{0\}\quad\text{且}\quad
\dim\mathcal S_1+\dim\mathcal S_2=\dim V\tag{3.24}
$$

此时称 $V$ 是二者的**直和 (direct sum)**：

$$
V=\mathcal S_1\oplus\mathcal S_2\tag{3.25}
$$

**求解**：$[S_1\ S_2]$ 此时非奇异，于是

$$
\begin{bmatrix}\alpha_1\\ \alpha_2\end{bmatrix}=[S_1\ S_2]^{-1}v\tag{3.27}
$$

### ⚠️ 正交补：对偶形式 vs 欧氏形式

**这是原书特别警告的一个坑，值得完整引用**：

> *"Orthogonal complements have gained a degree of **notoriety** in the robotics literature,
> owing to the **mistaken practice of treating 6D vectors as if they were Euclidean**,
> and using the Euclidean form of orthogonal complement instead of the dual form."*
> （原书引 Duffy (1990) 讨论此事。）

**两种正交性的区别**：

|                             | 欧氏正交                         | 对偶正交                                 |
| --------------------------- | -------------------------------- | ---------------------------------------- |
| 基于                        | 欧氏内积                         | $U$ 与 $V=U^{*}$ 之间的标量积            |
| $Y_1\perp Y_2$ 何时可能成立 | $Y_1,Y_2$ 在**同一个**向量空间   | $Y_1,Y_2$ 在**不同**的（互为对偶的）空间 |
| 适用于空间向量吗            | ❌**不适用**（$M^6$ 不是欧氏的） | ✅                                       |

> 🔑 **具体到本书**：约束力子空间 $\mathcal S^{\perp}\subseteq F^6$ 是 $\mathcal S\subseteq M^6$
> 的**对偶**正交补，二者**不在同一个空间**。
> 若错误地把 $M^6$ 当欧氏空间、在 $M^6$ 内部求"正交补"，得到的东西**没有物理意义**，
> 而且**不是坐标变换下的不变量**。
>
> 有些作者为了和这个错误划清界限，用"natural orthogonal complement"或
> "reciprocal complement"来称呼对偶形式。

**在动力学中的作用**（一句话）：

$$
\boxed{\ \text{运动被限制在 }\mathcal S\subseteq M^n
\ \Longrightarrow\ \text{约束力属于 }\mathcal S^{\perp}\subseteq F^n\ }
$$

### 例 3.1：用惯性做力的分解（重要，是 ABA 的种子）

给定运动子空间 $\mathcal S\subseteq M^6$ 和正定惯性 $I:M^6\to F^6$，定义

$$
\mathcal T_a=I\mathcal S,\qquad \mathcal T_c=\mathcal S^{\perp}
$$

**结论**：$\mathcal T_a\oplus\mathcal T_c=F^6$。

**证明（原书的，很漂亮）**：维数显然满足；
反证不交：若 $0\ne t\in\mathcal T_a\cap\mathcal T_c$，
则 $t=Is$（某个 $s\ne0\in\mathcal S$），又 $s^{\mathsf T}t=0$，
于是 $s^{\mathsf T}Is=0$——**与 $I$ 正定矛盾** ∎

于是任意 $f\in F^6$ 可唯一分解为 $f=f_a+f_c$：

$$
\boxed{\ f_a=IS(S^{\mathsf T}IS)^{-1}S^{\mathsf T}f\ },\qquad f_c=f-f_a\tag{3.28}
$$

**物理解释（原书原话）**：若 $I$ 是一个被约束在 $\mathcal S$ 中运动的刚体的惯性，
$f$ 是施加的力，则

- $f_a$ = **产生加速度**的那一部分
- $f_c$ = **被约束力抵消**的那一部分（约束力就是 $-f_c$）

若刚体静止且只受 $f$ 与约束力，则

$$
a=I^{-1}f_a=\underbrace{S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}}_{\text{表观逆惯性}}f
$$

> 🔑 **认出这个表达式**：$S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}$ 就是 §3.6 的
> **表观逆惯性 $\Phi$**（式 3.54），而它的"对偶版本"
> $I-IS(S^{\mathsf T}IS)^{-1}S^{\mathsf T}I$ 正是**第 7 章 ABA 的 Schur 补
> $I^A-UD^{-1}U^{\mathsf T}$**。
>
> **本章这个不起眼的例子，是全书最重要算法的数学核心。**

原书还提示：式 3.28 **必须与 $S$ 的具体选取无关**，
把 $S$ 换成 $SA$ 代入即可看到 $A$ 全部约掉。（`verify_ch03.py` 验证了这点。）

---

## 3.4 约束的分类 (Classification of Constraints)

原书图 3.1 的分类树：

```text
                    运动学约束 kinematic constraints
                              │
              ┌───────────────┴───────────────┐
        等式约束 equality              不等式约束 inequality
          phi(...) = 0                   phi(...) >= 0
              │                        <- 第 11 章（碰撞、弹跳、失去接触）
      ┌───────┴────────┐
  完整约束 holonomic   非完整约束 nonholonomic
   （位置变量的约束）      phi(q_1..q_n, qd_1..qd_n, t) = 0
   典型来自「滑动接触」     典型来自「滚动接触」
      │
  ┌───┴────┐
定常 scleronomic      非定常 rheonomic
 phi(q_1..q_n) = 0    phi(q_1..q_n, t) = 0
```

### 等式 vs 不等式

- **等式**：源自两刚体间的**永久物理接触**
- **不等式**：物体既可接触也可分离 ⟹ 描述**碰撞、弹跳、失去接触**（**第 11 章**）

### 完整 vs 非完整

**非完整约束**关于速度变量是**线性**的：

$$
\phi_{nh}(q,\dot q,t)=\phi^0_{nh}+\sum_{i=1}^{n}\phi^i_{nh}\,\dot q_i\tag{3.29}
$$

对比完整约束函数的时间导数：

$$
\frac{d}{dt}\phi_h(q,t)=\frac{\partial\phi_h}{\partial t}+\sum_{i=1}^{n}\frac{\partial\phi_h}{\partial q_i}\dot q_i\tag{3.30}
$$

> 🔑 **两式代数形式完全相同！原书由此给出一个极其重要的结论**：
>
> $$
> \boxed{\ \text{在速度层和加速度层，完整约束与非完整约束}\textbf{没有区别}\ }
> $$
>
> **这就是为什么本书的算法基本不用操心这个分类**——
> 算法工作在加速度层。区别只在于：**式 3.29 不可积**
> （$\phi_{nh}$ 不是任何函数的导数，否则它就是一个用导数表达的完整约束）。

**后果**：含非完整约束的系统，**位置自由度多于速度自由度**，
因而需要比速度变量更多的位置变量 ⟹ 必须用式 3.5 而非式 3.1。

### 定常 vs 非定常

- **定常 (scleronomic)**：只依赖位置变量。来自两刚体固连表面之间的**滑动接触**。
  **最简单也最常见。**
- **非定常 (rheonomic)**：还依赖时间。可看作"某些滑动自由度被强制按时间的
  规定函数变化"的定常约束。

> 💡 **非定常约束的作用是引入"运动学激励 (kinematic excitation)"，即规定运动。**
> 例：把一个转动关节的角度设成 $\theta(t)=\sin t$。

### 例子（原书图 3.2）

| 图  | 系统                                               | 类型         | 自由度                       |
| --- | -------------------------------------------------- | ------------ | ---------------------------- |
| (a) | 圆柱关节 (cylindrical joint)                       | 定常         | 2（滑动 + 转动）             |
| (b) | 移动关节 (prismatic joint)                         | 定常         | 1                            |
| (c) | 三个圆盘被一根**无质量、零直径、不可伸长的绳**围住 | **多体约束** | 满足$l_1+l_2+l_3+2\pi r=p$   |
| (d) | 球在平面上**纯滚动**                               | **非完整**   | 瞬时 3 自由度，位置 5 自由度 |

**图 3.2(c) 的价值**：原书用它说明"运动学约束**几乎总是**两个刚体之间的关系，
但**可以**构造涉及两个以上刚体、且**不能**化归为一组两体约束的约束"。

**图 3.2(d) 的价值**：这是非完整约束的**典型特征**——
球瞬时只有 3 个运动自由度（两个方向滚动 + 绕接触点自旋），
但通过适当的操作序列，**可以到达平面上任意点、任意姿态**，
所以有 **5 个位置自由度**。若球可以滑动，那就变成定常约束了。

**关于建模选择**（原书最后一段，很实用）：
图 3.2(a)(b) 的形状决定了物体**拉不开**；
但 (c) 的绳可能**松弛**，(d) 的球可能**离开**平面。
**如果你对将要作用的力了解得足够多、确信这些事不会发生**，
就可以按等式约束建模；否则必须按不等式约束建模（第 11 章）。

---

## 3.5 关节约束 (Joint Constraints)

### 定义与术语

> **原书对关节的定义：关节 = 两个刚体之间的任意运动学约束。**
> 因此关节可以施加 0 到 6 个约束。

两个刚体称为**前驱 (predecessor)** 和**后继 (successor)**，
关节"从前驱连到后继"。**关节速度**定义为后继相对前驱的速度：

$$
v_J=v_s-v_p\tag{3.31}
$$

### 关节约束的一般形式

$$
\boxed{\ v_J=S(q,t)\,\dot q+\sigma(q,t)\ }\tag{3.32}
$$

| 符号     | 名称                                                       |
| -------- | ---------------------------------------------------------- |
| $S$      | **运动子空间矩阵 (motion subspace matrix)**，$6\times n_f$ |
| $\sigma$ | **偏置速度 (bias velocity)**                               |

**$\sigma$ 是什么**：$\dot q=0$ 时 $v_J$ 的取值。
**只有当关节显含时间时才非零**（即非定常约束、或含时的非完整约束）。

不含时则简化为

$$
v_J=S(q)\,\dot q\tag{3.33}
$$

**维数关系**：

$$
\dim(\mathcal S)=n_f,\qquad n_c=6-n_f,\qquad \text{约束力}\in\mathcal S^{\perp}\subseteq F^6
$$

> 💡 **原书脚注里的一个实用细节**：若 $n_f=0$（固定关节），$S$ 是 $6\times0$ 矩阵。
> 这类矩阵遵循矩阵代数的常规，外加一条：$m\times0$ 与 $0\times n$ 之积是 $m\times n$ 零矩阵。
> **所以固定关节不需要特判**，公式自动退化。

### ⭐ 三个子空间：$\mathcal S$、$\mathcal T$、$\mathcal T_a$

关节传递的力 $f_J$（从前驱传给后继）可分解为**主动力**与**约束力**：

$$
f_J=T_a\,\tau+T\,\lambda\tag{3.34}
$$

| 矩阵  | 张成                                         | 含义                                         |
| ----- | -------------------------------------------- | -------------------------------------------- |
| $S$   | $\mathcal S\subseteq M^6$                    | 允许的**运动**方向                           |
| $T$   | $\mathcal T=\mathcal S^{\perp}\subseteq F^6$ | **约束力**方向（来自轴承等）                 |
| $T_a$ | $\mathcal T_a\subseteq F^6$                  | **主动力**方向（来自驱动器、弹簧、阻尼器等） |

满足

$$
T_a^{\mathsf T}S=\mathbf 1\tag{3.35}
$$

$$
T^{\mathsf T}S=0\tag{3.36}
$$

且 $\mathcal T\oplus\mathcal T_a=F^6$（$\mathcal T_a$ 可以是**任何**满足这个直和条件的子空间）。

> ⚠️ **我此前的笔记漏了 $T_a$。** 只讲 $S$ 和 $T$ 是不完整的：
> $T=\mathcal S^{\perp}$ 由 $S$ **唯一确定**，但 $T_a$ **不唯一**，
> 需要额外选择——例 3.2 说明了这个选择的物理意义。

由此立即得到两条**核心公式**：

$$
\boxed{\ S^{\mathsf T}f_J=\tau\ }\tag{3.37}
$$

$$
T^{\mathsf T}v_J=T^{\mathsf T}\sigma\tag{3.38}
$$

- 式 3.37：**从关节传递的空间力得到广义力**——
  **这就是第 5 章 RNEA 内推那一行 `τ_i = S_iᵀ f_i` 的出处**。
  注意约束力 $T\lambda$ 被 $S^{\mathsf T}$ 自动消掉（因 $T^{\mathsf T}S=0$）。
- 式 3.38：速度约束的隐式形式。

**式 3.35 的来历（功率平衡，原书给的）**：
关节交付的总功率是 $f_J\cdot v_J$，来自两个独立功率源：$\sigma$ 和 $\tau$。
归因于 $\tau$ 的那部分是 $f_J\cdot S\dot q$，又必须等于 $\tau\cdot\dot q$，于是

$$
\tau\cdot\dot q=f_J\cdot S\dot q=\tau^{\mathsf T}T_a^{\mathsf T}S\dot q\tag{3.39}
$$

对**所有** $\tau$ 和 $\dot q$ 成立 ⟹ $T_a^{\mathsf T}S=\mathbf 1$ ∎

### 加速度约束

对式 3.32 求导（注意用 §2.10 的运动坐标系求导法则）：

$$
a_J=S\ddot q+\mathring S\dot q+\mathring\sigma+v_s\times(S\dot q+\sigma)
=\boxed{S\ddot q+c_J+v_s\times v_J}\tag{3.40}
$$

$$
c_J=\mathring S\dot q+\mathring\sigma\quad(\text{实际上 } c_J=\mathring v_J)\tag{3.41}
$$

$\mathring S$、$\mathring\sigma$ 是在**随后继刚体运动的坐标系**中的表观导数：

$$
\mathring S=\frac{\partial S}{\partial t}+\sum_{i=1}^{n_p}\frac{\partial S}{\partial q_i}\dot q_i
\tag{3.42}
$$

$$
\mathring\sigma=\frac{\partial\sigma}{\partial t}+\sum_{i=1}^{n_p}\frac{\partial\sigma}{\partial q_i}\dot q_i\tag{3.43}
$$

> 🔑 **原书紧接着的这段提醒极其重要，能省掉大量无谓的复杂度**：
>
> *"readers should bear in mind that the relatively complicated equations above apply to
> a general case that is **almost never needed in practice**. This is partly because
> **most common joint types have the property that $\mathring S=0$ and $\sigma=0$
> (implying $c_J=0$)**, and partly because on the rare occasions when one encounters
> an explicit time dependency in a joint, it is almost always possible to implement
> the time-dependent term via **the hybrid dynamics algorithms of Chapter 9**
> rather than using $\sigma$."*
>
> **翻译成实践建议**：
>
> - 常见关节（转动、移动、螺旋）：$\mathring S=0$、$\sigma=0$ ⟹ $c_J=0$，直接省掉
> - 需要规定运动时：**用第 9 章的混合动力学，别用 $\sigma$**

隐式加速度约束（左乘 $T^{\mathsf T}$ 即得，比求导容易）：

$$
T^{\mathsf T}a_J=T^{\mathsf T}(c_J+v_s\times v_J)\tag{3.44}
$$

### 例 3.2：主动力子空间的物理意义（齿轮驱动的转动关节）

**设置**（原书图 3.3）：刚体 $B$ 经转动关节连到固定基座，
电机 $M$ 通过**小齿轮 + 齿轮**驱动该关节。$z$ 轴沿转轴。

$$
S=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix},
\qquad
T=\begin{bmatrix}
1&0&0&0&0\\ 0&1&0&0&0\\ 0&0&0&0&0\\ 0&0&1&0&0\\ 0&0&0&1&0\\ 0&0&0&0&1
\end{bmatrix}
$$

主动力在两齿轮接触处传递给 $B$，是 $y$ 方向、大小 $f$ 的线力。
接触点坐标 $(-r,0,0)$（$r$ = 齿轮半径），故主动空间力的 Plücker 坐标是

$$
\mathbf f=[0\ \ 0\ \ -rf\ \ 0\ \ f\ \ 0]^{\mathsf T}
$$

由 $T_a^{\mathsf T}S=1$ 解得

$$
T_a=\mathbf f/(\mathbf f^{\mathsf T}S)=[0\ \ 0\ \ 1\ \ 0\ \ -1/r\ \ 0]^{\mathsf T}
$$

> 🔑 **原书的评注是本例的全部价值**：
>
> $T_a$ 和 $\mathrm{range}(T_a)$ **都依赖 $r$**。这个依赖
> **对运动方程没有任何影响**，但**影响约束力的值**。
> 物理上，约束力来自关节的**轴承**。所以：
>
>
> | 任务                     | $r$ 是否相关 |
> | ------------------------ | ------------ |
> | **仿真系统运动**         | ❌ 无关      |
> | **为关节选择合适的轴承** | ✅ 必须考虑  |
>
> **这解释了 $T_a$ 为什么"不唯一却有意义"**：不同的 $T_a$ 对应
> 同一个运动、但不同的内力分配。做仿真可以随便选，做机械设计不能。

---

## 3.6 受约束刚体的动力学 (Dynamics of a Constrained Rigid Body)

**问题设置**（原书图 3.4）：刚体经一个运动子空间为 $\mathcal S$ 的关节连到固定基座。
受三个力：施加力 $f$、约束力 $f_c$、重力 $f_g$。

把不感兴趣的项收进偏置力 $p=v\times^{*}Iv-f_g$（**这正是 §2.14 式 2.69 给的自由度**）：

$$
f+f_c=Ia+p\tag{3.45}
$$

约束条件：

$$
v\in\mathcal S,\qquad f_c\in\mathcal S^{\perp}\tag{3.46}
$$

**目标**：求"加速度作为施加力的函数"。$f_c$ 未知，必须在求解过程中**算出或消掉**。

> 🔑 **原书用三种方法解同一个问题。这三种方法后面会长成三类不同的算法，
> 所以值得都掌握。**

### 法 1：用 $S$ 消掉 $f_c$ → 表观逆惯性

引入 $6\times n_f$ 矩阵 $S$ 张成 $\mathcal S$：

$$
v=S\dot q\tag{3.47}
$$

$$
 S^{\mathsf T}f_c=0\tag{3.48}
$$

$$
a=S\ddot q+\dot S\dot q\tag{3.49}
$$

式 3.48 意味着**左乘 $S^{\mathsf T}$ 就能消掉 $f_c$**：

$$
S^{\mathsf T}f=S^{\mathsf T}(Ia+p)\tag{3.50}
$$

代入 3.49 解出

$$
\ddot q=(S^{\mathsf T}IS)^{-1}S^{\mathsf T}(f-I\dot S\dot q-p)\tag{3.51}
$$

> 💡 $S^{\mathsf T}IS$ 是 $n_f\times n_f$ **对称正定**矩阵，**保证可逆**。
> （这条保证在第 7 章会变成"ABA 的 $D_i>0$，无需选主元"。）

再代回 3.49：

$$
a=\Phi f+b\tag{3.53}
$$

$$
\boxed{\ \Phi=S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}\ }\tag{3.54}
$$

$$
b=\dot S\dot q-\Phi(I\dot S\dot q+p)\tag{3.55}
$$

**$\Phi$ = 受约束刚体的表观逆惯性 (apparent inverse inertia)**，
$b$ = 偏置加速度（$f=0$ 时的加速度）。

**$\Phi$ 的性质**（原书列出，`verify_ch03.py` 全部验证）：

$$
\Phi=\Phi^{\mathsf T},\quad \Phi\succeq0,\quad
\mathrm{range}(\Phi)=\mathcal S,\quad \mathrm{null}(\Phi)=\mathcal S^{\perp},\quad
\mathrm{rank}(\Phi)=n_f
$$

> 🔑 **回看 §2.15**：那里说"受约束刚体没有 $I$，只有 $\Phi$"。
> 这里给出了 $\Phi$ 的显式构造。**$\Phi$ 奇异正是"某些方向推不动"的数学表达**。

### 法 2：用 $T$ 引入 $\lambda$ → 鞍点系统

引入 $6\times n_c$ 矩阵 $T$ 张成 $\mathcal S^{\perp}$：

$$
T^{\mathsf T}v=0\tag{3.56}
$$

$$
f_c=T\lambda\tag{3.57}
$$

$$
 T^{\mathsf T}a+\dot T^{\mathsf T}v=0\tag{3.58}
$$

联立得

$$
\begin{bmatrix}I&T\\ T^{\mathsf T}&0\end{bmatrix}
\begin{bmatrix}a\\ -\lambda\end{bmatrix}
=\begin{bmatrix}f-p\\ -\dot T^{\mathsf T}v\end{bmatrix}\tag{3.60}
$$

**原书对系数矩阵的评价**：*"symmetric and nonsingular, but **not positive definite**."*

> ⚠️ **不正定 ⟹ 不能用 Cholesky**，必须用 LDL$^{\mathsf T}$ 或 QR。
> 这个结论在第 8 章闭环系统里会再次出现。

### 法 3：广义坐标形式

由功率相等定义广义力：$\tau^{\mathsf T}\dot q=f^{\mathsf T}v=f^{\mathsf T}S\dot q$ 对所有 $\dot q$ 成立，故

$$
\tau=S^{\mathsf T}f\tag{3.61}
$$

代入 3.51 得标准形式

$$
H\ddot q+C=\tau\tag{3.62}
$$

$$
\boxed{\ H=S^{\mathsf T}IS\ }\tag{3.63}
$$

$$
\boxed{\ C=S^{\mathsf T}(I\dot S\dot q+p)\ }\tag{3.64}
$$

> 💡 **式 3.63 是"$H$ 从空间惯性投影而来"的最小例子**。
> 第 6 章的 $H_{ij}=S_i^{\mathsf T}I^{c}_i\,{}^{i}X_jS_j$ 就是它在树上的推广。

### ⭐ 法 3 的信息损失（原书的重要评注）

> *"Equation 3.61 implies a **small loss of information** in going from $f$ to $\tau$.
> **Infinitely many different values of $f$ map to the same value of $\tau$**,
> and therefore produce the same acceleration, but each produces a different value of $f_c$."*

**后果**：只知道 $\tau$ 时，可以算出 $f+f_c$ 的和，但**算不出各自的值**。

| 目标                                       | 丢失的信息是否重要 |
| ------------------------------------------ | ------------------ |
| **运动仿真**                               | ❌ 无关紧要        |
| **机械系统设计**（如计算关节轴承的动载荷） | ✅**相关**         |

> 🔑 这与例 3.2 的 $T_a$ 讨论是同一件事的两面：
> **广义坐标形式把内力信息压缩掉了。** 做仿真是优点（维数低），
> 做设计是缺点（要另外把内力找回来）。

---

## 3.7 多体系统的动力学 (Dynamics of a Multibody System)

本节用 **§3.2 的方法 1**（全收集 + 全施加约束）走一遍完整流程，
把 §3.5/3.6 的关节材料和 §3.2 的广义约束连起来。

**设置**：固定基座为 body 0；$N_B$ 个运动刚体（编号 $1..N_B$）；
$N_J$ 个关节（编号 $1..N_J$）。对每个关节 $j$，$p(j)$、$s(j)$ 是其前驱/后继的刚体编号。
**这 $2N_J$ 个数定义了系统的连通性。** 本节对连通性**不作任何限制**（可以有回路）。

### 六个步骤

#### 步骤 1：刚体的运动方程

$$
f=Ia+p\tag{3.66}
$$

$f,a,p$ 是 $6N_B$ 维向量，$I$ 是 $6N_B\times6N_B$ **块对角**矩阵（式 3.65）。
$p_i=v_i\times^{*}I_iv_i-f_{gi}$。

#### 步骤 2：由刚体运动得到关节运动

$$
v_J=P^{\mathsf T}v\tag{3.67}
$$

$$
a_J=P^{\mathsf T}a\tag{3.68}
$$

$$
P_{ij}=\begin{cases}
\mathbf 1_{6\times6}&i=s(j)\quad(\text{后继})\\
-\mathbf 1_{6\times6}&i=p(j)\quad(\text{前驱})\\
\mathbf 0_{6\times6}&\text{否则}
\end{cases}\tag{3.69}
$$

**$P$ 的结构**：$6N_B\times6N_J$，按 $N_B\times N_J$ 的 $6\times6$ 块组织。
**行 = 刚体，列 = 关节。**

- 每**列**至多两个非零块：后继行是 $+\mathbf 1$，前驱行是 $-\mathbf 1$。
  连到固定基座的关节只有**一个**非零块。
- 每**行** $i$：以 body $i$ 为后继的关节处是 $+\mathbf 1$，以它为前驱的处是 $-\mathbf 1$。

> 💡 **$P$ 就是图论里的关联矩阵 (incidence matrix)**，只是每个元素换成 $6\times6$ 块。
> **系统的全部拓扑信息都在 $P$ 里。**

#### 步骤 3：由关节力得到刚体力

$$
f=P\,f_J\tag{3.70}
$$

即：作用在 body $i$ 上的力 = （以它为后继的关节力之和）−（以它为前驱的关节力之和）。

> 🔑 **注意步骤 2 用 $P^{\mathsf T}$、步骤 3 用 $P$——这是对偶性的又一次体现**：
> 同一个矩阵，作用在运动上要转置，作用在力上不转置。

#### 步骤 4：运动约束

逐关节的加速度约束 $T_j^{\mathsf T}a_{Jj}+\dot T_j^{\mathsf T}v_{Jj}=0$ 收集起来：

$$
T^{\mathsf T}a_J+\dot T^{\mathsf T}v_J=0\tag{3.71}
$$

结合 3.67、3.68：

$$
\boxed{\ T^{\mathsf T}P^{\mathsf T}a+\dot T^{\mathsf T}P^{\mathsf T}v=0\ }\tag{3.72}
$$

#### 步骤 5：约束力

$$
f_J=T_a\tau+T\lambda\tag{3.73}
$$

（$T_a$、$T$ 是块对角矩阵。）

#### 步骤 6：最终方程

代入即得

$$
\boxed{\ \begin{bmatrix}I&PT\\ T^{\mathsf T}P^{\mathsf T}&0\end{bmatrix}
\begin{bmatrix}a\\ -\lambda\end{bmatrix}
=\begin{bmatrix}PT_a\tau-p\\ -\dot T^{\mathsf T}P^{\mathsf T}v\end{bmatrix}\ }\tag{3.74}
$$

**与 §3.2 的对应**（原书明确列出）：

| 式 3.74 中的                        | 对应 §3.2 中的   |
| ----------------------------------- | ----------------- |
| $T^{\mathsf T}P^{\mathsf T}$        | $K$               |
| $-\dot T^{\mathsf T}P^{\mathsf T}v$ | $k$               |
| $I,\ a,\ p$                         | $H,\ \ddot q,\ C$ |
| $PT_a\tau$                          | $\tau$            |
| $PT\lambda$                         | $\tau_c$          |
| 式 3.74 整体                        | **式 3.17**       |

> ✅ `verify_ch03.py` 用三体分支树完整实现了式 3.74，
> 解出的 $\ddot q$ 与 ABA 一致到机器精度（$10^{-16}$）。
> **这条对拍验证了本章的整个框架。**

### 原书的"讨论"——为什么这不是一个好算法

原书自己指出式 3.74 作为算法基础的**若干未解决问题**，很值得看：

1. **坐标系问题**：速度和加速度变量显然是各刚体的 Plücker 坐标，**但在哪个坐标系里？**
2. **没提位置变量**，也没说怎么算 $I$、$T$、$T_a\tau$、$\dot T^{\mathsf T}v$。
3. **怎么解**：系数矩阵**大**，但也**稀疏**——用稀疏分解会大幅提效。
4. **可能奇异**：若式 3.72 的约束不全线性无关，系数矩阵奇异，
   需要能处理奇异矩阵的求解过程。
5. **约束漂移**：式 3.72 只在**加速度层**指定约束，
   仿真中的舍入和截断误差会在位置和速度约束上累积（→ **§8.3**）。

> 🔑 **这五条正是后续各章要解决的问题清单**：
> 第 4 章解决 1、2（建模与坐标系），
> 第 5–7 章解决 3（用递推代替解大方程），
> 第 8 章解决 4、5（闭环、过约束、约束稳定化）。

---

## 本章要点回顾

1. **式 3.1 是全书的锚**：$H\ddot q+C=\tau$；$C$ 是"产生零加速度所需的 $\tau$"。
2. **构造运动方程只有两种操作**（收集、加约束），**顺序不同 = 算法不同**（四种方法表）。
3. **约束的隐式/显式二象性**：$K$（不许往哪走）vs $G$（只许往哪走），$KG=0$。
4. **Jourdain 虚功率原理**给出 $\tau_c=K^{\mathsf T}\lambda$ 和 $G^{\mathsf T}\tau_c=0$。
5. **子空间才是物理，矩阵不是**：$S\to SA$ 不改变任何东西。
6. **正交补必须用对偶形式**，不能把 $M^6$ 当欧氏空间（Duffy 1990）。
7. **关节有三个子空间**：$S$（运动）、$T=\mathcal S^{\perp}$（约束力）、$T_a$（主动力，不唯一）。
8. **$S^{\mathsf T}f_J=\tau$**（式 3.37）——RNEA 内推的理论依据。
9. **常见关节 $\mathring S=0$、$\sigma=0$ ⟹ $c_J=0$**，别被一般公式吓到。
10. **受约束刚体三种解法**：$\Phi=S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}$、鞍点系统、广义坐标。
11. **例 3.1 的 $IS(S^{\mathsf T}IS)^{-1}S^{\mathsf T}$ 是 ABA Schur 补的种子。**
12. **广义坐标形式会丢失内力信息**——仿真无所谓，机械设计要注意。

---

## 易错点

1. **$C$ 是向量、含重力**（与多数教材的 $C$ 矩阵不同）。
2. **把 $M^6$ 当欧氏空间求正交补**——原书专门警告的经典错误。
3. **忘了 $T_a$**：只有 $S$ 和 $T$ 是不完整的模型。
4. **被式 3.40–3.43 的一般形式吓退**：常见关节 $c_J=0$。
5. **在过约束系统里以为 $\lambda$ 唯一**：$\ddot q$ 唯一但 $\lambda$ 未必。
6. **对式 3.60 / 3.74 的鞍点矩阵用 Cholesky**——它不正定。

---

## 与其他章的联系

- ← 第 2 章：$\hat f=I\hat a+\hat v\times^{*}I\hat v$、式 2.69 的偏置力自由度、$\Phi$
- → 第 4 章：解决"讨论"中的坐标系与建模问题；$P$ 的结构变成 parent array
- → 第 5 章：式 3.37 $S^{\mathsf T}f_J=\tau$ = RNEA 内推
- → 第 6 章：式 3.63 $H=S^{\mathsf T}IS$ 的树上推广
- → 第 7 章：**例 3.1 / 式 3.54 是 ABA 的数学核心**；§3.2 方法 3 是 ABA 的结构
- → 第 8 章：§3.2 方法 2 = 闭环的标准做法；式 3.17 = 闭环的 KKT
- → 第 9 章：$\sigma$ 的替代方案是混合动力学
- → 第 11 章：§3.4 的不等式约束分支

---

## ✍️ 我的理解

<!-- 建议：用自己的话说清「四种构造方法」各对应哪一章的算法 -->

## ❓ 疑问与待办

- [ ]  自己推一遍式 3.20（投影法），确认 $G^{\mathsf T}$ 左乘就消掉了 $\tau_c$
- [ ]  验证式 3.28 在 $S\to SA$ 下不变（把 $A$ 消掉的那几步）
- [ ]  例 3.2 换一个齿轮半径 $r$，确认运动方程不变但约束力变
- [ ]  读完第 7 章后回来看例 3.1，确认 ABA 的 Schur 补就是它
- [ ]  找一个过约束系统的实例（如平面四杆按空间机构建模），观察 $\lambda$ 的不定性

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.39–64 逐节核对 -->
