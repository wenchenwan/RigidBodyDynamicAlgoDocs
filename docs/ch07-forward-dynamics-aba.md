# 第 7 章 正动力学·传播法 (Forward Dynamics — Propagation Methods) —— ABA

> **原书 pp. 119–140**，共 5 节。★ 三大核心算法之一，也是全书最精妙的一章。

> **一句话概括**：不构造 $H$，用**三趟递推**直接求出 $\ddot q$。
> 铰接体算法 (ABA)，复杂度 $O(N_B)$——**这是求解正动力学的理论下界**。

## 传播法的思想（原书开篇，值得逐字读）

> 正动力学问题给了我们**两组未知量**：关节加速度和关节约束力。
> 通常**无法在任何单个刚体上局部解出**其中任何一个；
> **但可以列出它们必须满足的方程**。
>
> 传播法的做法是：**局部计算这些方程的系数，并把系数传播给相邻刚体**，
> 直到最终到达某个点，在那里可以**局部地**解出动力学。
> 有了这一个刚体上的解，就能解相邻刚体，如此下去直到整个问题解完。

> 🔑 **原书自己给出的类比**：
> *"In principle, the process is **akin to solving a linear equation by first
> triangularizing the coefficient matrix, and then solving it by back-substitution**."*
>
> **传播法 = 三角化 + 回代。** §7.4 会把这件事做实（真的写出复合矩阵做高斯消元）。

**复杂度**：$O(N_B)$，**这是求解正动力学问题的理论最小值**。
代价是算法比惯性矩阵法**更复杂**。

## 原书节次对照

| 节 | 原书标题 | 页 | 内容 |
|---|---|---|---|
| 7.1 | Articulated-Body Inertia | 119 | **铰接体惯性的定义与性质**、多柄、例 7.1 |
| 7.2 | Calculating Articulated-Body Inertias | 123 | **两种策略：投影法与装配法** |
| 7.3 | The Articulated-Body Algorithm | 128 | **ABA 与表 7.1** |
| 7.4 | Alternative Assembly Formulae | 131 | **高斯消元的系统化推导**与变体 |
| 7.5 | Multiple Handles | 136 | 多柄装配、分治并行算法 |

---

## 7.1 铰接体惯性 (Articulated-Body Inertia)

### 定义

$$
\boxed{\ \textbf{铰接体惯性}=\text{一个刚体}\textbf{作为刚体系统的一部分时}\text{表现出来的惯性}\ }
$$

**对比（原书图 7.1）**：

| | 方程 | 系数 |
|---|---|---|
| (a) 孤立刚体 $B$ | $f=Ia+p$ (7.1) | 刚体惯性、刚体偏置力 |
| (b) $B$ 经关节连着 $B'$ | $f=I^Aa+p^A$ (7.2) | **铰接体惯性、铰接体偏置力** |

> **两式代数形式相同，但系数不同。**

### 两个关键术语

$$
\boxed{
\begin{aligned}
&\textbf{柄 (handle)}：I^A\text{ 与 }p^A\text{ 所指的那个刚体}\\
&\textbf{铰接体 (articulated body)}：\text{其动力学效应被计入 }I^A,p^A\text{ 的那个（子）系统}
\end{aligned}}
$$

### 与刚体惯性共有的性质

- 对称、**正定**矩阵
- $M^6\to F^6$ 的映射
- **遵循相同的坐标变换规则**

### ⚠️ 与刚体惯性的关键差别

> 原书原话：*"However, articulated-body inertias are mappings from
> **acceleration to force**, not **velocity to momentum**, so expressions like
> **$I^Av$ and $\tfrac12v^{\mathsf T}I^Av$ do not make sense**."*

$$
\boxed{\ I^A\text{ 没有「动量」和「动能」的含义！}\ }
$$

> 🔑 **这是最容易犯的概念错误。** $I$ 既能把 $a$ 映成力、也能把 $v$ 映成动量；
> 而 $I^A$ **只有前一种含义**。
> 看到 $\tfrac12v^{\mathsf T}I^Av$ 就该警觉——**那不是任何东西的动能**。

### 其他性质

**依赖关系**：$I^A$ 只依赖成员刚体的惯性和成员关节施加的约束，
因此是**关节位置变量的函数**，**与速度变量和各种力项无关**。

> 💡 **这一点与广义惯性矩阵 $H$ 完全类似**。
> 速度项、外力项（$f$ 之外的）等**只出现在偏置力里**。

**参数个数**（与第 2 章 §2.13-C 呼应）：

| | 3D | 2D |
|---|---|---|
| 刚体惯性 | 10 | 4 |
| **铰接体惯性** | **21** | **6** |

21 和 6 分别是对称 $6\times6$ 和 $3\times3$ 矩阵的独立元个数
——**铰接体惯性用满了所有自由度**。

**表观质量的方向性**：

> *"In general, the **apparent mass** of an articulated-body inertia
> **varies with the direction of applied force**, and it **does not have a
> centre of apparent mass**."*

**退化情形**：铰接体可以只含一个刚体，此时 $I^A=I$、$p^A=p$。
**这为递推计算提供了起点。**

### 加法（原书图 7.2）

两个铰接体 $A_1$、$A_2$，把它们的柄 $B_1$、$B_2$ **粘在一起**成为新刚体 $B$，则

$$
\boxed{\ I^A=I^A_1+I^A_2,\qquad p^A=p^A_1+p^A_2\ }
$$

> 💡 **注意与"用关节连接"的区别**：**粘死**才能直接相加；
> **用关节连**要用 §7.2.2 的装配公式（有 Schur 补）。

### 逆惯性（柄的自由度 < 6 时）

式 7.2 假定柄有 6 个自由度。否则必须用更一般的形式：

$$
a=\Phi^Af+b^A\tag{7.3}
$$

> **原书：这种形式的方程 *always exists*。**

| 情形 | $\Phi^A$ | $I^A$ |
|---|---|---|
| 柄有 6 自由度 | $\Phi^A=(I^A)^{-1}$，$b^A=-\Phi^Ap^A$ | 存在 |
| 柄自由度 < 6 | **奇异**，$\mathrm{rank}(\Phi^A)=$ 柄的运动自由度数 | **不存在** |

$\Phi^A$ 是对称**半正定**矩阵，变换规则同刚体逆惯性。

### 多柄 (Multiple Handles)

一个铰接体可以有**多个柄**。$h$ 个柄时式 7.3 推广为

$$
\begin{bmatrix}a_1\\ \vdots\\ a_h\end{bmatrix}
=\begin{bmatrix}\Phi^A_1&\cdots&\Phi^A_{1h}\\ \vdots&\ddots&\vdots\\ \Phi^A_{h1}&\cdots&\Phi^A_h\end{bmatrix}
\begin{bmatrix}f_1\\ \vdots\\ f_h\end{bmatrix}
+\begin{bmatrix}b^A_1\\ \vdots\\ b^A_h\end{bmatrix}\tag{7.4}
$$

**$\Phi^A_{ij}$ 是交叉耦合项**：描述"在柄 $j$ 上施力，在柄 $i$ 处引起的加速度"。

系数矩阵对称半正定，**秩等于这组柄总共拥有的独立运动自由度数**。
秩为 $6h$ 时式 7.4 可逆，得到式 7.2 的多柄版本。

### ⭐ 例 7.1：带槽的盒子 + 圆柱（把"方向性"讲透）

**设置**（原书图 7.3）：盒子上开一条槽，圆柱嵌在槽里。
圆柱相对盒子有 **2 个自由度**：绕 $x$ 轴转动、沿 $y$ 方向平移。接触**无摩擦**。
两者质心都在原点。盒子质量 $m_1$、绕质心转动惯量 $I_1$（各向同性）；
圆柱 $m_2$、$I_2$。

**把盒子当柄**，逐个方向做思想实验：

| 施加什么 | 谁在动 | 表观惯量 |
|---|---|---|
| 沿 $x$ 或 $z$ 的力 | 两者**一起**动 | $m_1+m_2$ |
| **沿 $y$ 的力** | **只有盒子**动 | $m_1$ |
| 绕 $y$ 或 $z$ 的力偶 | 两者**一起**转 | $I_1+I_2$ |
| **绕 $x$ 的力偶** | **只有盒子**转 | $I_1$ |

$$
I^A=\begin{bmatrix}
I_1&&&&&\\ &I_1{+}I_2&&&&\\ &&I_1{+}I_2&&&\\
&&&m_1{+}m_2&&\\ &&&&\mathbf{m_1}&\\ &&&&&m_1{+}m_2
\end{bmatrix}
$$

> 🔑 **看第 5 个对角元**：$m_1$ 而不是 $m_1+m_2$。
> **"盒子在 $x$、$z$ 方向显得比在 $y$ 方向更重。"**
> 这就是"表观质量随方向变化"的具体含义——
> 沿关节允许的方向推它很轻（圆柱会自己滑开），沿被约束的方向推它很重。

> 💡 **原书的收尾很重要**：本例中柄**恰好**有一个"表观质量中心"（在原点）——
> 任何沿过该点的直线施加的纯力都只产生纯线加速度。
> **但一般情况下这样的点并不存在。**

---

## 7.2 计算铰接体惯性

**原书给出两种基本策略**：

| 策略 | 做法 | 主要用途 |
|---|---|---|
| **投影法 (projection)** | 把**整个铰接体**的运动方程**投影到柄的运动空间**上 | **操作空间惯性**（机器人控制，Khatib 1987/1995） |
| **装配法 (assembly)** | 从组成部分**逐步装配**，沿途算出一串铰接体惯性 | **$O(n)$ 动力学算法** |

### 7.2.1 投影法

铰接体也是刚体系统，有广义坐标下的运动方程 $H\ddot q+C=\tau$（式 7.5），
其中 $C$ 涵盖除"作用在柄上的力"之外的一切力。

设 $J$ 是把 $\dot q$ 映到柄的空间速度的雅可比（$v=J\dot q$），则

$$
\tau=J^{\mathsf T}f\tag{7.6}
\qquad\qquad
a=J\ddot q+\dot J\dot q\tag{7.7}
$$

联立得

$$
a=JH^{-1}(J^{\mathsf T}f-C)+\dot J\dot q=\Phi^Af+b^A\tag{7.8}
$$

$$
\boxed{\ \Phi^A=JH^{-1}J^{\mathsf T}\ }\tag{7.9}
\qquad
b^A=\dot J\dot q-JH^{-1}C\tag{7.10}
$$

柄有完整 6 自由度时：

$$
\boxed{\ I^A=(JH^{-1}J^{\mathsf T})^{-1}\ }\tag{7.12}
\qquad
p^A=-I^Ab^A\tag{7.13}
$$

> 🔑 **式 7.12 就是操作空间惯性矩阵 $\Lambda$！**
> 我在第 8、9、11 章会反复见到 $(JH^{-1}J^{\mathsf T})^{-1}$——
> 现在知道它的正式身份了：**它是"以末端为柄"的铰接体惯性**。

**从式 7.12 能读出的性质**（原书指出）：

- $H$ 对称正定 ⟹ **$I^A$ 也对称正定**
- $I^A$ **显然不依赖速度或力项**

**$\Phi^A$ 的性质**（式 7.14）：

$$
\mathrm{rank}(\Phi^A)=\mathrm{rank}(J),\quad
\mathrm{range}(\Phi^A)=\mathrm{range}(J),\quad
\mathrm{null}(\Phi^A)=\mathrm{null}(J^{\mathsf T})
$$

**多柄推广**：把 $a,f,J$ 定义成各柄的堆叠向量/矩阵，式 7.6–7.14 **原样适用**。

### 7.2.2 装配法

**本节只考虑一个特例**（原书明说），即 ABA 所需的那个：
**每个铰接体都是浮动运动学树**（与固定基座无连接的运动学树）。它的特殊性质是：

$$
\boxed{
\begin{aligned}
&1.\ \text{每个铰接体}\textbf{恰有一个柄}\\
&2.\ \text{每个柄有}\textbf{完整 6 自由度}\\
&3.\ \text{每次装配操作 = 用}\textbf{一个关节}\text{把一个柄连到另一个柄}
\end{aligned}}
$$

（更一般的情形在 §7.4、§7.5。）

#### 推导（原书图 7.4）

装配前：$f_1=I^A_1a_1+p^A_1$（7.15），$f_2=I^A_2a_2+p^A_2$（7.16）。
装配后 $B_1$ 作为新柄：$f=I^Aa_1+p^A$（7.17）。

新关节引入未知力 $f_J$（从 $B_1$ 传到 $B_2$）：

$$
f_1=f-f_J,\qquad f_2=f_J
$$

$$
a_2-a_1=S\ddot q+c,\qquad S^{\mathsf T}f_J=\tau
$$

**解出 $\ddot q$**：

$$
\tau=S^{\mathsf T}f_2=S^{\mathsf T}\big(I^A_2(a_1+c+S\ddot q)+p^A_2\big)
$$

$$
\ddot q=(S^{\mathsf T}I^A_2S)^{-1}\big(\tau-S^{\mathsf T}(I^A_2(a_1+c)+p^A_2)\big)\tag{7.18}
$$

> 💡 **$S^{\mathsf T}I^A_2S$ 正定因而可逆**——这与第 3 章式 3.51 的
> $S^{\mathsf T}IS$ 是同一件事，也是 ABA **无需选主元**的保证。

**代回**得到装配公式：

$$
\boxed{\ I^A=I^A_1+I^A_2-I^A_2S(S^{\mathsf T}I^A_2S)^{-1}S^{\mathsf T}I^A_2\ }\tag{7.19}
$$

$$
p^A=p^A_1+p^A_2+I^A_2c+I^A_2S(S^{\mathsf T}I^A_2S)^{-1}\big(\tau-S^{\mathsf T}(I^A_2c+p^A_2)\big)\tag{7.20}
$$

> 🎯 **认出式 7.19 的结构**：$I^A_2-I^A_2S(S^{\mathsf T}I^A_2S)^{-1}S^{\mathsf T}I^A_2$
> **正是第 3 章例 3.1 的 $f_a=IS(S^{\mathsf T}IS)^{-1}S^{\mathsf T}f$ 的对偶**！
> 例 3.1 挑出"产生加速度的那部分力"，这里挑出"被关节消化掉的那部分惯性"。
> **第 3 章那个不起眼的例子，正是 ABA 的数学核心。**

#### 子树版本（实现用的形式）

定义 $A_1\cdots A_{N_B}$，其中 $A_i$ 含以 body $i$ 为根的子树、$i$ 为柄。则

$$
I^A_i=I_i+\sum_{j\in\mu(i)}I^a_j\tag{7.21}
\qquad\qquad
p^A_i=p_i+\sum_{j\in\mu(i)}p^a_j\tag{7.22}
$$

$$
I^a_j=I^A_j-I^A_jS_j(S_j^{\mathsf T}I^A_jS_j)^{-1}S_j^{\mathsf T}I^A_j\tag{7.23}
$$

$$
p^a_j=p^A_j+I^a_jc_j+I^A_jS_j(S_j^{\mathsf T}I^A_jS_j)^{-1}(\tau_j-S_j^{\mathsf T}p^A_j)\tag{7.24}
$$

> **这些方程由式 7.19、7.20 对 body $i$ 的每个子节点各应用一次得到**：
> 每次应用时，$A_1$ 是"body $i$ 加上已处理完的子树"，$A_2$ 是下一棵要加进来的子树。

#### ⭐ $I^a_j$ 与 $p^a_j$ 的三重含义

**含义 1（原书原话）**：$I^a_j$、$p^a_j$ 是**把 $I^A_j$、$p^A_j$ 跨过关节 $j$ 传播的结果**。

**含义 2**：它们可以看作

$$
\boxed{\ \text{经关节 }j\text{ 连到 body }j\text{ 的一个}\textbf{无质量柄}\text{的表观惯性与偏置力}\ }
$$

于是式 7.21、7.22 的作用就是**把所有这些无质量柄粘到 body $i$ 上**（回到图 7.2 的加法）。

**含义 3（最有洞察力的一条，式 7.25）**：

$$
\boxed{\ f_j=I^A_ja_j+p^A_j=I^a_ja_{\lambda(j)}+p^a_j\ }\tag{7.25}
$$

> 🔑 **$I^A_j$ 把关节力表达成"**子**的加速度"的函数；
> $I^a_j$ 把**同一个力**表达成"**父**的加速度"的函数。**
>
> 这就是"跨过关节传播"的确切含义：换了自变量。
> 内推时父节点只关心"用我的加速度怎么表达子传来的力"，
> 所以传上去的必须是 $I^a$ 而不是 $I^A$。

#### $I^a_jS_j=0$：一个可利用的性质

原书指出 $I^a_j$ 在提高效率上有**两个作用**：

1. 它让式 7.20 中两个含 $c$ 的项**并成一项**（就是式 7.24 的 $I^a_jc_j$）；
2. **$I^a_jS_j=0$**。

> 💡 **第 2 条的用处**：视 $S_j$ 的取值，$I^a_j$ 的**一行一列会是零**。
> 例如转动关节 $S_j=[0\,0\,1\,0\,0\,0]^{\mathsf T}$ ⟹ $I^a_j$ 的**第 3 行和第 3 列全零**。
> 利用这一点可以适度降低计算铰接体惯性的代价。
>
> **物理解释**：关节允许的方向上，子树"推不动"任何东西——
> 力沿 $S_j$ 施加时子树会自己让开。

---

## 7.3 铰接体算法 (The Articulated-Body Algorithm)

### 核心洞察

考虑 $B_1$（图 7.5b）。若 $f_1$ 是经关节 1 传递的力：

$$
f_1=I^A_1a_1+p^A_1\tag{7.26}
$$

$$
a_1=a_0+c_1+S_1\ddot q_1\tag{7.27}
\qquad\qquad
S_1^{\mathsf T}f_1=\tau_1\tag{7.28}
$$

三式联立解出

$$
\ddot q_1=(S_1^{\mathsf T}I^A_1S_1)^{-1}\big(\tau_1-S_1^{\mathsf T}I^A_1(a_0+c_1)-S_1^{\mathsf T}p^A_1\big)\tag{7.29}
$$

> 🔑 **原书紧接着的这段是全章的关键**：
>
> *"Consider what we have just achieved. **Simply knowing $I^A_1$ and $p^A_1$ has
> allowed us to calculate $\ddot q_1$ directly, without needing to know the values of
> any of the other acceleration variables.** This is possible because Eq. 7.26
> **already accounts for the dynamic effect of every body in $A_1$** on the
> acceleration response of $B_1$."*
>
> **$I^A$ 的全部价值就在这一句：它把"整棵子树的动力学效应"压缩成一个 $6\times6$ 矩阵，
> 于是关节 1 的加速度可以"就地"解出，不必先知道任何别的加速度。**

有了 $\ddot q_1$ 代回 7.27 得 $a_1$，就能对 $B_1$ 的子节点重复这个过程。一般地：

$$
\ddot q_i=(S_i^{\mathsf T}I^A_iS_i)^{-1}\big(\tau_i-S_i^{\mathsf T}I^A_i(a_{\lambda(i)}+c_i)-S_i^{\mathsf T}p^A_i\big)\tag{7.30}
$$

$$
a_i=a_{\lambda(i)}+c_i+S_i\ddot q_i\tag{7.31}
$$

### 三趟结构

$$
\boxed{
\begin{aligned}
&\textbf{趟 1（外推，根→叶）}：\text{算速度与偏置项 }v_i,c_i,p_i\\
&\textbf{趟 2（内推，叶→根）}：\text{算铰接体惯性与偏置力 }I^A_i,p^A_i\\
&\textbf{趟 3（外推，根→叶）}：\text{算加速度 }\ddot q_i,a_i
\end{aligned}}
$$

**趟 1** —— 式 7.32–7.36：

$$
v_{Ji}=S_i\dot q_i,\quad c_{Ji}=\mathring S_i\dot q_i,\quad
v_i=v_{\lambda(i)}+v_{Ji},\quad c_i=c_{Ji}+v_i\times v_{Ji},\quad
p_i=v_i\times^{*}I_iv_i-f^x_i
$$

> 与 RNEA 第一趟**高度相似**（原书明说）。

**趟 2** —— 式 7.37–7.40（即 7.21–7.24）。
⚠️ **原书的实现提醒**：这意味着 $I^a_j$、$p^a_j$ **不是对每个 $j$ 都算**，
只对满足 $\lambda(j)\ne0$ 的 $j$ 算。**若 body $i$ 无子节点，则 $I^A_i=I_i$、$p^A_i=p_i$。**

**趟 3** —— 式 7.41、7.42，且

$$
a_0=-a_g
$$

> 原书：*"By initializing the base acceleration to $-a_g$, rather than 0, we can
> simulate the effect of a uniform gravitational field with a fictitious upward
> acceleration. **This is more efficient than treating gravity as an external force.**"*

### 公共子表达式（式 7.43–7.50）

原书指出式 7.39–7.42 有大量公共子表达式，定义中间量：

$$
U_i=I^A_iS_i\tag{7.43}
\qquad
D_i=S_i^{\mathsf T}U_i\tag{7.44}
\qquad
u_i=\tau_i-S_i^{\mathsf T}p^A_i\tag{7.45}
\qquad
a'_i=a_{\lambda(i)}+c_i\tag{7.46}
$$

于是化简为

$$
\boxed{
\begin{aligned}
I^a_j&=I^A_j-U_jD_j^{-1}U_j^{\mathsf T}&&(7.47)\\
p^a_j&=p^A_j+I^a_jc_j+U_jD_j^{-1}u_j&&(7.48)\\
\ddot q_i&=D_i^{-1}(u_i-U_i^{\mathsf T}a'_i)&&(7.49)\\
a_i&=a'_i+S_i\ddot q_i&&(7.50)
\end{aligned}}
$$

⚠️ **式 7.48 里是 $I^a_j$（Schur 补之后的），不是 $I^A_j$**——这是最常见的 ABA bug。

### ⭐ 表 7.1：完整伪代码（原书）

```
────────────────────────────────────────────────────────────────
# ═══ 趟 1：外推 ═══
v₀ = 0
for i = 1 to N_B do
    [X_J, S_i, vJ, cJ] = jcalc(jtype(i), q_i, q̇_i)
    ⁱX_λ(i) = X_J · X_T(i)
    if λ(i) ≠ 0 then
        ⁱX₀ = ⁱX_λ(i) · λ⁽ⁱ⁾X₀
    end
    v_i   = ⁱX_λ(i) · v_λ(i) + vJ
    c_i   = cJ + v_i × vJ
    I^A_i = I_i
    p^A_i = v_i ×* I_i · v_i − ⁱX₀* · f^x_i
end

# ═══ 趟 2：内推 ═══
for i = N_B to 1 do
    U_i = I^A_i · S_i
    D_i = S_iᵀ · U_i
    u_i = τ_i − S_iᵀ · p^A_i
    if λ(i) ≠ 0 then
        I^a = I^A_i − U_i·D_i⁻¹·U_iᵀ
        p^a = p^A_i + I^a·c_i + U_i·D_i⁻¹·u_i
        I^A_λ(i) = I^A_λ(i) + λ⁽ⁱ⁾X_i* · I^a · ⁱX_λ(i)
        p^A_λ(i) = p^A_λ(i) + λ⁽ⁱ⁾X_i* · p^a
    end
end

# ═══ 趟 3：外推 ═══
a₀ = −a_g
for i = 1 to N_B do
    a' = ⁱX_λ(i) · a_λ(i) + c_i
    q̈_i = D_i⁻¹ · (u_i − U_iᵀ · a')
    a_i = a' + S_i · q̈_i
end
────────────────────────────────────────────────────────────────
```

**原书对伪代码的说明**：

- 局部变量：趟 1 的 `X_J, vJ, cJ`；趟 2 的 `I^a, p^a`；趟 3 的 `a'`
- 含时关节 ⟹ 给 `jcalc` 加第四个参数（时间）
- 速度变量不是位置变量导数 ⟹ 把 $\dot q_i,\ddot q_i$ 换成 $\alpha_i,\dot\alpha_i$
- **无外力时不需要算 ${}^{i}X_0$**
- 伪代码假定 $X$ 一算出来，$X^{*}$、$X^{-1}$、$(X^{*})^{-1}$ **立即可用**
  （实现见**附录 A**）
- 式 7.37、7.38 已从 $\mu(i)$ 改写成 $\lambda(i)$：
  **$I^A_i$、$p^A_i$ 在第一趟循环末尾初始化为 $I_i$、$p_i$**，
  其余计算在第二趟的 `if` 里做（同 §5.3 的技巧）

### 为什么必须是三趟

| 趟 | 方向 | 算什么 | 为什么不能合并 |
|---|---|---|---|
| 1 | 外推 | $v_i,c_i$，$I^A/p^A$ 初值 | 速度只依赖父节点 ⟹ **必须外推** |
| 2 | 内推 | $I^A,p^A,U,D,u$ | 铰接体惯性依赖**整个子树** ⟹ **必须内推** |
| 3 | 外推 | $\ddot q_i,a_i$ | 加速度依赖父节点的加速度 ⟹ **必须外推** |

方向是"外-内-外"，交替两次，**所以三趟是下界**。

> 💡 **对比 RNEA 只需两趟**：因为那里 $\ddot q$ 是已知输入。
> **ABA 比 RNEA 多的那一趟，本质上就是"解方程"。**

---

## 7.4 ⭐ 替代装配公式 —— 高斯消元的系统化推导

> **这一节把"ABA = 三角化 + 回代"从类比变成了字面事实。**

### 问题的复合矩阵形式

§7.2.2 解的问题是：给定五个方程（7.51–7.55）

$$
f_1=I^A_1a_1+p^A_1,\quad f_2=I^A_2a_2+p^A_2,\quad f=f_1+f_2,
\quad a_2=a_1+c+S\ddot q,\quad S^{\mathsf T}f_2=\tau
$$

求目标方程 $f=I^Aa_1+p^A$（7.56）的系数。

**原书把这五个方程组装成一个复合方程**：

$$
\begin{bmatrix}
1&-1&-1&0&0\\
0&1&0&0&0\\
0&0&1&-I^A_2&0\\
0&0&0&1&-S\\
0&0&S^{\mathsf T}&0&0
\end{bmatrix}
\begin{bmatrix}f\\ f_1\\ f_2\\ a_2\\ \ddot q\end{bmatrix}
=
\begin{bmatrix}0\\ p^A_1\\ p^A_2\\ c\\ \tau\end{bmatrix}
+
\begin{bmatrix}0\\ I^A_1\\ 0\\ 1\\ 0\end{bmatrix}a_1
$$

### 组装的五条规则（原书）

$$
\boxed{
\begin{aligned}
&1.\ \text{目标方程}\textbf{右端的未知量}(a_1)\ \to\ \text{复合方程右端的未知量}\\
&2.\ \text{目标方程}\textbf{左端的未知量}(f)\ \to\ \text{复合未知向量的}\textbf{顶部}\\
&3.\ \text{维数不是 6 的未知向量}(\ddot q)\ \to\ \text{复合未知向量的}\textbf{底部}\\
&4.\ \text{维数不是 6 的给定方程}(S^{\mathsf T}f_2=\tau)\ \to\ \text{复合方程的}\textbf{最后一行}\\
&5.\ \text{其余方程与未知量任意排列，使 }5\times5\text{ 系数矩阵}\textbf{尽可能接近上三角}
\end{aligned}}
$$

### 求解 = 高斯消元 + 回代

> *"Having constructed this equation, the original problem can be solved using
> a process of **Gaussian elimination and back-substitution**."*

高斯消元阶段结束时，**第五行**读作

$$
S^{\mathsf T}I^A_2S\,\ddot q=\tau-S^{\mathsf T}p^A_2-S^{\mathsf T}I^A_2c-S^{\mathsf T}I^A_2a_1
$$

即式 7.18。然后**回代**：先用 $a_1$ 表达 $a_2$，……最终得到用 $a_1$ 表达 $f$，
此时 **$a_1$ 的系数就是 $I^A$，其余项就是 $p^A$**。

> 🔑 **这正式坐实了第 6 章 §6.6 引用的 Ascher et al. (1997)**：
> **ABA 就是对某个置换后的系数矩阵做高斯消元。**
> 而 $I^A-UD^{-1}U^{\mathsf T}$ 就是消元产生的 **Schur 补**，$D_i$ 就是**主元**。
> $I^A$ 正定 ⟹ $D_i>0$ ⟹ **ABA 天然不需要选主元**。

### 两个变体（原书表 7.2）

想让这套推导适用于**不同的给定方程**，有两个明显的改动：

1. 把式 7.51、7.52、7.56 换成它们的**逆**（用 $\Phi^A$、$b^A$）；
2. 把式 7.54、7.55 换成**等价的隐式关节约束**：
   $T^{\mathsf T}a_2=T^{\mathsf T}(a_1+c)$ 与 $f_2=T_a\tau+T\lambda$。

| 版本 | 改动 | 适用条件 |
|---|---|---|
| **标准** | 无 | **两个柄都必须有完整 6 自由度** |
| **变体 2** | 只做改动 2 | 同上 |
| **变体 1** | 两个改动都做 | **只需 $T^{\mathsf T}(\Phi^A_1+\Phi^A_2)T$ 满秩** |

**变体 1 条件的物理含义**（原书）：

> *"Physically, this is equivalent to requiring that the joint **does not impose a
> redundant constraint**; that is, it does not duplicate a constraint that was
> already there."*

**充分条件**：两个柄中**至少有一个**具有完整 6 自由度。

> 🔑 **变体 1 的价值**：基于它的算法**能为任意运动学树、以及少数闭环系统
> 计算铰接体逆惯性**，而标准算法**只限于浮动运动学树**。
>
> ⚠️ **代价**（原书脚注）：两个变体都**不会顺便算出 $\ddot q$**，
> 但可以用 $\ddot q=T_a^{\mathsf T}(a_2-a_1-c)$ 补算。
> 基于变体 1 的算法若把主动关节力 $T_a\tau$ 换成一对等效外力，会更高效。

**一个副产品**（原书特别提到）：这个过程会产生**一些可能有用的矩阵恒等式**——
变体给出的 $I^A$ 公式必须等于标准公式，且都必须等于 $\Phi^A$ 公式的逆。

---

## 7.5 多柄 (Multiple Handles)

**动机**：装配运动学树只需单柄。但若想装配**一般刚体系统**，
或想在**装配顺序上有更大灵活性**，就必须能装配多柄铰接体。

> 🔑 **原书举的成果**：Featherstone (1999a,b) 的**分治算法**
> 能装配运动学树和闭环系统，使其动力学可在
> $O(N_B)$ 个处理器的并行计算机上以 **$O(\log N_B)$ 时间**算出。
>
> 这就是第 3 章 §3.2 "方法 4"（子系统两两合并）的落地。

**为什么这里必须用逆惯性**（原书的理由很重要）：

> *"We are using **inverse inertias** here because it is possible—and indeed
> **quite likely**—that the two handles within an articulated body **do not have a
> full 12 degrees of motion freedom between them**, even if they do have
> 6 degrees of freedom individually."*

此时系数矩阵**奇异**，$I^A$ 不存在，只能用 $\Phi^A$。

**假设**：柄 2、3 中至少一个有 6 自由度 ⟹ $\Phi_2$ 与 $\Phi_3$ 至少一个正定
⟹ **$\Phi_2+\Phi_3$ 正定**（这是求解的关键）。

**装配约束**：

$$
T^{\mathsf T}(a_3-a_2)=T^{\mathsf T}c\tag{7.60}
\qquad\qquad
f_3=T\lambda=-f_2\tag{7.61}
$$

⚠️ **式 7.61 的一个建模约定**：假定关节上的主动力**已经计入**式 7.57、7.58 的偏置加速度——
即把 $T_a\tau$ 换成了作用在柄 3 和柄 2 上的一对**等值反向外力**。

代入消元得

$$
T^{\mathsf T}(\Phi_2+\Phi_3)T\,\lambda=T^{\mathsf T}(\Phi_{21}f_1-\Phi_{34}f_4+c-b_3+b_2)
$$

---

## ABA vs CRBA

| 维度 | ABA | CRBA + 分解 |
|---|---|---|
| 复杂度 | $O(N_B)$（**理论下界**） | $O(nd^2)$ ~ $O(n^3)$ |
| 常数因子 | **大** | 小 |
| 能否得到 $H$ | ❌ | ✅ |
| 利用分支稀疏 | 自然的树递推 | ✅ 显式利用，$\rho^2$ 效应 |
| 内存 | $O(n)$ | $O(n^2)$ |
| 数值稳定性 | 好（$D_i>0$，**无需选主元**） | 好（重排 Cholesky） |
| 实现难度 | 中 | 低 |
| 加约束/接触 | 困难 | 容易（进 KKT） |
| 适用范围 | **仅浮动运动学树**（标准版） | 更广 |

**交叉点**：见第 10 章。原书第 6 章开篇已说明：
$O(nd^2)$ 算法在**刚体较少或分支足够多**时可匹敌甚至超过 $O(n)$ 算法。

---

## 易错点与陷阱

1. **把 $I^A$ 当成能算动量/动能的东西**。$I^Av$、$\tfrac12v^{\mathsf T}I^Av$ **无意义**（§7.1）。
2. **$p^a$ 里用了 $I^A$ 而不是 $I^a$**（式 7.48）——最常见的 ABA bug。
   症状：静态正确，一动就错。
3. **趟 3 的 $a'$ 忘了加 $c_i$**。
4. **$c_i$（趟 1 算）、$U_i,D_i,u_i$（趟 2 算）忘了存下来给趟 3 用**。
5. **$X$ / $X^{*}$ 用混**：$I^a$ 用三明治 $X^{*}I^aX$，$p^a$ 是力向量用 $X^{*}$。
6. **$D_i$ 接近零**：理论上 $I^A$ 正定 ⟹ $D_i>0$。若 $D_i\approx0$ 说明模型有问题
   （零惯量连杆、$S$ 退化）。**加断言**。
7. **重力技巧的位置**：$a_0=-a_g$ 在**趟 3** 开头，不是趟 1。
8. **对非浮动树或闭环用标准公式**：标准装配公式**要求两个柄都有 6 自由度**（§7.4）。

## 调试建议

- **RNEA 互验（最有效）**：`τ' = RNEA(q, q̇, ABA(q, q̇, τ))`，应有 $\tau'=\tau$。
  两算法完全独立，几乎不可能同时错成一致。
- **与 CRBA 路线对拍**：$\ddot q_{\text{ABA}}$ 应等于 $H^{-1}(\tau-C)$。
- **$I^aS=0$ 检验**：转动关节时 $I^a$ 的对应行列应为零。
- **$I^A$ 对称正定**、$D_i>0$。
- **自由落体**：单刚体 + 自由关节 + $\tau=0$ ⟹ $a=a_g$。
- **静止**：$\dot q=0$、$\tau=g(q)$ ⟹ $\ddot q=0$。
- **能量守恒**：无驱动无耗散时总能量漂移缓慢且有界。

全部实现在 `code/verify_all.py`。

## 与其他章的联系

- ← 第 2 章：全部代数工具；$I$ 的 10 参数 vs $I^A$ 的 21 参数（§2.13-C）
- ← 第 3 章：**例 3.1 / 式 3.54 是式 7.19 的数学核心**；§3.2 方法 3 = 装配法
- ← 第 5 章：趟 1、趟 3 与 RNEA 外推同构
- ↔ 第 6 章：另一条 FD 路线；§6.6 的 Ascher et al. (1997) 与本章 §7.4 互相印证
- → 第 8 章：变体 1 可用于部分闭环系统
- → 第 9 章：**铰接体混合动力学**（§9.2）直接改造本章算法
- → 第 10 章：与 CRBA 的代价交叉点
- → 第 11 章：$\Lambda=(JH^{-1}J^{\mathsf T})^{-1}$ 就是式 7.12 的投影法结果
- → 附录 A：$X$、$X^{*}$、$X^{-1}$ 的高效实现

---

## ✍️ 我的理解

<!-- 建议：用自己的话说清 I^A 与 I^c 的区别、以及式 7.25 的「换自变量」含义 -->

## ❓ 疑问与待办

- [ ] 实现 ABA，与 RNEA 互验
- [ ] 自己推一遍式 7.19、7.20（从 7.15–7.18 出发）
- [ ] **照着 §7.4 的五条规则写出复合矩阵，手工做一遍高斯消元**，
      确认最终得到式 7.19、7.20
- [ ] 验证 $I^a_jS_j=0$，看转动关节时零掉的是哪一行列
- [ ] 用例 7.1 的盒子+圆柱手算 $I^A$，与代码对拍
- [ ] 计算某个末端的 $\Lambda=(JH^{-1}J^{\mathsf T})^{-1}$，与投影法式 7.12 对照

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.119–140 逐节核对。
     此前版本遗漏了 §7.1 的多柄与例 7.1、§7.2.1 投影法（即操作空间惯性）、
     §7.4 的高斯消元系统化推导与两个变体、§7.5 多柄装配，现已补上。 -->
