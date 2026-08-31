# 第 5 章 逆动力学 (Inverse Dynamics) —— RNEA

> **原书 pp. 89–100**，共 5 节。★ 三大核心算法之一。

> **一句话概括**：给定 $q,\dot q,\ddot q$ 求 $\tau$。
> 递推牛顿-欧拉算法 (RNEA) 用**一趟外推 + 一趟内推**，以 $O(n)$ 代价解决。
> 本章还系统解释了**为什么递推算法比非递推算法快几个数量级**。

## 原书节次对照

| 节 | 原书标题 | 页 | 内容 |
|---|---|---|---|
| 5.1 | Algorithm Complexity | 89 | 运算次数、big-O、$O(N_B^p)=O(n^p)$ |
| 5.2 | Recurrence Relations | 90 | **递推关系为什么高效**（本章的理论核心） |
| 5.3 | The Recursive Newton-Euler Algorithm | 92 | 三步法、体坐标版、**表 5.1 伪代码** |
| 5.4 | The Original Version | 97 | Luh et al. (1980) 的 3D 版本与本版的对应 |
| 5.5 | Additional Notes | 99 | 历史与文献 |

> ⚠️ **我此前的笔记完全漏掉了 §5.1 和 §5.2。** §5.2 是全书"递推思想"的正式论述，
> 解释了整本书为什么长这样，**比 RNEA 本身更值得读**。

---

## 5.1 算法复杂度 (Algorithm Complexity)

### 计算代价的定义

$$
\boxed{\ \textbf{计算代价 (computational cost)}=\textbf{运算次数 (operations count)}\ }
$$

对动力学算法这种实数运算，**惯例是只数浮点运算**——
计算数组下标、递增循环变量之类的**整数运算被忽略**。

### big-O 的精确定义

算法具有 $O(f(n))$ 复杂度，当且仅当存在正常数 $C_{\min},C_{\max},n_{\min}$，
使得对所有 $n>n_{\min}$，计算代价**至少** $C_{\min}f(n)$、**至多** $C_{\max}f(n)$。

> 💡 **原书的一个有用视角**：$O(f(n))$ 可以看作**集合**——
> 所有随 $n\to\infty$ 按 $f(n)$ 比例增长的函数的集合。
> 于是 "算法 X 是 $O(f(n))$" 意思是"X 的代价函数**属于**这个集合"，
> 而 $O(f(n))=O(g(n))$ 是**集合相等**。

**取哪个函数当 $O$ 的参数**：集合中**最简单**的那个。
若精确代价是多项式 $a_0+a_1n+\cdots+a_pn^p$，则记作 $O(n^p)$。
一般只有**增长最快的那一项**才重要；若 $f$ 收敛到常数则记 $O(1)$。

### ⚠️ 复杂度依赖于"问题规模怎么度量"

> 原书举的例子很好：矩阵求逆的复杂度是 $O(n^3)$，
> **这依赖于 $n$ 是矩阵的维数**。
> 若 $n$ 改为矩阵的**元素个数**，同一个算法的复杂度就变成 $O(n^{1.5})$。

一般地，若 $m$ 和 $n$ 是两种规模度量且 $m=g(n)$，则 $O(f(m))=O(f(g(n)))$。

### 对刚体系统：$N_B$ 与 $n$ 可以混用

运动学树有两种合理的规模度量：**刚体数 $N_B$** 和**自由度数 $n$**。
若每个关节至少有 1 个自由度，则

$$
N_B\le n\le 6N_B
\quad\Longrightarrow\quad
\boxed{\ O(N_B^{\,p})=O(n^{p})\ }\tag{5.1}
$$

> 🔑 **所以引用复杂度时不必区分 $N_B$ 和 $n$。** 本书全程混用。

---

## 5.2 ⭐ 递推关系 (Recurrence Relations)

> **这一节解释了整本书为什么长这样。原书原话**：
>
> *"Modern dynamics algorithms are **recursive**, and they are **orders of magnitude
> faster** than their non-recursive predecessors. ... **It is the use of recurrence
> relations that accounts for their efficiency**."*

### 什么是递推关系

**递推关系 (recurrence relation)** = 用**前面的**元素定义序列**下一个**元素的公式。
它与一组初值一起，给出序列的递归定义。

**例（原书给的）**：Fibonacci 数列 $0,1,1,2,3,5,8,13,\dots$
由 $F_i=F_{i-1}+F_{i-2}$ 加初值 $F_0=0,F_1=1$ 递归定义。

**在刚体系统中**：按 §4.1.2 编号后，刚体速度序列 $v_0,v_1,v_2,\dots$
由递推关系 $v_i=v_{\lambda(i)}+S_i\dot q_i$ 加初值 $v_0$ 递归定义。

### 定量对比：速度计算

设无分支树、1-DoF 关节（$\lambda(i)=i-1$）。

**递推形式**：

$$
v_i=v_{i-1}+s_i\dot q_i\qquad(v_0=0)\tag{5.2}
$$

**非递推（闭式）形式**：

$$
v_i=\sum_{j=1}^{i}s_j\dot q_j\tag{5.3}
$$

设 $m$ = 一次"标量×向量"的代价，$a$ = 一次向量加法的代价：

| | 单次代价 | 算前 $n$ 个的总代价 | 复杂度 |
|---|---|---|---|
| **式 5.2（递推）** | $m+a$ | $n(m+a)$ | $O(n)$ ✅ |
| **式 5.3（闭式）** | $i\,m+(i-1)a$ | $\tfrac12\big(n(n{+}1)m+n(n{-}1)a\big)$ | $O(n^2)$ ❌ |

**低效的来源一目了然**（原书把它摊开写）：

$$
\begin{aligned}
v_1&=s_1\dot q_1\\
v_2&=s_1\dot q_1+s_2\dot q_2\\
&\ \vdots\\
v_n&=s_1\dot q_1+s_2\dot q_2+\cdots+s_n\dot q_n
\end{aligned}
$$

> 乘积 $s_1\dot q_1$ 被算了 **$n$ 遍**，而它只需要算 **1 遍**。
> $s_2\dot q_2$ 被算了 $n-1$ 遍，和 $s_1\dot q_1+s_2\dot q_2$ 也算了 $n-1$ 遍……
>
> 🔑 **递推关系提供了一种系统性的办法，避免在计算一串相关量时做不必要的重复。
> 这就是它高效的秘密。**

### 更极端的对比：加速度计算

**递推**（式 5.2 的时间导数）：

$$
a_i=a_{i-1}+s_i\ddot q_i+v_i\times s_i\dot q_i\tag{5.4}
$$

**闭式**：

$$
a_i=\sum_{j=1}^{i}s_j\ddot q_j+\sum_{j=1}^{i}\sum_{k=1}^{j-1}(s_k\dot q_k)\times(s_j\dot q_j)\tag{5.5}
$$

$$
\boxed{\ \text{式 5.4：}O(n)\qquad\text{式 5.5：}\mathbf{O(n^3)}\ }
$$

> 💡 **注意这里的差距比速度那里更大**（$n$ vs $n^3$，而不是 $n$ vs $n^2$）：
> **计算越复杂，递推的优势越大。**

### 闭式方程与历史教训

**闭式 (closed form)** = 直接用最基本的组成成分（关节轴向量、关节速度/加速度变量）
表达的方程。刚体系统的运动方程写成闭式是：

$$
\tau_i=\sum_{j=1}^{n}H_{ij}\ddot q_i+\sum_{j=1}^{n}\sum_{k=1}^{n}C_{ijk}\dot q_j\dot q_k+g_i\tag{5.6}
$$

（$C_{ijk}$ 是科氏力与离心力的系数——**这就是 Christoffel 符号**。）

> 🔑 **原书给出的历史数字，值得记住**：
>
> - 早期**非递推**算法就是用式 5.6 这类方程计算的，复杂度典型为 $\mathbf{O(n^4)}$
> - 递推算法的复杂度可低至 $\mathbf{O(n)}$
> - **RNEA 刚出现时，对一个 6 自由度机械臂，比它的非递推前身
>   Uicker/Kahn 算法快了大约 100 倍**（Hollerbach, 1980）
>
> *"Prior to the invention of efficient recursive algorithms, closed-form equations
> were the standard way of expressing and evaluating equations of motion."*

> 💡 **对今天的实践含义**：如果你在某本教材里看到 2 连杆机械臂那种
> 手推出来的、几页纸长的闭式运动方程——那是式 5.6 的形式。
> 它对**理解**有帮助，但**永远不要**把它作为计算的基础。

---

## 5.3 递推牛顿-欧拉算法 (The Recursive Newton-Euler Algorithm)

### 三个步骤（原书的组织方式）

$$
\boxed{
\begin{aligned}
&1.\ \text{计算树中每个刚体的}\textbf{速度和加速度}\\
&2.\ \text{计算}\textbf{产生这些加速度所需的力}\\
&3.\ \text{由作用在刚体上的力，计算}\textbf{关节间传递的力}
\end{aligned}}
$$

### 步骤 1：运动学（外推）

$$
v_i=v_{\lambda(i)}+S_i\dot q_i\qquad(v_0=0)\tag{5.7}
$$

$$
a_i=a_{\lambda(i)}+S_i\ddot q_i+\dot S_i\dot q_i\qquad(a_0=0)\tag{5.8}
$$

（式 5.8 就是 5.7 的时间导数。）

### 步骤 2：每个刚体所需的净力（局部）

$$
f^B_i=I_ia_i+v_i\times^{*}I_iv_i\tag{5.9}
$$

### 步骤 3：关节力（内推）

参考原书图 5.1：$f_i$ 是经关节 $i$ 从 $\lambda(i)$ 传给 $i$ 的力，
$f^x_i$ 是作用在 body $i$ 上的**外力**（力场、物理接触等，**作为算法的输入**）。

body $i$ 上的净力：

$$
f^B_i=f_i+f^x_i-\sum_{j\in\mu(i)}f_j
$$

整理成递推关系：

$$
\boxed{\ f_i=f^B_i-f^x_i+\sum_{j\in\mu(i)}f_j\ }\tag{5.10}
$$

> 💡 **原书特别指出**：式 5.10 **不需要初值**，因为无子节点的刚体 $\mu(i)=\varnothing$。

最后投影出广义力（即第 3 章的式 3.37）：

$$
\tau_i=S_i^{\mathsf T}f_i\tag{5.11}
$$

### 重力技巧（原书的原话）

> *"External forces can be used to model a variety of environmental influences,
> **including gravitational forces**. However, a **uniform** gravitational field can be
> modelled **more efficiently** as a **fictitious acceleration of the base**."*

做法：把式 5.8 的初值换成 $a_0=-a_g$。

⚠️ **副作用（原书明确写出）**：

$$
\boxed{\ \text{用了这个技巧后，算出的 }a_i\text{ 和 }f^B_i\text{ 不再是真实的加速度和净力}\ }
$$

它们分别偏移了重力加速度和重力力向量。**body $i$ 上的重力是 $I_ia_g$**。

> 💡 **记住 $f_g=I\,a_g$ 这个式子**：它说明"重力空间力 = 空间惯性 × 重力空间加速度"。
> 展开验证：$I[0;g]=[mc\times g;\ mg]$——力是 $mg$、关于原点的矩是 $c\times(mg)$ ✓

### 算法特征（原书 "Algorithm Features"）

**数据流向**：

| 方程 | 传播方向 |
|---|---|
| 5.7、5.8 | 根 → 叶（**外推**） |
| 5.10 | 叶 → 根（**内推**） |
| 5.9、5.11 | **局部**（每个刚体自己算） |

$$
\boxed{\ \text{RNEA 是}\textbf{两趟}\text{算法：一趟外推算 }v,a\text{；一趟内推算关节力}\ }
$$

> **刚体力 $f^B_i$ 的计算可以放在第一趟或第二趟，原书说"前者通常是更好的选择"。**

### RNEA 算出的不只是 $\tau$

> *"This algorithm calculates **more than just** $\tau_i$. Other quantities, like $v_i$ and
> $f_i$, can sometimes be useful. For example, one of the tasks of a machine designer
> is to **select joint bearings** of sufficient strength to withstand the dynamic
> reaction forces during operation. These forces can be calculated from $f_i$."*

> 🔑 **呼应第 3 章**：那里说广义坐标形式会**丢失内力信息**。
> RNEA 在内推过程中**天然算出了 $f_i$**——所以它既能做仿真也能做机械设计。
> 很多实现只返回 $\tau$ 而丢掉 $f_i$，**如果你要做受力分析，记得把它留出来**。

### 通过置零输入得到特化版本

| 置零的输入 | 保留 | 得到 |
|---|---|---|
| $f^x$、$\dot q$、$\ddot q$ | 虚拟基座加速度 | **重力补偿项**（机器人控制用） |
| 重力、$\ddot q$、$f^x$ | $\dot q$ | **科氏力与离心力项** |

**重力补偿的精简版**（原书式 5.12、5.13）：

$$
f_i=-I_ia_g+\sum_{j\in\mu(i)}f_j\tag{5.12}
\qquad\qquad
\tau_i=S_i^{\mathsf T}f_i\tag{5.13}
$$

> 💡 要算**真实的重力项**（而非补偿项），把 $-I_ia_g$ 换成 $+I_ia_g$。
> 注意这个精简版**完全不含运动学**——不需要算 $v_i$、$a_i$，快得多。

### 体坐标系版本

> *"In practice, the algorithm works best if the calculations pertaining to body $i$
> are performed in **body $i$ coordinates**."*

$$
v_i={}^{i}X_{\lambda(i)}v_{\lambda(i)}+S_i\dot q_i\qquad(v_0=0)\tag{5.14}
$$

$$
a_i={}^{i}X_{\lambda(i)}a_{\lambda(i)}+S_i\ddot q_i+\mathring S_i\dot q_i+v_i\times S_i\dot q_i
\qquad(a_0=-a_g)\tag{5.15}
$$

**相对基本版做了两处改动**（原书列出）：

1. 插入坐标变换，把 $v_{\lambda(i)}$、$a_{\lambda(i)}$ 从 $\lambda(i)$ 系变到 $i$ 系；
2. **$\dot S_i$ 被替换成 $\mathring S_i+v_i\times S_i$**（$\mathring S_i$ 是体坐标系中的表观导数）。

> 🔑 **第 2 点是关键**：这正是第 2 章 §2.10 那条提示的落地——
> "在与动系瞬时重合的**静**系中列写方程"，于是绝对导数被拆成
> "表观导数 + 叉乘项"，而叉乘项被显式写了出来。

利用 `jcalc` 提供的量：

$$
v_{Ji}=S_i\dot q_i\tag{5.16}
\qquad
c_{Ji}=\mathring S_i\dot q_i\tag{5.17}
$$

$$
v_i={}^{i}X_{\lambda(i)}v_{\lambda(i)}+v_{Ji}\qquad(v_0=0)\tag{5.18}
$$

$$
a_i={}^{i}X_{\lambda(i)}a_{\lambda(i)}+S_i\ddot q_i+c_{Ji}+v_i\times v_{Ji}\qquad(a_0=-a_g)\tag{5.19}
$$

$$
f_i=f^B_i-{}^{i}X_0^{*}f^x_i+\sum_{j\in\mu(i)}{}^{i}X_j^{*}f_j\tag{5.20}
$$

> ⚠️ **式 5.20 假设外力来自系统外部，因而在绝对（body 0）坐标系中给出**，
> 所以要做 ${}^{i}X_0^{*}$ 变换。

### ⭐ 表 5.1：完整伪代码（原书）

```
────────────────────────────────────────────────────────────────
v₀ = 0
a₀ = −a_g

for i = 1 to N_B do
    [X_J, S_i, vJ, cJ] = jcalc(jtype(i), q_i, q̇_i)
    ⁱX_λ(i) = X_J · X_T(i)
    if λ(i) ≠ 0 then
        ⁱX₀ = ⁱX_λ(i) · λ⁽ⁱ⁾X₀
    end
    v_i = ⁱX_λ(i) · v_λ(i) + vJ
    a_i = ⁱX_λ(i) · a_λ(i) + S_i·q̈_i + cJ + v_i × vJ
    f_i = I_i·a_i + v_i ×* I_i·v_i − ⁱX₀* · f^x_i
end

for i = N_B to 1 do
    τ_i = S_iᵀ · f_i
    if λ(i) ≠ 0 then
        f_λ(i) = f_λ(i) + λ⁽ⁱ⁾X_i* · f_i
    end
end
────────────────────────────────────────────────────────────────
```

**原书对伪代码的四点说明**：

1. `X_J`、`vJ`、`cJ` 是**局部变量**，每次迭代取新值。
2. 若树中有"速度变量不是位置变量导数"的关节，
   把 $\dot q_i$、$\ddot q_i$ 换成 $\alpha_i$、$\dot\alpha_i$ 即可。
3. 若有**显含时间**的关节，`jcalc` 会自动用式 3.32、3.41 算 $v_J$、$c_J$；
   伪代码里**唯一的改动是给 `jcalc` 多传一个时间参数**。
4. **${}^{i}X_0$ 唯一的用途是变换外力**。没有外力时，
   计算这些变换的 `if` 语句**可以整个省掉**。

### 从 $\mu(i)$ 到 $\lambda(i)$：内推的实现技巧

**式 5.20 是用 $\mu(i)$（子节点集）写的**，直译成伪代码是：

```
for i = N_B to 1 do
    f_i = f^B_i − ⁱX₀* f^x_i
    for each j in μ(i) do
        f_i = f_i + ⁱX_j* f_j
    end
end
```

**但同样的计算可以换个顺序做，结果相同**：

```
for i = 1 to N_B do
    f_i = f^B_i − ⁱX₀* f^x_i            ← 并入第一趟循环的最后一行
end
for i = N_B to 1 do
    if λ(i) ≠ 0 then
        f_λ(i) = f_λ(i) + λ⁽ⁱ⁾X_i* f_i   ← 「累加到父节点」
    end
end
```

> 🔑 **这就是为什么最终伪代码里看不到 $\mu(i)$**——
> 它被"每个节点把自己的力累加给父节点"这个等价写法替换掉了。
> **好处**：不需要显式维护子节点列表，只需要 `parent` 数组。
>
> ⚠️ **实现陷阱**：必须是 `+=` 而不是 `=`。分支节点有多个子节点，
> 写成赋值会丢掉除最后一个之外的所有分支。

---

## 5.4 原始版本 (The Original Version)

RNEA 最早由 **Luh et al. (1980a)** 提出（后收入 Brady et al. 1982），
用 **3D 向量**表述，表面上与 §5.3 差别很大。

> 🔑 **原书的判断**：*"the two algorithms are **actually almost the same**,
> the **only significant difference** being that one uses **spatial accelerations**
> while the other uses **classical accelerations**."*

### 3D 版本的一个缺点

> *"One of the drawbacks of using 3D vectors to express rigid-body dynamics is that
> the equations are **sensitive to joint type**."*

Luh et al. 的算法**同时针对转动和移动关节**，因此**充满**这样的式子：

```
variable = { 一个表达式    若关节是转动
           { 另一个表达式  若关节是移动
```

> 💡 **对比空间向量版本**：关节类型的差异被**完全封装在 `jcalc` 里**，
> 算法主体一行都不用改。这是空间记法在**软件工程**上的收益，
> 而不只是"少写点代数"。

### 两个版本的对应（原书表 5.2 摘要）

$$
\hat v_i=\begin{bmatrix}\omega_i\\ v_i\end{bmatrix},\qquad
\hat a_i=\begin{bmatrix}\dot\omega_i\\ \dot v_i-\omega_i\times v_i\end{bmatrix},\qquad
\hat f^B_i=\begin{bmatrix}N_i+c_i\times F_i\\ F_i\end{bmatrix},\qquad
\hat f_i=\begin{bmatrix}n_i\\ f_i\end{bmatrix}
$$

$$
{}^{i}X_{i-1}=\begin{bmatrix}{}^{i}E_{i-1}&0\\ -{}^{i}E_{i-1}\,{}^{i-1}r_i\times&{}^{i}E_{i-1}\end{bmatrix},
\qquad
\hat I_i=\begin{bmatrix}I^{cm}_i-m_ic_i\times c_i\times&m_ic_i\times\\ -m_ic_i\times&m_i\mathbf 1\end{bmatrix}
$$

**只有两处 3D 量与空间向量的分量不对应**（原书明确指出）：

| 不对应的量 | 原因 | 严重性 |
|---|---|---|
| **线加速度** | 空间 vs 经典加速度 | **非平凡**——它影响刚体间传播的量 |
| **净力矩 $N_i$** | $N_i$ 在**质心**处算，而非体坐标系原点 | 平凡——用 $N_i$ 算 $n_i$ 时立即修正 |

### ⭐ 一个有趣的性质：3D 版本不需要算线速度

> *"One interesting feature of the original algorithm is that **the linear velocity is
> never used in any subsequent expression, and therefore need not be calculated**.
> This is a special property that arises from using **classical accelerations**."*

**结果**：省掉线速度计算后，3D 版本比空间版本**略快**。

> 但原书立刻补了一句：*"the difference is only **a few percent**, and is easily
> **dwarfed** by other efficiency-related matters, such as the choice of programming
> language, the computer hardware, and so on."*
>
> 🔑 **这是很诚实的一段话**：不要为了几个百分点放弃空间记法带来的
> 可读性、可维护性和通用性。

### 两个版本的坐标系解释（原书脚注 1）

这条脚注把整件事讲透了：

| 版本 | 在什么坐标系中列写方程 |
|---|---|
| **表 5.1（空间）** | 与体坐标系**当前瞬时重合的静止**坐标系 |
| **表 5.2（3D）** | 与体坐标系当前瞬时重合、**不转动但具有与体坐标系原点相同线速度**的坐标系 |

> 🔑 **差别就在"跟不跟着平移"**：
> 空间版的观察者完全不动；3D 版的观察者跟着原点平移（但不跟着转）。
> 这正是 §2.11 说的"经典加速度是空间速度在带纯线速度 $\dot r$ 的坐标系中的表观导数"。

---

## 5.5 补充说明与历史 (Additional Notes)

**RNEA 的谱系**：

| 年份 | 工作 | 意义 |
|---|---|---|
| 1976 | Stepanenko & Vukobratovic | 递推牛顿-欧拉的**早期形式** |
| 1979 | Orin et al. | 同上 |
| **1980** | **Luh, Walker & Paul** | **通常被视为 RNEA 本身** |
| 1980 | Hollerbach | 递推 Lagrange 算法；**对比了各种算法的复杂度**，指出 Uicker/Kahn (1971) 是 $O(n^4)$ 而递推算法是 $O(n)$ |
| 1988/89 | Balafoutis et al.；He & Goldenberg | **已发表的最快版本之二** |

**当时的研究动机**（原书说明）：开发**快到足以用于机器人运动实时控制**的逆动力学算法。
早期例子是 Luh et al. (1980b)，但**操作空间 (operational-space) 表述最终更流行**（Khatib, 1987）。

**当时计算机慢，因而有两条支线**：

- **并行计算实现**（Lathrop, 1985）
- **符号优化技术**（Murray & Neuman, 1984；Neuman & Murray, 1987）→ **第 10 章 §10.4**

**一个新方向**（原书特别提到）：

> *"the development of efficient algorithms to calculate the **partial derivative of
> inverse dynamics with respect to a design parameter**."* （Fang & Pollard, 2003）

> 💡 **这条今天特别重要**：可微仿真、基于梯度的最优控制、
> 强化学习中的解析梯度，全都需要动力学的导数。
> 现代库（Pinocchio 的 `computeRNEADerivatives` 等）正是沿这条线发展的。

---

## 实现要点与易错点

1. **加速度交叉项是 $v_i\times v_{Ji}$**，用**已合成的** $v_i$，不是 $v_{\lambda(i)}$；
   **且没有系数 2**（第 2 章 §2.11）。
2. **重力符号**：$a_0=-a_g$（负号）。症状：重力补偿力矩方向全反。
3. **$\times$ 与 $\times^{*}$**：`a_i` 那行用 $\times$，`f_i` 那行用 $\times^{*}$。
4. **内推必须是 `+=`**（见 §5.3 末的实现技巧）。
5. **外力坐标系**：$f^x$ 在世界系给出，必须做 ${}^{i}X_0^{*}$ 变换。
6. **记住 $a_i$、$f^B_i$ 在用重力技巧后不是真值**，需要真值时加回 $a_g$、$I_ia_g$。
7. **常见关节 $c_J=0$**（第 3、4 章），可以省掉。

## 调试建议

- **RNEA ↔ ABA 互验**：$\mathrm{ID}(q,\dot q,\mathrm{FD}(q,\dot q,\tau))=\tau$
  ——**最有效的一条**，两算法完全独立。
- **与 $H\ddot q+C$ 对拍**：用 CRBA 独立算 $H$、用 RNEA 算 $C$。
- **重力技巧对拍**：分别用"$a_0=-a_g$"和"逐体施加外力 $f_g=I_ia_g$（变换到世界系）"，
  结果应完全一致。（`code/verify_all.py` 已验证。）
- **静态检验**：$\dot q=\ddot q=0$ 时输出应等于重力力矩。
- **单刚体检验**：$N_B=1$ + 自由关节，退化成牛顿-欧拉方程，可解析验证。

## 与其他章的联系

- ← 第 2 章：式 2.55 = 式 5.4；空间 vs 经典加速度 = §5.4 的全部差别
- ← 第 3 章：式 3.37 = 式 5.11
- ← 第 4 章：`model`、`jcalc`、$\lambda(i)<i$、$\mu(i)$
- → 第 6 章：$C=\mathrm{ID}(q,\dot q,0)$；以及"$n$ 次 RNEA 求 $H$"
- → 第 7 章：ABA 的第 1、3 趟与本章外推同构
- → 第 8 章：闭环系统的逆动力学（§8.12）
- → 第 9 章：混合动力学中已知加速度的那部分关节
- → 第 10 章：§10.4 符号优化；精确代价

---

## ✍️ 我的理解

<!-- 建议：不看笔记默写表 5.1 的两个循环 -->

## ❓ 疑问与待办

- [ ] 用 §5.2 的方法，自己数一遍式 5.4 与 5.5 的运算次数，确认 $O(n)$ vs $O(n^3)$
- [ ] 实现 RNEA，用 2 连杆平面臂的解析解验证
- [ ] 实现精简版式 5.12/5.13（只算重力补偿），与完整 RNEA 对拍
- [ ] 从 RNEA 中把 $f_i$ 提取出来，看看关节反力有多大（轴承选型的视角）
- [ ] 读第 10 章后回来补 RNEA 的精确运算次数

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.89–100 逐节核对。
     此前版本遗漏了 §5.1 复杂度 与 §5.2 递推关系两整节，现已补上。 -->
