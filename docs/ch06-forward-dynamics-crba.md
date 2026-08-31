# 第 6 章 正动力学·惯性矩阵法 (Forward Dynamics — Inertia Matrix Methods)

> **原书 pp. 101–118**，共 6 节。★ 三大核心算法之一（CRBA）。

> **一句话概括**：给定 $q,\dot q,\tau$ 求 $\ddot q$。
> 路线是**构造 $H$ 再解方程**。构造 $H$ 的最快算法是**复合刚体算法 (CRBA)**；
> 解方程时利用**分支诱导稀疏**做 $LTL$/$LTDL$ 分解。

## 两条正动力学路线（原书开篇）

对运动学树求正动力学，主要有两条路：

| 路线 | 做法 | 复杂度 | 本书章节 |
|---|---|---|---|
| **1. 惯性矩阵法** | 为**整个系统**建立运动方程，解出加速度变量 | $O(nd^2)$（最坏 $O(n^3)$） | **第 6 章** |
| **2. 传播法** | 把约束从一个刚体传播到下一个，使加速度**一个关节一个关节**地算出来 | $O(n)$ | 第 7 章 |

> 🔑 **原书对"$O(n^3)$ 算法"这个叫法的纠正**：
> *"Such algorithms have $O(n^3)$ complexity in the worst case, and are sometimes
> referred to collectively as **$O(n^3)$ algorithms**. However, it would be
> **more accurate** to describe them as **$O(nd^2)$ algorithms**, where $d$ is the depth
> of the tree."*

**权衡**（原书明确给出）：

- $O(n)$ 算法适合**刚体数很多**的系统；
- 但 $O(nd^2)$ 算法在**刚体较少**、或**分支足够多以致深度很小**的树上，
  可以**匹敌甚至略微超过** $O(n)$ 算法；
- $O(nd^2)$ 算法还是**闭环动力学算法的重要组成部分**（第 8 章）。

## 原书节次对照

| 节 | 原书标题 | 页 | 内容 |
|---|---|---|---|
| 6.1 | The Joint-Space Inertia Matrix | 102 | $H$ 的定义、$\mathrm{ID}_\delta$ 与表 6.1 |
| 6.2 | The Composite-Rigid-Body Algorithm | 104 | **CRBA 的推导与表 6.2** |
| 6.3 | A Physical Interpretation | 108 | **第二种推导：算法是怎么被发现的** |
| 6.4 | Branch-Induced Sparsity | 110 | 分支诱导稀疏、密度 $\rho$、无填充 |
| 6.5 | Sparse Factorization Algorithms | 112 | **表 6.3 $LTL$/$LTDL$、表 6.4 展开父数组、表 6.5 回代** |
| 6.6 | Additional Notes | 117 | 历史与两种方法的深层联系 |

公式与算法已数值验证：`python3 code/verify_crba_2link.py`、`code/verify_all.py`

---

## 6.1 关节空间惯性矩阵 (The Joint-Space Inertia Matrix)

$$
\tau=H(q)\,\ddot q+C(q,\dot q,f^x)\tag{6.1}
$$

> ⚠️ **注意 $C$ 的完整依赖是 $C(q,\dot q,f^x)$** —— 它还依赖**外力**。
> 原书说：*"The most useful piece of information they convey is the fact that
> **$H$ depends only on the position variables**."*

**若系统静止且除 $\tau$ 外无其他力**，则 $C=0$，方程简化为 $\tau=H\ddot q$。

### 三步法

$$
\boxed{
\begin{aligned}
&1.\ \text{计算 }C&&O(n)\\
&2.\ \text{计算 }H&&O(nd)\\
&3.\ \text{解 }H\ddot q=\tau-C&&O(nd^2)
\end{aligned}}
$$

- $d\le n$ ⟹ 步骤 2、3 的**最坏情形**是 $O(n^2)$、$O(n^3)$
- **另一个极端**：若 $d$ 有固定上界（实践中确实会出现），则**三步全是 $O(n)$**

### $H$ 的分块结构

$H$ 是**对称正定**的 $n\times n$ 矩阵，组织成 $N_B\times N_B$ 的子矩阵数组
（树形系统 $N_B=N_J$）。若关节 $i$ 的自由度是 $n_i$，则 $H_{ij}$ 是
占据块行 $i$、块列 $j$ 的 $n_i\times n_j$ 子矩阵。

$$
\boxed{\ \text{物理上，}H_{ij}\text{ 把}\textbf{关节 }j\textbf{ 的加速度}\text{关联到}\textbf{关节 }i\textbf{ 的力}\ }
$$

全 1-DoF 时 $n=N_B$，$H_{ij}$ 就是第 $i$ 行第 $j$ 列的标量。

### 用 ID 得到 $C$ 与 $H$

$$
C=\mathrm{ID}(model,q,\dot q,0,f^x)\tag{6.2}
$$

$$
H\ddot q=\mathrm{ID}(model,q,\dot q,\ddot q,f^x)-C\tag{6.3}
$$

> 🔑 **一次 ID 调用得到 $C$；两次 ID 调用之差得到 $H$ 与任意向量的乘积。**

### ⭐ 差分逆动力学 $\mathrm{ID}_\delta$ 与表 6.1

定义

$$
\mathrm{ID}_\delta(model,q,\ddot q)=\mathrm{ID}(model,q,\dot q,\ddot q,f^x)-\mathrm{ID}(model,q,\dot q,0,f^x)\tag{6.4}
$$

则 $H$ 的第 $\alpha$ 列 $=\mathrm{ID}_\delta(model,q,\delta_\alpha)$。
**逐列重复即得整个 $H$** —— 这就是 **Walker & Orin (1982) 的方法 1 和 2**。

> 🔑 **原书的关键观察**：
> *"There is a **great deal of cancellation** of terms on the right-hand side of Eq. 6.4.
> In particular, **every term depending on either $\dot q$ or $f^x$ cancels out**,
> which is why these vectors have been **omitted from the argument list** of $\mathrm{ID}_\delta$."*

于是可以写一个**大幅简化**的 RNEA 专门算 $\mathrm{ID}_\delta$（**原书表 6.1**）：

```
────────────────────────────────────────────────
a₀ = 0                          ← 注意：不是 −a_g！
for i = 1 to N_B do
    a_i = ⁱX_λ(i) · a_λ(i) + S_i·q̈_i
    f_i = I_i · a_i             ← 没有 v×*Iv 项！
end
for i = N_B to 1 do
    τ_i = S_iᵀ · f_i
    if λ(i) ≠ 0 then
        f_λ(i) = f_λ(i) + λ⁽ⁱ⁾X_i* · f_i
    end
end
────────────────────────────────────────────────
```

**与完整 RNEA 的差别**：**没有速度项、没有重力、没有外力**。
伪代码还省略了 $S_i$ 和 ${}^{i}X_{\lambda(i)}$ 的计算，
因为它们是"算 $C$ 时的副产品"。

> 💡 **这比我此前笔记里说的"跑 $n$ 次完整 RNEA"高效得多**。
> 但原书紧接着说：*"However, it is **not the most efficient way** to calculate $H$.
> The fastest algorithm for that job is the composite-rigid-body algorithm."*

---

## 6.2 复合刚体算法 (The Composite-Rigid-Body Algorithm)

### 推导：从动能出发（五步）

**第 1 步**：树的动能是各刚体动能之和

$$
T=\frac12\sum_{k=1}^{N_B}v_k^{\mathsf T}I_kv_k\tag{6.5}
$$

**第 2 步**：$v_k$ 关于 $\dot q$ 是线性的

$$
v_k=\sum_{i\in\kappa(k)}S_i\dot q_i\tag{6.6}
$$

**第 3 步**：代入

$$
T=\frac12\sum_{k=1}^{N_B}\sum_{i\in\kappa(k)}\sum_{j\in\kappa(k)}\dot q_i^{\mathsf T}S_i^{\mathsf T}I_kS_j\dot q_j\tag{6.7}
$$

**第 4 步**：**这是对"关节 $i$ 与关节 $j$ 都支撑刚体 $k$"的所有 $(i,j,k)$ 三元组求和**，
因此可以重排（用第 4 章的恒等式 4.3）：

$$
T=\frac12\sum_{i=1}^{N_B}\sum_{j=1}^{N_B}\sum_{k\in\nu(i)\cap\nu(j)}\dot q_i^{\mathsf T}S_i^{\mathsf T}I_kS_j\dot q_j\tag{6.8}
$$

其中

$$
\nu(i)\cap\nu(j)=\begin{cases}
\nu(i)&i\in\nu(j)\\
\nu(j)&j\in\nu(i)\\
\varnothing&\text{否则}
\end{cases}\tag{6.9}
$$

**第 5 步**：与 $T=\tfrac12\dot q^{\mathsf T}H\dot q$（式 6.10）对比即得

$$
\boxed{\ H_{ij}=\sum_{k\in\nu(i)\cap\nu(j)}S_i^{\mathsf T}I_kS_j\ }\tag{6.11}
$$

> 🎯 **式 6.11 是本章最重要的中间结论。** 它同时回答三个问题：
> **谁参与求和**（既在 $i$ 子树、又在 $j$ 子树里的刚体）、
> **何时为零**（交集为空）、**为什么会出现复合刚体惯性**（下一步）。

### 复合刚体惯性

$$
I^c_i=\sum_{j\in\nu(i)}I_j\tag{6.12}
\qquad\qquad
I^c_i=I_i+\sum_{j\in\mu(i)}I^c_j\tag{6.13}
$$

> **$I^c_i$ = 以 body $i$ 为根的子树、当作单个复合刚体时的惯性。**
> （**算法的名字就是从这来的。**）依据是第 2 章的"空间惯性可加"。

### ⭐ 结果：式 6.14

结合 6.11、6.12、6.9：

$$
\boxed{\ H_{ij}=\begin{cases}
S_i^{\mathsf T}I^{c}_{i}S_j & i\in\nu(j)\\[3pt]
S_i^{\mathsf T}I^{c}_{j}S_j & j\in\nu(i)\\[3pt]
0 & \text{否则}
\end{cases}\ }\tag{6.14}
$$

### ⚠️ 下标条件极易记反

**读法**：第一种情形 $i\in\nu(j)$ 意思是"**body $i$ 在关节 $j$ 的子树里**"，
即 **$j$ 是 $i$ 的祖先**；此时用的是 **$I^c_i$——后代 $i$ 的**复合刚体惯性。

**为什么必须是后代的 $I^c$**：回到式 6.9。$j$ 是祖先 ⟹ $\nu(i)\subseteq\nu(j)$ ⟹
交集取的是**小的那个** $\nu(i)$。

$$
\boxed{\ \textbf{一句话记忆：交集取小的，所以 }I^c\textbf{ 取后代的}\ }
$$

**无分支树的特例**（原书式 6.15、6.16）把这一点变得**显而易见**：

$$
I^c_i=I_i+I^c_{i+1}\tag{6.15}
\qquad\qquad
H_{ij}=S_i^{\mathsf T}I^{c}_{\max(i,j)}S_j\tag{6.16}
$$

> 🔑 **$\max(i,j)$！** 无分支树中编号大的就是后代，所以式 6.16 直接写着
> "取**编号较大**（= 更靠叶端）的那个 $I^c$"。**这是最好的记忆锚点。**

**记反的后果**：数值上会得到一个**对称、量纲正确、看起来完全合理但错误**的 $H$。
在 `code/verify_crba_2link.py` 的 2R 算例中，正确值 $H_{12}=0.15363$，记反得 $0.10755$。
**单看结果发现不了**，必须靠对拍（见"调试建议"）。

### 体坐标系版本

$$
I^c_i=I_i+\sum_{j\in\mu(i)}{}^{i}X_j^{*}I^c_j\,{}^{j}X_i\tag{6.17}
$$

$$
H_{ij}=\begin{cases}
S_i^{\mathsf T}I^{c}_{i}\,{}^{i}X_jS_j & i\in\nu(j)\\
S_i^{\mathsf T}\,{}^{i}X_j^{*}I^{c}_{j}S_j & j\in\nu(i)\\
0&\text{否则}
\end{cases}\tag{6.18}
$$

**为了高效计算，引入 ${}^{j}F_i$**：$I^c_iS_i$ 在 body $j$ 坐标系中的值。

$$
{}^{j}F_i={}^{j}X_i^{*}\,I^c_iS_i
$$

$$
{}^{\lambda(j)}F_i={}^{\lambda(j)}X_j^{*}\,{}^{j}F_i\qquad({}^{i}F_i=I^c_iS_i)\tag{6.19}
$$

$$
H_{ij}=\begin{cases}{}^{j}F_i^{\mathsf T}S_j& i\in\nu(j)\\ H_{ji}^{\mathsf T}& j\in\nu(i)\\ 0&\text{否则}\end{cases}\tag{6.20}
$$

> **需要为每一对满足 $j\in\kappa(i)$ 的 $(i,j)$ 准备一个 ${}^{j}F_i$**，
> 递推方式是：$i$ 取遍 $1..N_B$，对每个 $i$，$j$ 依次取
> $i,\lambda(i),\lambda(\lambda(i)),\dots$ 直到基座。

### ⭐ 表 6.2：完整伪代码（原书）

```
────────────────────────────────────────────────────────────────
H = 0
for i = 1 to N_B do
    I^c_i = I_i
end

for i = N_B to 1 do
    if λ(i) ≠ 0 then
        I^c_λ(i) = I^c_λ(i) + λ⁽ⁱ⁾X_i* · I^c_i · ⁱX_λ(i)
    end

    F    = I^c_i · S_i
    H_ii = S_iᵀ · F

    j = i
    while λ(j) ≠ 0 do
        F    = λ⁽ʲ⁾X_j* · F
        j    = λ(j)
        H_ij = Fᵀ · S_j
        H_ji = H_ijᵀ
    end
end
────────────────────────────────────────────────────────────────
```

**原书对伪代码的说明**：

- 开头两个准备步骤：把 $H$ 置零、把每个 $I^c_i$ 设为 $I_i$。
  **置零那步是提醒读者：后面的代码只会初始化 $H$ 中的非零子矩阵**
  （是否真的需要置零取决于 $H$ 怎么用）。
- `F` 是局部变量，初值 $I^c_iS_i$（即方程中的 ${}^{i}F_i$）。
- `while` 循环把 `F` 依次变换到路径上每个 body $j$ 的坐标系（即每个 $j\in\kappa(i)$）。
- **设置 $H_{ji}$ 只在需要访问上三角时才必要。**
- 伪代码**不含**计算 $S_i$ 和 ${}^{i}X_{\lambda(i)}$ 的代码，
  因为它们被假定是"用 RNEA 算 $C$ 时的副产品"。

**从 $\mu(i)$ 到 $\lambda(i)$**（与 §5.3 同样的技巧）：
式 6.17 直译需要遍历子节点，但可以改写成"每个节点把自己的 $I^c$ 累加给父节点"，
于是只需要 `parent` 数组。

### 复杂度

> **若 $d$ 是树的深度，则 `while` 循环对任何 $i$ 都不会执行超过 $d-1$ 次。**
> 因此 CRBA 的计算代价**不会超过 $O(nd)$**，
> 且 **$H$ 中的非零元数目也至多是 $O(nd)$**。

---

## 6.3 一个物理解释 (A Physical Interpretation)

> 原书用**第二种推导**重讲一遍 CRBA，*"the purpose being to show some of the
> **physical insights that played a role in the discovery** of this algorithm."*

### 思想实验：令 $\ddot q=\delta_\alpha$

系统静止、除 $\tau$ 外无外力 ⟹ $\tau=H\ddot q$。
令 $\ddot q=\delta_\alpha$，则 **$\tau$ 就是 $H$ 的第 $\alpha$ 列**。

$\ddot q$ 的第 $\alpha$ 个变量属于某个关节，称为**关节 $i$**；
$S_i$ 中对应的那一列称为 $s_\alpha$（一个空间运动向量）。

**树的运动状态**：关节 $i$ 上的加速度是 $s_\alpha$，其余关节加速度全为零。于是

$$
a_j=\begin{cases}s_\alpha& j\in\nu(i)\\ 0&\text{否则}\end{cases}\tag{6.21}
$$

> 🔑 **关键：其余关节被"冻结"，所以整个子树 $\nu(i)$ 像一个刚体那样以 $s_\alpha$ 运动。**

**无速度项**，故

$$
f^B_j=I_ja_j=\begin{cases}I_js_\alpha& j\in\nu(i)\\ 0&\text{否则}\end{cases}\tag{6.22}
$$

### 核心论证

> *"Now, **joint $j$ is the only connection between the subtree $\nu(j)$ and the rest of
> the system**. As there are no forces acting on the bodies other than those transmitted
> across the joints, it follows that **$f_j$ must be the net force acting on subtree $\nu(j)$**."*

$$
f_j=\sum_{k\in\nu(j)}f^B_k\tag{6.23}
$$

联立即得

$$
\boxed{\ f_j=\begin{cases}
I^c_js_\alpha& j\in\nu(i)\\[3pt]
I^c_is_\alpha& j\in\kappa(i)\\[3pt]
0&\text{否则}
\end{cases}\ }\tag{6.24}
$$

$$
\tau_j=S_j^{\mathsf T}f_j=\begin{cases}
S_j^{\mathsf T}I^c_js_\alpha& j\in\nu(i)\\
S_j^{\mathsf T}I^c_is_\alpha& j\in\kappa(i)\\
0&\text{否则}
\end{cases}\tag{6.25}
$$

### ⭐ 导致算法被发现的两个关键事实

> 原书原话：*"The **second case in Eq. 6.24** reveals the **two key facts** that led to
> the discovery of the composite-rigid-body algorithm:"*

$$
\boxed{
\begin{aligned}
&1.\ f_i\text{ 是一个}\textbf{关节运动向量}\text{与一个}\textbf{复合刚体惯性}\text{的乘积}\\
&2.\ \boldsymbol{f_j=f_i}\ \textbf{对所有 }j\in\kappa(i)\ \textbf{成立}
\end{aligned}}
$$

> 🔑 **第 2 条是伪代码那个 `while` 循环的全部理由**：
> 从关节 $i$ 一路到根，**传递的空间力是同一个** $f_i=I^c_is_\alpha$！
> 所以代码只需要把 `F` **做坐标变换**往上搬，
> **完全不需要重新计算**。
>
> **这也解释了 CRBA 为什么比"逐列跑 $\mathrm{ID}_\delta$"快**：
> 后者每列都要重新走一遍整棵树，前者一次算好 $F$ 然后沿路径投影。

**算法因此由四步组成**（原书总结）：

1. 用式 6.17 计算复合刚体惯性；
2. 在 body $i$ 坐标系中计算 $f_i=I^c_is_\alpha$；
3. 把 $f_i$ 变换到 $\kappa(i)$ 中每个刚体的坐标系；
4. 计算第 $\alpha$ 列主对角线以上的 $H$ 元素，并因对称性复制到第 $\alpha$ 行。

> 💡 **§6.2 与 §6.3 的唯一差别**：后者**逐列**处理，前者**按块**处理
> （一块 = 属于同一个关节的所有列）。若关节 $i$ 有 $n_i$ 个自由度，
> 则式 6.19 中的 ${}^{i}F_i$ 就是一个 $6\times n_i$ 矩阵，
> 各列是关节 $i$ 的每个 $\alpha$ 所对应的 $f_i=I^c_is_\alpha$。

---

## 6.4 分支诱导稀疏 (Branch-Induced Sparsity)

### 结论

由式 6.14：

$$
\boxed{\ H_{ij}=0\quad\text{当 }i\text{ 既不是 }j\text{ 的祖先、也不是后代、也不等于 }j\ }
$$

**这只可能发生在 $i$ 与 $j$ 位于树的不同分支上**，故称**分支诱导稀疏**。

### 原书图 6.2 的五个例子

| 例 | 拓扑 | 稀疏模式 |
|---|---|---|
| **(a)** | 每个刚体**直接**连到基座 | **极端情形**：每个非对角子矩阵都是零 |
| **(b)** | 无分支树 | **没有**分支诱导稀疏 |
| (c) | 典型情形 | 一般稀疏模式 |
| (d)(e) | 同一个图的**两种编号** | 得到的矩阵互为**对称置换** |

> 🔑 **一个很有用的一般结论**：
> *"More generally, $H$ will always be **block diagonal** (or a symmetric permutation
> thereof) **if the base has more than one child**."*
>
> 人形机器人若把躯干接在基座上、四肢从躯干分出，则 $H$ 不是块对角；
> 但若四条腿各自直接接到基座（如某些并联平台的建模），$H$ 就是块对角的。

> 💡 **关于编号的选择**（原书）：*"In principle, these two matrices have the
> **same underlying sparsity pattern**, and **the choice of numbering does not
> affect our ability to exploit it**."* —— 编号怎么选都不影响能否利用稀疏性。

### ⭐ 密度 $\rho$ 的经验法则

原书给了一条**很实用**的估计规则。设 $\rho$ = $H$ 中**非分支诱导零**的元素比例（$0<\rho\le1$），则

$$
\boxed{
\begin{aligned}
\text{计算 }H\text{ 的代价}&\approx\rho\times(\text{同尺寸稠密 }H\text{ 的代价})\\
\text{分解 }H\text{ 的代价}&\approx\rho^{2}\times(\text{同尺寸稠密矩阵的代价})
\end{aligned}}
$$

> *"It is **not unusual** to encounter densities of around **0.5**, and densities
> **close to zero** are possible."*

**总体效应**（原书结论）：

$$
\boxed{\ \text{惯性矩阵法在}\textbf{有分支}\text{的树上，比在}\textbf{同规模无分支}\text{的树上更高效}\ }
$$

> 💡 **这是很反直觉的一句话**：一般人会以为"分支让问题变复杂"，
> 实际上**分支让 $H$ 变稀疏，反而更快**。

### $LTL$ / $LTDL$ 与无填充性质

**CRBA 自动利用了稀疏性**（它只计算非零子矩阵）。
但若用**标准**分解算法去分解结果，它会把矩阵当稠密处理，做 $O(n^3)$ 次运算。

**解法**：分解成 $L^{\mathsf T}L$ 或 $L^{\mathsf T}DL$，并让分解算法**跳过分支诱导零**。

**在稀疏矩阵文献中的说法**（原书给出，很重要）：

| 本书记法 | 稀疏矩阵文献中的名字 |
|---|---|
| $H=L^{\mathsf T}L$ | **重排的 Cholesky 分解 (reordered Cholesky factorization)** |
| $H=L^{\mathsf T}DL$ | **重排的 $LDL^{\mathsf T}$ 分解** |

> 🔑 **推论**：$LTL$ 和 $LTDL$ 与标准 Cholesky、$LDL^{\mathsf T}$
> **具有相同的数值性质**。不必担心它们"不标准"。

**特殊性质——无填充 (no fill-in)**：

> *"the factorization proceeds **without fill-in**. In other words, **every
> branch-induced zero element in the matrix remains zero throughout the
> factorization process**."* （证明见 Featherstone (2005)。）

具有此性质的分解称为**最优的 (optimal)**（Duff et al., 1986）。**两个直接后果**：

1. $H$ 的稀疏模式**保留在因子中**，因此 $L$ 也是稀疏的
   ⟹ **回代过程也能利用稀疏性**；
2. 分解过程中**可以完全忽略**分支诱导零。

---

## 6.5 稀疏分解算法 (Sparse Factorization Algorithms)

### 输入要求

- $H$：$n\times n$，**对称正定**
- $\lambda$：$n$ 元整数数组，满足 $0\le\lambda(i)<i$
- **$H$ 的稀疏模式**：

$$
\boxed{\ \text{在 }H\text{ 的每一行 }k\text{ 上，主对角线}\textbf{以下}\text{的非零元
只出现在列 }\lambda(k),\lambda(\lambda(k)),\dots\ }
$$

> 💡 **若 $\lambda(k)=k-1$ 对所有 $k$ 成立**，算法把每个元素都当非零，
> 于是**退化成稠密的 $LTL$/$LTDL$ 分解**。
> 这是一个很好的性质：**同一份代码同时处理稀疏和稠密情形。**

### ⭐ 表 6.3：分解算法（原书）

```
────────────── LTL 分解 ──────────────      ────────────── LTDL 分解 ──────────────
for k = n to 1 do                          for k = n to 1 do
    H_kk = sqrt(H_kk)                          i = λ(k)
    i = λ(k)                                   while i ≠ 0 do
    while i ≠ 0 do                                 a = H_ki / H_kk
        H_ki = H_ki / H_kk                         j = i
        i = λ(i)                                   while j ≠ 0 do
    end                                                H_ij = H_ij − a·H_kj
    i = λ(k)                                           j = λ(j)
    while i ≠ 0 do                                 end
        j = i                                      H_ki = a
        while j ≠ 0 do                             i = λ(i)
            H_ij = H_ij − H_ki·H_kj            end
            j = λ(j)                       end
        end
        i = λ(i)
    end
end
```

**两个算法的共同点**（原书说明）：

- **就地 (in situ)** 操作给定矩阵；
- **从不访问主对角线以上的元素**；
- 外层循环从 $n$ 倒着访问每一行。任一阶段，第 $k+1..n$ 行已完成。

**三个计算区域**（原书图 6.4）：

| 区域 | 范围 | $LTL$ 做什么 | $LTDL$ 做什么 |
|---|---|---|---|
| **区域 1** | 元素 $H_{kk}$ | $H_{kk}=\sqrt{H_{kk}}$ | 什么都不做 |
| **区域 2** | 第 $k$ 行的第 $1..k-1$ 个元素 | $H_{ki}=H_{ki}/H_{kk}$ | $H'_{ki}=H_{ki}/H_{kk}$，最后写回 |
| **区域 3** | 第 $1..k-1$ 行的三角区 | $H_{ij}=H_{ij}-H_{ki}H_{kj}$ | $H_{ij}=H_{ij}-H'_{ki}H_{kj}$ |

**两者的差别**：$LTL$ 把区域 2 和 3 的计算放在**两个独立循环**里；
$LTDL$ **交错**这两个计算，从而**任何时刻只需记住一个 $H'_{ki}$ 值**（存在局部变量 `a` 中）。

**返回值**：
- $LTL$：$L$ 返回在 $H$ 的下三角
- $LTDL$：$D$ 在主对角线上，$L$ 的非对角元在其下方
  （因为 $L$ 是**单位**下三角，对角元已知为 1，不需要返回）

> 🔑 **稀疏性在哪里被利用**：内层循环**只遍历 $\lambda(k),\lambda(\lambda(k)),\dots$**。
> *"In effect, the algorithms **know where the zeros are, and simply skip over them**."*

**如何选择**：$LTDL$ **不需要开平方根**，数值上更稳、也更快（见下面的代价表）。
除非有特殊理由，**优先用 $LTDL$**。

### 表 6.4：展开的父数组

**问题**：$\lambda$ 有 $N_B$ 个元素，但分解算法需要 $n$ 个元素的数组。
$n=N_B$ 时无需处理；否则必须构造**展开的父数组 $\lambda'$**。

**做法**（原书）：把原连通图中每个**多自由度关节替换成单自由度关节链**，
得到每条弧代表一个关节变量的**展开图**，重新编号使弧 $i$ 代表 $\ddot q$ 的第 $i$ 个变量。

```
for i = 1 to n do
    λ'(i) = i − 1          ← 巧妙的初始化：一次正确设好 n − N_B 个元素
end
map(0) = 0
for i = 1 to N_B do
    map(i) = map(i − 1) + n_i
end
for i = 1 to N_B do
    λ'(map(i − 1) + 1) = map(λ(i))
end
```

**原书对这段代码的解释**：

- 第一步把 $\lambda'$ 初始化成 $(0,1,2,\dots,n-1)$。
  **这是个巧妙的技巧**：它一次就正确设好了 $n-N_B$ 个元素，只剩 $N_B$ 个要从 $\lambda$ 算。
- 第二个循环构造映射数组 $map(i)=\sum_{j=1}^{i}n_j$
  （关节 $i$ 的变量在 $\ddot q$ 中占据 $map(i-1)+1$ 到 $map(i)$）。
- 第三个循环**把两个映射合二为一**：$\lambda(i)$ 的值被转换成 $map(\lambda(i))$
  （正确的值），并插入到 $\lambda'$ 的 $map(i-1)+1$ 位置（正确的地方）。

**原书的例子**：关节 2 有 3 个自由度、其余各 1 个。
原父数组 $(0,1,1,2,2,3)$ ⟹ 展开父数组 $(0,1,2,3,1,4,4,5)$。

> 💡 **建议**：$\lambda'$ 算一次就**存进系统模型**，供后续分解算法反复使用。

### 表 6.5：乘法与回代算法

$L$ 继承了 $H$ 的稀疏模式，所以**回代也能用同样的技巧**。
原书表 6.5 给出 $Lx$、$L^{\mathsf T}x$、$L^{-1}x$、$L^{-\mathsf T}x$、$Hx$ 五个算法，例如：

```
x = L⁻¹x                        x = L⁻ᵀx
for i = 1 to n do               for i = n to 1 do
    j = λ(i)                        x_i = x_i / L_ii
    while j ≠ 0 do                  j = λ(i)
        x_i = x_i − L_ij·x_j        while j ≠ 0 do
        j = λ(j)                        x_j = x_j − L_ij·x_i
    end                                 j = λ(j)
    x_i = x_i / L_ii                end
end                             end
```

**原书的实现细节**：

- 若 $L$ 是**单位**下三角（$LTDL$ 的输出），把 $L_{ii}$ 换成 1 即可；
- $L^{-1}x$ 和 $L^{-\mathsf T}x$ 都**就地**工作；
- $Lx$ 允许 $y$ 与 $x$ 是同一向量（就地），但 $L^{\mathsf T}x$ 与 $Hx$ **要求 $y\ne x$**。

### 表 6.6：计算代价

设 $m,a,d,\sqrt{\ }$ 分别是乘、加、除、开方的代价：

| 运算 | $LTL$ | $LTDL$ |
|---|---|---|
| 分解 | $n\sqrt{\ }+D_1d+D_2(m+a)$ | $D_1d+D_2(m+a)$ |
| 回代 | $2nd+2D_1(m+a)$ | $nd+2D_1(m+a)$ |
| $Lx,\ L^{\mathsf T}x$ | $nm+D_1(m+a)$ | $D_1(m+a)$ |
| $L^{-1}x,\ L^{-\mathsf T}x$ | $nd+D_1(m+a)$ | $D_1(m+a)$ |
| $Hx$ | \multicolumn{2}{c}{$nm+2D_1(m+a)$} |

其中

$$
D_1=\sum_{k=1}^{n}(d_k-1)\tag{6.26}
\qquad\qquad
D_2=\sum_{k=1}^{n}\frac{d_k(d_k-1)}{2}\tag{6.27}
$$

$d_k=|\kappa(k)|$ = body $k$ 到基座路径上的关节数，可看作 **body $k$ 在树中的深度**。

> 💡 **$D_1$、$D_2$ 的含义**：分别是处理第 $k$ 行时区域 2、区域 3 计算的**总次数**。
> **$D_1$ 也是主对角线以下非零元的总数**，因此
> $$\boxed{\ H\text{ 中的非零元总数}=n+2D_1\ }$$
> （`code/verify_all.py` 验证了这条计数公式。）

**无分支诱导稀疏时** $D_1,D_2$ 取最大值：

$$
D_1=\frac{n^2-n}{2},\qquad D_2=\frac{n^3-n}{6}\tag{6.28}
$$

代入代价表得到的数字**与标准 Cholesky、$LDL^{\mathsf T}$ 完全相同** —— 又一次印证"重排"这个说法。

**上界**（$d$ = 树深度）：

$$
D_1\le n(d-1),\qquad D_2\le\frac{nd(d-1)}{2}\tag{6.29}
$$

$$
\boxed{\ \text{因此分解过程}\textbf{最坏}\text{是 }O(nd^2)\ }
$$

---

## 完整算例：2R 平面机械臂

**能手算到底的例子胜过十遍推导。**

### 模型

| 量 | 连杆 1 | 连杆 2 |
|---|---|---|
| 质量 | $m_1$ | $m_2$ |
| 质心距本关节 | $r_1$ | $r_2$ |
| 绕质心转动惯量（$z$） | $I_{1z}$ | $I_{2z}$ |
| 连杆长（关节 1→2） | $l_1$ | — |

$\lambda=(0,1)$，$S_1=S_2=[0\ 0\ 1\ 0\ 0\ 0]^{\mathsf T}$，${}^{2}X_1=X_{\mathrm{rot}z}(q_2)\,\mathrm{xlt}(l_1,0,0)$。

### 第 1 步：$H_{22}=S_2^{\mathsf T}I^c_2S_2$

$I^c_2=I_2$（叶节点）。$S_2=e_3$ ⟹ **取 $I_2$ 的 $(3,3)$ 元**：

$$
\boxed{H_{22}=I_{2z}+m_2r_2^{2}}
$$

即**连杆 2 绕关节 2 的转动惯量**（平行轴定理）。

### 第 2 步：$H_{21}=S_2^{\mathsf T}I^{c}_{2}\,{}^{2}X_{1}\,S_1$

（$i=2$ 后代、$j=1$ 祖先，**用后代的 $I^c_2$**。）

$$
{}^{2}X_{1}S_1=[0\ \ 0\ \ 1\ \ l_1\sin q_2\ \ l_1\cos q_2\ \ 0]^{\mathsf T}
$$

**几何检查**：角分量仍是绕 $z$ 的单位转动（转轴平行）；
线分量是长 $l_1$ 的杆绕关节 1 转动在关节 2 处产生的线速度，大小 $l_1$、方向随 $q_2$ 旋转 ✓

$I_2$ 的第 3 行是 $(0,\ 0,\ I_{2z}+m_2r_2^2,\ 0,\ m_2r_2,\ 0)$，点乘得

$$
\boxed{H_{21}=H_{12}=I_{2z}+m_2r_2^{2}+m_2l_1r_2\cos q_2}
$$

> 💡 **$\cos q_2$ 从哪来**：${}^{2}X_1S_1$ 的线分量与 $I_2$ 第 3 行里质心偏移项 $m_2r_2$ 的点乘。
> **臂伸直时 ($q_2=0$) 耦合最强，折成 90° 时耦合项消失** —— 符合物理直觉。

### 第 3 步：$H_{11}=S_1^{\mathsf T}I^{c}_{1}S_1$，$I^{c}_{1}=I_1+{}^{1}X_2^{*}I_2\,{}^{2}X_1$

$$
\boxed{H_{11}=I_{1z}+m_1r_1^{2}+I_{2z}+m_2\!\left(l_1^{2}+r_2^{2}+2l_1r_2\cos q_2\right)}
$$

**这正是任何机器人学教材里 2R 机械臂的经典结果。**
$m_2(l_1^2+r_2^2+2l_1r_2\cos q_2)$ 就是余弦定理——连杆 2 质心到关节 1 的距离平方。

### 自查点

- ✅ 对称
- ✅ $H_{22}$ 是常数（关节 2 下游没有别的关节）
- ✅ 只依赖 $q_2$ 不依赖 $q_1$（$q_1$ 只是整体旋转）
- ✅ $q_2=0$（伸直）时 $H_{11}$ 最大 —— **臂伸直时最"重"**
- ✅ 正定：$H_{11}H_{22}-H_{12}^2>0$

`code/verify_crba_2link.py` 逐项验证了以上全部结论。

---

## 6.6 补充说明与历史 (Additional Notes)

### 早期工作

| 年份 | 工作 | 特点 |
|---|---|---|
| 1965/67 | Uicker | DH 的 $4\times4$ 运动学 + $4\times4$ **伪惯性矩阵** + Lagrange 方法 |
| 1965 | Hooker & Margulies | 向量式（牛顿-欧拉）+ **增广刚体**，用于自由漂浮的卫星树 |
| 1975 | Paul | 早期工作的**优秀综述** |
| 1977 | Orlandea et al. | 程序 **ADAMS** 的**稀疏矩阵表述**——"一个特别成功的早期想法" |
| 1977 | Wittenburg | **最早的计算刚体动力学教科书** |

### CRBA 的谱系

- **首次出现**：**Walker & Orin (1982) 的方法 3**
  （方法 1、2 就是 §6.1 的 $\mathrm{ID}_\delta$ 逐列法）
- 更高效的版本与变体：Featherstone (1984, 1987)、Balafoutis & Patel (1989, 1991)、
  Lilly & Orin (1991)、Lilly (1993)、McMillan & Orin (1998)
- ⚠️ **但分支诱导稀疏这一现象直到 Featherstone (2005) 才被认识到**

> 💡 **这是个有意思的历史细节**：CRBA 1982 年就有了，
> 但"它在分支树上会自动变快"这件事，**过了 23 年才被发现**。

### ⭐ 两类方法的深层联系（原书最有价值的一段）

> *"Superficially, it would seem that the inertia-matrix methods and propagation
> methods are **very different**. However, **two interesting connections** have been found."*

**联系 1（Rodriguez et al., 1987/1991）**：借鉴 **Kalman 滤波理论**，
给出关节空间惯性矩阵的**两种分解**：
- 一种**蕴含 RNEA**；
- 另一种的**逆蕴含 ABA**。

**联系 2（Ascher et al., 1997）**：从一个形如式 3.18 的方程出发，证明

$$
\boxed{
\begin{aligned}
&\text{关节空间惯性矩阵}\ \leftarrow\ \text{对系数矩阵的}\textbf{一种置换}\text{做高斯消元}\\
&\text{铰接体算法 (ABA)}\ \leftarrow\ \text{对}\textbf{另一种置换}\text{做高斯消元}
\end{aligned}}
$$

> 🔑 **这正式证实了"ABA 本质上是树上的高斯消元"这个解读**（见第 7 章）。
> CRBA 路线和 ABA 不是两种不同的思想，**而是同一个消元过程的两种排序**。

### 第三种策略

> 严格说，惯性矩阵法和传播法只是**三种**可能策略中的两种。
> 第三种是：把各刚体的运动方程与约束方程**组装成一个大矩阵方程**（如第 3 章 §3.7）。
> **这个策略对闭环系统非常有用，但用在运动学树上竞争不过前两种。**
> 例子见 Baraff (1996)。

---

## 什么时候用 CRBA 而不是 ABA

| 场景 | 推荐 | 理由 |
|---|---|---|
| $n$ 较小（典型 6/7 轴机械臂） | **CRBA** | 常数因子小 |
| 分支多且浅（人形、四足、灵巧手） | **CRBA** | $\rho$ 小，$O(nd^2)$ 极快 |
| $n$ 很大的**长链** | **ABA** | $O(n)$ 优势显现，且分支稀疏性用不上 |
| **需要 $H$ 本身** | **CRBA** | ABA 根本不构造 $H$ |
| 有接触约束、闭环 | **CRBA** | 原书：$O(nd^2)$ 算法是闭环算法的重要组成部分 |

**什么时候你会需要 $H$ 本身**：操作空间控制（$\Lambda=(JH^{-1}J^{\mathsf T})^{-1}$）、
阻抗控制、零空间投影、接触力求解（KKT 矩阵里要 $H$）、
模型分析（特征值、条件数）、最优控制/MPC。

> 🔑 **这是 CRBA 在实践中比 ABA 更常见的真正原因——很多应用要的不只是 $\ddot q$。**

---

## 易错点与陷阱

1. **$H_{ij}$ 的下标条件记反**（最隐蔽）。
   式 6.14 要求 **$i\in\nu(j)$（$j$ 是祖先）时用 $I^c_i$（后代的）**。
   记反会得到**对称、量纲正确、数值看似合理但错误**的 $H$。
   记忆锚点：**无分支树的式 6.16 写着 $I^c_{\max(i,j)}$**。
2. **$I^c$ 的累加顺序**：必须内推；先累加到父节点，用的是**已完整**的 $I^c_i$。
3. **$F$ 的变换方向**：$F$ 是**力向量**，用 $X^{*}$。
4. **只填了下三角忘了对称**（是否需要取决于 $H$ 怎么用，见 §6.2 的说明）。
5. **用了标准稠密 Cholesky**：会产生填充，稀疏优势全失。
6. **忘了 $C$**：$H\ddot q=\tau$ 是错的，必须 $\tau-C$。
7. **$H$ 数值不对称**导致 Cholesky 失败：只算下三角再镜像，或强制 $\tfrac12(H+H^{\mathsf T})$。
8. **多自由度关节时忘了展开父数组**（表 6.4）。

## 调试建议

- **CRBA vs $n$ 次 $\mathrm{ID}_\delta$**（表 6.1）：应逐元素相等。
- **CRBA vs $\sum_kJ_k^{\mathsf T}I_kJ_k$**：检验变换与 $S$ 的用法。
- **CRBA vs $\partial^2T/\partial\dot q_i\partial\dot q_j$**（数值二阶差分）：
  这一路**完全不碰 $I^c$、$X^{*}$ 和下标条件**，
  是**唯一能可靠检出"下标写反"**的对拍方式。
- **$LTL$ 检验**：$L^{\mathsf T}L$ 应还原 $H$；且 $L$ 的稀疏模式应与 $H$ 相同（无填充）。
- **非零元计数**：$H$ 的非零元数应等于 $n+2D_1$（式 6.26）。
- **与 ABA 对拍**：两条 FD 路线的 $\ddot q$ 应一致。

以上全部实现在 `code/verify_all.py` 与 `code/verify_crba_2link.py` 中。

## 与其他章的联系

- ← 第 2 章：惯性可加性（式 6.12/6.13 的依据）、$X^{*}IX$
- ← 第 3 章：$H$ 的对称正定；式 3.63 $H=S^{\mathsf T}IS$ 是式 6.14 的最小情形
- ← 第 4 章：$\kappa,\nu,\mu$、恒等式 4.3（式 6.7→6.8 的重排依据）、式 4.6
- ← 第 5 章：$C=\mathrm{ID}(q,\dot q,0)$；表 6.1 是 RNEA 的简化版
- ↔ 第 7 章：另一条 FD 路线；§6.6 揭示二者是同一消元的两种排序
- → 第 8 章：$O(nd^2)$ 算法是闭环算法的组成部分
- → 第 10 章：§10.3.2 有详细的代价分析

---

## ✍️ 我的理解

<!-- 建议：说清 §6.3 的「两个关键事实」，尤其 f_j = f_i 对所有 j∈κ(i) -->

## ❓ 疑问与待办

- [ ] 实现表 6.1 的 $\mathrm{ID}_\delta$，逐列构造 $H$，与 CRBA 对拍
- [ ] 照着 §"完整算例"手算一遍 2R 的 $H_{22}\to H_{21}\to H_{11}$
- [ ] 自己推一遍式 6.7 → 6.8 的重排（用第 4 章恒等式 4.3）
- [ ] 实现表 6.3 的 $LTDL$，验证无填充性质
- [ ] 用第 10 章 §10.3.2 补上精确代价
- [ ] 读 Ascher et al. (1997)，看清 CRBA 与 ABA 是同一消元的两种置换

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.101–118 逐节核对。
     此前版本对 H_ij 下标条件的注释写反，已按式 6.14 修正并确认。 -->
