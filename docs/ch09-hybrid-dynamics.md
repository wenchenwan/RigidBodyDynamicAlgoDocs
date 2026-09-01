# 第 9 章 混合动力学与其他专题 (Hybrid Dynamics and Other Topics)

> **原书 pp. 171–194**，共 7 节。

> **一句话概括**：四个专题——**混合动力学**（部分关节已知力、部分已知加速度）、
> **浮动基**、**齿轮**、**动力学等价**。

## 原书节次对照

| 节 | 标题 | 页 |
|---|---|---|
| 9.1 | Hybrid Dynamics | 171 |
| 9.2 | Articulated-Body Hybrid Dynamics | 176 |
| 9.3 | Floating Bases | 179 |
| 9.4 | Floating-Base Forward Dynamics | 181 |
| 9.5 | Floating-Base Inverse Dynamics | 183 |
| 9.6 | Gears | 186 |
| 9.7 | **Dynamic Equivalence** | 189 |

> ⚠️ **重要更正**：我此前的笔记把"冲量动力学与碰撞"放在本章——**这是错的**。
> 原书把接触与碰撞放在**独立的第 11 章**。本章只有上面 7 节。

---

## 9.1 混合动力学 (Hybrid Dynamics)

### 问题定义

每个关节 $i$ 有一个力变量 $\tau_i$ 和一个加速度变量 $\ddot q_i$：

| 问题 | 给定 | 求 |
|---|---|---|
| **正动力学 FD** | **每个**关节的 $\tau_i$ | 全部加速度 |
| **逆动力学 ID** | **每个**关节的 $\ddot q_i$ | 全部力 |
| **混合动力学** | 每个关节给 $\tau_i$ **或** $\ddot q_i$ **之一** | 未知的加速度与力 |

> **混合动力学 = 在一部分关节上做正动力学、在其余关节上做逆动力学。
> 它把 FD 和 ID 都作为特例包含在内。**

### ⭐ 三类用途（原书讲得很实在）

**用途 1：引入规定运动 (prescribed motions)**

- 研究机构或传动装置在**输入端规定运动**条件下的动力学
- 复摆对给定**运动学激励**的响应
- 动画角色在空中的**翻滚行为**

**用途 2：为欠驱动系统设计轨迹**

> 例：为动画角色设计一串肢体动作，使其**起跳、翻筋斗、双脚落地**。

**用途 3：⭐ 简化仿真——去掉不重要的/高频的动力学**

> 原书这段非常有价值：
>
> 机器人关节的运动通常由**高增益反馈控制系统**控制。仿真整个机器人时，
> 仿真器必须把驱动器和控制系统的动力学也算进去。
> **但如果某个控制系统可以视为"完美"**（受控关节的实际运动≈指令运动），
> 那就可以**丢掉这个控制系统和它的驱动器，直接规定受影响关节的运动**。
>
> **若被去掉的部件正是原系统中最高频动力学的来源，
> 那么仿真器现在可以用更长的积分步长。**
> （例见 Hu et al. 2005。）

> 🔑 **这是混合动力学最被低估的用途**：它不是"另一种动力学问题"，
> 而是一个**降低仿真刚性、加大步长**的工程手段。

### 方程

定义 $fd$ = **正动力学关节**集合（力已知、加速度未知）。
$Q$ 是把 FD 变量排到前面的置换矩阵：

$$
\begin{bmatrix}H_{11}&H_{12}\\ H_{21}&H_{22}\end{bmatrix}
\begin{bmatrix}\ddot q_1\\ \ddot q_2\end{bmatrix}
=\begin{bmatrix}\tau_1\\ \tau_2\end{bmatrix}-\begin{bmatrix}C_1\\ C_2\end{bmatrix}\tag{9.1}
$$

未知量是 $\ddot q_1$ 和 $\tau_2$。把它们并到左边：

$$
\boxed{\ \begin{bmatrix}H_{11}&0\\ H_{21}&-\mathbf 1\end{bmatrix}
\begin{bmatrix}\ddot q_1\\ \tau_2\end{bmatrix}
=\begin{bmatrix}\tau_1\\ 0\end{bmatrix}-\begin{bmatrix}C'_1\\ C'_2\end{bmatrix}\ }\tag{9.2}
$$

$$
\begin{bmatrix}C'_1\\ C'_2\end{bmatrix}
=\begin{bmatrix}C_1+H_{12}\ddot q_2\\ C_2+H_{22}\ddot q_2\end{bmatrix}=QC'
$$

### ⭐ $C'$ 的物理含义与算法

> *"Physically, $C'$ is **the force required to impart zero acceleration to each
> forward-dynamics joint and the given acceleration to each inverse-dynamics joint**."*

因此可以**用一次逆动力学调用算出**：

$$
\boxed{\ C'=\mathrm{ID}\!\left(q,\dot q,\ Q^{\mathsf T}\begin{bmatrix}0\\ \ddot q_2\end{bmatrix}\right)\ }\tag{9.3}
$$

> 🔑 **这一步很漂亮**：一次 RNEA 就同时处理了"FD 关节置零"和"ID 关节用给定值"。

### 四步算法

$$
\boxed{
\begin{aligned}
&1.\ \text{用式 9.3 算 }C'\\
&2.\ \text{算 }H_{11}\\
&3.\ \text{解 }H_{11}\ddot q_1=\tau_1-C'_1\\
&4.\ \text{用 }\tau=C'+\mathrm{ID}_\delta\!\left(Q^{\mathsf T}[\ddot q_1;\ 0]\right)\text{ 算 }\tau_2
\end{aligned}}
$$

（步骤 4 中的 $\mathrm{ID}_\delta$ 就是第 6 章表 6.1 那个乘 $H$ 的简化算法。）

### 步骤 2 的两种做法

**简单做法**：用 CRBA 算出整个 $H$，取子矩阵 $H_{11}$。

**高效做法（原书表 9.1）**：修改 CRBA，**只做得到 $H_{11}$ 所必需的最少计算**。
需要一个新量：

$$
\nu(fd)=\bigcup_{i\in fd}\nu(i)\tag{9.5}
$$

即**被至少一个正动力学关节支撑的刚体集合**。它可以由 $fd$ 和 $\lambda$ 算出：

```
nu(fd) = fd
for i = 1 to N_B do
    if lam(i) in nu(fd) then  nu(fd) = nu(fd) union {i}
end
```

**修改后的 CRBA 与原版结构相同，只是做了子集的计算**：
$I^c_i$ 只对 $i\in\nu(fd)$ 计算，$H_{ij}$ 只对 $i,j\in fd$ 计算。

> 💡 **原书的工程建议**：$\nu(fd)$ 的计算代价很小，
> 每次调用混合动力学例程时重算也无妨；
> 但因为它只依赖 $fd$ 和 $\lambda$，**一次仿真过程中通常不变**，
> 所以也可以在开始时算一次存起来。
> **实现上可以用 $N_B+1$ 个布尔值的数组表示。**

### ⭐ 步骤 3：$H_{11}$ 也有分支诱导稀疏

> *"$H_{11}$ **inherits a portion of whatever branch-induced sparsity is present in $H$**."*

因此可以用第 6 章 §6.5 的稀疏分解。前提是拿到描述 $H_{11}$ 稀疏模式的**修改父数组 $\lambda'$**。

**怎么构造 $\lambda'$**（原书图 9.1，思路极清晰）：

$$
\boxed{\ \text{从原连通图出发，把每一对被}\textbf{逆动力学关节}\text{连接的节点}\textbf{合并}\ }
$$

合并后节点变少，需要重新编号。**原书的例子**：关节 1、2 是 ID 关节
⟹ 节点 1、2 与节点 0 合并，节点 3–7 重编号为 1–5。

```
map(0) = 0 ;  j = 1
for i = 1 to n do
    if i in fd then
        map(i) = j
        lam'(j)  = map(lam(i))
        j = j + 1
    else
        map(i) = map(lam(i))
    end
end
```

> ⚠️ 若修改图中有多自由度关节，$\lambda'$ **还需要按第 6 章 §6.5 的表 6.4 展开**。

> 💡 **合并的直觉**：ID 关节的加速度是给定的，
> 所以它连接的两个刚体在"未知量"的意义上是**刚性连在一起**的——
> 合并节点正是这件事的图论表达。

---

## 9.2 铰接体混合动力学

> **ABA 很容易改造成混合动力学**：在**正动力学关节**上用第 7 章的常规铰接体方程，
> 在**逆动力学关节**上用下面的方程。

### 推导

$$
f_i=I^A_ia_i+p^A_i
\qquad
f_i=I^a_ia_{\lambda(i)}+p^a_i
\qquad
a_i=a_{\lambda(i)}+\dot S_i\dot q_i+S_i\ddot q_i\tag{9.6, 9.7, 9.8}
$$

> **ABA 的主要步骤是由 $I^A_i,p^A_i$ 算出 $I^a_i,p^a_i$。
> 若关节 $i$ 是 FD 关节，直接的困难是 $\ddot q_i$ 未知，所以第一步要先求出 $\ddot q_i$。
> 但对 ID 关节，这个困难根本不存在。**

把 9.8 代入 9.6 并与 9.7 比较，立刻得到

$$
\boxed{\ I^a_i=I^A_i\ }
\qquad\qquad
\boxed{\ p^a_i=p^A_i+I^A_i(\dot S_i\dot q_i+S_i\ddot q_i)\ }\tag{9.9, 9.10}
$$

> 🔑 **对比 FD 关节的式 7.47、7.48**：
> - FD：$I^a=I^A-UD^{-1}U^{\mathsf T}$（**Schur 补，消元**）
> - ID：$I^a=I^A$（**惯性原样传递，不消元**）
>
> **统一视角**：$fd$ 类关节做**消元**（未知量），
> 非 $fd$ 类做**代入**（已知量，效应并入偏置力）。
> 这正是解线性方程组时"消元 vs 代入"的选择，只是在树上做。
>
> **RNEA 是全代入的极端，ABA 是全消元的极端，混合动力学是两者的插值。**

### 表 9.2 / 9.3 的另外两处改动

除了式 9.9、9.10，ID 关节的方程还有两处改动：

1. **$S_i\ddot q_i$ 被加进了 $c_i$**（趟 1），因此**趟 3 里不再出现 $a'_i$**；
2. **趟 3 新增一个方程**算 $\tau_i=S_i^{\mathsf T}f_i=S_i^{\mathsf T}(I^A_ia_i+p^A_i)$。

**伪代码结构**（表 9.3，由表 7.1 改造）：

```
趟 1: 对每个 i
        v_i = X[i,lam] v_lam + vJ
        if i in fd: c_i = cJ + v_i x vJ
        else:       c_i = cJ + v_i x vJ + S_i qdd_i     <- 差别在这
        I^A_i = I_i ;  p^A_i = v_i x* I_i v_i - X[i,0]^* f^x_i

趟 2: for i = N_B to 1
        if i in fd:                                     <- 消元
            U_i, D_i, u_i ;  I^a = I^A_i - U D^-1 U^T ; p^a = p^A_i + I^a c_i + U D^-1 u
        else:                                           <- 代入
            I^a = I^A_i ;                               p^a = p^A_i + I^a c_i
        累加 I^a, p^a 到父节点

趟 3: for i = 1 to N_B
        if i in fd:  a' = X[i,lam] a_lam + c_i ; qdd_i = D^-1 (u_i - U^T a') ; a_i = a' + S_i qdd_i
        else:        a_i = X[i,lam] a_lam + c_i ; tau_i = S_i^T (I^A_i a_i + p^A_i)
```

---

## 9.3 浮动基 (Floating Bases)

### 转换成固定基系统（三步）

$$
\boxed{
\begin{aligned}
&1.\ \text{所有刚体与关节编号}\textbf{加 1}\text{（浮动基成为 body 1，最小关节号为 2）}\\
&2.\ \text{加一个固定基座（新 body 0）和一个连接二者的}\textbf{6-DoF 关节}\text{（新 joint 1）}\\
&3.\ N_B,N_J\ \textbf{各加 1}，n\ \textbf{加 6}
\end{aligned}}
$$

做完这三步，**固定基的所有技术和算法都可以原样使用**。
但浮动基系统是重要的特例，**值得专门处理**。

### 关键选择：$S_1=\mathbf 1_{6\times6}$

关节 1 有 6 个自由度，$S_1$ 必须是 $6\times6$ 满秩矩阵。
**任意满秩矩阵都可以**——两种选择只是两种速度变量的选择，不是本质变化。
所以取**最方便的**：**在浮动基坐标系中取单位阵**。

于是

$$
v_1=\dot q_1,\qquad a_1=\ddot q_1,\qquad \tau_1=f_1
$$

> 即 **$\dot q_1,\ddot q_1,\tau_1$ 就是空间向量 $v_1,a_1,f_1$ 的 Plücker 坐标**。
>
> ⚠️ 原书脚注：${}^{1}S_1=\mathbf 1_{6\times6}$，但**在任何其他坐标系中
> ${}^{i}S_1={}^{i}X_1\mathbf 1_{6\times6}$ 就不是单位阵了**。
> 又因 $S_1$ 固连于 body 1，有 $\mathring S_1=0$、$\dot S_1=v_1\times S_1$。

**$f_1$ 与 $f^x_1$ 含义本质相同**（都是作用在浮动基上的外力），
总外力实际是二者之和。**我们总有自由决定怎么在 $f^x_1$ 与 $f_1$ 之间分配。**

### 浮动基运动方程

由 $H_{1i}=S_1^{\mathsf T}I^c_iS_i=I^c_iS_i$ 和 $H_{11}=I^c_1$，且 $C_1=p^c_1$
（**整个浮动基系统作为复合刚体的空间偏置力**），得到

$$
\boxed{\ \begin{bmatrix}I^c_0&F\\ F^{\mathsf T}&H\end{bmatrix}
\begin{bmatrix}a_0\\ \ddot q\end{bmatrix}
+\begin{bmatrix}p^c_0\\ C\end{bmatrix}
=\begin{bmatrix}0\\ \tau\end{bmatrix}\ }\tag{9.13}
$$

$$
F=\begin{bmatrix}F_1&F_2&\cdots&F_{N_B}\end{bmatrix},\qquad F_i=I^c_iS_i\tag{9.12}
$$

> 🔑 **$F$ 的物理含义**（原书）：$6\times n$ 矩阵，各列是
> **"支持每个关节变量上的单位加速度所需的、作用在浮动基上的空间力"**。
> 与第 6 章 §6.3 的 $f_i=I^c_is_\alpha$ 完全一致。

**这是 $n+6$ 个方程、$n+6$ 个未知量**，对应浮动基系统的 $n+6$ 个自由度。

**完整状态描述**：

| 描述 | 需要 |
|---|---|
| 加速度 | $a_0$ **和** $\ddot q$ |
| 速度 | $v_0$ **和** $\dot q$ |
| 位置 | ${}^{0}X_{ref}$ **和** $q$ |

### ⚠️ 例 9.2：消去 $a_0$ 会毁掉稀疏性

对式 9.13 消去 $a_0$（第一行乘 $F^{\mathsf T}(I^c_0)^{-1}$ 减去第二行）：

$$
H^{fl}\ddot q+C^{fl}=\tau\tag{9.14}
$$

$$
H^{fl}=H-F^{\mathsf T}(I^c_0)^{-1}F
\qquad
C^{fl}=C-F^{\mathsf T}(I^c_0)^{-1}p^c_0\tag{9.15, 9.16}
$$

> ⚠️ **原书的警告**：
> *"this equation is **not necessarily a good choice for computational purposes**
> because the subtraction of $F^{\mathsf T}(I^c_0)^{-1}F$ from $H$ in Eq. 9.15
> **erases any branch-induced sparsity** that may have been present in $H$.
> Thus, $H^{fl}$ is a **dense** matrix; and **it is perfectly possible for $H^{fl}$
> to have more nonzero elements than $H$, despite being smaller**."*

> 🔑 **这条很重要**：$H^{fl}$ 形式漂亮（直接给出 $\ddot q$ 与 $\tau$ 的关系），
> 但**计算上是个陷阱**——人形机器人的 $H$ 本来 70% 是零，消去 $a_0$ 后全填满。
> **实践中应该直接解式 9.13 的鞍点系统，而不是构造 $H^{fl}$。**

---

## 9.4 浮动基正动力学

> 若已有支持 6-DoF 关节的固定基软件，**直接就能算浮动基动力学**。
> 若需要专用的浮动基软件，**ABA 和 CRBA 都能改造**。

**浮动基 ABA（表 9.4）与表 7.1 的四处差别**：

1. **$v_0$ 现在是输入**（不再恒为零）
2. **趟 2 被扩展**，要算出 $I^A_0$ 和 $p^A_0$
3. $a_0$ 由下式解出：
   $$I^A_0a_0+p^A_0=0\tag{9.17}$$
4. **重力加速度在最后才加到 $a_0$ 上**

> ⚠️ **一个容易忽略的实现细节**：这个实现要求 $f^x_i$ 和 $a_g$
> 在**浮动基坐标系**中给出（而非参考坐标系），
> 这样就不需要知道 ${}^{0}X_{ref}$。因此记作 ${}^{0}f^x_i$、${}^{0}a_g$。
>
> **注意 ${}^{0}a_g$ 不是常量**——它随浮动基转动而变化：
> $${}^{0}a_g={}^{0}X_{ref}\,{}^{ref}a_g$$

> 🔑 **式 9.17 的物理含义**：浮动基上没有外部约束，
> 所以"整个铰接体（= 整个系统）受到的净力为零"⟹ $I^A_0a_0+p^A_0=0$。
> 这就是自由漂浮系统加速度的来源。

---

## 9.5 浮动基逆动力学

### 推导

引入**相对加速度** $a^r_i$（相对浮动基的加速度），则

$$
f_i=\sum_{j\in\nu(i)}I_ja_0+\sum_{j\in\nu(i)}\big(I_ja^r_j+v_j\times^{*}I_jv_j-f^x_j\big)\tag{9.19}
$$

$a_0$ 的系数正是复合刚体惯性 $I^c_i$，其余当作偏置力：

$$
f_i=I^c_ia_0+p^c_i
\qquad
p^c_i=\sum_{j\in\nu(i)}\big(I_ja^r_j+v_j\times^{*}I_jv_j-f^x_j\big)\tag{9.20, 9.21}
$$

> **$p^c_i$ 的含义**：*"the force that would be required to support the motion of
> all the bodies in $\nu(i)$ **if $a_0$ happened to be zero**."*

递推形式：

$$
p^c_i=p_i+\sum_{j\in\mu(i)}p^c_j
\qquad
p_i=I_ia^r_i+v_i\times^{*}I_iv_i-f^x_i\tag{9.22, 9.23}
$$

**$i=0$ 时式 9.20 可以直接解出 $a_0$**：

$$
\boxed{\ a_0=-(I^c_0)^{-1}p^c_0\ }\tag{9.24}
$$

然后

$$
\tau_i=S_i^{\mathsf T}f_i=S_i^{\mathsf T}(I^c_ia_0+p^c_i)\tag{9.25}
$$

**三趟算法**（表 9.6）：趟 1 外推算 $v_i,a^r_i,p_i$ 并初始化 $I^c_i$；
趟 2 内推累加 $I^c$、$p^c$；趟 3 解出 $a_0$ 再外推算 $\tau_i$。

**重力处理**：把 $a^r_0$ 初始化为 $-{}^{0}a_g$ 而非零。
这使重力项出现在 $p^c_i$ 中，进而使正确的重力加速度出现在 $a_0$ 中。
**或者**：趟 1 开头把 $a^r_0$ 置零，趟 3 末尾把 $a_g$ 加到 $a_0$ 上。

### ⭐ 例 9.3：动量守恒的简洁形式

系统总空间动量 $=\sum_{i=0}^{N_B}I_iv_i$，**无外力时守恒**。原书证明它也等于

$$
\boxed{\ \text{momentum}=I^c_0v_0+F\dot q\ }
$$

**证明**（用第 4 章式 4.3 的求和换序）：

$$
\sum_{i=0}^{N_B}I_iv_i=\sum_iI_i\Big(v_0+\sum_{j\in\kappa(i)}S_j\dot q_j\Big)
=\sum_iI_iv_0+\sum_{j=1}^{N_B}\sum_{i\in\nu(j)}I_iS_j\dot q_j
=I^c_0v_0+\sum_jI^c_jS_j\dot q_j=I^c_0v_0+F\dot q
$$

> 🔑 **这个式子是浮动基控制的基石**：
> $I^c_0$ 是整个系统的**复合惯性**（关于浮动基），$F\dot q$ 是关节运动的贡献。
> 现代人形机器人控制里的 **centroidal momentum matrix** 就是它的变体。
> 而它的推导只用了第 4 章那条不起眼的恒等式 4.3。

---

## 9.6 齿轮 (Gears)

> ⚠️ **重要更正**：我此前的笔记说"最简处理是在 $H$ 的对角元上加 $r^2I_{rotor}$"。
> **原书不是这么做的**——它把转子建模成**独立刚体**，
> 并用 **§8.11 的闭环函数技术**处理齿轮约束。

### 齿轮约束的本质

> **从数学角度看，一对齿轮就是一个涉及两个关节变量的代数方程。**
> 典型情况下一个变量是另一个的常数倍，**但非线性关系也是可能的**。
> 还可能涉及**两个以上**变量——例如差速器通常有两个输出、一或两个输入，
> 因而总共涉及**三或四个**变量。

**处理方法**：**与 §8.11 处理运动学回路完全相同**——
选一组独立关节变量 $y$，定义把 $y$ 映到 $q$ 的函数 $\gamma$。

### 原书的算例（图 9.2）：齿轮驱动的两连杆机械臂

**系统组成**：固定基座、大臂、小臂、两个电机、两个大齿轮、两个小齿轮。
每个电机含**定子**和**转子**。定子分别固定在基座和大臂上；转子各自固连一个小齿轮；
两个大齿轮分别固连在大臂和小臂上。

**四个运动刚体**：

| 刚体 | 组成 |
|---|---|
| body 1 | 大臂 + 第一个大齿轮 + 第二个电机的**定子** |
| body 2 | 第一个小齿轮 + 第一个电机的**转子** |
| body 3 | 小臂 + 第二个大齿轮 |
| body 4 | 第二个小齿轮 + 第二个电机的**转子** |

**关节 1、3** = 肩关节和肘关节；**关节 2、4** = 电机内部让转子相对定子转动的轴承。

**齿轮约束**：

$$
q_2=\rho_1q_1,\qquad q_4=\rho_2q_3
$$

**独立变量**取 $y=[q_1\ q_3]^{\mathsf T}$，则

$$
\gamma(y)=\begin{bmatrix}y_1\\ \rho_1y_1\\ y_2\\ \rho_2y_2\end{bmatrix}
\qquad
G=\frac{\partial\gamma}{\partial y}=\begin{bmatrix}1&0\\ \rho_1&0\\ 0&1\\ 0&\rho_2\end{bmatrix}
\qquad
g=\dot G\dot y=0
$$

> 💡 **$g$ 只在系统中含非线性齿轮对时才非零。**

然后用第 3 章式 3.21 / 第 8 章式 8.45 的投影法即可：$G^{\mathsf T}HG\,\ddot y=G^{\mathsf T}(\tau-C-Hg)$。

> 🔑 **这种建模方式比"加 $\rho^2I_{rotor}$ 到对角元"更完整**：
> - 它自动包含了**转子的陀螺效应**（转子的完整空间惯性被建模了）；
> - 它能处理**非线性齿轮**和**多输入/输出**（差速器）；
> - "加对角元"实际上是这个模型在"转子惯量各向同性、忽略陀螺效应"
>   假设下的**退化结果**。
>
> **实践建议**：多数应用中"加对角元"够用且快得多；
> 但高速转子（如反作用飞轮、CMG）**必须**用完整建模。

---

## 9.7 动力学等价 (Dynamic Equivalence)

> **两个不同的刚体系统完全可能有相同的运动方程**，称为**动力学等价**。
> 重要特例：两个系统只在**惯性参数**上不同。

### 两个后果

**后果 1（负面）**：

$$
\boxed{\ \text{刚体系统的惯性参数}\textbf{无法}\text{从其动力学行为的测量中辨识出来}\ }
$$

只能辨识一个**可观测子集**，称为**基惯性参数 (base inertia parameters)**
（Khalil & Dombre, 2002）。

> 💡 **这解释了机器人参数辨识的一个基本困难**：
> 你可以完美拟合所有实验数据，却仍然得不到真实的 $m$、$c$、$\bar I$——
> 因为有无穷多组参数给出同一个运动方程。

**后果 2（正面，本节主题）**：

$$
\boxed{\ \text{可以}\textbf{调整}\text{系统模型的惯性参数而不改变其行为}\ }
$$

> **若调整后有更多参数变成零，则计算动力学的代价可以降低。**

### 条件

考虑两刚体系统（图 9.3），运动方程是式 9.13 在两体情形的应用：

$$
\begin{bmatrix}I_1+X_J^{\mathsf T}I_2X_J&X_J^{\mathsf T}I_2S\\ S^{\mathsf T}I_2X_J&S^{\mathsf T}I_2S\end{bmatrix}
\begin{bmatrix}a_1\\ \ddot q\end{bmatrix}
+\begin{bmatrix}C_f\\ C_\tau\end{bmatrix}
=\begin{bmatrix}f\\ \tau\end{bmatrix}\tag{9.29}
$$

**把 $I_\Delta$ 加到 $I_1$、从 $I_2$ 减掉**。要使式 9.29 的系数不变，$I_\Delta$ 必须满足

$$
\boxed{\ I_\Delta S\equiv0,\qquad
I_\Delta\equiv X_J^{\mathsf T}I_\Delta X_J,\qquad
S^{\mathsf T}v_1\times^{*}I_\Delta v_1\equiv0\ }\tag{9.30}
$$

**这些条件必须对每个关节位形、每个 $v_1\in M^6$ 都成立**，
因此**只依赖关节类型**——不需要知道系统的任何其他信息就能找出合适的 $I_\Delta$。

原书给出了转动、移动等关节类型的 $I_\Delta$ 表。例如**转动关节**：

$$
I_\Delta=\begin{bmatrix}
I&0&0&0&-mc_z&0\\ 0&I&0&mc_z&0&0\\ 0&0&0&0&0&0\\
0&mc_z&0&m&0&0\\ -mc_z&0&0&0&m&0\\ 0&0&0&0&0&m
\end{bmatrix}
$$

**移动关节**：$I_\Delta=\begin{bmatrix}\bar I&0\\ 0&0\end{bmatrix}$（任意对称 $\bar I$）。

> ⚠️ **原书脚注的一个前提**：推出这组条件需要假设关节的运动子空间
> $\mathrm{range}(S)$ 在 $B_1$ 或 $B_2$ 坐标系之一中是常量。
> **这个假设对所有关节类型并不都成立。**

> 🔑 **实用价值**：给定一个 URDF 模型，可以在**不改变任何动力学行为**的前提下
> 把某些惯性参数"搬"到相邻连杆上，使更多参数变成零，
> 从而让**符号化代码生成**（§10.4）产生更短的代码。

---

## 本章要点回顾

1. **混合动力学 = 部分关节 FD、部分 ID**，FD 与 ID 都是它的特例。
2. **最被低估的用途**：把完美控制的关节改成规定运动，
   **去掉高频动力学从而加大积分步长**。
3. **$C'$ 是"给 FD 关节零加速度、给 ID 关节指定加速度所需的力"**，一次 ID 调用即得。
4. **$H_{11}$ 继承部分分支稀疏性**；$\lambda'$ 由"合并被 ID 关节连接的节点"得到。
5. **ABA 混合版**：FD 关节**消元**（Schur 补），ID 关节**代入**（$I^a=I^A$）。
6. **浮动基取 $S_1=\mathbf 1$** ⟹ $\dot q_1,\ddot q_1,\tau_1$ 就是 $v_1,a_1,f_1$。
7. **式 9.13 的 $F_i=I^c_iS_i$** = "支持关节 $i$ 单位加速度所需的、作用在浮动基上的力"。
8. ⚠️ **消去 $a_0$ 得到的 $H^{fl}$ 是稠密的**，可能比 $H$ 的非零元还多——**不要这么算**。
9. **浮动基 ABA**：$v_0$ 是输入；$a_0$ 由 $I^A_0a_0+p^A_0=0$ 解出；${}^{0}a_g$ 不是常量。
10. **动量 $=I^c_0v_0+F\dot q$**（例 9.3），只用第 4 章恒等式 4.3 就能证。
11. **齿轮用 §8.11 的闭环函数技术处理**，转子作为独立刚体建模。
12. **惯性参数无法从动力学行为辨识**（只能辨识基惯性参数）；
    反过来可以调整参数以降低计算代价。

---

## 易错点

1. **把接触/碰撞当成本章内容**——那是**第 11 章**。
2. **构造 $H^{fl}$ 来做浮动基仿真**——毁掉稀疏性，见例 9.2。
3. **浮动基的四元数不归一化**——会造成整个子树的尺度误差（第 4 章 §4.5）。
4. **忘了 ${}^{0}a_g$ 随基座姿态变化**（不是常量）。
5. **混合动力学里同一关节既给 $\tau$ 又给 $\ddot q$**（过约束）或都不给（欠定）。
6. **齿轮只加对角项而忽略陀螺效应**——高速转子场合会出错。

## 与其他章的联系

- ← 第 4 章：浮动基的建模（§4.1）、恒等式 4.3（例 9.3 用到）
- ← 第 5 章：$C'=\mathrm{ID}(\dots)$
- ← 第 6 章：表 6.1 的 $\mathrm{ID}_\delta$、$\lambda'$ 的展开、稀疏分解
- ← 第 7 章：式 7.25 是 §9.2 的起点
- ← 第 8 章：**§9.6 齿轮直接复用 §8.11 的闭环函数技术**
- → 第 10 章：§9.7 的参数调整服务于 §10.4 的符号化简
- → 第 11 章：接触与碰撞

---

## ✍️ 我的理解

<!-- 建议：说清「消元 vs 代入」这个统一视角 -->

## ❓ 疑问与待办

- [ ] 实现混合动力学 ABA（表 9.3），用"全 fd"退化成 ABA、"全非 fd"退化成 RNEA 来验证
- [ ] 实现表 9.1 的修改版 CRBA，与"算完整 $H$ 再取子块"对拍
- [ ] 实现浮动基 ABA，验证无外力时 $I^c_0v_0+F\dot q$ 守恒
- [ ] 构造一个人形模型，比较 $H$ 与 $H^{fl}$ 的非零元个数（验证例 9.2 的警告）
- [ ] 用图 9.2 的齿轮机械臂算例，对比"完整建模"与"加 $\rho^2I_{rotor}$"的差异
- [ ] 用式 9.30 找一组 $I_\Delta$，验证调整后运动方程确实不变

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.171–194 逐节核对。
     此前版本的错误已修正：(1) 把接触/碰撞误放在本章（实为第 11 章）；
     (2) 齿轮的处理方式（原书用闭环函数技术并把转子建模成独立刚体，
     而非简单在对角元加 ρ²I_rotor）；(3) 完全遗漏 §9.7 动力学等价。 -->
