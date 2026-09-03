# 推导补充

> 书上跳过、或者自己想通后要记下来的推导。
> **建议**：每读完一章，挑一两个关键结论自己重推一遍写在这里。
> 能推出来才算真懂。

**记号**：D1–D4 属第 2 章，空间向量与同名 3D 向量并存，按原书戴帽子
（$\hat v=[\omega;\ v_O]$、$\hat f=[n_O;\ f]$）；D5、D6 属第 6、7 章，按原书去掉帽子（$v_i,\ f_i,\ p^A_i$）。

---

## 目录

- [D1. $X^{*}=X^{-\mathsf T}$](#d1)
- [D2. 空间惯性矩阵的两行](#d2)
- [D3. 单刚体运动方程](#d3)
- [D4. 空间加速度与经典加速度的关系](#d4)
- [D5. ABA 的递推关系](#d5)
- [D6. $H_{ij}$ 的闭式表达与稀疏性](#d6)

---

<a id="d1"></a>
## D1. $X^{*}=X^{-\mathsf T}$（对偶变换）

**要证**：力向量的坐标变换矩阵是运动向量变换矩阵的**逆的转置**（原书 §2.1 对偶坐标、§2.8）。

**证明**：功率 $P=\hat f^{\mathsf T}\hat v$ 是标量物理量，与坐标系无关。

设 ${}^{B}\hat v=X\,{}^{A}\hat v$，${}^{B}\hat f=X^{*}\,{}^{A}\hat f$，则

$$
{}^{B}\hat f^{\mathsf T}\,{}^{B}\hat v
=\left(X^{*}\,{}^{A}\hat f\right)^{\mathsf T}\left(X\,{}^{A}\hat v\right)
={}^{A}\hat f^{\mathsf T}\left(X^{*\mathsf T}X\right){}^{A}\hat v
$$

要对**任意** $\hat f,\hat v$ 都等于 ${}^{A}\hat f^{\mathsf T}\,{}^{A}\hat v$，
必须 $X^{*\mathsf T}X=\mathbf 1_6$，即

$$
X^{*}=X^{-\mathsf T}\qquad\blacksquare
$$

**同理可得** $\hat v\times^{*}=-(\hat v\times)^{\mathsf T}$（原书 §2.9）：
若 $\hat m\in M^6$、$\hat f\in F^6$ 都固连在速度为 $\hat v$ 的刚体上，标量 $\hat f^{\mathsf T}\hat m$ 不随时间变化，
用式 2.29、2.30 求导：

$$
0=\frac{d}{dt}\big(\hat f^{\mathsf T}\hat m\big)
=(\hat v\times^{*}\hat f)^{\mathsf T}\hat m+\hat f^{\mathsf T}(\hat v\times\hat m)
=\hat f^{\mathsf T}\big[(\hat v\times^{*})^{\mathsf T}+\hat v\times\big]\hat m
$$

对任意 $\hat f,\hat m$ 成立 ⟹ $(\hat v\times^{*})^{\mathsf T}=-\hat v\times$。$\blacksquare$

**要点**：这两条恒等式都不是定义，是**能量守恒逼出来的必然结果**。
`code/verify_all.py` 第 2 章一节对它们做了数值核对。

---

<a id="d2"></a>
## D2. 空间惯性矩阵的两行

**要证**：

$$
I_O=\begin{bmatrix}\bar I+m\,c\times c\times^{\mathsf T} & m\,c\times\\ m\,c\times^{\mathsf T} & m\,\mathbf{1}_3\end{bmatrix}
$$

**推导**：由 $\hat h=I\hat v$，$\hat v=[\omega;\ v_O]$，$\hat h=[h_\omega;\ h_v]$。

**下面一行（线动量）**：质心速度 $v_C=v_O+\omega\times c$，

$$
h_v=m\,v_C=m(v_O+\omega\times c)=m\,v_O-m\,c\times\omega
$$

（用了 $\omega\times c=-c\times\omega$。）
写成矩阵行：$\begin{bmatrix}-m\,c\times & m\,\mathbf{1}\end{bmatrix}=\begin{bmatrix}m\,c\times^{\mathsf T} & m\,\mathbf{1}\end{bmatrix}$ ✓

**上面一行（关于 $O$ 的角动量）**：绕质心的角动量 $\bar I\omega$，
加上质心处线动量对 $O$ 的矩：

$$
h_\omega=\bar I\omega+c\times(m\,v_C)
=\bar I\omega+m\,c\times(v_O+\omega\times c)
$$

$$
=\bar I\omega+m\,c\times v_O-m\,c\times c\times\omega
=\left(\bar I-m\,c\times c\times\right)\omega+m\,c\times v_O
$$

再用 $c\times^{\mathsf T}=-c\times$，即 $-m\,c\times c\times=m\,c\times c\times^{\mathsf T}$：

$$
h_\omega=\left(\bar I+m\,c\times c\times^{\mathsf T}\right)\omega+m\,c\times v_O
\qquad\blacksquare
$$

**注**：上左块 $\bar I-m\,c\times c\times = \bar I+m\left((c\cdot c)\mathbf{1}-cc^{\mathsf T}\right)$
正是**平行轴定理**（把绕质心的惯量搬到 $O$ 点）。

---

<a id="d3"></a>
## D3. 单刚体运动方程

**要证**（原书式 2.68）：$\hat f=I\hat a+\hat v\times^{*}I\hat v$

**推导**：牛顿第二定律的空间形式是 $\hat f=\dfrac{d\hat h}{dt}$，$\hat h=I\hat v$。

在**随体坐标系**中，$I$ 是常量，所以表观导数

$$
\frac{\mathring d\hat h}{dt}=I\,\frac{\mathring d\hat v}{dt}=I\hat a
$$

$\hat h$ 是**力向量**，所以用 $\times^{*}$ 版本的求导公式（§2.10）：

$$
\hat f=\frac{d\hat h}{dt}=\frac{\mathring d\hat h}{dt}+\hat v\times^{*}\hat h
=I\hat a+\hat v\times^{*}I\hat v\qquad\blacksquare
$$

**另一条路**：$\hat f=\frac{d}{dt}(I\hat v)=I\hat a+\dot I\hat v$，
配合 $\dot I=\hat v\times^{*}I-I\,\hat v\times$（式 2.65）和 $\hat v\times\hat v=0$，
得到同样结果。

**按 $3\times3$ 分块展开**：在质心系（$c=0$）展开，恰好得到牛顿方程 $f=m\,\ddot c$（式 2.70）
加欧拉方程 $n_C=\bar I_C\dot\omega+\omega\times\bar I_C\omega$（式 2.71），见 ch02 例 2.6；
`code/verify_ch02.py` 逐项验证。关键一步是用 D4 的式 2.48 把空间加速度的线分量换成
经典加速度 $\ddot c$——离心项正是在这一步被吸收的。

---

<a id="d4"></a>
## D4. 空间加速度与经典加速度

**要证**（原书式 2.47、2.48）：设 $O'$ 是此刻恰好经过 $O$ 的**物质点**，$r(t)$ 是它的位置，则

$$
\dot v_O=\ddot r-\omega\times\dot r
\qquad\qquad
\hat a_O=\begin{bmatrix}\dot\omega\\ \ddot r-\omega\times\dot r\end{bmatrix}\tag{2.47, 2.48}
$$

即 **经典加速度的线分量 $=$ 空间加速度的线分量 $+\ \omega\times v_O$**；
并解释为什么空间加速度的合成公式里既没有离心项、也没有科氏项的系数 2。

### 1. 式 2.47 的推导

$v_O$ 的定义是"此刻位于 $O$ 的那个物质点的速度"。刚体的速度场是
$v_P=v_O+\omega\times\overrightarrow{OP}$（原书 §2.2），所以物质点 $O'$ 在**任何**时刻都满足

$$
\dot r=v_O+\omega\times r
$$

（$r$ 从固定点 $O$ 量起；这条式子对所有 $t$ 成立，不只是 $r=0$ 的那一刻。）
两边对时间求导：

$$
\ddot r=\dot v_O+\dot\omega\times r+\omega\times\dot r
$$

取 $O'$ 恰好经过 $O$ 的那一刻，$r=0$：

$$
\ddot r=\dot v_O+\omega\times\dot r
\quad\Longrightarrow\quad
\dot v_O=\ddot r-\omega\times\dot r\qquad\blacksquare
$$

又因为此刻 $\dot r=v_O$，所以 $\ddot r=\dot v_O+\omega\times v_O$：
**物质点的加速度 = 速度场在 $O$ 点的时间变化率 + $\omega\times v_O$**。
第一项是空间加速度的线分量（欧拉描述：站在固定点 $O$ 看"流过来的速度"怎么变），
第二项是物质点顺着弯曲的流线走时额外"感受到"的加速度（拉格朗日描述）。
这正是原书用式 2.50 说的：经典加速度 $\hat a'_O$ 是把坐标原点**跟着 $O'$ 平移**时看到的表观导数。

### 2. 用"悖论"检验

刚体绕过 $O$ 的固定轴以 $\omega$ 匀速转动（原书图 2.6）。
$v_O=0$、$\dot\omega=0$ ⟹ $\hat a_O=0$。但对刚体上任一点 $P$（$\overrightarrow{OP}=r_P$）：

$$
\ddot r_P=\omega\times(\omega\times r_P)=-\omega^2r_P\ne0
$$

换到 $P$ 处看更清楚：速度场在固定点 $P$ 的取值 $v_P=\omega\times r_P$ 是**常量**，
所以空间加速度在 $P$ 处的线分量 $\dot v_P=0$；而由第 1 节的结论，
经过 $P$ 的物质点的加速度 $=\dot v_P+\omega\times v_P=\omega\times(\omega\times r_P)=-\omega^2r_P$ ✓。
`code/verify_ch02.py` §2.11 用 $\omega=2.1$、$r_P=(0.5,0,0)$ 算出了这两个数：$0$ 与 $-2.205$。

### 3. 为什么加速度合成里没有离心项和系数 2

经典的相对运动公式（点 $P$ 在以 $\omega$ 转动的坐标系里以相对速度 $\dot\rho$ 运动）有五项：

$$
\ddot r_P=\ddot r_{O}+\dot\omega\times\rho
+\underbrace{\omega\times(\omega\times\rho)}_{\text{离心}}
+\underbrace{2\,\omega\times\dot\rho}_{\text{科氏}}
+\ddot\rho
$$

系数 2 来自两处各贡献一个 $\omega\times\dot\rho$：
一处是对牵连项 $\omega\times\rho$ 求导，另一处是把相对导数换成绝对导数。

空间向量的合成只有一个交叉项。由 $\hat v_i=\hat v_{i-1}+\hat s_i\dot q_i$（原书式 2.16）直接求导：

$$
\hat a_i=\hat a_{i-1}+\dot{\hat s}_i\dot q_i+\hat s_i\ddot q_i,
\qquad
\dot{\hat s}_i=\hat v_i\times\hat s_i
$$

（$\hat s_i$ 固连于 body $i$，故用式 2.29；它同时固连于 body $i-1$，
而 $\hat v_{i-1}\times\hat s_i=(\hat v_i-\hat s_i\dot q_i)\times\hat s_i=\hat v_i\times\hat s_i$，两种写法一致。）于是

$$
\hat a_i=\hat a_{i-1}+\hat s_i\ddot q_i+\hat v_i\times\hat s_i\dot q_i\tag{2.55}
$$

同一件事的一般形式：$\hat v_2=\hat v_1+\hat v_{\text{rel}}$，其中 $\hat v_{\text{rel}}$ 在随 body 1 运动的坐标系里表达，
用 §2.10 的求导公式 $\dot{\hat v}_{\text{rel}}=\mathring{\hat v}_{\text{rel}}+\hat v_1\times\hat v_{\text{rel}}$，得
$\hat a_2=\hat a_1+\hat a_{\text{rel}}+\hat v_1\times\hat v_{\text{rel}}$。

**那些项去哪儿了？** 没有消失，而是被两次"吸收"：

1. **定义层面**：式 2.48 把 $-\omega\times\dot r$ 吃进了空间加速度的线分量，
   所以 $\hat a$ 本身不含向心加速度——这就是第 2 节里 $\hat a_O=0$ 而各点都在加速的原因；
2. **动力学层面**：运动方程 $\hat f=I\hat a+\hat v\times^{*}I\hat v$ 的第二项把它们还回来。
   例 2.6（见 D3）在质心系展开时，$\hat a$ 线分量里的 $-\omega\times v_C$
   与 $\hat v\times^{*}I\hat v$ 里的 $m\,\omega\times v_C$ 恰好相消，重新得到 $f=m\ddot c$；
   而 $\omega\times\bar I_C\omega$ 就是欧拉方程里的陀螺项。

**所以**：空间加速度把"离心 / 科氏"这些**依赖参考点的项**从运动学里剥离出来，
集中放进一个与参考点无关的项 $\hat v\times^{*}I\hat v$ 里。运动学递推因此只剩一个交叉项，
这就是 RNEA 只有几行的原因（式 2.56 的双重求和显示了不用递推时这些项会有 $O(i^2)$ 个）。

---

<a id="d5"></a>
## D5. ABA 的递推关系

**要证**（原书表 7.1 趟 2 的更新式，即式 7.19、7.20 / 7.23、7.24 的 $U,D,u$ 写法）：
以 body $i$ 为柄的铰接体，其惯性 $I^A_i$、偏置力 $p^A_i$ 在消掉关节加速度 $\ddot q_i$ 之后，对父节点表现为

$$
U_i=I^A_iS_i,\qquad D_i=S_i^{\mathsf T}U_i,\qquad u_i=\tau_i-S_i^{\mathsf T}p^A_i
\tag{7.43, 7.44, 7.45}
$$

$$
I^a_i=I^A_i-U_iD_i^{-1}U_i^{\mathsf T}
\qquad\qquad
p^a_i=p^A_i+I^a_ic_i+U_iD_i^{-1}u_i
$$

并且 $I^A_{\lambda(i)}\leftarrow I^A_{\lambda(i)}+I^a_i$、$p^A_{\lambda(i)}\leftarrow p^A_{\lambda(i)}+p^a_i$（坐标变换略）。

**设定**（原书 §7.2.2 的特例：浮动运动学树，每个铰接体恰有一个柄）。
把所有量都写在 body $i$ 的坐标系里，省略 $X$：

- 铰接体的定义：通过关节 $i$ 传给子树 $i$ 的力满足 $f_i=I^A_ia_i+p^A_i$（式 7.15、7.16 的形式）；
- 关节运动学：$a_i=a_{\lambda(i)}+c_i+S_i\ddot q_i$（式 7.31），记 $a'_i=a_{\lambda(i)}+c_i$（式 7.46）；
- 关节力平衡：$\tau_i=S_i^{\mathsf T}f_i$（约束力被 $S_i^{\mathsf T}$ 消掉，第 3 章式 3.37）。

### 1. 解出 $\ddot q_i$

$$
\tau_i=S_i^{\mathsf T}\big(I^A_i(a'_i+S_i\ddot q_i)+p^A_i\big)
=U_i^{\mathsf T}a'_i+D_i\ddot q_i+S_i^{\mathsf T}p^A_i
$$

$$
\Longrightarrow\quad
\ddot q_i=D_i^{-1}\big(u_i-U_i^{\mathsf T}a'_i\big)\tag{7.30}
$$

这就是趟 3 的回代公式。$D_i=S_i^{\mathsf T}I^A_iS_i$ 正定（$I^A_i$ 正定、$S_i$ 列满秩），
所以**永远可逆、无需选主元**——与第 3 章式 3.51 的 $S^{\mathsf T}IS$ 是同一件事。

### 2. 回代，消掉 $\ddot q_i$

$$
f_i=I^A_ia'_i+U_i\ddot q_i+p^A_i
=I^A_ia'_i+U_iD_i^{-1}u_i-U_iD_i^{-1}U_i^{\mathsf T}a'_i+p^A_i
=\underbrace{\big(I^A_i-U_iD_i^{-1}U_i^{\mathsf T}\big)}_{I^a_i}a'_i+p^A_i+U_iD_i^{-1}u_i
$$

再把 $a'_i=a_{\lambda(i)}+c_i$ 代回：

$$
f_i=I^a_i\,a_{\lambda(i)}+\underbrace{\big(p^A_i+I^a_ic_i+U_iD_i^{-1}u_i\big)}_{p^a_i}\qquad\blacksquare
$$

⚠️ **$p^a_i$ 里乘 $c_i$ 的是 $I^a_i$ 而不是 $I^A_i$**——因为 $c_i$ 是 $a'_i$ 的一部分，
它经历了与 $a_{\lambda(i)}$ 相同的消元。原书式 7.24 写成
$p^a_i=p^A_i+I^a_ic_i+I^A_iS_i(S_i^{\mathsf T}I^A_iS_i)^{-1}(\tau_i-S_i^{\mathsf T}p^A_i)$，
最后一项就是 $U_iD_i^{-1}u_i$。若把 $I^a_ic_i$ 展开成 $I^A_ic_i-U_iD_i^{-1}U_i^{\mathsf T}c_i$，
也可以写成 $p^a_i=p^A_i+I^A_ic_i+U_iD_i^{-1}(u_i-U_i^{\mathsf T}c_i)$。

### 3. 并到父节点：装配

$f_i$ 现在是 $a_{\lambda(i)}$ 的线性函数。父体 $\lambda$ 自己的方程是 $I_\lambda a_\lambda+p_\lambda$
（$p_\lambda=v_\lambda\times^{*}I_\lambda v_\lambda-f^x_\lambda$，即 $p^A$ 的初值），
它还要通过各子关节把 $f_j$ 传给每棵子树 $j\in\mu(\lambda)$，所以经关节 $\lambda$ 传入的总力为

$$
f_\lambda=I_\lambda a_\lambda+p_\lambda+\sum_{j\in\mu(\lambda)}\big(I^a_ja_\lambda+p^a_j\big)
=\Big(I_\lambda+\sum_j I^a_j\Big)a_\lambda+\Big(p_\lambda+\sum_j p^a_j\Big)
$$

这就是式 7.21、7.22：$I^A_\lambda=I_\lambda+\sum_{j}{}^{\lambda}X^{*}_jI^a_j\,{}^{j}X_\lambda$，
$p^A_\lambda=p_\lambda+\sum_j{}^{\lambda}X^{*}_jp^a_j$。
由于 $I^a_j$、$p^a_j$ 只依赖子树 $j$ 内部的量，趟 2 从叶到根一遍就能算完——**这是 $O(n)$ 的来源**。

### 4. 认出 Schur 补

第 1、2 步其实是在对下面这个 $(6+n_i)$ 阶分块方程做一步高斯消元：

$$
\begin{bmatrix}I^A_i&U_i\\ U_i^{\mathsf T}&D_i\end{bmatrix}
\begin{bmatrix}a'_i\\ \ddot q_i\end{bmatrix}
=\begin{bmatrix}f_i-p^A_i\\ \tau_i-S_i^{\mathsf T}p^A_i\end{bmatrix}
$$

（第一行是 $f_i=I^A_i(a'_i+S_i\ddot q_i)+p^A_i$，第二行是它左乘 $S_i^{\mathsf T}$。）
用第二行消去 $\ddot q_i$，第一行的系数就变成 Schur 补 $I^A_i-U_iD_i^{-1}U_i^{\mathsf T}=I^a_i$。
ABA 的趟 2 = 沿树从叶到根依次消元，趟 3 = 回代——这就是第 7 章 §7.4
"ABA 是运动方程的高斯消元"这一解释的单步版本。

### 5. 两条一致性检查

$$
I^a_iS_i=U_i-U_iD_i^{-1}D_i=0
\qquad\qquad
S_i^{\mathsf T}p^a_i=S_i^{\mathsf T}p^A_i+0+D_iD_i^{-1}u_i=\tau_i
$$

第一条说：消元后的铰接体沿关节自由方向**没有惯性**（$I^a_i$ 奇异，秩 $\le 6-n_i$）——
父体沿这些方向推它不费力，力全被关节"放掉"了；
第二条说：无论 $a_{\lambda(i)}$ 是什么，$S_i^{\mathsf T}f_i=S_i^{\mathsf T}(I^a_ia_{\lambda(i)}+p^a_i)=\tau_i$ 都自动成立，
关节力平衡已经内置在 $p^a_i$ 里。

**数值验证**：`code/verify_all.py` 第 7 章一节（ABA ↔ RNEA 互验、ABA vs $H^{-1}(\tau-C)$、
CRBA+LTL 回代 == ABA）。

---

<a id="d6"></a>
## D6. $H_{ij}$ 的闭式表达与稀疏性

**要证**：

$$
H_{ij}=\sum_{k\in\nu(i)\cap\nu(j)}\left({}^{k}X_{i}S_i\right)^{\mathsf T}I_k\left({}^{k}X_{j}S_j\right)
$$

由此立刻得到稀疏性，以及 $j$ 为 $i$ 的祖先时的 $H_{ij}=S_i^{\mathsf T}I^{c}_{i}\,{}^{i}X_{j}S_j$。

**推导**：见 [`docs/ch06-forward-dynamics-crba.md` §二](../docs/ch06-forward-dynamics-crba.md)，
那里有完整的五步版本。核心是三点：

1. $T=\sum_k\tfrac12 v_k^{\mathsf T}I_k v_k$，
   $v_k=\sum_{l\in\kappa(k)}{}^{k}X_{l}S_l\dot q_l$ 关于 $\dot q$ 线性；
2. 于是 $\partial v_k/\partial\dot q_i={}^{k}X_{i}S_i$（$i\in\kappa(k)$ 时）否则为零，
   而 **$i\in\kappa(k)\iff k\in\nu(i)$**；
3. 二阶偏导要求两个因子同时非零 ⟹ 求和范围是 $\nu(i)\cap\nu(j)$。

**稀疏性**：树上 $\nu(i)\cap\nu(j)\ne\varnothing$ 当且仅当 $i,j$ 在同一条根到叶的路径上，
即有祖先-后代关系。否则 $H_{ij}=0$。$\blacksquare$

**$I^c$ 为什么取后代的**：$j$ 是祖先 ⟹ $\nu(i)\subseteq\nu(j)$ ⟹ 交集 $=\nu(i)$。
**交集取小的，所以 $I^c$ 取后代的。**

**数值验证**：[`code/verify_crba_2link.py`](../code/verify_crba_2link.py)
