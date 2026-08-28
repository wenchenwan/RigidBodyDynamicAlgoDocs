# 推导补充

> 书上跳过、或者自己想通后要记下来的推导。
> **建议**：每读完一章，挑一两个关键结论自己重推一遍写在这里。
> 能推出来才算真懂。

---

## 目录

- [D1. $X^{*}=X^{-\mathsf T}$](#d1)
- [D2. 空间惯性矩阵的两行](#d2)
- [D3. 单刚体运动方程](#d3)
- [D4. 空间加速度与经典加速度的关系](#d4)〔待补〕
- [D5. ABA 的递推关系](#d5)〔待补〕
- [D6. $H$ 的稀疏性](#d6)〔待补〕

---

<a id="d1"></a>
## D1. $X^{*}=X^{-\mathsf T}$（对偶变换）

**要证**：力向量的坐标变换矩阵是运动向量变换矩阵的**逆的转置**。

**证明**：功率 $P=\mathbf{f}^{\mathsf T}\mathbf{v}$ 是标量物理量，与坐标系无关。

设 ${}^{B}\mathbf{v}=X\,{}^{A}\mathbf{v}$，${}^{B}\mathbf{f}=X^{*}\,{}^{A}\mathbf{f}$，则

$$
{}^{B}\mathbf{f}^{\mathsf T}\,{}^{B}\mathbf{v}
=\left(X^{*}{}^{A}\mathbf{f}\right)^{\mathsf T}\left(X\,{}^{A}\mathbf{v}\right)
={}^{A}\mathbf{f}^{\mathsf T}\left(X^{*\mathsf T}X\right){}^{A}\mathbf{v}
$$

要对**任意** $\mathbf{f},\mathbf{v}$ 都等于 ${}^{A}\mathbf{f}^{\mathsf T}\,{}^{A}\mathbf{v}$，
必须 $X^{*\mathsf T}X=\mathbf{1}_6$，即

$$
X^{*}=X^{-\mathsf T}\qquad\blacksquare
$$

**同理可得** $\mathbf{v}\times^{*}=-(\mathbf{v}\times)^{\mathsf T}$：
对 $\frac{d}{dt}(\mathbf{f}^{\mathsf T}\mathbf{v})$ 用两个求导公式展开，
两个叉乘项必须相消。

**要点**：这两条恒等式都不是定义，是**能量守恒逼出来的必然结果**。

---

<a id="d2"></a>
## D2. 空间惯性矩阵的两行

**要证**：

$$
I_O=\begin{bmatrix}\bar I+m\,c\times c\times^{\mathsf T} & m\,c\times\\ m\,c\times^{\mathsf T} & m\,\mathbf{1}_3\end{bmatrix}
$$

**推导**：由 $\mathbf{h}=I\mathbf{v}$，$\mathbf{v}=[\omega;\ v_O]$，$\mathbf{h}=[h_\omega;\ h_v]$。

**下面一行（线动量）**：质心速度 $v_C=v_O+\omega\times c$，

$$
h_v=m\,v_C=m(v_O+\omega\times c)=m\,v_O-m\,c\times\omega
$$

（用了 $\omega\times c=-c\times\omega$。）
写成矩阵行：$\begin{bmatrix}-m\,c\times & m\,\mathbf{1}\end{bmatrix}
=\begin{bmatrix}m\,c\times^{\mathsf T} & m\,\mathbf{1}\end{bmatrix}$ ✓

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

**要证**：$\mathbf{f}=I\mathbf{a}+\mathbf{v}\times^{*}I\mathbf{v}$

**推导**：牛顿第二定律的空间形式是 $\mathbf{f}=\dfrac{d\mathbf{h}}{dt}$，$\mathbf{h}=I\mathbf{v}$。

在**随体坐标系**中，$I$ 是常量，所以表观导数

$$
\frac{\mathring d\mathbf{h}}{dt}=I\,\frac{\mathring d\mathbf{v}}{dt}=I\mathbf{a}
$$

$\mathbf{h}$ 是**力向量**，所以用 $\times^{*}$ 版本的求导公式：

$$
\mathbf{f}=\frac{d\mathbf{h}}{dt}=\frac{\mathring d\mathbf{h}}{dt}+\mathbf{v}\times^{*}\mathbf{h}
=I\mathbf{a}+\mathbf{v}\times^{*}I\mathbf{v}\qquad\blacksquare
$$

**另一条路**：$\mathbf{f}=\frac{d}{dt}(I\mathbf{v})=I\mathbf{a}+\dot I\mathbf{v}$，
配合 $\dot I=\mathbf{v}\times^{*}I-I\,\mathbf{v}\times$ 和 $\mathbf{v}\times\mathbf{v}=0$，
得到同样结果。

**练习**：把这条式子按 $3\times3$ 分块展开，验证它就是牛顿方程
$f=m\,a_C$ 加欧拉方程 $n_C=\bar I\dot\omega+\omega\times\bar I\omega$。
〔待做〕

---

<a id="d4"></a>
## D4. 空间加速度与经典加速度〔待补〕

**目标**：证明 $a_P^{\text{classical}}=\dot v_O+\omega\times v_O$，
并说明为什么加速度合成里没有系数 2。

<!-- 读第 2 章时补 -->

---

<a id="d5"></a>
## D5. ABA 的递推关系〔待补〕

**目标**：推出

$$
I^{a}=I^{A}-UD^{-1}U^{\mathsf T},\qquad
\mathbf{p}^{a}=\mathbf{p}^{A}+I^{a}\mathbf{c}+UD^{-1}u
$$

**思路提示**：从子节点的 $\mathbf{f}=I^{A}\mathbf{a}+\mathbf{p}^{A}$ 出发，
代入 $\mathbf{a}=\mathbf{a}_{\lambda}+\mathbf{c}+S\ddot q$，
再用 $\tau=S^{\mathsf T}\mathbf{f}$ 解出 $\ddot q$ 并回代消元。
最后会看到 Schur 补的结构。

<!-- 读第 7 章时补 -->

---

<a id="d6"></a>
## D6. $H$ 的稀疏性〔待补〕

**目标**：证明 $H_{ij}\ne0\iff i\in\nu(j)$ 或 $j\in\nu(i)$。

**思路提示**：$H_{ij}=\partial^2 T/\partial\dot q_i\partial\dot q_j$，
而 $T=\sum_k\tfrac12\mathbf{v}_k^{\mathsf T}I_k\mathbf{v}_k$，
$\mathbf{v}_k$ 只依赖 $\kappa(k)$ 上的关节速度。

<!-- 读第 6 章时补 -->
