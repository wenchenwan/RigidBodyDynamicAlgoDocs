# 空间向量代数公式速查

> 一页纸版本。打印出来贴在桌上的那种。**已按原书全文核对**，附原书公式编号。
> 详细解释见 [`docs/ch02-spatial-vector-algebra.md`](../docs/ch02-spatial-vector-algebra.md)。

## 1. 基本对象

$$
\hat v=\begin{bmatrix}\omega\\ v_O\end{bmatrix}\in M^6,
\qquad
\hat f=\begin{bmatrix}n_O\\ f\end{bmatrix}\in F^6
$$

⚠️ 注意排列：速度是 **(角, 线)**，力是 **(矩, 力)**。

**功率（唯一的天然配对）**：

$$
P=\hat f^{\mathsf T}\hat v=n_O\cdot\omega+f\cdot v_O
$$

**动能**：$\;T=\tfrac12\,\hat v^{\mathsf T}I\,\hat v\;$，
系统层面 $\;T=\tfrac12\dot q^{\mathsf T}H\dot q$

## 2. 坐标变换

设 $E={}^{B}R_{A}$（旋转），$r$ = $B$ 原点在 $A$ 系中的位置。

$$
{}^{B}X_{A}=\begin{bmatrix}E & 0\\ -E\,r\times & E\end{bmatrix}
\qquad
{}^{B}X_{A}^{*}=\begin{bmatrix}E & -E\,r\times\\ 0 & E\end{bmatrix}
$$

$$
\boxed{X^{*}=X^{-\mathsf T}}\qquad
{}^{B}X_{A}^{*}=\left({}^{A}X_{B}\right)^{\mathsf T}\qquad
\left({}^{B}X_{A}\right)^{-1}={}^{A}X_{B}
$$

$$
{}^{C}X_{A}={}^{C}X_{B}\ {}^{B}X_{A}\qquad (XY)^{*}=X^{*}Y^{*}
$$

**惯性的变换（三明治）**：

$$
{}^{B}I={}^{B}X_{A}^{*}\ {}^{A}I\ {}^{A}X_{B}=\left({}^{A}X_{B}\right)^{\mathsf T}{}^{A}I\ {}^{A}X_{B}
$$

## 3. 叉乘

$$
\hat v\times=\begin{bmatrix}\omega\times & 0\\ v_O\times & \omega\times\end{bmatrix}
\qquad
\hat v\times^{*}=\begin{bmatrix}\omega\times & v_O\times\\ 0 & \omega\times\end{bmatrix}
$$

$$
\boxed{\hat v\times^{*}=-\left(\hat v\times\right)^{\mathsf T}}
\qquad
\hat v\times\hat v=0
\qquad
X(\hat v\times)X^{-1}=(X\hat v)\times
$$

## 4. 求导（坐标系以 $\hat v$ 运动）

$$
\frac{d\hat m}{dt}=\frac{\mathring{d}\hat m}{dt}+\hat v\times\hat m
\quad(\hat m\in M^6)
\qquad
\frac{d\hat g}{dt}=\frac{\mathring{d}\hat g}{dt}+\hat v\times^{*}\hat g
\quad(\hat g\in F^6)
$$

（$\mathring d$ = 表观导数，把基当常量）

## 5. 速度与加速度的合成

$$
\hat v_2=\hat v_1+\hat v_{\text{rel}}
\qquad
\boxed{\hat a_2=\hat a_1+\hat a_{\text{rel}}+\hat v_1\times\hat v_{\text{rel}}}
$$

⚠️ **没有系数 2**，没有单独的离心项和科氏项——它们被空间加速度的定义吸收了。

**空间加速度 ↔ 经典加速度**：

$$
a_{P}^{\text{classical}}=\underbrace{\dot v_O}_{\hat a\text{ 的线分量}}+\ \omega\times v_O
$$

## 6. 空间惯性

$$
I_O=\begin{bmatrix}
\bar I+m\,c\times c\times^{\mathsf T} & m\,c\times\\[2pt]
m\,c\times^{\mathsf T} & m\,\mathbf{1}_3
\end{bmatrix}
$$

（$c\times^{\mathsf T}=-c\times$；上左块 $=\bar I-m\,c\times c\times$）

**性质**：对称、正定、**可加**（同坐标系下 $I=I_1+I_2$）、只有 10 个独立参数。

**动量**：$\;\hat h=I\hat v$

## 7. 单刚体运动方程

$$
\boxed{\ \hat f=I\,\hat a+\hat v\times^{*}I\,\hat v\ }
$$

等价：$\;\dot I=\hat v\times^{*}I-I\,\hat v\times\;$，
$\;\hat f=\frac{d}{dt}(I\hat v)=I\hat a+\dot I\hat v$

## 8. 系统运动方程

$$
\boxed{\ H(q)\,\ddot q+C(q,\dot q)=\tau\ }
$$

- $H$：对称、正定、只依赖 $q$、**分支诱导稀疏**（$H_{ij}\ne0\iff$ 祖先-后代关系）
- $C$ 是**向量**（含科氏+离心+重力），$C=\mathrm{ID}(q,\dot q,\mathbf 0)$

**关节约束**：

$$
v_{\text{rel}}=S_i\dot q_i,
\qquad S_i^{\mathsf T}T_i=0,
\qquad \boxed{\tau_i=S_i^{\mathsf T}f_i}
$$

**闭环**：$H\ddot q+C=\tau+K^{\mathsf T}\lambda$，$K\ddot q=k$

**冲量**：$H\,\Delta\dot q=\iota$（$q$ 不变，$C$ 消失）

## 9. 常见关节的 $S$（在关节坐标系中，绕/沿 $z$ 轴）

| 关节 | $S$ |
|---|---|
| revolute | $[0\ 0\ 1\ 0\ 0\ 0]^{\mathsf T}$ |
| prismatic | $[0\ 0\ 0\ 0\ 0\ 1]^{\mathsf T}$ |
| helical (螺距 $h$) | $[0\ 0\ 1\ 0\ 0\ h]^{\mathsf T}$ |
| free / floating | $\mathbf{1}_6$ |

## 10. 受约束刚体（第 3 章 §3.6）

$$
\Phi=S(S^{\mathsf T}IS)^{-1}S^{\mathsf T}\ (3.54)
\qquad
H=S^{\mathsf T}IS\ (3.63)
\qquad
C=S^{\mathsf T}(I\dot S\dot q+p)\ (3.64)
$$

$\mathrm{range}(\Phi)=\mathcal S$，$\mathrm{null}(\Phi)=\mathcal S^{\perp}$，$\mathrm{rank}(\Phi)=n_f$。

**约束的两种描述**（式 3.11）：

|  | 位置 | 速度 | 加速度 |
|---|---|---|---|
| 隐式 | $\phi(q)=0$ | $K\dot q=0$ | $K\ddot q=k$ |
| 显式 | $q=\gamma(y)$ | $\dot q=G\dot y$ | $\ddot q=G\ddot y+g$ |

$KG=0$，$Kg=k$（式 3.12）；$\tau_c=K^{\mathsf T}\lambda$（3.15）；$G^{\mathsf T}\tau_c=0$（3.16）。

**KKT**（式 3.17）：$\begin{bmatrix}H&K^{\mathsf T}\\ K&0\end{bmatrix}\begin{bmatrix}\ddot q\\ -\lambda\end{bmatrix}=\begin{bmatrix}\tau-C\\ k\end{bmatrix}$
（**对称但不正定**）

**投影法**（式 3.20/3.21）：$H_G\ddot y+C_G=u$，$H_G=G^{\mathsf T}HG$，$C_G=G^{\mathsf T}(C+Hg)$，$u=G^{\mathsf T}\tau$

## 11. 三个惯性（第 6、7 章）

| | 子树状态 | 关系式 | 参数 |
|---|---|---|---|
| $I$ | 只有自己 | $f=I a+v\times^{*}Iv$ | **10** |
| $I^c$ | **焊死** | $I^c_i=I_i+\sum_{j\in\mu(i)}I^c_j$ (6.13) | 10 |
| $I^A$ | **自由** | $f=I^Aa+p^A$ | **21** |

$$
H_{ij}=S_i^{\mathsf T}I^{c}_{i}S_j\ \ (i\in\nu(j)，\text{即 }j\text{ 是祖先})\ (6.14)
\qquad
\text{无分支：}H_{ij}=S_i^{\mathsf T}I^{c}_{\max(i,j)}S_j\ (6.16)
$$

$$
I^a_j=I^A_j-U_jD_j^{-1}U_j^{\mathsf T}\ (7.47)
\qquad
p^a_j=p^A_j+I^a_jc_j+U_jD_j^{-1}u_j\ (7.48)
$$

$$
U_i=I^A_iS_i,\quad D_i=S_i^{\mathsf T}U_i,\quad u_i=\tau_i-S_i^{\mathsf T}p^A_i,\quad a'_i=a_{\lambda(i)}+c_i
$$

$$
f_j=I^A_ja_j+p^A_j=I^a_ja_{\lambda(j)}+p^a_j\ (7.25)
\qquad
I^A=(JH^{-1}J^{\mathsf T})^{-1}\ (7.12)=\Lambda
$$

## 12. 闭环与接触（第 8、11 章）

$$
K_{lj}=\epsilon_{lj}T_k^{\mathsf T}S_j\ (8.20)
\qquad
A=KH^{-1}K^{\mathsf T}\ (8.38)
\qquad
\text{mobility}=n-r\ (8.52)
$$

Baumgarte：$\ddot e+2\alpha\dot e+\beta^2e=\text{noise}$，$\alpha=\beta=1/T_{stab}$

$$
\zeta=n\cdot v\ (11.2)
\qquad
\dot\zeta=n\cdot a+\dot n\cdot v\ (11.3)
\qquad
\boxed{\dot\zeta\ge0,\ \lambda\ge0,\ \dot\zeta\lambda=0}\ (11.4)
$$

$$
M=T^{\mathsf T}H^{-1}T\ (11.35)
\qquad
t_i=(J_{sc(i)}-J_{pc(i)})^{\mathsf T}n_i=\Big(\tfrac{\partial\phi_i}{\partial q}\Big)^{\mathsf T}
$$

$$
\iota=I\Delta v\ (11.55)
\qquad
\iota=I^A\Delta v\ (11.56)
\qquad
u=H\Delta\dot q\ (11.57)
$$

$$
\lambda=\frac{-(1+e)\,n\cdot(v_2-v_1)}{n\cdot(I_1^{-1}+I_2^{-1})\,n}\ (11.65)
$$

## 13. 代价与复杂度（第 10 章表 10.1）

| 算法 | 乘法 | 加法 |
|---|---|---|
| RNEA | $93n-108$ | $81n-100$ |
| CRBA | $10n^2+22n-32$ | $6n^2+37n-43$ |
| F&S | $\tfrac16n^3+\tfrac32n^2-\tfrac23n$ | $\tfrac16n^3+n^2-\tfrac76n$ |
| ABA | $224n-259$ | $205n-248$ |

$$
\boxed{\ n\le8：O(n^3)\text{ 路线更快}\qquad n\ge9：\text{ABA 更快}\ }
$$

**真实复杂度**：CRBA 是 $O(nd)$，完整 FD 是 $O(nd^2)$（$d$ = 树深度）。
**$\kappa(H)$ 最坏 $\approx4N_B^4$** ⟹ 长链仿真要谨慎，生成树选**深度最小**的。

## 14. 快速自查（写代码前过一遍）

- [ ] 这个量属于 $M^6$ 还是 $F^6$？→ 决定用 $X$ 还是 $X^*$、用 $\times$ 还是 $\times^*$
- [ ] 两个要相加的空间向量在**同一个坐标系**里吗？
- [ ] 分量排列是 (角, 线) 吗？
- [ ] 加速度合成项是 $v_i\times v_J$，且**没有 2**？
- [ ] 重力用的是 $a_0=-a_g$（负号）？
- [ ] $H_{ij}$ 的下标：**$i$ 后代、$j$ 祖先、$I^c$ 取后代**（记忆锚点：式 6.16 的 $\max$）？
- [ ] ABA 的 $p^a$ 里用的是 $I^a$ 不是 $I^A$？
- [ ] 闭环：没忘 $\tau^a$？鞍点矩阵没用 Cholesky？
- [ ] 接触：$n$ 指向外侧 / 由前驱指向后继？
