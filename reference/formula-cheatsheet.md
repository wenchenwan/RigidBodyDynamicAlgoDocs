# 空间向量代数公式速查

> 一页纸版本。打印出来贴在桌上的那种。
> 详细解释见 [`docs/ch02-spatial-vector-algebra.md`](../docs/ch02-spatial-vector-algebra.md)。

## 1. 基本对象

$$
\mathbf{v}=\begin{bmatrix}\omega\\ v_O\end{bmatrix}\in M^6,
\qquad
\mathbf{f}=\begin{bmatrix}n_O\\ f\end{bmatrix}\in F^6
$$

⚠️ 注意排列：速度是 **(角, 线)**，力是 **(矩, 力)**。

**功率（唯一的天然配对）**：

$$
P=\mathbf{f}^{\mathsf T}\mathbf{v}=n_O\cdot\omega+f\cdot v_O
$$

**动能**：$\;T=\tfrac12\,\mathbf{v}^{\mathsf T}I\,\mathbf{v}\;$，
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
\mathbf{v}\times=\begin{bmatrix}\omega\times & 0\\ v_O\times & \omega\times\end{bmatrix}
\qquad
\mathbf{v}\times^{*}=\begin{bmatrix}\omega\times & v_O\times\\ 0 & \omega\times\end{bmatrix}
$$

$$
\boxed{\mathbf{v}\times^{*}=-\left(\mathbf{v}\times\right)^{\mathsf T}}
\qquad
\mathbf{v}\times\mathbf{v}=0
\qquad
X(\mathbf{v}\times)X^{-1}=(X\mathbf{v})\times
$$

## 4. 求导（坐标系以 $\mathbf{v}$ 运动）

$$
\frac{d\mathbf{m}}{dt}=\frac{\mathring{d}\mathbf{m}}{dt}+\mathbf{v}\times\mathbf{m}
\quad(\mathbf{m}\in M^6)
\qquad
\frac{d\mathbf{g}}{dt}=\frac{\mathring{d}\mathbf{g}}{dt}+\mathbf{v}\times^{*}\mathbf{g}
\quad(\mathbf{g}\in F^6)
$$

（$\mathring d$ = 表观导数，把基当常量）

## 5. 速度与加速度的合成

$$
\mathbf{v}_2=\mathbf{v}_1+\mathbf{v}_{\text{rel}}
\qquad
\boxed{\mathbf{a}_2=\mathbf{a}_1+\mathbf{a}_{\text{rel}}+\mathbf{v}_1\times\mathbf{v}_{\text{rel}}}
$$

⚠️ **没有系数 2**，没有单独的离心项和科氏项——它们被空间加速度的定义吸收了。

**空间加速度 ↔ 经典加速度**：

$$
a_{P}^{\text{classical}}=\underbrace{\dot v_O}_{\mathbf{a}\text{ 的线分量}}+\ \omega\times v_O
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

**动量**：$\;\mathbf{h}=I\mathbf{v}$

## 7. 单刚体运动方程

$$
\boxed{\ \mathbf{f}=I\,\mathbf{a}+\mathbf{v}\times^{*}I\,\mathbf{v}\ }
$$

等价：$\;\dot I=\mathbf{v}\times^{*}I-I\,\mathbf{v}\times\;$，
$\;\mathbf{f}=\frac{d}{dt}(I\mathbf{v})=I\mathbf{a}+\dot I\mathbf{v}$

## 8. 系统运动方程

$$
\boxed{\ H(q)\,\ddot q+C(q,\dot q)=\tau\ }
$$

- $H$：对称、正定、只依赖 $q$、**分支诱导稀疏**（$H_{ij}\ne0\iff$ 祖先-后代关系）
- $C$ 是**向量**（含科氏+离心+重力），$C=\mathrm{ID}(q,\dot q,\mathbf 0)$

**关节约束**：

$$
\mathbf{v}_{\text{rel}}=S_i\dot q_i,
\qquad S_i^{\mathsf T}T_i=0,
\qquad \boxed{\tau_i=S_i^{\mathsf T}\mathbf{f}_i}
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

## 10. 快速自查（写代码前过一遍）

- [ ] 这个量属于 $M^6$ 还是 $F^6$？→ 决定用 $X$ 还是 $X^*$、用 $\times$ 还是 $\times^*$
- [ ] 两个要相加的空间向量在**同一个坐标系**里吗？
- [ ] 分量排列是 (角, 线) 吗？
- [ ] 加速度合成项是 $\mathbf{v}_i\times\mathbf{v}_J$，且**没有 2**？
- [ ] 重力用的是 $\mathbf{a}_0=-\mathbf{a}_g$（负号）？
