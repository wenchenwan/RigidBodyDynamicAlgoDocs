# 第 2 章 空间向量代数 (Spatial Vector Algebra)

> **一句话概括**：建立一套 6 维向量代数，把刚体的角运动/线运动、力矩/力
> 统一在两个互为对偶的向量空间 $M^6$（运动）与 $F^6$（力）中，
> 并给出这套代数的全部运算规则。

> ★ **这是全书最重要的一章。** 后面 8 章全部建立在它之上。
> 建议投入总阅读时间的 30%，读到能不看书默写出坐标变换和叉乘的矩阵形式为止。

## 本章定位

第 1 章说"我们要用 6D 向量"，这一章就把这句话兑现成完整的数学工具：
定义对象 → 定义坐标 → 定义运算（加、数乘、标量积、坐标变换、叉乘、求导）
→ 用它写出单刚体的运动方程。学完这一章，你手上就有了一套能直接写算法的语言。

## 前置依赖

- 3D 向量代数、叉乘、反对称矩阵 $a\times$
- 旋转矩阵与坐标系变换
- 刚体运动学：角速度、刚体上任意点的速度公式 $v_P = v_O + \omega\times p$

## 核心概念

### 2.1 为什么是两个空间：$M^6$ 与 $F^6$

**关键洞察**：运动量和力量**不属于同一个向量空间**。

- 运动空间 $M^6$：装空间速度 $\mathbf{v}$、空间加速度 $\mathbf{a}$、关节运动子空间的基向量
- 力空间 $F^6$：装空间力 $\mathbf{f}$、空间动量 $\mathbf{h}$

两者是**对偶空间** (dual spaces)：$F^6 = (M^6)^*$。它们之间唯一的天然运算是
**配对 (pairing)**，其物理意义是**功率**：

$$
P = \mathbf{f}\cdot\mathbf{v} = \mathbf{f}^{\mathsf T}\mathbf{v}
$$

**为什么这件事这么重要**：因为 $M^6$ 上没有天然的内积。
"两个空间速度的点积"没有物理意义，也不是坐标变换下的不变量。
一旦承认这一点，下面几件初学者觉得"多余"的事就全都自然了：

| 现象 | 原因 |
|---|---|
| 力向量的坐标变换用 $X^*$ 而不是 $X$ | $X^*$ 是 $X$ 在对偶空间上的诱导变换 |
| 存在两个叉乘算子 $\times$ 和 $\times^*$ | 一个作用在 $M^6$ 上，一个作用在 $F^6$ 上 |
| 空间惯性 $I$ 不是"对称矩阵"那么简单 | 它是映射 $I: M^6 \to F^6$，其对称性是关于对偶配对而言的 |

> 💡 **记忆口诀**：带 `*` 的运算 = 作用在力上的。不带 `*` = 作用在运动上的。

### 2.2 空间速度 (Spatial Velocity)

刚体的运动状态由一个 6D 向量完全描述：

$$
\mathbf{v} = \begin{bmatrix}\omega \\ v_O\end{bmatrix} \in M^6
$$

- $\omega$：刚体角速度（与参考点无关）
- $v_O$：**刚体上此刻与原点 $O$ 重合的那个物质点**的线速度

⚠️ **第一个大坑**：$v_O$ **不是**"原点的速度"。原点是空间中的固定点，不会动。
$v_O$ 是把刚体想象成无限延展的刚性体，取其上恰好位于 $O$ 处的那个点的速度。
理解这一点，后面所有"为什么加速度不含离心项"的问题都会迎刃而解。

刚体上任意点 $P$（位置向量 $p$）的速度：

$$
v_P = v_O + \omega\times p
$$

**螺旋解释**：任何刚体瞬时运动都等价于绕某条轴的旋转 + 沿该轴的平移
（Chasles 定理），即一个**螺旋 (screw)**。螺距 (pitch)：

$$
h = \frac{\omega\cdot v_O}{\omega\cdot\omega}
$$

- $h = 0$：纯转动 → $\mathbf{v}$ 是**线向量 (line vector)**
- $\omega = 0$：纯平移 → $\mathbf{v}$ 是**自由向量 (free vector)**

### 2.3 空间力 (Spatial Force)

作用在刚体上的合力系同样由一个 6D 向量描述：

$$
\mathbf{f} = \begin{bmatrix}n_O \\ f\end{bmatrix} \in F^6
$$

- $n_O$：关于原点 $O$ 的合力矩
- $f$：合力（与参考点无关）

⚠️ **第二个坑：分量排列是反的**。空间速度是 $[\omega;\ v_O]$——
**与参考点无关的量在上**；空间力是 $[n_O;\ f]$——**与参考点相关的量在上**。
这不是笔误，是为了让对偶配对 $\mathbf{f}^{\mathsf T}\mathbf{v}$ 直接给出功率：

$$
\mathbf{f}^{\mathsf T}\mathbf{v} = n_O\cdot\omega + f\cdot v_O
$$

（力矩配角速度、力配线速度，正好是功率的两项。）

关于另一点 $P$ 的力矩：$n_P = n_O - p\times f$。

### 2.4 Plücker 坐标 (Plücker Coordinates)

上面写出的 6 个数字是在一组特定基下的坐标。这组基叫 **Plücker 基**：

**$M^6$ 的基**：$\{\mathbf{d}_{Ox},\mathbf{d}_{Oy},\mathbf{d}_{Oz},\ \mathbf{d}_{x},\mathbf{d}_{y},\mathbf{d}_{z}\}$

- $\mathbf{d}_{Ox}$：绕过 $O$ 的 $x$ 轴的单位角速度
- $\mathbf{d}_{x}$：沿 $x$ 方向的单位平移速度

**$F^6$ 的基**：$\{\mathbf{e}_{Ox},\mathbf{e}_{Oy},\mathbf{e}_{Oz},\ \mathbf{e}_{x},\mathbf{e}_{y},\mathbf{e}_{z}\}$

- $\mathbf{e}_{Ox}$：绕 $x$ 轴的单位力偶
- $\mathbf{e}_{x}$：沿过 $O$ 的 $x$ 轴的单位力

这两组基互为**对偶基 (reciprocal basis)**：$\mathbf{e}_i\cdot\mathbf{d}_j = \delta_{ij}$。
正因为选了对偶基，配对才能简单地写成转置乘法 $\mathbf{f}^{\mathsf T}\mathbf{v}$，
而不需要插入一个度量矩阵。

> **一个坐标系 = 一个 Plücker 坐标系**。原书里说"在 body $i$ 的坐标系中表示"，
> 意思就是用以 body $i$ 的原点和轴向定义的那组 Plücker 基。

### 2.5 坐标变换 (Coordinate Transforms)

设坐标系 $B$ 相对坐标系 $A$：旋转矩阵 $E = {}^{B}R_{A}$（把 $A$ 系分量转成 $B$ 系分量），
$B$ 原点在 $A$ 系中的位置为 $r$。则：

**运动向量的变换**

$$
{}^{B}X_{A} = \begin{bmatrix} E & 0 \\ -E\,r\times & E \end{bmatrix},
\qquad
{}^{B}\mathbf{v} = {}^{B}X_{A}\,{}^{A}\mathbf{v}
$$

**力向量的变换**

$$
{}^{B}X_{A}^{*} = \begin{bmatrix} E & -E\,r\times \\ 0 & E \end{bmatrix},
\qquad
{}^{B}\mathbf{f} = {}^{B}X_{A}^{*}\,{}^{A}\mathbf{f}
$$

**两者的关系**（务必记住，考试必考级别的重点）：

$$
X^{*} = X^{-\mathsf T} = (X^{-1})^{\mathsf T}
\qquad\Longleftrightarrow\qquad
{}^{B}X_{A}^{*} = \left({}^{A}X_{B}\right)^{\mathsf T}
$$

**为什么必须这样**：功率是物理量，不能随坐标系改变。

$$
{}^{B}\mathbf{f}^{\mathsf T}\,{}^{B}\mathbf{v}
= \left(X^{*}\,{}^{A}\mathbf{f}\right)^{\mathsf T}\left(X\,{}^{A}\mathbf{v}\right)
= {}^{A}\mathbf{f}^{\mathsf T}\left(X^{*\mathsf T}X\right){}^{A}\mathbf{v}
$$

要让它等于 ${}^{A}\mathbf{f}^{\mathsf T}\,{}^{A}\mathbf{v}$，就必须
$X^{*\mathsf T}X = \mathbf{1}$，即 $X^{*} = X^{-\mathsf T}$。**对偶性直接由能量守恒逼出来。**

**性质**

- $X$ 一般**不是**正交矩阵（只有纯旋转时才是），所以 $X^{-1}\ne X^{\mathsf T}$
- 复合：${}^{C}X_{A} = {}^{C}X_{B}\,{}^{B}X_{A}$（和普通变换一样按顺序串）
- $\left({}^{B}X_{A}\right)^{-1} = {}^{A}X_{B}$，且 $(XY)^* = X^*Y^*$

### 2.6 空间叉乘 (Spatial Cross Products)

对 $\mathbf{v} = [\omega;\ v_O]$，定义两个算子：

**作用在运动向量上**（$M^6\to M^6$）：

$$
\mathbf{v}\times = \begin{bmatrix}\omega\times & 0 \\ v_O\times & \omega\times\end{bmatrix}
$$

**作用在力向量上**（$F^6\to F^6$）：

$$
\mathbf{v}\times^{*} = -\left(\mathbf{v}\times\right)^{\mathsf T}
= \begin{bmatrix}\omega\times & v_O\times \\ 0 & \omega\times\end{bmatrix}
$$

**关键恒等式**：

$$
\mathbf{v}\times^{*} = -\left(\mathbf{v}\times\right)^{\mathsf T}
$$

同样可由功率守恒推出：$\frac{d}{dt}(\mathbf{f}^{\mathsf T}\mathbf{v})$ 展开时
两个叉乘项必须相消。

**性质**

- $\mathbf{v}\times\mathbf{v} = 0$（同一向量自叉为零，RNEA 里会用到）
- $\mathbf{v}\times$ 满足 Jacobi 恒等式，$(M^6, \times)$ 构成李代数
  （即 $\mathfrak{se}(3)$，$\mathbf{v}\times$ 就是 $\mathrm{ad}_{\mathbf v}$，
  $\mathbf{v}\times^*$ 就是 $\mathrm{ad}^*_{\mathbf v}$）
- 变换律：$X(\mathbf{v}\times)X^{-1} = (X\mathbf{v})\times$

### 2.7 在运动坐标系中求导（全章最容易错的地方）

设坐标系本身以空间速度 $\mathbf{v}_{\text{frame}}$ 运动。对任意向量：

$$
\underbrace{\frac{d\mathbf{m}}{dt}}_{\text{绝对导数}}
= \underbrace{\frac{\mathring{d}\mathbf{m}}{dt}}_{\text{表观导数}} + \mathbf{v}_{\text{frame}}\times\mathbf{m}
\qquad (\mathbf{m}\in M^6)
$$

$$
\frac{d\mathbf{g}}{dt}
= \frac{\mathring{d}\mathbf{g}}{dt} + \mathbf{v}_{\text{frame}}\times^{*}\mathbf{g}
\qquad (\mathbf{g}\in F^6)
$$

其中"表观导数 (apparent derivative)"= 只对分量求导、把基当成常量。
这是 3D 里 $\frac{d\mathbf{u}}{dt} = \frac{\mathring{d}\mathbf{u}}{dt}+\omega\times\mathbf{u}$
的 6D 推广，**注意运动量用 $\times$、力量用 $\times^*$**。

### 2.8 空间加速度 (Spatial Acceleration) —— 全书最反直觉的概念

定义极其简单：

$$
\mathbf{a} = \dot{\mathbf{v}} = \begin{bmatrix}\dot\omega \\ \dot v_O\end{bmatrix}
$$

⚠️ **第三个大坑，也是最大的坑**：$\dot v_O$ **不是**任何物质点的加速度。

设 $P$ 是刚体上此刻与原点重合的物质点，它的**经典加速度**（真实物理加速度）为：

$$
a_{P}^{\text{classical}} = \dot v_O + \omega\times v_O
$$

即：

$$
\boxed{\ \text{经典加速度} = \text{空间加速度的线分量} + \omega\times v_O\ }
$$

**为什么要用这么反直觉的定义**？因为它让**加速度的合成变成（几乎）线性**：

$$
\mathbf{v}_2 = \mathbf{v}_1 + \mathbf{v}_{\text{rel}}
\quad\Longrightarrow\quad
\mathbf{a}_2 = \mathbf{a}_1 + \mathbf{a}_{\text{rel}} + \mathbf{v}_1\times\mathbf{v}_{\text{rel}}
$$

对比经典记法：那里会冒出**离心加速度** $\omega\times(\omega\times r)$ 和
**科氏加速度** $2\omega\times v$ 两组项，而且随着连杆串联层层嵌套，展开式爆炸增长。
空间记法下只剩一个 $\mathbf{v}_1\times\mathbf{v}_{\text{rel}}$，而且**没有系数 2**。

> 💡 **这就是 RNEA 只有几行的根本原因。** 那个"消失的因子 2"和"消失的离心项"
> 并没有消失，它们被吸收进了空间加速度的定义里，最终会在
> $\mathbf{v}\times^{*}I\mathbf{v}$ 那一项重新出现。

**什么时候必须换回经典加速度**：当你要和真实的加速度计读数、
或者要计算某个具体点的轨迹时。日常算动力学时**一律用空间加速度**。

### 2.9 空间动量 (Spatial Momentum) 与空间惯性 (Spatial Inertia)

**空间动量**：

$$
\mathbf{h} = \begin{bmatrix}h_\omega \\ h_v\end{bmatrix} = I\mathbf{v} \in F^6
$$

- $h_v = m\,v_C$：线动量（$v_C$ 为质心速度）
- $h_\omega$：关于 $O$ 的角动量

**空间刚体惯性**：设质量 $m$、质心位置 $c$（相对 $O$）、绕质心的转动惯量 $\bar I$，

$$
I_O = \begin{bmatrix}
\bar I + m\,c\times c\times^{\mathsf T} & m\,c\times \\[2pt]
m\,c\times^{\mathsf T} & m\,\mathbf{1}_3
\end{bmatrix}
$$

注意 $c\times^{\mathsf T} = -c\times$，所以左下块也可写成 $-m\,c\times$，
左上块也可写成 $\bar I - m\,c\times c\times$（这是平行轴定理的 6D 形式）。

**推导线索**（自己动手验证一遍，收益很大）：
- 线动量 $h_v = m(v_O + \omega\times c) = m\,v_O - m\,c\times\omega$ → 给出下面一行
- 角动量 $h_\omega = \bar I\omega + c\times(m v_C)$，展开即得上面一行

**性质**

- 对称：$I = I^{\mathsf T}$（这里的对称性指关于对偶配对，物理上对应
  动能 $T = \tfrac12\mathbf{v}^{\mathsf T}I\mathbf{v}$ 良定义）
- 正定（$m>0$ 且 $\bar I$ 正定时）
- 可加：同一坐标系下，刚性连接的两个刚体 $I = I_1 + I_2$
  （**复合刚体惯性的基础，第 6 章 CRBA 直接用**）
- 只有 **10 个独立参数**（$m$、$m c$ 共 3 个、$\bar I$ 共 6 个），
  虽然矩阵有 36 个元素、对称后 21 个。实现时存 10 个数即可。

**坐标变换**：

$$
{}^{B}I = {}^{B}X_{A}^{*}\ {}^{A}I\ {}^{A}X_{B}
= \left({}^{A}X_{B}\right)^{\mathsf T}\ {}^{A}I\ {}^{A}X_{B}
$$

**怎么记**：$I$ 是 $M^6\to F^6$ 的映射，所以右边要先把运动量转过去（$X$），
左边再把结果的力量转回来（$X^*$）。三明治结构，两侧的变换方向相反。

### 2.10 单刚体运动方程

$$
\boxed{\ \mathbf{f} = I\mathbf{a} + \mathbf{v}\times^{*}I\,\mathbf{v}\ }
$$

**推导**（干净利落，值得记住）：

在随体坐标系中 $I$ 是常量，所以表观导数 $\mathring{\dot{\mathbf h}} = I\mathbf{a}$。
动量是力向量，用 $\times^*$ 版本的求导公式：

$$
\mathbf{f} = \frac{d\mathbf{h}}{dt}
= \frac{\mathring d(I\mathbf{v})}{dt} + \mathbf{v}\times^{*}(I\mathbf{v})
= I\mathbf{a} + \mathbf{v}\times^{*}I\mathbf{v}
$$

**物理意义**：
- $I\mathbf{a}$：牛顿-欧拉方程的"质量×加速度"部分
- $\mathbf{v}\times^{*}I\mathbf{v}$：**偏置力 (bias force)**，
  一次性打包了所有离心力和陀螺力矩项（在 3D 记法里就是欧拉方程的
  $\omega\times(\bar I\omega)$ 加上离心项）

**这一条式子 = 牛顿方程 + 欧拉方程。** 展开它就能看到熟悉的 3D 形式。

**等价形式**：$\dot I = \mathbf{v}\times^{*}I - I\,\mathbf{v}\times$，
于是 $\mathbf{f} = \frac{d}{dt}(I\mathbf{v}) = I\mathbf{a} + \dot I\mathbf{v}$
（用到 $\mathbf{v}\times\mathbf{v}=0$）。

**动能**：$T = \tfrac12\,\mathbf{v}^{\mathsf T}I\,\mathbf{v}$

## 计算实现要点

**不要真的用 6×6 矩阵做变换。** 这是新手最常见的性能问题。

| 对象 | 朴素存法 | 推荐存法 | 说明 |
|---|---|---|---|
| 变换 $X$ | 6×6 = 36 个数 | $(E, r)$ = 12 个数 | 用 $E,r$ 直接算，避免 $6\times6$ 矩阵乘（216 乘法 → 约 30 乘法） |
| 空间惯性 $I$ | 6×6 = 36 个数 | $(m, c, \bar I)$ = 10 个数 | 变换公式有专门的紧凑形式 |
| 叉乘 $\mathbf{v}\times$ | 构造 6×6 再乘 | 直接按块写公式 | $\mathbf{v}\times\mathbf{m}$ 只要 2 次 3D 叉乘 + 1 次加法 |

对应 `spatial_v2` 里的函数：`plux`（构造/分解 Plücker 变换）、
`rbi`（构造刚体惯性）、`crm`/`crf`（$\times$ 与 $\times^*$）、`Xrot*`/`xlt`（基本变换）。

〔待核对〕函数名以你下载的版本为准。

## 易错点与陷阱汇总

1. **$v_O$ 不是原点的速度**，是刚体上瞬时与原点重合的物质点的速度。
2. **分量顺序**：$\mathbf{v}=[\omega; v_O]$ 但 $\mathbf{f}=[n_O; f]$，角在上、矩在上。
   跨库对照时尤其小心：部分文献/库用 (线, 角) 排列，两者**不能混用**。
3. **$X^{*}\ne X$，$X^{-1}\ne X^{\mathsf T}$**。$X$ 只有在纯旋转时才正交。
4. **空间加速度 ≠ 物质点的加速度**，差一个 $\omega\times v_O$。
5. **$\times$ 与 $\times^{*}$ 不能混用**。判断依据只有一条：被作用的对象是运动量还是力量。
6. **$I$ 的变换是三明治 $X^{*}IX$**，不是相似变换 $XIX^{-1}$，也不是 $X^{\mathsf T}IX$
   （注意 $X^*$ 用的是**逆**的转置）。
7. **不同刚体的空间向量相加前必须先变换到同一坐标系**。
   6 个数字看起来能加，但物理上无意义。这是调试时最常见的 bug。

## 与其他章的联系

- ← 第 1 章：兑现"用空间向量"的承诺
- → 第 3 章：把单刚体的 $\mathbf{f}=I\mathbf{a}+\mathbf{v}\times^*I\mathbf{v}$ 推广到多体系统
- → 第 5 章 RNEA：外推用 $\mathbf{v},\mathbf{a}$ 的合成公式，内推用 $X^*$ 传递力
- → 第 6 章 CRBA：直接用"空间惯性可加"这条性质
- → 第 7 章 ABA：把 $I$ 推广成铰接体惯性 $I^A$

---

## ✍️ 我的理解

<!-- 建议自测：不看书默写出 X、X*、v×、v×*、I 的矩阵形式，以及运动方程 -->

## ❓ 疑问与待办

- [ ] 手推一遍空间惯性矩阵 $I_O$ 的两行，验证与书上一致
- [ ] 手推一遍 $X^{*} = X^{-\mathsf T}$
- [ ] 用具体数值验证：空间加速度线分量 + $\omega\times v_O$ = 经典加速度
- [ ] 核对本章小节编号与原书是否一致

## 📌 与原文的出入

<!-- 发现不符时记录在此并修正正文 -->
