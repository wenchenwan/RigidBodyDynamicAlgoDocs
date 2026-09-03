# 第 2 章 空间向量代数 (Spatial Vector Algebra)

> **原书 pp. 7–38**，共 17 节。全书最重要的一章，后面 9 章全部建立在它之上。

> **一句话概括**：建立两个互为对偶的 6 维向量空间 $M^6$（运动）与 $F^6$（力），
> 并给出这套代数的全部运算规则，最后写出单刚体运动方程 $\hat f = I\hat a + \hat v\times^{*}I\hat v$。

## 📌 你在 PDF 上标注的三处困惑

本章笔记针对这三处写了专门小节，可直接跳读：

| 位置        | 你的批注                                   | 本笔记对应小节                                                               |
| ----------- | ------------------------------------------ | ---------------------------------------------------------------------------- |
| p.23 §2.5  | 「这里的理解仍然有困惑」                   | [§2.5 线向量与自由向量](#25-线向量与自由向量-line-vectors-and-free-vectors) |
| p.40 §2.13 | 「$g_i$ 的实际含义是什么，需要考虑？？？」 | [§2.13-B 惯性的并矢表示](#213-b-惯性的并矢表示g_i-到底是什么)               |
| p.42 §2.13 | 「why???」（为什么 10 个参数）             | [§2.13-C 为什么恰好是 10 个参数](#213-c-为什么恰好是-10-个参数)             |

所有结论都有可运行验证：`python3 code/verify_ch02.py`

## 原书节次对照

| 节   | 原书标题                      | 页 | 本笔记                               |
| ---- | ----------------------------- | -- | ------------------------------------ |
| 2.1  | Mathematical Preliminaries    | 7  | 数学预备：四个向量空间、对偶基、并矢 |
| 2.2  | Spatial Velocity              | 10 | 空间速度                             |
| 2.3  | Spatial Force                 | 13 | 空间力                               |
| 2.4  | Plücker Notation             | 15 | Plücker 记法                        |
| 2.5  | Line Vectors and Free Vectors | 16 | **线向量与自由向量** ⚠️            |
| 2.6  | Scalar Product                | 17 | 标量积与对偶性                       |
| 2.7  | Using Spatial Vectors         | 18 | 使用规则清单                         |
| 2.8  | Coordinate Transforms         | 20 | 坐标变换                             |
| 2.9  | Spatial Cross Products        | 23 | 空间叉乘                             |
| 2.10 | Differentiation               | 25 | 微分                                 |
| 2.11 | Acceleration                  | 28 | 空间加速度 vs 经典加速度             |
| 2.12 | Momentum                      | 31 | 空间动量                             |
| 2.13 | Inertia                       | 32 | **空间惯性** ⚠️                    |
| 2.14 | Equation of Motion            | 35 | 运动方程                             |
| 2.15 | Inverse Inertia               | 36 | 逆惯性                               |
| 2.16 | Planar Vectors                | 37 | 平面向量                             |
| 2.17 | Further Reading               | 38 | 与螺旋理论/李代数的关系              |

---

## 2.1 数学预备 (Mathematical Preliminaries)

### 四个向量空间

原书明确区分四类向量，这个区分是理解后面一切的前提：

| 记号           | 名称                          | 特点                                     |
| -------------- | ----------------------------- | ---------------------------------------- |
| $\mathbb{R}^n$ | 坐标向量 (coordinate vectors) | $n$ 元实数组，即 $n\times1$ 列矩阵       |
| $\mathbb{E}^n$ | 欧氏向量 (Euclidean vectors)  | **定义了欧氏内积**，因而有"长度"和"夹角" |
| $M^n$          | 空间运动向量                  | 速度、加速度、运动自由度方向             |
| $F^n$          | 空间力向量                    | 力、冲量、动量、约束力方向               |

> 🔑 **本章最关键的一句话**：**空间向量不是欧氏向量。**
> $M^6$ 上**没有**内积，$F^6$ 上也没有。$\hat m\cdot\hat m$ 和 $\hat f\cdot\hat f$ **是没有定义的**。
> 唯一定义的是 $M^6$ 与 $F^6$ **之间**的标量积。
>
> 后面所有"为什么要分两个空间/两种变换/两种叉乘"的问题，答案都在这一句里。

**术语区分**：**抽象向量 (abstract vector)** 是被表示的那个几何/物理对象；
**坐标向量**是它在某组基下的一列数。二者需要区分时，坐标向量加下划线（$\underline v$）。

### 帽子与下划线记号

- **帽子**：空间向量与 3D 向量同时出现、且重名时，空间向量戴帽子——$\hat v,\ \hat f$。
  本笔记在容易混淆处保留帽子，其余场合省略。
- **下划线**：坐标向量与抽象向量需区分时用。

### 对偶基与互易条件 (Dual Basis / Reciprocity Condition)

设 $D=\{d_1,\dots,d_n\}\subset U$、$E=\{e_1,\dots,e_n\}\subset V$，若

$$
d_i\cdot e_j=\begin{cases}1 & i=j\\ 0 & \text{否则}\end{cases}
$$

则称 $D,E$ 构成**对偶基**，二者互为对方的**互易基 (reciprocal)**，记 $E=D^{*}$。

**要点**：

- 一组对偶基包含**两套基向量**，覆盖**两个**向量空间；
- 给定其中一套，另一套由互易条件**唯一确定**；
- **正交归一基是对偶基的特例**——当 $U=V=\mathbb{E}^n$ 且 $D=E$ 时。
  这解释了为什么在普通欧氏空间里我们从来不需要操心"对偶"：那里两套基恰好重合了。

### 对偶坐标 (Dual Coordinates)

用对偶基定义的坐标系叫对偶坐标系。它的**核心性质**：

$$
u\cdot v=\underline u^{\mathsf T}\underline v
$$

即抽象向量的标量积 = 坐标向量的普通转置乘法，**中间不需要插入度量矩阵**。
分量也有干净的表达：$u_i=e_i\cdot u$，$v_i=d_i\cdot v$。

**代价**：坐标变换需要**两个**矩阵。若 $X$ 变换 $\underline u$，则变换 $\underline v$ 的矩阵记作 $X^{*}$，且

$$
\boxed{X^{*}=X^{-\mathsf T}}\tag{原书 2.13}
$$

> 📝 这正是你在 p.25 批注的「对运动向量和力向量分别进行坐标变换的关系」。
> 它的来源就是上面这条：要让 $\underline u^{\mathsf T}\underline v$ 在变换后不变，
> 必须 $(X^{*}\underline v)^{\mathsf T}(X\underline u)=\underline v^{\mathsf T}\underline u$，
> 即 $X^{*\mathsf T}X=\mathbf 1$。**这是纯代数结论，与物理无关**；
> 物理上它对应功率守恒（见 §2.6）。

### 算子 $a\cdot$ 与 $a\times$

原书把 $a\cdot$ 和 $a\times$ 看作**算子**：

- $a\cdot$ 把 $b$ 映到标量 $a\cdot b$；若 $a$ 是坐标向量则 $a\cdot = a^{\mathsf T}$
- $a\times$ 把 $b$ 映到 $a\times b$；是一个方阵

约定：$a\times b\times c$ 读作 $(a\times)(b\times)c$；$E\,a\times$ 读作 $E(a\times)$。
其他结合方式必须加括号。

### 并矢与并矢张量 (Dyads and Dyadics)

**这一小节是 §2.13 的伏笔，直接关系到你标注的 $g_i$ 问题，务必读懂。**

**并矢 (dyad)**：形如 $a\,b\cdot$ 的表达式，是一个**线性算子**：

$$
(a\,b\cdot)\,c = a\,(b\cdot c)
$$

即"先用 $b\cdot$ 把 $c$ 压成一个标量，再用这个标量去缩放 $a$"。
结果永远是 $a$ 的倍数，所以**并矢的秩恒为 1**。

坐标形式：$a\,b\cdot \;\longleftrightarrow\; \underline a\,\underline b^{\mathsf T}$（外积，秩 1 矩阵）。

**并矢的类型是"跨空间"的**：可以 $a\in U$、$b\in V$、$c\in W$，
只要 $V$ 与 $W$ 之间定义了标量积。此时这个并矢是 $W\to U$ 的映射。

**并矢张量 (dyadic)**：一般的线性映射，可写成 $r$ 个线性无关并矢之和：

$$
L=\sum_{i=1}^{r}a_i\,b_i\cdot \qquad\longleftrightarrow\qquad
L=\sum_{i=1}^{r}\underline a_i\underline b_i^{\mathsf T}\ (\text{秩 } r \text{ 的 } m\times n \text{ 矩阵})
$$

> 💡 **翻译成线性代数**：并矢 = 秩 1 矩阵（外积）；并矢张量 = 秩 $r$ 矩阵写成 $r$ 个外积之和。
> 这就是"矩阵的秩 1 分解"。**空间惯性会用这个语言来定义**（§2.13）。

---

## 2.2 空间速度 (Spatial Velocity)

### 构造过程

取刚体 $B$ 和空间中一个**固定点** $O$（位置任意）。刚体的运动由两个 3D 向量描述：

- $\omega$：角速度
- $v_O$：**刚体上此刻与 $O$ 重合的那个物质点**的线速度

⚠️ **第一个坑**：$v_O$ **不是**"原点的速度"（原点固定不动）。
原书 §2.11 说得更透彻：$v_O$ 是**物质点流过 $O$ 的"流量"的度量**（flow of points through $O$）。
理解这一点，后面空间加速度的反直觉就消失了。

刚体上任意点 $P$ 的速度：

$$
v_P=v_O+\omega\times\overrightarrow{OP}\tag{2.1}
$$

**关键**：右端两项各自依赖 $O$ 的位置，但**和不依赖**——依赖性相消了。

### 速度场的视角（原书例 2.1）

式 2.1 给空间中每个点配了一个向量，**它定义了一个向量场**——刚体的速度场：

$$
V(P)=v_O+\omega\times\overrightarrow{OP}\tag{2.5}
$$

> 💡 **这个视角很有用**：空间速度 $\hat v$ 就是"整个速度场"这个对象本身，
> 而不是某个点的速度。于是 $\hat v_{\text{sum}}=\hat v_1+\hat v_2$
> 对应 $V_{\text{sum}}(P)=V_1(P)+V_2(P)$ 对**所有** $P$ 成立。
> 这解释了空间向量为什么能像普通向量那样相加。

### Plücker 基

在 $O$ 处建笛卡尔坐标系 $Oxyz$，定义 $M^6$ 上的基（原书图 2.1b）：

$$
D_O=\{d_{Ox},d_{Oy},d_{Oz},\ d_x,d_y,d_z\}\subset M^6\tag{2.2}
$$

- $d_{Ox},d_{Oy},d_{Oz}$：绕直线 $Ox,Oy,Oz$ 的**单位转动**
- $d_x,d_y,d_z$：沿 $x,y,z$ 方向的**单位平移**

于是刚体运动 = 六个基本运动之和：

$$
\hat v=\omega_x d_{Ox}+\omega_y d_{Oy}+\omega_z d_{Oz}+v_{Ox}d_x+v_{Oy}d_y+v_{Oz}d_z\tag{2.3}
$$

$$
\hat v_O=\begin{bmatrix}\omega\\ v_O\end{bmatrix}\in\mathbb{R}^6\tag{2.4}
$$

**每一项都依赖坐标系，但整体不变**（原书例 2.2 证明了这点）。

### 例 2.2 的要点：基向量本身会随参考点改变

取 $P=(r,0,0)$，则

$$
d_{Px}=d_{Ox},\qquad d_{Py}=d_{Oy}+r\,d_z,\qquad d_{Pz}=d_{Oz}-r\,d_y
$$

**怎么理解 $d_{Py}=d_{Oy}+rd_z$**：想象刚体以单位角速度绕 $Py$ 转动。
此时刚体上位于 $O$ 的物质点会有一个大小为 $r$、沿 $z$ 方向的线速度。
所以"绕 $Py$ 的单位转动" = "绕 $Oy$ 的单位转动" + "沿 $z$ 的大小为 $r$ 的平移"。

> 💡 **这个例子的价值**：它说明 Plücker 基向量**不是抽象符号**，
> 每一个都是一个具体的刚体运动，可以画出来、可以叠加。
> 换参考点时坐标变了、基也变了，两者的变化恰好抵消。

---

## 2.3 空间力 (Spatial Force)

完全平行的构造。最一般的力系 = 过 $O$ 的线力 $f$ + 力偶 $n_O$（关于 $O$ 的总力矩）。

$$
n_P=n_O+f\times\overrightarrow{OP}\tag{2.6}
$$

⚠️ **注意 2.6 与 2.1 的符号差异**：$v_P=v_O+\omega\times\overrightarrow{OP}$，
而 $n_P=n_O+f\times\overrightarrow{OP}$。两式形式完全一致
（都是"第二个 3D 量 + 第一个 3D 量 × 位移"），只是角色对调。

$F^6$ 上的 Plücker 基（原书图 2.3b）：

$$
E_O=\{e_x,e_y,e_z,\ e_{Ox},e_{Oy},e_{Oz}\}\subset F^6\tag{2.7}
$$

- $e_x,e_y,e_z$：沿 $x,y,z$ 的**单位力偶**
- $e_{Ox},e_{Oy},e_{Oz}$：沿直线 $Ox,Oy,Oz$ 的**单位线力**

$$
\hat f_O=\begin{bmatrix}n_O\\ f\end{bmatrix}\tag{2.9}
$$

---

## 2.4 Plücker 记法 (Plücker Notation)

**一般记法**（原书 §2.4）：对 $M^6$ 与 $F^6$ 的一般元素

$$
\hat m_O=\begin{bmatrix}m\\ m_O\end{bmatrix},
\qquad
\hat f_O=\begin{bmatrix}f_O\\ f\end{bmatrix}
$$

**排列约定**：本书一律 **角在前、线在后 (angular-before-linear)**。

⚠️ 原书明确提醒：**文献中另一种排列（线在前）也完全合法**，
数学上没有区别，但**软件是按某一种写死的**，混用必然出错。

| 约定                           | 排列     | 典型使用者            |
| ------------------------------ | -------- | --------------------- |
| **本书 / `spatial_v2` / RBDL** | (角, 线) | $[\omega;v]$、$[n;f]$ |
| 部分螺旋理论文献 / 部分库      | (线, 角) | $[v;\omega]$、$[f;n]$ |

> 这是跨资料/跨库移植代码的**头号错误来源**。第一件事永远是确认排列顺序。

---

## 2.5 线向量与自由向量 (Line Vectors and Free Vectors)

> 📌 **你在 p.23 标注了「这里的理解仍然有困惑」。本节写得最详细。**

### A. 两个定义

| 类型                       | 定义                                        | 由几个数确定 | $M^6$ 中的例子 | $F^6$ 中的例子 |
| -------------------------- | ------------------------------------------- | ------------ | -------------- | -------------- |
| **线向量 (line vector)**   | 由一条**有向直线** + 一个**大小**刻画       | 5            | 纯转动         | 沿一条线的纯力 |
| **自由向量 (free vector)** | 由一个**方向** + 一个**大小**刻画（无位置） | 3            | 纯平移         | 纯力偶         |

**为什么线向量是 5 个数**：直线在 3D 中有 4 个自由度（方向 2 + 位置 2），
再加一个大小 = 5。而 Plücker 坐标有 6 个数，所以**线向量满足一个约束**——就是下面的判据。

**为什么自由向量是 3 个数**：只有方向和大小，"作用在哪条线上"没有意义
（力偶可以平移到任何地方，纯平移也一样）。

### B. 判据

设 $\hat s$ 是任意空间向量（运动或力），$s$ 与 $s_O$ 是它的两个 3D 分量。

> **约定**：$s$ 是**方向部分**（运动向量取 $\omega$，力向量取 $f$），
> $s_O$ 是**关于 $O$ 的矩部分**（运动取 $v_O$，力取 $n_O$）。
> ⚠️ 注意这与 §2.4 的书写顺序不同：运动向量是 $[s;s_O]$，力向量是 $[s_O;s]$。
> 判据说的是 3D 分量本身，与它们在 6D 列向量里排第几无关。

$$
\boxed{s=0\ \Longrightarrow\ \hat s\ \text{是自由向量}}
$$

$$
\boxed{s\cdot s_O=0\ \Longrightarrow\ \hat s\ \text{是线向量}}
$$

线向量的**直线**是满足 $\overrightarrow{OP}\times s=s_O$ 的所有点 $P$，方向由 $s$ 给出。

**判据从哪来**：设刚体绕过 $P$、方向 $s$ 的直线做单位转动。
则 $\omega=s$，而 $O$ 处物质点的速度是

$$
v_O=\omega\times\overrightarrow{PO}=s\times(-\overrightarrow{OP})=\overrightarrow{OP}\times s
$$

所以 $s_O=\overrightarrow{OP}\times s$，从而 $s\cdot s_O=s\cdot(\overrightarrow{OP}\times s)=0$ ✓
（一个向量点乘"含自己的叉积"必为零）。

### C. ⭐ 分解定理——困惑的核心

原书给了**两条**分解，很容易混在一起。它们回答的是**不同的问题**。

#### 分解一：指定直线过哪个点 → 分解唯一

> 原书原话：*"Any spatial vector can be expressed as the sum of a line vector and a
> free vector. If the line vector must pass through a given point, then the expression is unique."*

**问题**：给定 $\hat s$ 和一个**你指定的点** $P$，把 $\hat s$ 拆成
"过 $P$ 的线向量" + "自由向量"。

**答案**：

$$
\underbrace{\begin{bmatrix}s\\ \overrightarrow{OP}\times s\end{bmatrix}}_{\text{过 }P\text{、方向 }s\text{ 的线向量}}
\ +\
\underbrace{\begin{bmatrix}0\\ s_O-\overrightarrow{OP}\times s\end{bmatrix}}_{\text{自由向量}}
\ =\ \begin{bmatrix}s\\ s_O\end{bmatrix}
$$

**为什么唯一**：线向量的方向部分必须等于 $s$（否则加起来第一块对不上），
而"方向 $s$ + 必须过 $P$"就把这条直线**完全钉死**了。
直线定了，线向量就定了，自由向量只能是差值。**没有任何自由度剩下。**

**物理翻译（运动向量）**：

> **任何刚体运动，都可以看成"绕一条过你指定的点 $P$ 的轴转动" + "一个平移"。
> 一旦你指定了 $P$，这个拆法就唯一了。**

这非常直观：你想描述一个门的运动，可以说"绕合页转"，也可以说
"绕房间中心转 + 一个平移"——**两种说法都对**，取决于你把轴放在哪。
指定了轴的位置，剩下的平移就被确定了。

#### 分解二：要求平移**平行于**转轴 → 这就是螺旋（Chasles 定理）

> 原书原话：*"Any spatial vector, other than a free vector, can be expressed uniquely as
> the sum of a line vector and a **parallel** free vector."*

**问题**：不指定点，而是要求**自由向量平行于 $s$**。

**答案**：

$$
\underbrace{\begin{bmatrix}s\\ s_O-h\,s\end{bmatrix}}_{\text{线向量}}
\ +\
\underbrace{\begin{bmatrix}0\\ h\,s\end{bmatrix}}_{\text{平行的自由向量}},
\qquad
\boxed{h=\frac{s\cdot s_O}{s\cdot s}}
$$

**验证是线向量**：$s\cdot(s_O-hs)=s\cdot s_O-h(s\cdot s)=0$ ✓

**这条轴在哪**：过点 $\dfrac{s\times s_O}{s\cdot s}$（轴上离 $O$ 最近的点）。

**物理翻译**：这就是 **Chasles 定理**——

> **任何刚体瞬时运动 = 绕某条特定轴的转动 + 沿同一条轴的平移。**
> 这条轴叫**螺旋轴 (screw axis)**，$h$ 叫**螺距 (pitch)**。

**两个特例**：

- $h=0$（即 $s\cdot s_O=0$）：纯转动，$\hat s$ 本身就是线向量
- $s=0$：纯平移，$\hat s$ 是自由向量，螺距无定义（这就是原书说"other than a free vector"的原因）

#### 两条分解的关系

|              | 分解一                                        | 分解二                   |
| ------------ | --------------------------------------------- | ------------------------ |
| 你指定什么   | **点** $P$                                    | **平行**这个条件         |
| 自由度       | 你可以任选$P$，得到**无穷多组**分解，每组唯一 | 分解**唯一**（无可选项） |
| 自由向量方向 | 一般**不**平行于 $s$                          | **平行**于 $s$           |
| 名字         | ——                                          | Chasles 定理 / 螺旋分解  |

> 🔑 **一句话说清**：分解一是"**你选轴的位置**，平移随之确定"；
> 分解二是"**让平移平行于轴**，轴的位置随之确定"。
> 分解二是分解一中取 $P=\dfrac{s\times s_O}{s\cdot s}$ 的那一个特例。

**验证**：`python3 code/verify_ch02.py` 的 §2.5 部分逐条验证了以上四个断言。

### D. 大小 (Magnitude)

空间向量**没有欧氏意义上的大小**。但对四类特殊向量——转动、平移、线力、力偶——
可以定义有意义的大小。**只能同类比较**：转动的大小可以和另一个转动比，
但不能和平移或力比。

⚠️ 因此 Plücker 基向量（都是单位向量）**不仅依赖坐标系的位置和姿态，还依赖单位制**：
单位转动是"每单位时间 1 弧度"，单位平移是"每单位时间 1 个长度单位"。

---

## 2.6 标量积与对偶性 (Scalar Product)

$$
\hat m\cdot\hat f=\hat f\cdot\hat m \quad(\text{有定义}),
\qquad
\hat m\cdot\hat m,\ \hat f\cdot\hat f \quad(\textbf{无定义})
$$

物理意义：若 $\hat f$ 作用在速度为 $\hat m$ 的刚体上，则 $\hat m\cdot\hat f$ 是**功率**。

**对偶关系的三个实际后果**（原书列得很清楚）：

1. 在 $M^6$ 与 $F^6$ 上使用**对偶坐标**；
2. 力与运动遵循**不同的坐标变换规则**；
3. 存在**两个叉乘算子**，一个作用于运动、一个作用于力（§2.9）。

$D_O$ 与 $E_O$ 恰好构成对偶基（把 $d_{Ox}$ 当 $d_1$、$e_x$ 当 $e_1$，依此类推），于是

$$
\hat m\cdot\hat f=\underline{\hat m}^{\mathsf T}\underline{\hat f}
=n_O\cdot\omega+f\cdot v_O\tag{2.12}
$$

> 💡 **现在能解释 §2.3 那个"排列反了"的疑惑了**：
> 速度是 $[\omega;v_O]$、力是 $[n_O;f]$，看似不对称，
> 其实是为了让 $D_O$ 与 $E_O$ **恰好互易**——
> $d_{Ox}$（绕 $Ox$ 的单位转动）与 $e_x$（沿 $x$ 的单位力偶）配对得 1，
> $d_x$（沿 $x$ 的单位平移）与 $e_{Ox}$（过 $Ox$ 的单位力）配对得 1。
> 这样标量积才能写成裸的转置乘法，力矩配角速度、力配线速度，正好是功率的两项。

---

## 2.7 使用空间向量的规则 (Using Spatial Vectors)

原书 §2.7 给了一份**规则清单**，把空间向量与物理现象对应起来。
这份清单是全书的"公理表"，值得抄下来贴桌上：

| 规则                                 | 内容                                                                                                                                                                    |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **usage 用途**                       | $M^6$ 装：速度、加速度、无穷小位移、运动自由与约束的**方向**；$F^6$ 装：力、动量、冲量、力自由与约束的**方向**                                                          |
| **uniqueness 唯一性**                | $M^6$ 的元素与"刚体所有可能运动"一一对应；$F^6$ 的元素与"所有可能的力系(wrench)"一一对应                                                                                |
| **relative velocity 相对速度**       | $\hat v_{\text{rel}}=\hat v_2-\hat v_1$                                                                                                                                 |
| **rigid connection 刚性连接**        | 刚性连接的两个物体速度**相同**                                                                                                                                          |
| **summation of forces 力的叠加**   | 作用于同一刚体的$\hat f_1,\hat f_2$ 等效于 $\hat f_1+\hat f_2$                                                                                                          |
| **action and reaction 作用反作用**   | $B_1$ 对 $B_2$ 施力 $\hat f$，则 $B_2$ 对 $B_1$ 施力 $-\hat f$（牛顿第三定律）                                                                                          |
| **scalar product 标量积**            | 功率$=\hat f\cdot\hat v$                                                                                                                                                |
| **scaling 缩放**                     | $\alpha\hat f$ 作用于 $\beta\hat v$ 交付 $\alpha\beta$ 倍功率                                                                                                           |
| **acceleration 加速度**              | 空间加速度是空间速度的时间导数，是**真正的向量**，与速度遵循**同样的求和规则**：$\hat v_{\text{rel}}=\hat v_2-\hat v_1\Rightarrow\hat a_{\text{rel}}=\hat a_2-\hat a_1$ |
| **summation of inertias 惯性叠加** | 刚性连接的复合体：$I=I_1+I_2$                                                                                                                                           |
| **momentum 动量**                    | $\hat h=I\hat v$                                                                                                                                                        |
| **equation of motion 运动方程**      | $\hat f=\dfrac{d}{dt}(I\hat v)=I\hat a+\hat v\times^{*}I\hat v$                                                                                                         |

> 🔑 **"acceleration 是真正的向量、与速度同样求和"这一条是全书效率的来源。**
> 经典记法里加速度**不**满足简单求和（要冒出离心项和科氏项），
> 空间记法里满足。第 5 章 RNEA 只有几行，根子就在这。

### 例 2.3：运动学链（原书 §2.7）

这是全书第一次出现递推公式，第 5 章的 RNEA 就是它的推广。

关节 $i$ 从 body $i-1$ 连到 body $i$。定义**关节速度**为"到"物体相对"从"物体的速度：

$$
v_{Ji}=v_i-v_{i-1}\tag{2.14}
$$

单自由度关节：

$$
v_{Ji}=s_i\dot q_i\tag{2.15}
$$

$s_i$ 是**关节轴向量 (axis vector)**，$\dot q_i$ 是速度变量。
转动关节 ⟹ $s_i$ 是单位转动向量；移动关节 ⟹ $s_i$ 是单位平移向量。

于是有三种等价写法：

$$
v_i=v_{i-1}+s_i\dot q_i\quad(v_0=0)\tag{2.16}
$$

$$
v_i=\sum_{j=1}^{i}s_j\dot q_j\tag{2.17}
$$

$$
v_i=\begin{bmatrix}s_1&s_2&\cdots&s_i&0&\cdots&0\end{bmatrix}\dot q=J_i\,\dot q\tag{2.18}
$$

$J_i$ 是 body $i$ 的**物体雅可比 (body Jacobian)**，$6\times N$。

> 💡 **原书特意指出**：从计算角度看 2.16 比 2.17 高效，因为它**复用了前一步的结果**。
> 这一句就是全书"递推公式 (recurrence relation)"思想的第一次亮相——
> 第 5 章 §5.2 会正式讲它。

⚠️ **记号提醒**：本章用**小写** $s_i$ 表示单自由度关节的轴向量。
第 3、4 章推广到多自由度关节后改用**大写** $S_i$（$6\times n_i$ 矩阵）。

---

## 2.8 坐标变换 (Coordinate Transforms)

拆成"旋转"和"平移"两步。

### 纯旋转（共原点）

$$
{}^{B}X_{A}=\begin{bmatrix}E&0\\0&E\end{bmatrix}
\qquad
{}^{B}X_{A}^{*}=\begin{bmatrix}E&0\\0&E\end{bmatrix}\tag{2.19, 2.20}
$$

> 💡 **纯旋转时两者完全相同**。这是因为旋转矩阵正交，$E^{-\mathsf T}=E$。
> 只有平移才让 $X$ 与 $X^{*}$ 分道扬镳。

### 纯平移（同姿态，$r=\overrightarrow{OP}$）

$$
{}^{P}X_{O}=\begin{bmatrix}\mathbf 1&0\\-r\times&\mathbf 1\end{bmatrix}
\qquad
{}^{P}X_{O}^{*}=\begin{bmatrix}\mathbf 1&-r\times\\0&\mathbf 1\end{bmatrix}\tag{2.21, 2.22}
$$

### 一般变换

$$
{}^{B}X_{A}=\begin{bmatrix}E&0\\0&E\end{bmatrix}\begin{bmatrix}\mathbf 1&0\\-r\times&\mathbf 1\end{bmatrix}
=\begin{bmatrix}E&0\\-E\,r\times&E\end{bmatrix}\tag{2.24}
$$

$$
{}^{B}X_{A}^{*}=\begin{bmatrix}E&-E\,r\times\\0&E\end{bmatrix}\tag{2.25}
$$

$$
{}^{A}X_{B}=\begin{bmatrix}E^{\mathsf T}&0\\ r\times E^{\mathsf T}&E^{\mathsf T}\end{bmatrix}
\qquad
{}^{A}X_{B}^{*}=\begin{bmatrix}E^{\mathsf T}&r\times E^{\mathsf T}\\0&E^{\mathsf T}\end{bmatrix}\tag{2.26, 2.27}
$$

**与齐次变换的对比**（原书特意给出）：

$$
{}^{B}T_{A}=\begin{bmatrix}E&-E\,r\\0&1\end{bmatrix}\tag{2.28}
$$

> 💡 **注意 $4\times4$ 齐次变换与 $6\times6$ Plücker 变换的关系**：
> 都由同一对 $(E,r)$ 决定，但**作用对象不同**——
> 齐次变换作用于**位置**（描述位形），Plücker 变换作用于**速度/力**（描述微分量与作用量）。
> 二者互补，不是替代。

**函数记法**（原书表 2.2）：$\mathrm{rot}(E)$、$\mathrm{xlt}(r)$，于是 $ {}^{B}X_{A}=\mathrm{rot}(E)\,\mathrm{xlt}(r)$。

**方向约定**（原书特别提醒）：$\mathrm{rotx}(\theta)$ 读作
"把坐标系 $A$ 绕 $x$ 轴转 $\theta$ 使之与 $B$ 重合"，其值是 ${}^{B}X_{A}$。
**有些作者用相反的约定**，跨资料时务必确认。

---

## 2.9 空间叉乘 (Spatial Cross Products)

### 定义

设 $\hat m\in M^6$、$\hat f\in F^6$ **固连在一个速度为 $\hat v$ 的刚体上**，
除此之外不发生变化。定义两个算子：

$$
\dot{\hat m}=\hat v\times\hat m\tag{2.29}
$$

$$
\dot{\hat f}=\hat v\times^{*}\hat f\tag{2.30}
$$

> 💡 **叉乘算子的本质是"微分算子"**。3D 中 $\dot r=\omega\times r$ 就是这个意思；
> 空间向量把它推广了，只是需要两个版本。

### 从基向量的导数推出矩阵（原书图 2.5）

> 📝 **你在 p.32 批注了「参考图 2.5a」——这里展开讲。**

原书用一张**乘法表**定义这两个算子：**表中每个元素 = 顶行那个基向量的时间导数，
假定它以左列那个基向量的速度运动**。

**例（图 2.5a，第一行）**：坐标系以 $d_{Ox}$ 的速度运动
（= 绕直线 $Ox$ 以单位角速度转动）。经过 $\delta t$：

- 直线 $Ox$ **不动** ⟹ $d_{Ox}(t+\delta t)=d_{Ox}(t)$ ⟹ $\dot d_{Ox}=0$
- 直线 $Oy$ 在 $y$–$z$ 平面内转了 $\delta t$ 弧度 ⟹
  $d_{Oy}(t+\delta t)=d_{Oy}(t)+\delta t\,d_{Oz}$ ⟹ $\dot d_{Oy}=d_{Oz}$

**例（图 2.5b，第四行）**：坐标系以 $d_x$ 的速度运动（= 沿 $x$ 单位平移）。

- $Ox$ 不受影响 ⟹ $\dot d_{Ox}=0$
- $Oy$ 沿 $x$ 平移了 $\delta t$ ⟹ $d_{Oy}(t+\delta t)=d_{Oy}(t)+\delta t\,d_z$ ⟹ $\dot d_{Oy}=d_z$

> 🔑 **为什么"平移会改变一个转动基向量"**：
> $d_{Oy}$ 是"绕直线 $Oy$ 的单位转动"。把这条直线沿 $x$ 平移一点，
> 得到的是"绕新直线的单位转动"，它等于"绕原直线的单位转动 + 一个沿 $z$ 的平移"
> （回忆 §2.2 例 2.2 的 $d_{Py}=d_{Oy}+rd_z$）。
> **所以平移确实会改变转动基向量。** 这正是 $\hat v\times$ 矩阵左下块有 $v_O\times$ 的原因。

### 矩阵形式

$$
\hat v\times=\begin{bmatrix}\omega\times&0\\ v_O\times&\omega\times\end{bmatrix}\tag{2.31}
$$

$$
\hat v\times^{*}=\begin{bmatrix}\omega\times&v_O\times\\0&\omega\times\end{bmatrix}=-(\hat v\times)^{\mathsf T}\tag{2.32}
$$

分量展开：

$$
\begin{bmatrix}\omega\\v_O\end{bmatrix}\times\begin{bmatrix}m\\m_O\end{bmatrix}
=\begin{bmatrix}\omega\times m\\ \omega\times m_O+v_O\times m\end{bmatrix}\tag{2.33}
$$

$$
\begin{bmatrix}\omega\\v_O\end{bmatrix}\times^{*}\begin{bmatrix}n_O\\f\end{bmatrix}
=\begin{bmatrix}\omega\times n_O+v_O\times f\\ \omega\times f\end{bmatrix}\tag{2.34}
$$

**其他性质**：$\hat v\times\hat v=0$；$X(\hat v\times)X^{-1}=(X\hat v)\times$；
$(M^6,\times)$ 构成李代数 $\mathfrak{se}(3)$。

---

## 2.10 微分 (Differentiation)

### 在运动坐标系中求导

$$
\dot{\hat m}=\mathring{\hat m}+\hat v_A\times\hat m\tag{2.43}
$$

$$
\dot{\hat f}=\mathring{\hat f}+\hat v_A\times^{*}\hat f\tag{2.44}
$$

- $\mathring{(\cdot)}$（**圆圈记号**）：**表观导数 (apparent rate of change)**，
  即坐标系中速度为 $\hat v_A$ 的观察者看到的变化率——只对分量求导，把基当常量
- $\dot{(\cdot)}$（**点记号**）：绝对导数

对照 3D 的熟悉公式：$\dot u=\mathring u+\omega_A\times u$（式 2.38）。

⚠️ **运动量用 $\times$，力量用 $\times^{*}$，不能混。**

### 原书的重要提示

> *"In formulating the dynamics algorithms that appear in subsequent chapters,
> we generally **avoid** performing differentiation in moving coordinates by formulating
> the algorithms in **stationary** coordinate systems that happen to coincide with
> the moving ones at the current instant."*

🔑 **这句话极其重要，能省掉大量困惑**：
后面所有算法里的 $\hat v_i$、$\hat a_i$，都是在"此刻恰好与 body $i$ 坐标系重合的**静止**坐标系"
中表示的。所以算法里**看不到** $\mathring{(\cdot)}$，直接就是分量求导。
这就是为什么 RNEA 的伪代码里没有任何"表观导数"的痕迹。

### 例 2.4：变换矩阵的导数

$$
\frac{d}{dt}\,{}^{B}X_{A}={}^{B}(\hat v_A-\hat v_B)\times\,{}^{B}X_{A}\tag{2.45}
$$

即：**变换矩阵的导数 = 两坐标系相对速度的叉乘算子 × 该变换矩阵**。

---

## 2.11 空间加速度 (Acceleration) —— 全章最反直觉的概念

### 定义

$$
\hat a_O=\frac{d}{dt}\begin{bmatrix}\omega\\v_O\end{bmatrix}=\begin{bmatrix}\dot\omega\\ \dot v_O\end{bmatrix}\tag{2.46}
$$

### 原书的"悖论"及其化解

> 考虑一个绕空间中固定直线**匀角速转动**的刚体（原书图 2.6）。
> 它的空间速度是常量，所以**空间加速度为零**；
> 然而刚体上几乎每个点都在做圆周运动，**都在加速**。

**原书给的两点化解**：

1. **空间加速度是整个刚体的属性**，不是某个物质点的属性；
2. **$v_O$ 不是某个特定物质点的速度，而是"物质点流过 $O$ 的流量"的度量。**
   因此 $\dot v_O$ 不是某点的加速度，而是**流量变化率**。
   流量恒定 ⟹ $\dot v_O=0$，哪怕每个点都在加速。

> 🔑 **第 2 点是理解本章的钥匙。** 把 $v_O$ 想成"流场在 $O$ 处的取值"，
> 一切就顺了：定常流场的时间导数为零，但流线是弯的，
> 跟着流走的粒子（物质点）当然有加速度。**这是欧拉描述 vs 拉格朗日描述的区别。**

### 与经典加速度的定量关系

引入 $O'$ = 此刻与 $O$ 重合的**物质点**，$r(t)$ = $O'$ 的位置。
则 $\dot r$ 是 $O'$ 的速度、$\ddot r$ 是它的加速度。原书推出：

$$
\dot v_O=\ddot r-\omega\times\dot r\tag{2.47}
$$

$$
\hat a_O=\begin{bmatrix}\dot\omega\\ \ddot r-\omega\times\dot r\end{bmatrix}\tag{2.48}
$$

**经典加速度 (classical acceleration)**——原书正式定义为一个 6D 向量：

$$
\hat a'_O=\begin{bmatrix}\dot\omega\\ \ddot r\end{bmatrix}\tag{2.49}
$$

两者的关系：

$$
\hat a_O=\hat a'_O+\begin{bmatrix}0\\ \dot r\end{bmatrix}\times\hat v_O\tag{2.50}
$$

$$
\boxed{\ \text{经典加速度的线分量}=\text{空间加速度的线分量}+\omega\times v_O\ }
$$

**原书对 $\hat a'$ 的解释**：经典加速度是空间速度在"原点位于 $O'$（而非 $O$）
的 Plücker 坐标系"中的**表观导数**——即该坐标系带有纯线速度 $\dot r$。

> 💡 **一句话区分**：
> **空间加速度**在**固定**坐标系里求导（原点不动）；
> **经典加速度**在**跟着物质点平移**的坐标系里求导（原点跟着 $O'$ 走）。

**数值例子**（`code/verify_ch02.py`）：刚体绕 $z$ 轴以 $\omega=2.1$ 匀速转动，
点 $P=(0.5,0,0)$：

```
空间加速度线分量 = [0, 0, 0]          <- 恒为零
P 点经典加速度    = [-2.205, 0, 0]    <- 正是向心加速度 -w^2 r = -2.1^2 * 0.5
```

### 为什么值得用这个反直觉的定义

**因为空间加速度是真正的向量，与速度遵循同样的求和与坐标变换规则**：

$$
\hat v_{\text{rel}}=\hat v_2-\hat v_1
\quad\xrightarrow{\ \text{直接求导}\ }\quad
\hat a_{\text{rel}}=\hat a_2-\hat a_1
$$

对运动学链（原书式 2.55，由 2.16 求导而来）：

$$
\boxed{\ a_i=a_{i-1}+s_i\ddot q_i+v_i\times s_i\dot q_i\ }\tag{2.55}
$$

⚠️ **注意交叉项只有一个 $v_i\times s_i\dot q_i$，且没有系数 2。**
经典记法里会冒出离心项 $\omega\times(\omega\times r)$ 和科氏项 $2\omega\times v$，
且随连杆串联层层嵌套、展开式爆炸增长。

> 🔑 **"消失"的离心项和那个因子 2 并没有真的消失**，
> 它们被**吸收进了空间加速度的定义**（式 2.48 的 $-\omega\times\dot r$），
> 最终会在运动方程的 $\hat v\times^{*}I\hat v$ 项里重新出现。
> **这就是 RNEA 只有几行的根本原因。**

式 2.56、2.57 给出另外两种写法：

$$
a_i=\sum_{j=1}^{i}s_j\ddot q_j+\sum_{j=1}^{i}\sum_{k=1}^{j-1}(s_k\times s_j)\,\dot q_j\dot q_k\tag{2.56}
$$

$$
a_i=J_i\ddot q+\dot J_i\dot q\tag{2.57}
$$

> 式 2.56 显式地把"速度乘积项"暴露出来了：双重求和那部分就是科氏/离心力的来源，
> 共 $O(i^2)$ 项——**这正是不用递推就会爆炸的证据**。

---

## 2.12 空间动量 (Momentum)

质心 $C$ 处：线动量 $h=m v_C$，固有角动量 $h_C=\bar I_C\omega$。

$$
\hat h_C=\begin{bmatrix}\bar I_C\,\omega\\ m\,v_C\end{bmatrix}\tag{2.59}
$$

关于任意点 $O$ 的动量矩 = 固有角动量 + 线动量对该点的矩：

$$
h_O=h_C+\overrightarrow{OC}\times h\tag{2.58}
$$

$$
\hat h_O=\begin{bmatrix}\mathbf 1&\overrightarrow{OC}\times\\0&\mathbf 1\end{bmatrix}\hat h_C\tag{2.60}
$$

> 💡 **注意式 2.58 与式 2.6（力矩的搬移）形式完全相同**——
> 原书明说"essentially the same formula"。
> 这不是巧合：**动量是力空间的元素**，遵循力向量的一切规则。
> 线动量甚至是个**线向量**，其作用线过质心。

**动量按力向量变换**：${}^{B}\hat h={}^{B}X_A^{*}\,{}^{A}\hat h$。

**动能**：$T=\tfrac12\hat v\cdot\hat h=\tfrac12\hat v^{\mathsf T}I\hat v$

---

## 2.13 空间惯性 (Inertia)

### A. 矩阵形式

$$
\hat h=I\hat v\tag{2.61}
$$

质心处（$c=0$）：

$$
I_C=\begin{bmatrix}\bar I_C&0\\0&m\mathbf 1\end{bmatrix}\tag{2.62}
$$

一般点 $O$（$c=\overrightarrow{OC}$）：

$$
\boxed{\ I_O=\begin{bmatrix}\bar I_C+m\,c\times c\times^{\mathsf T}&m\,c\times\\
m\,c\times^{\mathsf T}&m\,\mathbf 1\end{bmatrix}\ }\tag{2.63}
$$

其中 $\bar I_C+m\,c\times c\times^{\mathsf T}$ 就是刚体**绕 $O$** 的转动惯量
（**平行轴定理的 6D 形式**；注意 $c\times^{\mathsf T}=-c\times$，
所以也可写成 $\bar I_C-m\,c\times c\times$）。

**性质**：$\bar I_C$ 对称 ⟹ $I_C,I_O$ 对称；$m>0$ 且 $\bar I_C$ 正定 ⟹ $I_C,I_O$ 正定。

### 2.13-B 惯性的并矢表示：$g_i$ 到底是什么

> 📌 **你在 p.40 批注：「$g_i$ 的实际含义是什么，需要考虑？？？」**

原书的说法：

$$
I=\sum_{i=1}^{6}g_i\,g_i\cdot\qquad(g_i\in F^6)\tag{2.64}
$$

#### 第一层：这在数学上说什么

回忆 §2.1：并矢 $f\,g\cdot$ 把 $\hat m$ 映到 $f(g\cdot\hat m)$，
坐标形式是外积 $\underline f\,\underline g^{\mathsf T}$（秩 1）。

若 $f=g$，并矢是**对称的**，矩阵也对称。

所以式 2.64 说的就是：

$$
\boxed{\ I=\sum_{i=1}^{6}\underline g_i\,\underline g_i^{\mathsf T}
=G\,G^{\mathsf T}\quad(G=[\underline g_1\ \cdots\ \underline g_6])\ }
$$

**这就是"对称正定矩阵可以分解成 $GG^{\mathsf T}$"**——线性代数的老结论。
构造 $G$ 至少有两条现成路子：

| 方法                                             | $g_i$ 是什么                |
| ------------------------------------------------ | --------------------------- |
| **Cholesky** $I=LL^{\mathsf T}$                  | $g_i$ = $L$ 的第 $i$ 列     |
| **特征分解** $I=\sum\lambda_iu_iu_i^{\mathsf T}$ | $g_i=\sqrt{\lambda_i}\,u_i$ |

#### 第二层：⚠️ $g_i$ 不唯一

**这是理解上的关键。** 对任意正交矩阵 $Q$，$G$ 与 $GQ$ 给出**同一个** $I$：

$$
(GQ)(GQ)^{\mathsf T}=GQQ^{\mathsf T}G^{\mathsf T}=GG^{\mathsf T}=I
$$

`verify_ch02.py` 里对同一个 $I$ 用两种方法各求一组：

```
Cholesky 的 g_1 = [ 1.141 -0.520  0.932 ...]
特征分解的 g_1 = [-0.295  0.164 -0.061 ...]      <- 完全不同，但都对
```

> 🔑 **所以"$g_i$ 是什么"这个问题，问的方式需要调整**：
> 不存在"那个 $g_i$"。存在的是**"一组 $g_i$"**，有无穷多组。
> **原书需要的只是它们的存在性，不是具体取值。**

#### 第三层：物理含义——它们是过质点的单位力

虽然不唯一，但有一组**特别有物理意义**的：把刚体看成**质点系**。

设质点 $k$ 质量 $m_k$、位置 $p_k$。该点线速度由空间速度给出：

$$
v_{P_k}=v_O+\omega\times p_k=\underbrace{\begin{bmatrix}-p_k\times&\mathbf 1\end{bmatrix}}_{J_k\ (3\times6)}\hat v
$$

动能：

$$
T=\sum_k\tfrac12 m_k\|v_{P_k}\|^2=\tfrac12\hat v^{\mathsf T}\Big(\underbrace{\sum_k m_kJ_k^{\mathsf T}J_k}_{=\,I}\Big)\hat v
$$

把每个 $m_kJ_k^{\mathsf T}J_k$ 按 $J_k^{\mathsf T}$ 的三列拆开。$J_k^{\mathsf T}$ 的第 $a$ 列是

$$
\begin{bmatrix}p_k\times e_a\\ e_a\end{bmatrix}
$$

**这恰好是"沿方向 $e_a$、过点 $p_k$ 的单位线力"的 Plücker 坐标！**
（力 $f=e_a$，关于 $O$ 的矩 $n_O=p_k\times e_a$。）

$$
\boxed{\ g=\sqrt{m_k}\times\big(\text{过质点 }p_k\text{、沿 }e_a\text{ 的单位力}\big)\ }
$$

> 🔑 **回答你的问题**：$g_i$ **是力向量**（$g_i\in F^6$，原书写明了）。
> 物理上可以取成"**过某个质点的单位线力，按 $\sqrt{质量}$ 缩放**"。
> 一个刚体拆成若干质点，每个质点贡献 3 个这样的 $g$；
> 由于秩最多 6，总可以压缩成 6 个。

**能量解释（最简洁的一个）**：

$$
T=\tfrac12\hat v^{\mathsf T}I\hat v=\tfrac12\sum_{i=1}^{6}(g_i\cdot\hat v)^2
$$

> **$g_i$ 是一组"把动能配成平方和"的力向量。**
> 每个 $g_i\cdot\hat v$ 是一个标量（量纲是 $\sqrt{能量}$），动能就是它们的平方和。

#### 第四层：原书为什么要引入这个表示——为了求 $\dot I$

**这才是式 2.64 的真正用途。** 因为 $g_i$ 可以取成**固连在刚体上**的向量，
它们的导数就是 $\dot g_i=\hat v\times^{*}g_i$（式 2.30）。于是

$$
\dot I=\sum_i\big(\dot g_i\,g_i\cdot+g_i\,\dot g_i\cdot\big)
=\hat v\times^{*}I-I\,\hat v\times\tag{2.65}
$$

$$
\boxed{\ \dot I=\hat v\times^{*}I-I\,\hat v\times\ }
$$

**这一步是通往运动方程的桥**（§2.14 立刻要用）。

> 🔑 **总结你的困惑**：$g_i$ 是一个**构造性工具**，
> 引入它是为了能对 $I$ 求导。个体不唯一、也不重要，
> 重要的是"$I$ 可以写成一堆固连于刚体的力向量的对称并矢之和"这个**事实**，
> 因为固连向量的导数有现成公式。
> **不必纠结某个具体的 $g_i$ 是哪个力——原书也不关心。**

`verify_ch02.py` 验证了：两种分解都成立且不同、质点系构造给出合法刚体惯性、
以及由并矢表示确实推出式 2.65。

### 2.13-C 为什么恰好是 10 个参数

> 📌 **你在 p.42 批注「why???」**

原书原文：*"Ten parameters are required to define a spatial rigid-body inertia.
More general kinds of inertia, such as the articulated-body inertias in Chapter 7,
require up to 21 parameters, which is the maximum number of independent values
in a symmetric 6×6 matrix."*

#### 正向数：物理参数

$$
\underbrace{m}_{1}\ +\ \underbrace{c}_{3}\ +\ \underbrace{\bar I_C}_{6\ (\text{对称 }3\times3)}\ =\ \boxed{10}
$$

#### 反向数：从 21 减出来（更能说明"为什么"）

一般的对称 $6\times6$ 矩阵有 $\frac{6\times7}{2}=21$ 个独立元。
但**刚体**惯性的三个分块各自被强约束：

| 分块           | 一般对称阵 | 刚体惯性                                                     | 省下   |
| -------------- | ---------- | ------------------------------------------------------------ | ------ |
| 左上$3\times3$ | 6（对称）  | $\bar I_C+mc\times c\times^{\mathsf T}$，仍是任意对称阵 → 6 | 0      |
| 右上$3\times3$ | 9（任意）  | $m\,c\times$，**必须反对称** → 3                            | **6**  |
| 右下$3\times3$ | 6（对称）  | $m\,\mathbf 1$，**必须是标量乘单位阵** → 1                  | **5**  |
| **合计**       | **21**     | **10**                                                       | **11** |

$$
21-11=10\quad\checkmark
$$

> 🔑 **一句话回答"why"**：因为刚体惯性的 $6\times6$ 矩阵**不是任意的对称阵**——
> 它的右下块必须是 $m\mathbf 1$（**5 个约束**），右上块必须反对称（**6 个约束**）。
> 这 11 个约束把 21 砍到 10。

#### 为什么铰接体惯性需要满 21 个

第 7 章的**铰接体惯性 $I^A$** 描述的是"一整棵子树、关节可自由活动"时的等效惯性。
它**不再是单个刚体的惯性**，因此上面两个结构约束都不成立：

```
一般（铰接体）惯性的右下块：
[[20.06,  3.29,  5.25],
 [ 3.29,  7.79, -2.32],      <- 不是 m*1
 [ 5.25, -2.32,  5.91]]
```

**物理直觉**：铰接体在不同方向上"表观质量"不同——
沿关节允许运动的方向推它很轻，沿被约束的方向推它很重。
单个刚体没有这种方向性（各向同性的 $m\mathbf 1$），所以参数少。

> 💡 这条对比是理解第 7 章的重要铺垫：
> **$I^A$ 比 $I$ "更一般"，代价是参数从 10 涨到 21。**

### D. 惯性的坐标变换

由并矢表示，$ {}^{B}g_i={}^{B}X_A^{*}\,{}^{A}g_i$，代入即得

$$
\boxed{\ {}^{B}I={}^{B}X_{A}^{*}\ {}^{A}I\ {}^{A}X_{B}\ }\tag{2.66}
$$

**怎么记**：$I$ 是 $M^6\to F^6$ 的映射，所以右边先把运动量转过去（$X$），
左边再把结果的力量转回来（$X^{*}$）。**三明治结构，两侧变换方向相反。**

### E. 四类空间并矢张量（原书表 2.5）

| 映射         | 例子               | 并矢形式    | 变换公式              | 类型              |
| ------------ | ------------------ | ----------- | --------------------- | ----------------- |
| $M^6\to F^6$ | 惯性$I$            | $f\,f\cdot$ | $X^{*}I\,X^{-1}$      | 合同 (congruence) |
| $F^6\to M^6$ | 逆惯性$\Phi$       | $m\,m\cdot$ | $X\Phi X^{*-1}$       | 合同              |
| $M^6\to M^6$ | $\hat v\times$     | $m\,f\cdot$ | $X(\cdot)X^{-1}$      | 相似 (similarity) |
| $F^6\to F^6$ | $\hat v\times^{*}$ | $f\,m\cdot$ | $X^{*}(\cdot)X^{*-1}$ | 相似              |

> 💡 **规则**：把并矢张量看成"从 From 空间映到 To 空间"，
> 则并矢的**第一个**向量属于 **To** 空间，**第二个**向量属于 **From 的对偶**空间。
>
> **相似变换**保秩和特征值，**不保**对称性与正定性；
> **合同变换**保秩、对称性、正定性，**不保**特征值。
> 所以"$I$ 正定"在换坐标系后仍然成立（合同），
> 而"$\hat v\times$ 的特征值"在换坐标系后不变（相似）。

### F. 其他性质

- **10 个参数**（见 2.13-C）
- **可加性**：复合刚体的惯性 = 各部分惯性之和（**第 6 章 CRBA 的直接依据**）
- **动能**：$T=\tfrac12\hat v^{\mathsf T}I\hat v$

---

## 2.14 运动方程 (Equation of Motion)

$$
\hat f=\frac{d}{dt}(I\hat v)=I\hat a+(\hat v\times^{*}I-I\hat v\times)\hat v
=I\hat a+\hat v\times^{*}I\hat v\tag{2.68}
$$

（用到 $\hat v\times\hat v=0$ 和式 2.65。）

$$
\boxed{\ \hat f=I\hat a+\hat v\times^{*}I\hat v\ }
$$

### 偏置力形式（重要，后面章节反复用）

原书常写成**非齐次线性方程**：

$$
\hat f=I\hat a+\hat p\tag{2.69}
$$

$\hat p$ 称为**偏置力 (bias force)**，定义为"**产生零加速度所需的力**"。

> 💡 **原书特别强调式 2.69 比 2.68 灵活在哪**：
> 式 2.68 里的 $\hat f$ **永远**是净力；而式 2.69 允许你把**已知的**力分量
> 从 $\hat f$ 挪到 $\hat p$ 里去。
>
> **例**：净力 = 未知力 $\hat f_u$ + 已知重力 $\hat f_g$。可以选择
>
> $$
> \hat f_u=I\hat a+\hat p,\qquad \hat p=\hat v\times^{*}I\hat v-\hat f_g
> $$
>
> **第 5 章的重力技巧、第 7 章 ABA 的 $\hat p^A$，用的都是这个自由度。**
> 认出这一点，后面 $\hat p^A$ 里为什么会混进外力和关节力就不奇怪了。

### 例 2.6：展开回牛顿 + 欧拉方程

在**质心系**（$c=0$）展开式 2.68：

$$
\begin{bmatrix}n_C\\f\end{bmatrix}
=\begin{bmatrix}\bar I_C&0\\0&m\mathbf 1\end{bmatrix}
\begin{bmatrix}\dot\omega\\ \ddot c-\omega\times v_C\end{bmatrix}
+\begin{bmatrix}\omega\times&v_C\times\\0&\omega\times\end{bmatrix}
\begin{bmatrix}\bar I_C&0\\0&m\mathbf 1\end{bmatrix}
\begin{bmatrix}\omega\\ v_C\end{bmatrix}
=\begin{bmatrix}\bar I_C\dot\omega+\omega\times\bar I_C\omega\\ m\,\ddot c\end{bmatrix}
$$

于是恰好恢复牛顿方程（左，式 2.70）与欧拉方程（右，式 2.71）：

$$
f=m\,\ddot c
\qquad\qquad
n_C=\bar I_C\dot\omega+\omega\times\bar I_C\omega\tag{2.70, 2.71}
$$

> 🔑 **一条 6D 方程 = 牛顿方程 + 欧拉方程。**
> 注意推导中用了式 2.48 把空间加速度换成 $\ddot c$（经典加速度）——
> **这正是"离心项被吸收"的具体位置**：
> $v_C\times m\mathbf 1 v_C=m\,v_C\times v_C=0$ 消掉一项，
> 而 $\hat a$ 的线分量 $\ddot c-\omega\times v_C$ 加上 $m\omega\times v_C$ 又还原成 $m\ddot c$。

`verify_ch02.py` 逐项验证了这个展开。

---

## 2.15 逆惯性 (Inverse Inertia)

$$
\hat a=\Phi\hat f+\hat b\tag{2.72}
$$

- $\Phi=I^{-1}$：**逆惯性**，$F^6\to M^6$ 的对称并矢张量
- $\hat b=-\Phi\hat p$：**偏置加速度 (bias acceleration)**，即 $\hat f=0$ 时的加速度

### 为什么需要它：可以描述受约束刚体

> 原书：*"The advantage of Eq. 2.72 over Eq. 2.69 is that it can be applied to a
> **constrained** rigid body."*

**关键性质**：

$$
\boxed{\ \mathrm{range}(\Phi)=\text{该刚体自由运动的子空间},\qquad
\mathrm{rank}(\Phi)=\text{该刚体的运动自由度数}\ }
$$

于是：

| 情形                         | $\Phi$     | $I$        |
| ---------------------------- | ---------- | ---------- |
| 无约束刚体（6 自由度）       | 满秩，可逆 | 存在       |
| **受约束刚体（< 6 自由度）** | **奇异**   | **不存在** |

> 🔑 **这就是引入 $\Phi$ 的理由**：受约束刚体（比如被固定在一根轴上的连杆）
> 根本**没有**惯性矩阵——因为某些方向上"再大的力也产生不了加速度"，
> 对应 $I$ 要取无穷大。但它的**逆惯性**完全良定义，只是秩不满。
> 用 $\hat a=\Phi\hat f+\hat b$ 就能统一处理有约束和无约束的情形。

公式（原书 2.73、2.74）就是 2.62、2.63 的逆：

$$
\Phi_C=\begin{bmatrix}\bar I_C^{-1}&0\\0&\tfrac1m\mathbf 1\end{bmatrix}
$$

---

## 2.16 平面向量 (Planar Vectors)

若所有刚体都平行于同一平面运动，可以用**平面向量代数**——空间向量代数的"瘦身版"。
平面刚体有 3 个自由度，所以平面向量属于 $M^3$ 与 $F^3$。

**推导方法**：从空间向量出发限制到 $x$–$y$ 平面，
**删掉第 1、2、6 个基向量及对应坐标**：

$$
\begin{bmatrix}\omega_x\\ \omega_y\\ \omega_z\\ v_{Ox}\\ v_{Oy}\\ v_{Oz}\end{bmatrix}
\Longrightarrow
\begin{bmatrix}\omega\\ v_{Ox}\\ v_{Oy}\end{bmatrix}
\qquad
\begin{bmatrix}n_{Ox}\\ n_{Oy}\\ n_{Oz}\\ f_x\\ f_y\\ f_z\end{bmatrix}
\Longrightarrow
\begin{bmatrix}n_O\\ f_x\\ f_y\end{bmatrix}
$$

（角分量的 $z$ 下标被省略。）

**叉乘**：

$$
\begin{bmatrix}\omega\\v_{Ox}\\v_{Oy}\end{bmatrix}\times
=\begin{bmatrix}0&0&0\\ v_{Oy}&0&-\omega\\ -v_{Ox}&\omega&0\end{bmatrix}\tag{2.75}
$$

**坐标变换**（原点移到 $(x,y)$ 后绕新原点转 $\theta$，$c=\cos\theta$、$s=\sin\theta$）：

$$
X=\begin{bmatrix}1&0&0\\ sx-cy&c&s\\ cx+sy&-s&c\end{bmatrix}\tag{2.77}
$$

**惯性**（$m$ 质量，$(c_x,c_y)$ 质心，$I_C$/$I_O$ 绕质心/原点的转动惯量）：

$$
I=\begin{bmatrix}I_C+m(c_x^2+c_y^2)&-mc_y&mc_x\\ -mc_y&m&0\\ mc_x&0&m\end{bmatrix}
=\begin{bmatrix}I_O&-mc_y&mc_x\\ -mc_y&m&0\\ mc_x&0&m\end{bmatrix}\tag{2.78}
$$

**平面刚体惯性需要 4 个参数**（$m,c_x,c_y,I_C$）——
对照空间的 10 个，用 2.13-C 的方法数一遍是很好的练习。

**继承的性质**：$\hat v\times^{*}=-(\hat v\times)^{\mathsf T}$、$X^{*}=X^{-\mathsf T}$ 都照样成立。

> 💡 **实用价值**：很多机构（平面连杆、轮式移动机器人的平面模型）本来就是平面的。
> 用 $3\times3$ 代替 $6\times6$，运算量降到约 1/4，且不会出现"数值上应为零的
> 面外分量因舍入而漂移"的问题。

---

## 2.17 延伸阅读：与其他 6D 理论的关系

原书 §2.17 明确给出了对应关系：

| 本书                                | 对应                                                |
| ----------------------------------- | --------------------------------------------------- |
| $M^6$                               | 李代数$\mathfrak{se}(3)$                            |
| $F^6$                               | 其对偶$\mathfrak{se}^{*}(3)$                        |
| $M^3,F^3$（平面）                   | $\mathfrak{se}(2),\mathfrak{se}^{*}(2)$             |
| 空间速度 / 空间力                   | 螺旋理论的**twist / wrench**                        |
| $\hat v\times$ / $\hat v\times^{*}$ | $\mathrm{ad}_{\hat v}$ / $\mathrm{ad}^{*}_{\hat v}$ |
| $X$ / $X^{*}$                       | $\mathrm{Ad}_g$ / $\mathrm{Ad}_g^{-\mathsf T}$      |

$\mathfrak{se}(3)$ 是特殊欧氏群 $SE(3)$ 对应的李代数
（形式上：李代数 = 李群在单位元处的切空间）。

> ⚠️ **原书指出的重要差别**：**motor 代数**（von Mises、Brand）
> 把所有量放在**同一个**向量空间里，而空间向量刻意分成两个。
> Featherstone 1987 年的旧书更接近 motor 代数，本书改了。
> **所以读老文献时要注意这个差异。**

---

## 易错点汇总

1. **$v_O$ 不是原点的速度**，而是物质点流过 $O$ 的流量度量（§2.11）。
2. **分量排列**：本书 (角, 线)；部分文献 (线, 角)。跨库移植的头号坑。
3. **$X^{*}\ne X$**，且 $X^{-1}\ne X^{\mathsf T}$（只有纯旋转时 $X^{*}=X$）。
4. **空间加速度 ≠ 物质点加速度**，差 $\omega\times v_O$；原书把后者正式命名为
   "经典加速度 $\hat a'$"并给了式 2.50。
5. **$\times$ 与 $\times^{*}$ 不能混**，判据只有一条：作用对象是运动还是力。
6. **$I$ 的变换是 $X^{*}IX$**（合同），不是相似变换。
7. **不同坐标系的空间向量不能直接相加**——6 个数看着能加，物理上无意义。
8. **$g_i$ 不唯一**（§2.13-B），不要试图找"那个" $g_i$。
9. **线向量的两种分解回答不同问题**（§2.5-C），不要混。
10. **受约束刚体没有 $I$，只有 $\Phi$**（§2.15）。

---

## 自测清单

读完本章应该能不看书回答：

- [ ]  为什么 $M^6$ 上没有内积？这导致了哪三个后果？
- [ ]  对偶基的互易条件是什么？正交归一基与它是什么关系？
- [ ]  写出 $X$、$X^{*}$、$\hat v\times$、$\hat v\times^{*}$、$I_O$ 的矩阵形式
- [ ]  为什么 $X^{*}=X^{-\mathsf T}$？（给出代数的和物理的两种理由）
- [ ]  线向量的判据是什么？为什么是 5 个数？
- [ ]  说清 §2.5 两种分解各自回答什么问题
- [ ]  一个匀速转动的刚体，空间加速度是多少？其上一点的经典加速度是多少？
- [ ]  式 2.55 的交叉项为什么没有系数 2？那个 2 去哪了？
- [ ]  $I=\sum g_ig_i\cdot$ 里 $g_i$ 属于哪个空间？唯一吗？引入它是为了做什么？
- [ ]  为什么刚体惯性是 10 个参数而铰接体惯性是 21 个？
- [ ]  $\Phi$ 奇异意味着什么？

---

## 与其他章的联系

- ← 第 1 章：兑现"用 6D 向量"的承诺
- → 第 3 章：把单刚体的 $\hat f=I\hat a+\hat v\times^{*}I\hat v$ 推广到多体系统；
  §2.15 的 $\Phi$ 在那里用于受约束刚体
- → 第 4 章：$s_i$ 推广成 $S_i$，$X_T$/$X_J$ 的拆分
- → 第 5 章：式 2.55 就是 RNEA 外推的核心；式 2.69 的偏置力自由度 = 重力技巧
- → 第 6 章：惯性可加性 = CRBA 的依据
- → 第 7 章：$I$ 的 10 参数 vs $I^A$ 的 21 参数
- → 附录 A：空间向量运算的高效实现

---

## ✍️ 我的理解

<!-- 建议自测：不看书默写 X、X*、v×、v×*、I 的矩阵形式与运动方程 -->

## ❓ 疑问与待办

- [ ]  手推式 2.63 的两行（`notes/derivations.md` D2 已给出，自己再推一遍）
- [ ]  手推式 2.65 $\dot I=\hat v\times^{*}I-I\hat v\times$（用并矢表示）
- [ ]  用 2.13-C 的数法，自己数出平面刚体惯性的 4 个参数
- [ ]  读附录 A，对照 §2.8 的紧凑实现
- [ ]  例 2.2 推广到一般位置的 $P$（原书说"easily extended"，自己做一遍）

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.7–38 逐节核对。发现不符请记录在此并修正正文 -->
