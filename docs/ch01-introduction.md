# 第 1 章 绪论 (Introduction)

> **原书 pp. 1–6**，共 5 节。篇幅很短，作用是"说明书"。

> **一句话概括**：交代本书要算什么（正/逆动力学）、用什么工具（6D 空间向量）、
> 用什么组织方式（基于模型的算法），以及全书怎么读。

## 原书节次对照

| 节 | 原书标题 | 页 |
|---|---|---|
| 1.1 | Dynamics Algorithms | 1 |
| 1.2 | Spatial Vectors | 3 |
| 1.3 | Units and Notation | 4 |
| 1.4 | Readers' Guide | 5 |
| 1.5 | Further Reading | 6 |

**读法**：快速通读。读完全书后回头再看一遍，会发现这章其实是全书摘要。

---

## 1.1 动力学算法 (Dynamics Algorithms)

### 两个（后来是三个）基本计算

原书开宗明义定义了两个核心计算：

| 计算 | 定义 | 主要用途 |
|---|---|---|
| **正动力学 (forward dynamics)** | 给定作用力，求系统的**加速度响应** | 主要用于**仿真** |
| **逆动力学 (inverse dynamics)** | 给定期望加速度，求**必须施加的力** | 运动控制、轨迹规划、机械设计，**以及作为正动力学计算的一个组成部分** |

第 9 章引入第三种：

| **混合动力学 (hybrid dynamics)** | 部分加速度、部分力已知，求其余 |

> 💡 **注意"as a component in a forward-dynamics calculation"这一句**：
> 逆动力学不只是"控制用的"，它还是第 6 章 FD 路线的第一步（用来算 $C$）。
> 这解释了为什么第 5 章要排在第 6、7 章之前。

### 标准形式（式 1.1）

$$
\boxed{\ \tau=H(q)\,\ddot q+C(q,\dot q)\ }\tag{1.1}
$$

| 符号 | 含义 |
|---|---|
| $q,\dot q,\ddot q$ | 位置、速度、加速度变量向量 |
| $\tau$ | 施加力向量 |
| $H$ | **惯性项矩阵**，写成 $H(q)$ 表示它是 $q$ 的函数 |
| $C$ | **力项向量**：科氏力、离心力、**重力**，以及 $\tau$ 之外作用于系统的**任何其他力** |

⚠️ **本书的 $C$ 是向量，不是矩阵**，而且**包含重力和其他外力**。
很多机器人学教材写成 $H\ddot q+C(q,\dot q)\dot q+g(q)=\tau$，那里的 $C$ 是矩阵。
**跨教材对照时这是头号混淆点。**

**术语**：$H$ 与 $C$ 是方程的**系数 (coefficients)**，$\tau$ 与 $\ddot q$ 是**变量 (variables)**。
典型情形是其中一个变量已知、另一个未知。

### 系统模型 (system model)

原书一个重要的观念澄清：

> 习惯上写 $H(q)$ 和 $C(q,\dot q)$，但**更准确的写法是 $H(model,q)$ 和 $C(model,q,\dot q)$**。

`model` 指描述某个具体刚体系统的一整套数据：刚体和关节的数目、连接方式、
每个部件的全部参数（惯性参数、几何参数等）。

**关键区分**：

| | 描述什么 | 例子 |
|---|---|---|
| **系统模型 (system model)** | **系统本身** | `model` 数据结构 |
| **数学模型 (mathematical model)** | 系统**行为的某个方面** | 式 1.1 |

### 基于模型的算法 (model-based algorithm)

把两个计算封装成函数：

$$
\ddot q=\mathrm{FD}(model,q,\dot q,\tau)\tag{1.2}
$$

$$
\tau=\mathrm{ID}(model,q,\dot q,\ddot q)\tag{1.3}
$$

显然 $\mathrm{FD}=H^{-1}(\tau-C)$、$\mathrm{ID}=H\ddot q+C$。

> **原书强调这两个式子的意义不在于结果，而在于它们清楚地显示了输入和输出——
> 尤其是 `model` 在两种情况下都是输入。**

**这样做的最大好处**（原书原话）：

> *"a single piece of computer code can be written, tested, documented, and so on,
> to calculate the dynamics of **any** rigid-body system in a broad class."*

一份代码，写一次、测一次、文档一次，适用于一整类机构。

### 两大系统类别

| 类别 | 定义 | 难度 |
|---|---|---|
| **运动学树 (kinematic tree)** | 不含运动学回路的刚体系统 | 较易 |
| **闭环系统 (closed-loop system)** | 不是运动学树的刚体系统 | **显著更难** |

（精确定义在第 4 章。）

**算法的二维分类**：按"做什么计算"× "适用什么系统类别"。
本书主体 = 运动学树和闭环系统的正/逆动力学算法。

---

## 1.2 空间向量 (Spatial Vectors)

### 动机：一条方程 vs 两条方程

3D 记法下，单刚体运动方程要写**两条**（式 1.4）：

$$
f=m\,a_C\qquad\text{和}\qquad n_C=I\dot\omega+\omega\times I\omega\tag{1.4}
$$

空间记法下写**一条**（式 1.5）：

$$
\hat f=I\hat a+\hat v\times^{*}I\hat v\tag{1.5}
$$

⚠️ **原书提醒**：式 1.5 与 1.4 有**符号冲突**（$f$、$I$、$a$ 重名），
解决办法是给空间量**戴帽子**（$\hat f$）。第 2 章会正式引入这个约定。

### 另一个例子：惯性叠加

两个刚体刚性连接成一体：

$$
I_{\text{new}}=I_1+I_2
$$

> **这一条替代了 3D 做法中的三条**：算新质量、算新质心、算绕新质心的转动惯量。
> **第 6 章 CRBA 的全部效率就建立在这一条上。**

### 收益的量化

原书给了两个数字：

- 正文：空间记法通常把**代数量减少至少 4 倍**
- 前言：典型是 **4 到 6 倍**的缩减

> *"With the barrier of algebra out of the way, the analyst is free to state a problem
> more succinctly, to solve it in fewer steps, and to arrive at a more compact solution."*

### 对写代码的好处

原书举了个例子：如果类型声明得当，式 1.5 可以直接写成

```cpp
f = I*a + v.cross(I*v);
```

编译器知道 `I` 是刚体惯性、`v` 是空间运动向量，
就能把 `I*v` 编译成调用空间运算库中专门优化的例程。

> 💡 **这正是 `code/spatial.py` 的设计思路**，也是 Pinocchio 等现代库
> 用 C++ 模板做类型分派的理由：**类型即语义**。

---

## 1.3 单位与记号 (Units and Notation)

- **角度用弧度**。除此之外方程与单位制无关，只要求所选单位制**自洽**；
  少数明确提到单位的地方用 SI。
- 记号大体遵循 ISO：**变量斜体**、**常量与函数正体**、**向量与矩阵粗斜体**；
  向量小写、矩阵大写。

**几个需要记住的细节**：

| 记号 | 含义 |
|---|---|
| $\mathbf 0$、$\mathbf 1$ | 零矩阵、单位矩阵（$\mathbf 0$ 也表示零向量） |
| $A^{-\mathsf T}$ | **逆的转置**，即 $(A^{-1})^{\mathsf T}$ |
| $a\times$ | 把 $b$ 映到 $a\times b$ 的**算子** |
| **数组 (array)** | 一个编号的列表；若 $\lambda$ 是数组，其元素 $i$ 写作 $\lambda(i)$ |
| 前置上标 ${}^{A}v$ | 标识**坐标系** |
| 帽子 $\hat v$ | 空间向量（与 3D 量重名时） |
| 下划线 $\underline v$ | 坐标向量（与抽象向量需区分时） |

> 💡 **$\lambda(i)$ 写成"数组的元素"而不是"函数值"**，这个措辞是有意的：
> 第 4 章的 parent array 在代码里就是一个数组。

原书 p.265 有符号表，本仓库的对应物是 [`reference/notation.md`](../reference/notation.md)。

---

## 1.4 读者指南 (Readers' Guide)

**原书自己给的三段划分**：

| 部分 | 章节 | 内容 |
|---|---|---|
| **预备 (preparation)** | 2, 3, 4 | 空间向量代数、运动方程的建立与分析、系统模型 |
| **主要算法 (main algorithms)** | 5, 6, 7, 8 | RNEA、CRBA、ABA、闭环 |
| **附加专题 (additional topics)** | 9, 10, 11 | 混合动力学与浮动基、精度与效率、接触与碰撞 |

**原书的几条建议**：

- *"Readers with enough background may find they can skip straight to Chapter 5."*
  —— 底子够的可以直接跳到第 5 章。
- 第 2、3 章是**最数学的两章**，而且**讲得比后续算法所需的最低限度更深**。
  （言下之意：这两章可以先读个够用，回头再深挖。）
- 第 4 章定义了后面大量要用的量。
- 第 5–8 章**按复杂度递增排列**（算法越来越复杂）。
- 第 9 章起假定读者掌握前面的基础，但**各自相对独立**。

**关于例子和伪代码**：

> *"every major algorithm is presented both as a set of equations and as a pseudocode
> program. In many cases, the two are side-by-side."*

例子在书中**分布不均**：有些是为了图解正文概念，有些则是"用例子讲更方便"的内容。
**所以例子不能跳。**

**关于源码**：本书不附带源码（原书理由：这类软件过时太快），
但作者网站上有大部分算法的实现。原书给的地址是
`http://users.rsise.anu.edu.au/~roy/spatial/`（并注明"is or was"），
找不到就搜 "Roy Featherstone"。

> 📝 这就是通称的 **`spatial_v2`** MATLAB 工具箱。
> 本仓库的 [`code/`](../code/) 是一份用 Python 重写的等价实现，
> 与书中伪代码逐行对应，可直接对照阅读。

---

## 1.5 延伸阅读 (Further Reading)

原书给的书单，按用途分类：

**入门教材**（面向本科生，例题和习题多）：
Amirouche (2006)、Moon (1998)、Shabana (2001)；Huston (1990) 更偏教程式。

**进阶**：Stejskal & Valášek (1996)、Roberson & Schwertassek (1988)、Wittenburg (1977)。

**面向游戏/VR**：Coutinho (2001)。

**机器人动力学**：Balafoutis & Patel (1991)、Lilly (1993)、Yamane (2004)，
以及作者自己的 **Featherstone (1987) 《Robot Dynamics Algorithms》**
——原书明说 *"has been superseded by the present volume"*（已被本书取代）。

**6D 向量的其他体系**：Ball (1900)（螺旋理论主要著作）、Brand (1953)、
von Mises (1924a,b)（motor 代数）、Murray et al. (1994)、Selig (1996)（李群/李代数）。

> ⚠️ 原书特别注明：*"the treatment of spatial vectors in Featherstone (1987)
> is a little different from that presented here."*
> **读 1987 年那本旧书时要注意记号差异**（旧书更接近 motor 代数，
> 所有量放在同一个空间里；本书刻意分成 $M^6$ 与 $F^6$ 两个空间）。

---

## 本章要点

1. **两个（三个）基本计算**：FD 用于仿真，ID 用于控制**且是 FD 的组成部分**。
2. **式 1.1 是全书的坐标原点**：$\tau=H\ddot q+C$，注意 $C$ 是**含重力的向量**。
3. **`model` 是算法的输入**——这就是"基于模型的算法"，一份代码通用于一类机构。
4. **两大系统类别**：运动学树（易）vs 闭环系统（难）。
5. **空间向量的收益**：代数量减少 4–6 倍，一条方程替代两条，惯性可加。
6. **符号约定**：帽子=空间向量，下划线=坐标向量，$A^{-\mathsf T}$=逆的转置，$\lambda(i)$=数组元素。

---

## 易错点

1. **$C$ 是向量且含重力**——与多数机器人学教材的 $C$ 矩阵不是一回事。
2. **不要把空间向量理解成"$4\times4$ 齐次变换的另一种写法"**。
   齐次变换描述**位形**，空间向量描述**速度、加速度、力、动量**。二者互补。
3. **"spatial" 不是"三维空间的"**，而是 Featherstone 沿用的术语，
   指"同时含角分量和线分量的 6D 量"。
4. **ID 与 FD 不是简单互逆**：ID 是纯递推 $O(n)$，FD 需要解方程或更复杂的递推。

---

## 与其他章的联系

- → 第 2 章：兑现"6D 向量"（§1.2 的展开）
- → 第 4 章：兑现"`model`"（§1.1 的展开）
- → 第 5 章：$\mathrm{ID}$ 的实现
- → 第 6、7 章：$\mathrm{FD}$ 的两条实现路线
- → 第 9 章：混合动力学
- → 第 10 章：效率的定量讨论

---

## ✍️ 我的理解

<!-- 读完后用自己的话复述 -->

## ❓ 疑问与待办

- [ ] 下载 `spatial_v2`，与本仓库 `code/` 对照
- [ ] 确认自己惯用的教材里 $C$ 的定义，避免与本书混用

## 📌 与原文的出入

<!-- 本笔记已按原书 pp.1–6 逐节核对 -->
