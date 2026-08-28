# 符号表 (Notation)

> 本书符号密度极高，建议读的时候把这一页放在手边。
> 遇到新符号随时补充，并注明首次出现的章节。

## 一、空间向量与算子

| 符号 | 名称 | 空间/类型 | 说明 |
|---|---|---|---|
| $M^6$ | 运动空间 (motion space) | — | 装速度、加速度 |
| $F^6$ | 力空间 (force space) | — | 装力、动量；$F^6=(M^6)^*$ |
| $\mathbf{v}$ | 空间速度 (spatial velocity) | $M^6$ | $[\omega;\ v_O]$ |
| $\mathbf{a}$ | 空间加速度 (spatial acceleration) | $M^6$ | $\dot{\mathbf v}$，**不是**物质点加速度 |
| $\mathbf{f}$ | 空间力 (spatial force) | $F^6$ | $[n_O;\ f]$ |
| $\mathbf{h}$ | 空间动量 (spatial momentum) | $F^6$ | $\mathbf{h}=I\mathbf{v}$ |
| $I$ | 空间刚体惯性 | $M^6\!\to\!F^6$ | $6\times6$，对称正定，10 个独立参数 |
| $X$ | 运动向量坐标变换 | $M^6\!\to\!M^6$ | ${}^{B}X_{A}$：$A$ 系 → $B$ 系 |
| $X^{*}$ | 力向量坐标变换 | $F^6\!\to\!F^6$ | $X^{*}=X^{-\mathsf T}$ |
| $\mathbf{v}\times$ | 空间叉乘（作用于运动） | $M^6\!\to\!M^6$ | $\mathrm{ad}_{\mathbf v}$ |
| $\mathbf{v}\times^{*}$ | 空间叉乘（作用于力） | $F^6\!\to\!F^6$ | $=-(\mathbf{v}\times)^{\mathsf T}$，$\mathrm{ad}^*_{\mathbf v}$ |

> 💡 **`*` 的含义统一**：带星号 = 力空间的对应物。

## 二、3D 量

| 符号 | 含义 |
|---|---|
| $\omega$ | 角速度 (angular velocity) |
| $v_O$ | 刚体上瞬时与原点 $O$ 重合的**物质点**的线速度 |
| $n_O$ | 关于原点 $O$ 的力矩 |
| $f$ | 合力 |
| $m$ | 质量 |
| $c$ | 质心 (centre of mass) 位置，相对参考点 |
| $\bar I$ | 绕**质心**的 $3\times3$ 转动惯量 |
| $E$ | $3\times3$ 旋转矩阵 |
| $r$ | 坐标系原点间的位置向量 |
| $a\times$ | 向量 $a$ 的 $3\times3$ 反对称矩阵，$a\times b = (a\times)b$ |
| $\mathbf{1}_3$ | $3\times3$ 单位阵 |

## 三、系统与拓扑

| 符号 | 含义 | 章节 |
|---|---|---|
| $n$ | 系统总自由度 | 第 3 章 |
| $N_B$ | 刚体数目（不含 body 0） | 第 4 章 |
| $n_i$ | 关节 $i$ 的自由度 | 第 3 章 |
| $d$ | 树深度 (depth) | 第 10 章 |
| $\lambda(i)$ | body $i$ 的**父节点** (parent) | 第 4 章 |
| $\mu(i)$ | body $i$ 的**子节点集合** (children) | 第 4 章 |
| $\nu(i)$ | 以 $i$ 为根的**子树** (subtree)，含 $i$ | 第 4 章 |
| $\kappa(i)$ | 从 $i$ 到根路径上的关节集合 (support) | 第 4 章 |
| body 0 | 固定基座 / 惯性系 | 第 4 章 |

**核心约定**：**正则编号** $\lambda(i) < i$ —— 父的编号小于子的编号。

## 四、关节与模型

| 符号 | 含义 | 维度 |
|---|---|---|
| $q,\ \dot q,\ \ddot q$ | 关节位置 / 速度 / 加速度 | $n$（球/自由关节 $\dim q\ne n$） |
| $\tau$ | 广义力 / 关节力矩 | $n$ |
| $S_i$ | 运动子空间矩阵 (motion subspace) | $6\times n_i$ |
| $T_i$ | 约束力子空间矩阵 | $6\times(6-n_i)$，$S^{\mathsf T}T=0$ |
| $X_{T}(i)$ | 树变换 (tree transform)，**常量** | 描述连杆几何 |
| $X_{J}(i)$ | 关节变换 (joint transform)，依赖 $q_i$ | 描述关节位形 |
| ${}^{i}X_{\lambda(i)}$ | 父 → 子的完整变换 | $=X_J X_T$ |
| $\mathbf{v}_{J}$ | 关节速度 | $=S_i\dot q_i$ |
| $\mathbf{c}_{J}$ | 关节偏置加速度 | $=\dot S_i\dot q_i$，常见关节为 0 |

## 五、动力学量

| 符号 | 含义 | 章节 |
|---|---|---|
| $H(q)$ | 关节空间惯性矩阵 (JSIM)，$n\times n$ 对称正定 | 第 3 章 |
| $C(q,\dot q)$ | 偏置力**向量**（科氏+离心+重力+外力） | 第 3 章 |
| $\mathbf{a}_g$ | 重力加速度对应的空间加速度 | 第 5 章 |
| $\mathbf{f}^{x}_i$ | 作用在 body $i$ 上的外力 | 第 5 章 |
| $I^{c}_i$ | 复合刚体惯性（子树**焊死**） | 第 6 章 |
| $I^{A}_i$ | 铰接体惯性（子树关节**自由**） | 第 7 章 |
| $\mathbf{p}^{A}_i$ | 铰接体偏置力 | 第 7 章 |
| $U_i, D_i, u_i$ | ABA 的中间量 | 第 7 章 |
| $\mathbf{c}_i$ | 偏置加速度 $=\mathbf{c}_J+\mathbf{v}_i\times\mathbf{v}_J$ | 第 7 章 |
| $L$ | $LTL$ 分解的下三角因子，$H=L^{\mathsf T}L$ | 第 6 章 |
| $K$ | 环雅可比 (loop Jacobian) | 第 8 章 |
| $\lambda$ | Lagrange 乘子（⚠️ 与 parent 的 $\lambda(i)$ 同字母，靠上下文区分） | 第 8 章 |
| $\Lambda$ | 操作空间惯性矩阵 $=(JH^{-1}J^{\mathsf T})^{-1}$ | 第 9 章 |
| $\iota$ | 广义冲量 | 第 3/9 章 |

⚠️ **符号冲突提醒**：$\lambda$ 在第 4 章是 parent 函数，在第 8 章是 Lagrange 乘子。
原书如此，靠上下文区分——前者永远带下标 $\lambda(i)$，后者不带。

## 六、上下标约定

| 写法 | 含义 |
|---|---|
| ${}^{B}X_{A}$ | 左上标 = 目标坐标系，右下标 = 源坐标系 |
| ${}^{A}\mathbf{v}$ | 左上标 = 该向量所在的坐标系 |
| $X^{*}$ | 星号 = 力空间版本 |
| $X^{\mathsf T}$ | 转置 |
| $I^{c},\ I^{A},\ I^{a}$ | 上标区分惯性的种类（复合 / 铰接体 / 传给父的等效值） |
| $\mathbf{f}^{x}$ | 上标 $x$ = external（外力） |

**记忆法**：变换的上下标像分数约分——
${}^{C}X_{B}\ {}^{B}X_{A} = {}^{C}X_{A}$，中间的 $B$ 消掉。

## 七、与其他文献/库的对照

⚠️ **头号跨资料错误来源：分量排列顺序**。

| 约定 | 排列 | 使用者 |
|---|---|---|
| **Featherstone（本书）** | **(角, 线)**：$[\omega;\ v]$、$[n;\ f]$ | 本书、`spatial_v2`、RBDL |
| 部分李群教材 / 部分库接口 | (线, 角)：$[v;\ \omega]$ | 部分 screw theory 文献、部分库 |

跨库移植代码时**第一件事就是确认排列顺序**，否则会得到"看起来差不多但就是不对"的结果。

| 本书 | 李群记法 | Pinocchio |
|---|---|---|
| $\mathbf{v}\in M^6$ | $\mathfrak{se}(3)$ | `Motion` |
| $\mathbf{f}\in F^6$ | $\mathfrak{se}(3)^*$ | `Force` |
| $X$ | $\mathrm{Ad}_g$ | `SE3::act` |
| $X^*$ | $\mathrm{Ad}_g^{-\mathsf T}$ | — |
| $\mathbf{v}\times$ | $\mathrm{ad}_{\mathbf v}$ | `Motion::cross` |
| $\mathbf{v}\times^*$ | $\mathrm{ad}^*_{\mathbf v}$ | `Motion::crossForce` |

---

## ✍️ 随读补充

<!-- 遇到本表没有的符号，记在这里 -->

| 符号 | 含义 | 出处 |
|---|---|---|
|  |  |  |
