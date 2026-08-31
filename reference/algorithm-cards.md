# 算法卡片速查

> 三大算法的伪代码集中版，便于对照实现。
> 详细解释见各章笔记。所有算法都假设**正则编号** $\lambda(i)<i$。

## 总览

| 算法 | 解决 | 输入 | 输出 | 趟数 | 复杂度 | 原书表 | 章节 |
|---|---|---|---|---|---|---|---|
| **RNEA** | 逆动力学 ID | $q,\dot q,\ddot q$ | $\tau$ | 2（外-内） | $O(n)$ | 表 5.1 | [第 5 章](../docs/ch05-inverse-dynamics-rnea.md) |
| **$\mathrm{ID}_\delta$** | 算 $H\ddot q$ | $q,\ddot q$ | $H\ddot q$ | 2 | $O(n)$ | 表 6.1 | [第 6 章](../docs/ch06-forward-dynamics-crba.md) |
| **CRBA** | 构造 $H$ | $q$ | $H$ | 1（内） | $O(nd)$ | 表 6.2 | [第 6 章](../docs/ch06-forward-dynamics-crba.md) |
| **$LTL$/$LTDL$** | 分解 $H$ | $H,\lambda$ | $L$（或 $L,D$） | — | $O(nd^2)$ | 表 6.3 | [第 6 章](../docs/ch06-forward-dynamics-crba.md) |
| **ABA** | 正动力学 FD | $q,\dot q,\tau$ | $\ddot q$ | 3（外-内-外） | $O(N_B)$ | 表 7.1 | [第 7 章](../docs/ch07-forward-dynamics-aba.md) |
| **$C-\tau^a$** | 闭环的偏置力 | $q,\dot q$ | $C-\tau^a$ | 3 | $O(n)$ | 表 8.3 | [第 8 章](../docs/ch08-closed-loop-systems.md) |
| **混合 ABA** | 混合动力学 | 部分 $\tau$、部分 $\ddot q$ | 其余 | 3 | $O(n)$ | 表 9.3 | [第 9 章](../docs/ch09-hybrid-dynamics.md) |
| **浮动基 ID** | 浮动基逆动力学 | $q,\dot q,\ddot q$ | $\tau,a_0$ | 3 | $O(n)$ | 表 9.6 | [第 9 章](../docs/ch09-hybrid-dynamics.md) |
| **接触仿真** | 含接触的仿真 | — | — | — | — | 表 11.1 | [第 11 章](../docs/ch11-contact-and-impact.md) |

**通用子程序**：`[X_J, S, c_J] = jcalc(jtype, q, q̇)`，$\mathbf{v}_J = S\dot q$

---

## 卡片 1：RNEA（逆动力学）

```
────────────────────────────────────────────────────────────
τ = ID(model, q, q̇, q̈, [f^x])                    O(n)
────────────────────────────────────────────────────────────
v₀ = 0 ;  a₀ = −a_g                        # 重力技巧

for i = 1 to N_B:                          # ── 外推 ──
    [X_J, S_i, v_J, c_J] = jcalc(jtype(i), q_i, q̇_i)
    ⁱX_λ(i) = X_J · X_T(i)
    v_i = ⁱX_λ(i)·v_λ(i) + v_J
    a_i = ⁱX_λ(i)·a_λ(i) + S_i·q̈_i + c_J + v_i × v_J
    f_i = I_i·a_i + v_i ×* I_i·v_i − ⁱX₀*·f^x_i

for i = N_B to 1:                          # ── 内推 ──
    τ_i = S_iᵀ·f_i
    if λ(i) ≠ 0:  f_λ(i) += λ⁽ⁱ⁾X_i*·f_i
────────────────────────────────────────────────────────────
```

**衍生用法**

| 目标 | 调用 |
|---|---|
| 偏置力 $C$ | `ID(q, q̇, 0)` |
| 重力项 $g(q)$ | `ID(q, 0, 0)` |
| $H$ 的第 $j$ 列 | `ID(q, 0, e_j) − ID(q, 0, 0)` |
| 计算力矩控制 | `ID(q, q̇, q̈_d + K_d ė + K_p e)` |

---

## 卡片 2：CRBA（构造关节空间惯性矩阵）

```
────────────────────────────────────────────────────────────
H = crba(model, q)                          O(n²) / O(nd)
────────────────────────────────────────────────────────────
for i = 1 to N_B:  I^c_i = I_i
H = 0

for i = N_B to 1:                          # ── 内推 ──
    if λ(i) ≠ 0:
        I^c_λ(i) += λ⁽ⁱ⁾X_i*·I^c_i·ⁱX_λ(i)

    F    = I^c_i·S_i
    H_ii = S_iᵀ·F

    j = i
    while λ(j) ≠ 0:
        F    = λ⁽ʲ⁾X_j*·F
        j    = λ(j)
        H_ij = Fᵀ·S_j
        H_ji = H_ijᵀ
────────────────────────────────────────────────────────────
```

**$H$ 元素的闭式表达**（`H_ij = F' * S_j` 那一行的数学形式）

$$
H_{ij}=S_i^{\mathsf T}\ I^{c}_{i}\ {}^{i}X_{j}\ S_j
\qquad\text{要求 }j\text{ 是 }i\text{ 的祖先，}I^{c}\text{ 取}\textbf{后代}\text{ }i\text{ 的}
$$

⚠️ 下标条件极易记反，且记反后 $H$ 仍然对称、量纲仍然正确。
完整推导、物理解释与算例见
[第 6 章](../docs/ch06-forward-dynamics-crba.md) 的「$H$ 的元素：完整解析」一节。

**配套：FD 完整流程**

```
C  = ID(q, q̇, 0)              # O(n)
H  = crba(q)                   # O(n²) / O(nd)
H  = Lᵀ L                      # LTL 分解，O(nd²)，不产生填充
Lᵀ y = τ − C                   # 前代
L  q̈ = y                       # 回代
```

---

## 卡片 3：ABA（正动力学）

```
────────────────────────────────────────────────────────────
q̈ = FDab(model, q, q̇, τ, [f^x])                  O(n)
────────────────────────────────────────────────────────────
# ═══ 趟 1：外推 ═══
v₀ = 0
for i = 1 to N_B:
    [X_J, S_i, v_J, c_J] = jcalc(jtype(i), q_i, q̇_i)
    ⁱX_λ(i) = X_J · X_T(i)
    v_i   = ⁱX_λ(i)·v_λ(i) + v_J
    c_i   = c_J + v_i × v_J
    I^A_i = I_i
    p^A_i = v_i ×* I_i·v_i − ⁱX₀*·f^x_i

# ═══ 趟 2：内推（消元） ═══
for i = N_B to 1:
    U_i = I^A_i·S_i
    D_i = S_iᵀ·U_i
    u_i = τ_i − S_iᵀ·p^A_i
    if λ(i) ≠ 0:
        I^a = I^A_i − U_i·D_i⁻¹·U_iᵀ          # Schur 补
        p^a = p^A_i + I^a·c_i + U_i·D_i⁻¹·u_i  # ⚠️ 用 I^a 不是 I^A
        I^A_λ(i) += λ⁽ⁱ⁾X_i*·I^a·ⁱX_λ(i)
        p^A_λ(i) += λ⁽ⁱ⁾X_i*·p^a

# ═══ 趟 3：外推（回代） ═══
a₀ = −a_g
for i = 1 to N_B:
    a'  = ⁱX_λ(i)·a_λ(i) + c_i
    q̈_i = D_i⁻¹·(u_i − U_iᵀ·a')
    a_i = a' + S_i·q̈_i
────────────────────────────────────────────────────────────
```

⚠️ 趟 1 算的 $\mathbf{c}_i$、趟 2 算的 $U_i,D_i,u_i$ 都要**存下来**给趟 3 用。

---

## 卡片 4：$\mathrm{ID}_\delta$（表 6.1）—— 计算 $H\ddot q$ 的简化 RNEA

```
────────────────────────────────────────────────────────────
τ = ID_δ(model, q, q̈)                            O(n)
「所有依赖 q̇ 与 f^x 的项都相消」⟹ 无速度项、无重力、无外力
────────────────────────────────────────────────────────────
a₀ = 0                              ← 不是 −a_g！
for i = 1 to N_B:
    a_i = ⁱX_λ(i)·a_λ(i) + S_i·q̈_i
    f_i = I_i·a_i                   ← 没有 v ×* I v 项
for i = N_B to 1:
    τ_i = S_iᵀ·f_i
    if λ(i) ≠ 0:  f_λ(i) += λ⁽ⁱ⁾X_i*·f_i
────────────────────────────────────────────────────────────
```

用途：逐列构造 $H$（Walker & Orin 1982 的方法 1、2）；第 9 章混合动力学的步骤 4。

---

## 卡片 5：$LTL$ / $LTDL$ 稀疏分解（表 6.3）

**前提**：$H$ 对称正定；$\lambda$ 满足 $0\le\lambda(i)<i$；
$H$ 第 $k$ 行主对角线以下的非零元只在列 $\lambda(k),\lambda(\lambda(k)),\dots$。

```
──────── LTL（H = LᵀL） ────────      ──────── LTDL（H = LᵀDL） ────────
for k = n to 1:                        for k = n to 1:
    H_kk = sqrt(H_kk)                      i = λ(k)
    i = λ(k)                               while i ≠ 0:
    while i ≠ 0:                               a = H_ki / H_kk
        H_ki = H_ki / H_kk                     j = i
        i = λ(i)                               while j ≠ 0:
    i = λ(k)                                       H_ij = H_ij − a·H_kj
    while i ≠ 0:                                   j = λ(j)
        j = i                                  H_ki = a
        while j ≠ 0:                           i = λ(i)
            H_ij = H_ij − H_ki·H_kj    end
            j = λ(j)
        i = λ(i)
end
```

**回代**（表 6.5）：

```
x = L⁻¹x                            x = L⁻ᵀx
for i = 1 to n:                     for i = n to 1:
    j = λ(i)                            x_i = x_i / L_ii
    while j ≠ 0:                        j = λ(i)
        x_i −= L_ij·x_j                 while j ≠ 0:
        j = λ(j)                            x_j −= L_ij·x_i
    x_i = x_i / L_ii                        j = λ(j)
```

（$L$ 为单位下三角时把 $L_{ii}$ 换成 1。）
**性质**：重排 Cholesky / $LDL^{\mathsf T}$，**无填充**，数值性质与标准版相同。
$LTDL$ **不需要开方**，通常更优。

**完整 FD 流程**：

```
C = ID(q, q̇, 0);   H = crba(q);   H = LᵀL
Lᵀ y = τ − C;      L q̈ = y
```

---

## 卡片 6：混合动力学 ABA（表 9.3）

```
趟 1: if i ∈ fd:  c_i = cJ + v_i × vJ
      else:       c_i = cJ + v_i × vJ + S_i q̈_i        ← 差别①

趟 2: if i ∈ fd:                                        ← 消元（Schur 补）
          U_i = I^A_i S_i ;  D_i = S_iᵀU_i ;  u_i = τ_i − S_iᵀp^A_i
          I^a = I^A_i − U_i D_i⁻¹ U_iᵀ
          p^a = p^A_i + I^a c_i + U_i D_i⁻¹ u_i
      else:                                             ← 代入
          I^a = I^A_i
          p^a = p^A_i + I^a c_i
      累加 I^a, p^a 到父节点

趟 3: if i ∈ fd:  a' = ⁱX_λ a_λ + c_i ; q̈_i = D_i⁻¹(u_i − U_iᵀa') ; a_i = a' + S_i q̈_i
      else:       a_i = ⁱX_λ a_λ + c_i ; τ_i = S_iᵀ(I^A_i a_i + p^A_i)   ← 差别②
```

> **统一视角**：RNEA = 全代入；ABA = 全消元；混合动力学 = 逐关节选择。

---

## 卡片 7：浮动基逆动力学（表 9.6）

```
趟 1（外推）: a^r₀ = −⁰a_g
    v_i   = ⁱX_λ v_λ + vJ
    a^r_i = ⁱX_λ a^r_λ + cJ + v_i × vJ + S_i q̈_i
    I^c_i = I_i
    p^c_i = I_i a^r_i + v_i ×* I_i v_i − ⁱX₀* ⁰f^x_i

趟 2（内推）: I^c_λ += λX_i* I^c_i ⁱX_λ ;  p^c_λ += λX_i* p^c_i

趟 3（外推）: ⁰a₀ = −(I^c₀)⁻¹ p^c₀
    ⁱa₀ = ⁱX_λ λa₀
    τ_i = S_iᵀ(I^c_i ⁱa₀ + p^c_i)
```

**浮动基运动方程**（式 9.13）：

$$
\begin{bmatrix}I^c_0&F\\ F^{\mathsf T}&H\end{bmatrix}
\begin{bmatrix}a_0\\ \ddot q\end{bmatrix}
+\begin{bmatrix}p^c_0\\ C\end{bmatrix}=\begin{bmatrix}0\\ \tau\end{bmatrix},
\qquad F_i=I^c_iS_i
$$

⚠️ **不要**构造 $H^{fl}=H-F^{\mathsf T}(I^c_0)^{-1}F$——它稠密，可能比 $H$ 的非零元还多。
**动量** $=I^c_0v_0+F\dot q$（例 9.3）。

---

## 卡片 8：闭环系统

```
1. 选生成树（应选深度最小的，见第 10 章 §10.2）；剩余边 = 闭环关节
2. 树部分用 CRBA / RNEA 得 H 与 C − τ^a（表 8.3）
3. 用 K_lj = ε_lj T_kᵀ S_j 算 K；用式 8.15 算 k；用式 8.22 算 k_stab
4. 解式 8.36：
      [ H   Kᵀ ] [ q̈  ]   [ τ − C + τ^a ]
      [ K   0  ] [ −λ ] = [  k + k_stab ]
```

**三种解法**：

| 方法 | 做法 |
|---|---|
| 1 | 直接解鞍点系统（**对称不定 ⟹ 用 LDLᵀ 不能用 Cholesky**） |
| 2 | 先解 $\lambda$：$A=KH^{-1}K^{\mathsf T}$，用 $LTL$ 的六步流程（$A=Y^{\mathsf T}Y$，从不显式求 $H^{-1}$） |
| 3 | 独立坐标：$G^{\mathsf T}HG\,\ddot y=G^{\mathsf T}(\tau-C+\tau^a-Hg)$ |

**Baumgarte**：$\alpha=\beta=1/T_{stab}$；工业机器人 $T_{stab}\approx0.1$ 合理。
**目的是稳定，不是精度。**

---

## 卡片 9：接触动力学（表 11.1）

**LCP**：$\dot\zeta=M\lambda+d$，$\dot\zeta\ge0$，$\lambda\ge0$，$\dot\zeta^{\mathsf T}\lambda=0$

$$
M=T^{\mathsf T}H^{-1}T,\qquad d=T^{\mathsf T}H^{-1}(\tau-C)+\dot T^{\mathsf T}\dot q,
\qquad t_i=(J_{sc(i)}-J_{pc(i)})^{\mathsf T}n_i
$$

**等价 QP**（实践中更常用，求解器更易得）：

$$
\min_\lambda\ \tfrac12\lambda^{\mathsf T}M\lambda+\lambda^{\mathsf T}d
\quad\text{s.t.}\quad\lambda\ge0
$$

**仿真流程（表 11.1）**：

```
1. 向前积分，把「活动接触」当等式约束，忽略其他接触
2. 监测两类事件：几何事件（接触产生/失去）、负接触力
3. 检测到事件 ⟹ 插值回到最早事件时刻
4. 接触产生 ⟹ 冲量动力学（§11.7），更新两个接触集合
5. 几何失接触 ⟹ 从两个集合中移除
6. 负接触力 ⟹ 解 LCP/QP，更新两个集合
7. 回到 1
```

> 🔑 **核心**：只在**出现负 $\lambda_i$ 或发生碰撞**时才解 LCP，其余时间用便宜的等式约束。

**四个状态**（图 11.4）：

| 状态 | 条件 | 处理 |
|---|---|---|
| 1 | $\phi>0$ | 忽略 |
| 2a | $\phi=0,\dot\phi<0$ | **冲量** |
| 2b | $\phi=0,\dot\phi>0$ | 同状态 1 |
| 3 | $\phi=0,\dot\phi=0$ | **接触力**（LCP） |

**两体碰撞**：$\lambda=\dfrac{-(1+e)\,n\cdot(v_2-v_1)}{n\cdot(I_1^{-1}+I_2^{-1})\,n}$

**冲量方程**：$\iota=I\Delta v$（刚体）、$\iota=I^A\Delta v$（铰接体柄）、$u=H\Delta\dot q$（系统）

---

## 实现检查清单

写完代码后逐条过：

- [ ] `parent[i] < i` 对所有 $i$ 成立（加断言）
- [ ] 变换顺序是 ${}^{i}X_{\lambda(i)} = X_J\,X_T$，没写反
- [ ] $I_i$ 存在 body $i$ 自己的坐标系里（因此是常量）
- [ ] 力向量用 $X^*$ 和 $\times^*$，运动向量用 $X$ 和 $\times$
- [ ] 加速度交叉项是 $\mathbf{v}_i\times\mathbf{v}_J$，**没有系数 2**
- [ ] 重力是 $\mathbf{a}_0=-\mathbf{a}_g$（**负号**）
- [ ] 内推时力是**累加** `+=`，不是赋值 `=`（分支节点会丢数据）
- [ ] ABA 的 $\mathbf{p}^a$ 里有 $I^{a}\mathbf{c}_i$ 项，且用的是 $I^a$ 不是 $I^A$
- [ ] CRBA 填了对称的另一半 `H_ji = H_ij'`
- [ ] FD 解的是 $H\ddot q = \tau - C$，没漏掉 $C$
- [ ] 非对称关节的极性正确（必要时用表 4.2 的极性反转）
- [ ] 多自由度关节时用了展开的父数组 $\lambda'$（表 6.4）
- [ ] 闭环：没忘 $\tau^a$；鞍点矩阵没用 Cholesky
- [ ] 浮动基：四元数每步归一化；${}^{0}a_g$ 随基座姿态更新（不是常量）
- [ ] 接触：$n$ 指向表面外侧 / 由前驱指向后继

## 验证清单

- [ ] **RNEA ↔ ABA 互验**：`ID(q, q̇, FD(q, q̇, τ)) == τ` ← **最重要的一条**
- [ ] **CRBA ↔ $n$ 次 RNEA**：两种方式算的 $H$ 逐元素相等
- [ ] $\|H - H^{\mathsf T}\|$ 在机器精度量级
- [ ] $H$ 能通过 Cholesky（正定）
- [ ] 自由落体：单刚体 + 自由关节 + $\tau=0$ ⟹ $\mathbf{a}=\mathbf{a}_g$
- [ ] 静止：$\dot q=0$、$\tau=g(q)$ ⟹ $\ddot q=0$
- [ ] 能量守恒：无驱动无耗散时总能量漂移缓慢且有界
- [ ] 结构性零：随机 $q$ 下，非祖先后代关系的 $H_{ij}$ 恒为 0
