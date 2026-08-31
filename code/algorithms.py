"""RNEA / CRBA / ABA 的直译实现。

刻意与书中伪代码逐行对应，不做任何性能优化——目的是可读与可验证。
配套 docs/ch05, ch06, ch07。
"""
import numpy as np
import spatial as sp


def _kinematics(model, q):
    """公共前置：各关节的 Xup = ⁱX_λ(i) 与 S。"""
    NB = model['NB']
    Xup, S = [None] * (NB + 1), [None] * (NB + 1)
    for i in range(1, NB + 1):
        XJ, Si = sp.jcalc(model['jtype'][i], q[i - 1])
        Xup[i] = XJ @ model['Xtree'][i]
        S[i] = Si
    return Xup, S


def a_gravity(model):
    """重力对应的空间加速度 a_g = [0; g]（线分量为重力加速度）。"""
    return np.concatenate([np.zeros(3), model['gravity']])


# ────────────────────────────── RNEA ──────────────────────────────

def rnea(model, q, qd, qdd, fext=None, gravity=True):
    """逆动力学 τ = ID(q, q̇, q̈)。fext[i] 为世界系下作用于 body i 的外力。"""
    NB = model['NB']
    Xup, S = _kinematics(model, q)
    v = [np.zeros(6)] * (NB + 1)
    a = [np.zeros(6)] * (NB + 1)
    f = [None] * (NB + 1)
    X0 = [None] * (NB + 1)                    # ⁱX_0

    a[0] = -a_gravity(model) if gravity else np.zeros(6)   # ← 重力技巧

    for i in range(1, NB + 1):                # 外推
        lam = model['parent'][i]
        X0[i] = Xup[i] if lam == 0 else Xup[i] @ X0[lam]
        vJ = S[i] * qd[i - 1]
        v[i] = Xup[i] @ v[lam] + vJ
        a[i] = Xup[i] @ a[lam] + S[i] * qdd[i - 1] + sp.crm(v[i]) @ vJ
        f[i] = model['I'][i] @ a[i] + sp.crf(v[i]) @ model['I'][i] @ v[i]
        if fext is not None and fext[i] is not None:
            f[i] = f[i] - sp.Xstar(X0[i]) @ fext[i]

    tau = np.zeros(NB)
    for i in range(NB, 0, -1):                # 内推
        tau[i - 1] = S[i] @ f[i]
        lam = model['parent'][i]
        if lam != 0:
            f[lam] = f[lam] + Xup[i].T @ f[i]      # Xup[i].T = λ⁽ⁱ⁾X*_i
    return tau


def bias_force(model, q, qd, **kw):
    """C(q, q̇) = ID(q, q̇, 0)。"""
    return rnea(model, q, qd, np.zeros(model['NB']), **kw)


def gravity_torque(model, q):
    """g(q) = ID(q, 0, 0)。"""
    n = model['NB']
    return rnea(model, q, np.zeros(n), np.zeros(n))


# ────────────────────────────── CRBA ──────────────────────────────

def crba(model, q):
    """关节空间惯性矩阵 H。"""
    NB = model['NB']
    Xup, S = _kinematics(model, q)
    Ic = [None] + [model['I'][i].copy() for i in range(1, NB + 1)]
    H = np.zeros((NB, NB))

    for i in range(NB, 0, -1):
        lam = model['parent'][i]
        if lam != 0:
            Ic[lam] = Ic[lam] + Xup[i].T @ Ic[i] @ Xup[i]
        F = Ic[i] @ S[i]
        H[i - 1, i - 1] = S[i] @ F
        j = i
        while model['parent'][j] != 0:
            F = Xup[j].T @ F                  # 力向量沿树上传：×X*
            j = model['parent'][j]
            H[i - 1, j - 1] = F @ S[j]
            H[j - 1, i - 1] = H[i - 1, j - 1]
    return H


def H_via_rnea(model, q):
    """用 n 次 RNEA 独立构造 H（第 5 章「用法 4」），作为 CRBA 的对照。"""
    n = model['NB']
    z = np.zeros(n)
    base = rnea(model, q, z, z, gravity=False)
    H = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n); e[j] = 1.0
        H[:, j] = rnea(model, q, z, e, gravity=False) - base
    return H


# ────────────────────────────── ABA ──────────────────────────────

def aba(model, q, qd, tau, fext=None):
    """正动力学 q̈ = FD(q, q̇, τ)，铰接体算法。"""
    NB = model['NB']
    Xup, S = _kinematics(model, q)
    v = [np.zeros(6)] * (NB + 1)
    c = [None] * (NB + 1)
    IA = [None] + [model['I'][i].copy() for i in range(1, NB + 1)]
    pA = [None] * (NB + 1)
    X0 = [None] * (NB + 1)

    # 趟 1：外推
    for i in range(1, NB + 1):
        lam = model['parent'][i]
        X0[i] = Xup[i] if lam == 0 else Xup[i] @ X0[lam]
        vJ = S[i] * qd[i - 1]
        v[i] = Xup[i] @ v[lam] + vJ
        c[i] = sp.crm(v[i]) @ vJ                       # c_J = 0（常量 S）
        pA[i] = sp.crf(v[i]) @ model['I'][i] @ v[i]
        if fext is not None and fext[i] is not None:
            pA[i] = pA[i] - sp.Xstar(X0[i]) @ fext[i]

    # 趟 2：内推（消元）
    U = [None] * (NB + 1); d = np.zeros(NB + 1); u = np.zeros(NB + 1)
    for i in range(NB, 0, -1):
        U[i] = IA[i] @ S[i]
        d[i] = S[i] @ U[i]
        u[i] = tau[i - 1] - S[i] @ pA[i]
        lam = model['parent'][i]
        if lam != 0:
            Ia = IA[i] - np.outer(U[i], U[i]) / d[i]        # Schur 补
            pa = pA[i] + Ia @ c[i] + U[i] * (u[i] / d[i])   # ⚠️ 用 Ia 不是 IA
            IA[lam] = IA[lam] + Xup[i].T @ Ia @ Xup[i]
            pA[lam] = pA[lam] + Xup[i].T @ pa

    # 趟 3：外推（回代）
    a = [np.zeros(6)] * (NB + 1)
    a[0] = -a_gravity(model)
    qdd = np.zeros(NB)
    for i in range(1, NB + 1):
        lam = model['parent'][i]
        ap = Xup[i] @ a[lam] + c[i]
        qdd[i - 1] = (u[i] - U[i] @ ap) / d[i]
        a[i] = ap + S[i] * qdd[i - 1]
    return qdd


def fd_crba(model, q, qd, tau):
    """另一条 FD 路线：解 H q̈ = τ − C。"""
    H = crba(model, q)
    C = bias_force(model, q, qd)
    return np.linalg.solve(H, tau - C)


# ─────────────────────── 能量 / 运动学量 ───────────────────────

def body_velocities(model, q, qd):
    NB = model['NB']
    Xup, S = _kinematics(model, q)
    v = [np.zeros(6)] * (NB + 1)
    for i in range(1, NB + 1):
        v[i] = Xup[i] @ v[model['parent'][i]] + S[i] * qd[i - 1]
    return v


def kinetic_energy(model, q, qd):
    v = body_velocities(model, q, qd)
    return 0.5 * sum(v[i] @ model['I'][i] @ v[i] for i in range(1, model['NB'] + 1))


def jacobians(model, q):
    """每个 body 的 6×n 雅可比 J_k（在 body k 坐标系中），v_k = J_k q̇。"""
    NB = model['NB']
    Xup, S = _kinematics(model, q)
    X0 = [None] * (NB + 1); X0[0] = np.eye(6)
    for i in range(1, NB + 1):
        X0[i] = Xup[i] @ X0[model['parent'][i]]
    J = [None] * (NB + 1)
    for k in range(1, NB + 1):
        Jk = np.zeros((6, NB))
        j = k
        while j != 0:
            Jk[:, j - 1] = X0[k] @ np.linalg.inv(X0[j]) @ S[j]
            j = model['parent'][j]
        J[k] = Jk
    return J


def H_via_jacobians(model, q):
    """H = Σ_k J_kᵀ I_k J_k（定义 C）。"""
    J = jacobians(model, q)
    return sum(J[k].T @ model['I'][k] @ J[k] for k in range(1, model['NB'] + 1))


def H_via_energy(model, q, h=1e-5):
    """H_ij = ∂²T/∂q̇_i∂q̇_j，中心二阶差分（定义 A）。

    这条路完全不碰 Ic / X* / 下标条件，是最独立的对拍方式。
    """
    n = model['NB']
    H = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            ea = np.zeros(n); ea[a] = h
            eb = np.zeros(n); eb[b] = h
            H[a, b] = (kinetic_energy(model, q, ea + eb)
                       - kinetic_energy(model, q, ea - eb)
                       - kinetic_energy(model, q, -ea + eb)
                       + kinetic_energy(model, q, -ea - eb)) / (4 * h * h)
    return H
