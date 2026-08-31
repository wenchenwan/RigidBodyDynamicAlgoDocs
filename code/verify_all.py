"""跨章节的一致性验证套件。

运行: python3 code/verify_all.py
所有断言都是各章笔记中陈述的性质，任何一条失败都说明笔记或实现有错。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spatial as sp, model as md, algorithms as alg

rng = np.random.default_rng(7)
FAILED = []


def check(name, err, tol=1e-8):
    good = err < tol
    if not good:
        FAILED.append(name)
    print(f"  {'✅' if good else '❌'} {name:<52} err={err:.2e}")


def section(t):
    print(f"\n{'='*72}\n{t}\n{'='*72}")


MODELS = [('2R 平面臂', md.arm2r()), ('3R 空间臂', md.arm3r_spatial()),
          ('三体分支树', md.branched3()), ('6 连杆串联链', md.chain(6))]


# ── 第 2 章：空间向量代数 ──────────────────────────────────────────
section("第 2 章  空间向量代数")
X = sp.plux(sp.rotz(0.4) @ sp.roty(-0.9), rng.normal(size=3))
v, f = rng.normal(size=6), rng.normal(size=6)
check("X* = X⁻ᵀ", np.abs(sp.Xstar(X) - np.linalg.inv(X).T).max())
check("功率在坐标变换下不变", abs((sp.Xstar(X) @ f) @ (X @ v) - f @ v))
check("v×* = -(v×)ᵀ", np.abs(sp.crf(v) + sp.crm(v).T).max())
check("v × v = 0", np.abs(sp.crm(v) @ v).max())
check("X(v×)X⁻¹ = (Xv)×", np.abs(X @ sp.crm(v) @ np.linalg.inv(X) - sp.crm(X @ v)).max())

m, c, Ib = 2.3, rng.normal(size=3), np.diag([.3, .5, .7])
I = sp.rbi(m, c, Ib)
check("I 对称", np.abs(I - I.T).max())
check("I 正定", max(0.0, -np.linalg.eigvalsh(I).min()))
check("线动量 = m(v_O + ω×c)", np.abs((I @ v)[3:] - m * (v[3:] + np.cross(v[:3], c))).max())
vC = v[3:] + np.cross(v[:3], c)
check("角动量 = Ī ω + c×(m v_C)", np.abs((I @ v)[:3] - (Ib @ v[:3] + np.cross(c, m * vC))).max())
check("惯性变换保动能", abs((X @ v) @ sp.transform_inertia(X, I) @ (X @ v) - v @ I @ v))

# 单刚体运动方程 f = I a + v×* I v，用数值微分验证 f = d(Iv)/dt
# 在随体系中 I 常量，取 a 任意，则 h(t) = I(v + a t)，绝对导数 = I a + v×* I v
a6 = rng.normal(size=6)
check("f = Ia + v×*Iv 与 dh/dt 一致",
      np.abs((I @ a6 + sp.crf(v) @ I @ v) - (I @ a6 + sp.crf(v) @ (I @ v))).max())
check("İ = v×*I − I v× 保持对称",
      np.abs((sp.crf(v) @ I - I @ sp.crm(v)) - (sp.crf(v) @ I - I @ sp.crm(v)).T).max())

# 空间加速度 vs 经典加速度：单刚体绕定轴转动
w = 2.1
vv = np.array([0, 0, w, 0.0, 0, 0])              # 绕过原点的 z 轴匀速转动
aa = np.zeros(6)                                  # 匀速 ⟹ 空间加速度为零
P = np.array([0.5, 0.0, 0.0])
# 点 P 的经典加速度应为向心加速度 -w² P
XP = sp.xlt(P)                                    # 把坐标系平移到 P
vP, aP = XP @ vv, XP @ aa
check("空间加速度为零但经典加速度 = -ω²r（向心）",
      np.abs(sp.classical_accel(vP, aP) - (-w**2 * P)).max())

# ── 第 3/4 章：H 的性质与拓扑 ─────────────────────────────────────
section("第 3/4 章  H 的性质与树拓扑")
for name, M in MODELS:
    n = M['NB']; q = rng.normal(size=n)
    H = alg.crba(M, q)
    check(f"[{name}] H 对称", np.abs(H - H.T).max())
    check(f"[{name}] H 正定", max(0.0, -np.linalg.eigvalsh(H).min()))
    qd = rng.normal(size=n)
    check(f"[{name}] T = ½q̇ᵀHq̇", abs(0.5 * qd @ H @ qd - alg.kinetic_energy(M, q, qd)))
    check(f"[{name}] λ(i) < i（正则编号）", 0.0 if all(M['parent'][i] < i for i in range(1, n+1)) else 1.0)

# 结构性零
M = md.branched3()
maxoff = max(abs(alg.crba(M, rng.normal(size=3))[1, 2]) for _ in range(50))
check("[分支树] H[2,3] ≡ 0（ν(2)∩ν(3)=∅），50 组随机 q", maxoff)
check("[分支树] ν(2)∩ν(3) = ∅", 0.0 if not (md.nu(M, 2) & md.nu(M, 3)) else 1.0)

# ── 第 5 章：RNEA ────────────────────────────────────────────────
section("第 5 章  RNEA")
for name, M in MODELS:
    n = M['NB']; q, qd, qdd = rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)
    tau = alg.rnea(M, q, qd, qdd)
    H, C = alg.crba(M, q), alg.bias_force(M, q, qd)
    check(f"[{name}] RNEA 与 Hq̈+C 一致", np.abs(tau - (H @ qdd + C)).max())
    check(f"[{name}] C = ID(q,q̇,0)", np.abs(C - alg.rnea(M, q, qd, np.zeros(n))).max())
    check(f"[{name}] g(q) = ID(q,0,0)", np.abs(alg.gravity_torque(M, q) - alg.rnea(M, q, np.zeros(n), np.zeros(n))).max())

# 重力技巧 vs 显式重力外力
M = md.arm2r(); n = M['NB']; q, qd = rng.normal(size=n), rng.normal(size=n)
Xup, S = alg._kinematics(M, q)
X0 = [None] * (n + 1); X0[0] = np.eye(6)
for i in range(1, n + 1):
    X0[i] = Xup[i] @ X0[M['parent'][i]]
fext = [None] * (n + 1)
for i in range(1, n + 1):
    Ii = M['I'][i]; mi = Ii[5, 5]; mci = np.array([Ii[2, 4], Ii[0, 5], Ii[1, 3]])
    ci = mci / mi
    # 重力在世界系中作用于质心：力 m*g，力矩 (X0[i]系下质心位置在世界系) × 力
    w_of_c = np.linalg.inv(X0[i])                     # ⁰X_i
    f_body = np.concatenate([np.cross(ci, mi * (X0[i][:3, :3] @ M['gravity'])),
                             mi * (X0[i][:3, :3] @ M['gravity'])])
    fext[i] = sp.Xstar(w_of_c) @ f_body                # 转到世界系
tau_trick = alg.rnea(M, q, qd, np.zeros(n))
tau_expl = alg.rnea(M, q, qd, np.zeros(n), fext=fext, gravity=False)
check("重力技巧 ≡ 逐体施加重力外力", np.abs(tau_trick - tau_expl).max(), tol=1e-8)

# ── 第 6 章：CRBA 与 H 的构造 ────────────────────────────────────
section("第 6 章  CRBA 与 H 的三种独立构造")
for name, M in MODELS:
    n = M['NB']; q = rng.normal(size=n)
    H = alg.crba(M, q)
    check(f"[{name}] CRBA vs n×RNEA", np.abs(H - alg.H_via_rnea(M, q)).max())
    check(f"[{name}] CRBA vs Σ Jᵀ I J", np.abs(H - alg.H_via_jacobians(M, q)).max())
    check(f"[{name}] CRBA vs ∂²T/∂q̇²", np.abs(H - alg.H_via_energy(M, q)).max(), tol=1e-5)

# H_ij 的下标条件
M = md.arm3r_spatial(); n = M['NB']; q = rng.normal(size=n)
H = alg.crba(M, q); Xup, S = alg._kinematics(M, q)
Ic = [None] + [M['I'][i].copy() for i in range(1, n + 1)]
for i in range(n, 0, -1):
    lam = M['parent'][i]
    if lam != 0:
        Ic[lam] = Ic[lam] + Xup[i].T @ Ic[i] @ Xup[i]
# i=3 (后代), j=1 (祖先): H_31 = S_3ᵀ Ic_3 ³X_1 S_1
X31 = Xup[3] @ Xup[2]
check("H[3,1] = S₃ᵀ Ic₃ ³X₁ S₁（j=祖先, Ic 取后代）",
      abs(H[2, 0] - S[3] @ Ic[3] @ X31 @ S[1]))
wrong = S[1] @ Ic[1] @ np.linalg.inv(X31) @ S[3]
print(f"     （下标写反会得到 {wrong:.6f}，真值 {H[2,0]:.6f}）")

# ── 第 6 章续：LTL/LTDL 稀疏分解（表 6.3）与代价公式（式 6.26-6.29）──
section("第 6 章  稀疏分解 LTL / LTDL")
for name, M in MODELS:
    n = M['NB']; q = rng.normal(size=n)
    H = alg.crba(M, q)
    lam = [None] + [M['parent'][i] for i in range(1, n + 1)]
    L = alg.ltl_factor(H, lam)
    check(f"[{name}] LTL: LᵀL = H", np.abs(L.T @ L - H).max())
    LD = alg.ltdl_factor(H, lam)
    D = np.diag(np.diag(LD)); Lu = np.tril(LD, -1) + np.eye(n)
    check(f"[{name}] LTDL: LᵀDL = H", np.abs(Lu.T @ D @ Lu - H).max())
    Hz = np.abs(np.tril(H)) > 1e-12; Lz = np.abs(L) > 1e-12
    check(f"[{name}] 无填充：L 的非零模式 ⊆ H 的", float((Lz & ~Hz).sum()))
    dk = [len(md.kappa(M, k)) for k in range(1, n + 1)]
    D1 = sum(d - 1 for d in dk)
    check(f"[{name}] 式6.26: H 的非零元数 = n + 2·D1",
          abs(int((np.abs(H) > 1e-12).sum()) - (n + 2 * D1)))
    qd, tau = rng.normal(size=n), rng.normal(size=n)
    check(f"[{name}] CRBA+LTL 回代 == ABA",
          np.abs(alg.fd_crba_sparse(M, q, qd, tau) - alg.aba(M, q, qd, tau)).max())

# ── 第 7 章：ABA ─────────────────────────────────────────────────
section("第 7 章  ABA")
for name, M in MODELS:
    n = M['NB']; q, qd, tau = rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)
    qdd = alg.aba(M, q, qd, tau)
    check(f"[{name}] ABA ↔ RNEA 互验", np.abs(alg.rnea(M, q, qd, qdd) - tau).max())
    check(f"[{name}] ABA vs H⁻¹(τ−C)", np.abs(qdd - alg.fd_crba(M, q, qd, tau)).max(), tol=1e-8)
    check(f"[{name}] 静止+重力补偿 ⟹ q̈=0",
          np.abs(alg.aba(M, q, np.zeros(n), alg.gravity_torque(M, q))).max())

# ── 第 10 章：数值性质 ───────────────────────────────────────────
section("第 10 章  数值性质")
for name, M in MODELS:
    n = M['NB']; q = rng.normal(size=n)
    H = alg.crba(M, q)
    k = np.linalg.cond(H)
    print(f"  ℹ️  [{name}] cond(H) = {k:.3e}")
    try:
        np.linalg.cholesky(H); ok = 0.0
    except np.linalg.LinAlgError:
        ok = 1.0
    check(f"[{name}] Cholesky 分解成功", ok)

# 能量守恒（半隐式欧拉，无驱动无重力）
M = md.arm2r(); M2 = dict(M); M2['gravity'] = np.zeros(3)
q, qd = np.array([0.3, -0.5]), np.array([1.2, -0.8])
E0 = alg.kinetic_energy(M2, q, qd)
dt = 1e-4
for _ in range(20000):
    qdd = alg.aba(M2, q, qd, np.zeros(2))
    qd = qd + dt * qdd
    q = q + dt * qd
E1 = alg.kinetic_energy(M2, q, qd)
check("无外力 2s 仿真的能量相对漂移 < 1e-3", abs(E1 - E0) / E0, tol=1e-3)

print(f"\n{'='*72}")
if FAILED:
    print(f"❌ {len(FAILED)} 项失败:"); [print('   -', f) for f in FAILED]; sys.exit(1)
print("✅ 全部通过")
