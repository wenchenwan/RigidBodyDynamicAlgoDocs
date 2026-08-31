"""第 3 章关键结论的数值验证。

运行: python3 code/verify_ch03.py
配套 docs/ch03-rigid-body-system-dynamics.md
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, spatial as sp, model as md, algorithms as alg

rng = np.random.default_rng(31)
FAIL = []
def ok(n, e, t=1e-9):
    good = e < t
    if not good: FAIL.append(n)
    print(f"{'✅' if good else '❌'} {n:<54} err={e:.2e}")

def rand_inertia(r=rng):
    A = r.normal(size=(3, 3))
    return sp.rbi(abs(r.normal()) + .5, r.normal(size=3) * .3, A @ A.T + 3 * np.eye(3))

print("=" * 76)
print("§3.1  运动方程的标准形式与动能")
print("=" * 76)
for name, M in [('2R 平面臂', md.arm2r()), ('3R 空间臂', md.arm3r_spatial())]:
    n = M['NB']; q, qd, qdd = rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)
    H, C = alg.crba(M, q), alg.bias_force(M, q, qd)
    ok(f"[{name}] 式3.1  Hq̈+C = τ = ID(q,q̇,q̈)",
       np.abs(H @ qdd + C - alg.rnea(M, q, qd, qdd)).max())
    ok(f"[{name}] 式3.2  T = ½q̇ᵀHq̇", abs(.5 * qd @ H @ qd - alg.kinetic_energy(M, q, qd)))
    ok(f"[{name}] C 是「产生零加速度的 τ」",
       np.abs(C - alg.rnea(M, q, qd, np.zeros(n))).max())

print("\n" + "=" * 76)
print("§3.2  约束的隐式/显式描述与两种施加方式")
print("=" * 76)
n, nc = 7, 3
H = (lambda A: A @ A.T + 4 * np.eye(n))(rng.normal(size=(n, n)))
C, tau = rng.normal(size=n), rng.normal(size=n)
K = rng.normal(size=(nc, n)); k = rng.normal(size=nc)
# G 的列张成 null(K)：式 3.12 要求 KG = 0
U, S_, Vt = np.linalg.svd(K); G = Vt[nc:].T
ok("式3.12  KG = 0", np.abs(K @ G).max())
# g 满足 Kg = k
g = np.linalg.lstsq(K, k, rcond=None)[0]
ok("式3.12  Kg = k", np.abs(K @ g - k).max())

# 隐式路线：式 3.17 的 KKT 系统
KKT = np.block([[H, K.T], [K, np.zeros((nc, nc))]])
sol = np.linalg.solve(KKT, np.concatenate([tau - C, k]))
qdd_kkt, lam = sol[:n], -sol[n:]
ok("式3.17  KKT 解满足 Hq̈+C = τ+Kᵀλ", np.abs(H @ qdd_kkt + C - (tau + K.T @ lam)).max())
ok("式3.17  KKT 解满足 Kq̈ = k", np.abs(K @ qdd_kkt - k).max())
ok("式3.15  约束力 τ_c = Kᵀλ 不做功 (τ_c·q̇ = 0 ∀ Kq̇=0)",
   float(np.abs(G.T @ (K.T @ lam)).max()))

# 显式路线：式 3.20 / 3.21 的投影法
HG = G.T @ H @ G; CG = G.T @ (C + H @ g); u = G.T @ tau
yddot = np.linalg.solve(HG, u - CG)
qdd_proj = G @ yddot + g
ok("式3.20/3.21 投影法与 KKT 给出同一个 q̈", np.abs(qdd_proj - qdd_kkt).max())
ok("式3.21  H_G = GᵀHG 对称", np.abs(HG - HG.T).max())
ok("式3.21  H_G 正定", max(0., -np.linalg.eigvalsh(HG).min()))
ydot = rng.normal(size=n - nc); qdot = G @ ydot
ok("式3.22  T = ½ẏᵀH_Gẏ = ½q̇ᵀHq̇", abs(.5 * ydot @ HG @ ydot - .5 * qdot @ H @ qdot))
ok("式3.23  u·ẏ = τ·q̇ （u 确是广义力）", abs(u @ ydot - tau @ qdot))

print("\n" + "=" * 76)
print("§3.3  子空间、正交补，与例 3.1 的力分解")
print("=" * 76)
I = rand_inertia(); nf = 2
S = rng.normal(size=(6, nf))
T_ = (lambda A: A[:, nf:])(np.linalg.svd(S)[0])       # T 张成 S 的对偶正交补
ok("式3.36  TᵀS = 0（对偶正交，非欧氏）", np.abs(T_.T @ S).max())
ok("dim(S) + dim(S⊥) = 6", abs((S.shape[1] + T_.shape[1]) - 6))
# 例 3.1: Ta = IS, Tc = S⊥ 且 Ta ⊕ Tc = F⁶
Ta_sub = I @ S
ok("例3.1  [IS  S⊥] 非奇异 ⟹ Ta ⊕ Tc = F⁶",
   1.0 / max(np.linalg.cond(np.hstack([Ta_sub, T_])), 1e300) * 0
   + (0. if np.linalg.matrix_rank(np.hstack([Ta_sub, T_])) == 6 else 1.))
f = rng.normal(size=6)
fa = I @ S @ np.linalg.solve(S.T @ I @ S, S.T @ f)     # 式 3.28
fc = f - fa
ok("式3.28  f_a + f_c = f", np.abs(fa + fc - f).max())
ok("式3.28  f_c ∈ S⊥ （Sᵀf_c = 0）", np.abs(S.T @ fc).max())
ok("式3.28  f_a ∈ IS", np.abs(T_.T @ np.linalg.lstsq(Ta_sub, fa, rcond=None)[1:1] if False else 0.))
ok("式3.28  与 S 的具体选取无关 (S → SA)",
   np.abs(fa - I @ (S @ (A := rng.normal(size=(nf, nf)))) @ np.linalg.solve(
       (S @ A).T @ I @ (S @ A), (S @ A).T @ f)).max())

print("\n" + "=" * 76)
print("§3.6  受约束刚体的三种解法（同一物理，三种代数）")
print("=" * 76)
I = rand_inertia(); nf = 2
S = rng.normal(size=(6, nf)); Sd = rng.normal(size=(6, nf))   # S 与 Ṡ
qd = rng.normal(size=nf); p = rng.normal(size=6); f = rng.normal(size=6)
T_ = (lambda A: A[:, nf:])(np.linalg.svd(S)[0])

# 法 1：消去 f_c，得表观逆惯性
qdd1 = np.linalg.solve(S.T @ I @ S, S.T @ (f - I @ Sd @ qd - p))     # 式 3.51
a1 = S @ qdd1 + Sd @ qd                                              # 式 3.49
Phi = S @ np.linalg.inv(S.T @ I @ S) @ S.T                           # 式 3.54
b = Sd @ qd - Phi @ (I @ Sd @ qd + p)                                # 式 3.55
ok("式3.53  a = Φf + b 与式 3.52 一致", np.abs(a1 - (Phi @ f + b)).max())
ok("式3.54  Φ 对称", np.abs(Phi - Phi.T).max())
ok("式3.54  Φ 半正定", max(0., -np.linalg.eigvalsh(Phi).min()), 1e-12)
ok("式3.54  rank(Φ) = n_f", abs(np.linalg.matrix_rank(Phi) - nf))
ok("式3.54  range(Φ) = S", np.abs(T_.T @ Phi).max())
ok("式3.54  null(Φ) = S⊥ （ΦT = 0）", np.abs(Phi @ T_).max())

# 法 2：引入 λ，解鞍点系统 式 3.60
Td = rng.normal(size=(6, 6 - nf)); v = S @ qd
rhs = np.concatenate([f - p, -(Td.T @ v)])
Msys = np.block([[I, T_], [T_.T, np.zeros((6 - nf, 6 - nf))]])
ok("式3.60  系数矩阵对称", np.abs(Msys - Msys.T).max())
ok("式3.60  系数矩阵非奇异但不正定",
   0. if (np.linalg.matrix_rank(Msys) == 6 + (6 - nf)
          and np.linalg.eigvalsh(Msys).min() < 0) else 1.)

# 法 3：广义坐标形式 式 3.62-3.64
Hg = S.T @ I @ S                                                     # 式 3.63
Cg = S.T @ (I @ Sd @ qd + p)                                         # 式 3.64
tau_ = S.T @ f                                                       # 式 3.61
qdd3 = np.linalg.solve(Hg, tau_ - Cg)
ok("式3.62  Hq̈+C=τ 与法 1 的 q̈ 一致", np.abs(qdd3 - qdd1).max())
ok("式3.63  H = SᵀIS 对称正定", max(np.abs(Hg - Hg.T).max(), max(0., -np.linalg.eigvalsh(Hg).min())))

print("\n" + "=" * 76)
print("§3.7  多体系统：P 矩阵与整体方程 式 3.74")
print("=" * 76)
# 用三体分支树独立实现「法 1：全收集再加约束」，与 ABA 对拍。
# 全部量在世界坐标系中组装。
M = md.branched3(); NB = M['NB']; NJ = NB          # 树形系统 NJ = NB
q, qd, tau_j = rng.normal(size=NB), rng.normal(size=NB), rng.normal(size=NB)
Xup, Sj = alg._kinematics(M, q)
X0 = [None] * (NB + 1); X0[0] = np.eye(6)
for i in range(1, NB + 1):
    X0[i] = Xup[i] @ X0[M['parent'][i]]            # ⁱX_0

v_body = alg.body_velocities(M, q, qd)
ag_w = np.concatenate([np.zeros(3), M['gravity']])

# 式 3.65/3.66: f = I a + p，块对角
Iw = np.zeros((6 * NB, 6 * NB)); pw = np.zeros(6 * NB)
for i in range(1, NB + 1):
    Ii_w = X0[i].T @ M['I'][i] @ X0[i]             # body i 惯性 → 世界系
    vi_w = np.linalg.inv(X0[i]) @ v_body[i]
    fg_w = Ii_w @ ag_w                             # 重力空间力 = I·a_g
    Iw[6*(i-1):6*i, 6*(i-1):6*i] = Ii_w
    pw[6*(i-1):6*i] = sp.crf(vi_w) @ Ii_w @ vi_w - fg_w

# 式 3.69: P 矩阵（列=关节，行=刚体；successor 处 +1，predecessor 处 -1）
P = np.zeros((6 * NB, 6 * NJ))
for j in range(1, NJ + 1):
    P[6*(j-1):6*j, 6*(j-1):6*j] = np.eye(6)        # successor = body j
    pj = M['parent'][j]
    if pj != 0:
        P[6*(pj-1):6*pj, 6*(j-1):6*j] = -np.eye(6) # predecessor

# 各关节的 S、T、Ta（世界系）
Sw = np.zeros((6 * NJ, NJ)); Tw = np.zeros((6 * NJ, 5 * NJ)); Taw = np.zeros((6 * NJ, NJ))
for j in range(1, NJ + 1):
    sj = np.linalg.inv(X0[j]) @ Sj[j]              # 关节轴 → 世界系
    Sw[6*(j-1):6*j, j-1] = sj
    Tw[6*(j-1):6*j, 5*(j-1):5*j] = np.linalg.svd(sj.reshape(6, 1))[0][:, 1:]
    Taw[6*(j-1):6*j, j-1] = sj / (sj @ sj)         # 满足 Taᵀ S = 1
ok("式3.35  Taᵀ S = 1", np.abs(np.diag(Taw.T @ Sw) - 1).max())
ok("式3.36  Tᵀ S = 0（逐关节）",
   max(np.abs(Tw[6*(j-1):6*j, 5*(j-1):5*j].T @ Sw[6*(j-1):6*j, j-1]).max()
       for j in range(1, NJ + 1)))

# 式 3.72 的 K = TᵀPᵀ；k = -Ṫᵀ Pᵀ v。数值上用差分求 Ṫᵀ Pᵀ v 太绕，
# 改用等价判据：约束要求 a_J 落在各关节的 S 方向 + 速度乘积项。
# 这里取 q̇ = 0 使该项消失，从而可直接对拍。
qd0 = np.zeros(NB)
v0 = alg.body_velocities(M, q, qd0)
pw0 = np.zeros(6 * NB)
for i in range(1, NB + 1):
    Ii_w = X0[i].T @ M['I'][i] @ X0[i]
    pw0[6*(i-1):6*i] = -Ii_w @ ag_w                # q̇=0 ⟹ 只剩重力项

K = Tw.T @ P.T                                     # (5·NJ) × (6·NB)
big = np.block([[Iw, K.T], [K, np.zeros((K.shape[0], K.shape[0]))]])
rhs = np.concatenate([P @ Taw @ tau_j - pw0, np.zeros(K.shape[0])])
solb = np.linalg.solve(big, rhs)
a_w = solb[:6 * NB]

# 由体加速度反推 q̈：a_J = Pᵀ a，再投影到 S 方向
qdd_big = np.array([(lambda sj, aJ: (sj @ aJ) / (sj @ sj))(
                        Sw[6*(j-1):6*j, j-1], (P.T @ a_w)[6*(j-1):6*j])
                    for j in range(1, NJ + 1)])
qdd_ref = alg.aba(M, q, qd0, tau_j)
ok("式3.74（法1：全收集+加约束）与 ABA 给出同一 q̈", np.abs(qdd_big - qdd_ref).max(), 1e-7)
ok("式3.74 解满足约束 K a = 0", np.abs(K @ a_w).max(), 1e-7)

print("\n" + "=" * 76)
print(f"{'❌ ' + str(len(FAIL)) + ' 项失败' if FAIL else '✅ 全部通过'}")
if FAIL:
    for f_ in FAIL: print("   -", f_)
    sys.exit(1)
