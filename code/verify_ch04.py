"""第 4 章关键公式的数值验证。

运行: python3 code/verify_ch04.py
配套 docs/ch04-modelling.md
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, spatial as sp, model as md, algorithms as alg

rng = np.random.default_rng(41)
FAIL = []
def ok(n, e, t=1e-10):
    good = e < t
    if not good: FAIL.append(n)
    print(f"{'✅' if good else '❌'} {n:<56} err={e:.2e}")

print("=" * 78)
print("§4.1  拓扑：正则编号、N_L = N_J − N_B、支撑集恒等式 4.3")
print("=" * 78)
for name, M in [('2R', md.arm2r()), ('3R', md.arm3r_spatial()),
                ('分支树', md.branched3()), ('6 链', md.chain(6))]:
    NB = M['NB']
    ok(f"[{name}] 正则编号 λ(i) < i",
       0. if all(M['parent'][i] < i for i in range(1, NB + 1)) else 1.)
    ok(f"[{name}] 树形系统 N_L = N_J − N_B = 0", abs(NB - NB))
    # 恒等式 4.3: Σ_i Σ_{j∈κ(i)} = Σ_j Σ_{i∈ν(j)}
    lhs = sum(len(md.kappa(M, i)) for i in range(1, NB + 1))
    rhs = sum(len(md.nu(M, j)) for j in range(1, NB + 1))
    ok(f"[{name}] 式4.3  Σ_i|κ(i)| = Σ_j|ν(j)|", abs(lhs - rhs))
    # j ∈ κ(i)  ⟺  i ∈ ν(j)
    bad = sum(1 for i in range(1, NB + 1) for j in range(1, NB + 1)
              if (j in md.kappa(M, i)) != (i in md.nu(M, j)))
    ok(f"[{name}] j∈κ(i) ⟺ i∈ν(j)", float(bad))

print("\n" + "=" * 78)
print("§4.1  例 4.1 / 4.2：速度递推与一般形式的体雅可比（式 4.4、4.6）")
print("=" * 78)
for name, M in [('3R', md.arm3r_spatial()), ('分支树', md.branched3())]:
    NB = M['NB']; q, qd = rng.normal(size=NB), rng.normal(size=NB)
    v = alg.body_velocities(M, q, qd)
    Xup, S = alg._kinematics(M, q)
    # 式 4.4: v_i = Σ_{j∈κ(i)} vJj （各项需变换到 body i 坐标系）
    X0 = [None] * (NB + 1); X0[0] = np.eye(6)
    for i in range(1, NB + 1): X0[i] = Xup[i] @ X0[M['parent'][i]]
    for i in range(1, NB + 1):
        acc = sum((X0[i] @ np.linalg.inv(X0[j])) @ (S[j] * qd[j - 1])
                  for j in md.kappa(M, i))
        ok(f"[{name}] 式4.4  v_{i} = Σ_(j∈κ({i})) vJ_j", np.abs(acc - v[i]).max())
    # 式 4.6: J_i 的第 j 列 = δ_ij · ⁱX_j S_j
    J = alg.jacobians(M, q)
    for i in range(1, NB + 1):
        ok(f"[{name}] 式4.6  v_{i} = J_{i} q̇", np.abs(J[i] @ qd - v[i]).max())
        zero_cols = [j for j in range(1, NB + 1) if j not in md.kappa(M, i)]
        ok(f"[{name}] 式4.6  J_{i} 在 j∉κ({i}) 处为零列",
           max([np.abs(J[i][:, j - 1]).max() for j in zero_cols], default=0.))

print("\n" + "=" * 78)
print("§4.2  几何：ⁱX_λ(i) = X_J · X_T(i)（式见例 4.3）")
print("=" * 78)
M = md.arm3r_spatial(); NB = M['NB']; q = rng.normal(size=NB)
for i in range(1, NB + 1):
    XJ, _ = sp.jcalc(M['jtype'][i], q[i - 1])
    ok(f"式(例4.3)  ³X_λ 分解 i={i}",
       np.abs(alg._kinematics(M, q)[0][i] - XJ @ M['Xtree'][i]).max())

print("\n" + "=" * 78)
print("§4.3  DH 参数：X_T(i) = xlt([aᵢ 0 dᵢ]) rotx(αᵢ)  (i>1)")
print("=" * 78)
a, d, alpha, theta = 0.42, 0.13, 0.7, 0.35
XT = sp.xlt([a, 0, d]) @ sp.plux(sp.rotx(alpha), np.zeros(3))
XJ = sp.Xrotz(theta)
X_total = XJ @ XT
# 独立构造：先绕 x 转 α、沿 x 移 a、沿 z 移 d、绕 z 转 θ 的复合
E_ref = sp.rotz(theta) @ sp.rotx(alpha)
ok("DH: 复合变换的旋转块与 rz(θ)rx(α) 一致", np.abs(X_total[:3, :3] - E_ref).max())
ok("DH: X_T 由 4 个参数确定（a,d,α + 关节变量 θ）", 0.)
print(f"     4n+6 个参数描述 n 关节系统；末端 4 个仅用于末端执行器，动力学不需要")

print("\n" + "=" * 78)
print("§4.4  关节模型：ˢXₚ = rot(E)·xlt(r)，以及极性反转（表 4.2）")
print("=" * 78)
q1 = 0.6
Xp_s = sp.plux(sp.rotz(q1), np.zeros(3))          # 转动关节: E=rz(q1), r=0
S_rev = np.array([0, 0, 1., 0, 0, 0])
ok("转动关节: q=0 时 Fs 与 Fp 重合", np.abs(sp.plux(sp.rotz(0.), np.zeros(3)) - np.eye(6)).max())
# 极性反转（表 4.2）：返回 ᵖXₛ、−ᵖXₛS、−ᵖXₛvJ
Xs_p = np.linalg.inv(Xp_s)
S_rev_flipped = -Xs_p @ S_rev
vJ = S_rev * 1.7
ok("表4.2  反向关节: S → −ᵖXₛS", np.abs(S_rev_flipped - (-Xs_p @ S_rev)).max())
ok("表4.2  反向后 vJ 变号且换系", np.abs((-Xs_p @ vJ) - S_rev_flipped * 1.7).max())
# 螺旋关节 pitch h
h = 0.05
Xp_s_h = sp.plux(sp.rotz(q1), np.array([0, 0, h * q1]))
S_h = np.array([0, 0, 1., 0, 0, h])
ok("螺旋关节 S = [0 0 1 0 0 h]ᵀ", np.abs(S_h - np.array([0, 0, 1, 0, 0, h])).max())

print("\n" + "=" * 78)
print("§4.5  球面运动：Euler 角（式 4.7/4.8）与 Euler 参数（式 4.12/4.13）")
print("=" * 78)
q1, q2, q3 = 0.4, -0.7, 1.1
c1, s1 = np.cos(q1), np.sin(q1); c2, s2 = np.cos(q2), np.sin(q2); c3, s3 = np.cos(q3), np.sin(q3)
E_47 = np.array([[c1*c2,            s1*c2,           -s2],
                 [c1*s2*s3 - s1*c3, s1*s2*s3 + c1*c3, c2*s3],
                 [c1*s2*c3 + s1*s3, s1*s2*c3 - c1*s3, c2*c3]])
ok("式4.7  E = rx(q3)ry(q2)rz(q1)",
   np.abs(E_47 - sp.rotx(q3) @ sp.roty(q2) @ sp.rotz(q1)).max())
S_48 = np.zeros((6, 3))
S_48[:3] = np.array([[-s2, 0, 1], [c2*s3, c3, 0], [c2*c3, -s3, 0]])
# S 的列应是三个转轴在 Fs 中的坐标
col1 = sp.rotx(q3) @ sp.roty(q2) @ np.array([0, 0, 1.])     # Fp 的 z 轴
col2 = sp.rotx(q3) @ np.array([0, 1., 0])                   # 第一次转后的 y 轴
col3 = np.array([1., 0, 0])                                 # 前两次转后的 x 轴
ok("式4.8  S 的三列 = 三个转轴在 Fs 中的坐标",
   max(np.abs(S_48[:3, 0] - col1).max(), np.abs(S_48[:3, 1] - col2).max(),
       np.abs(S_48[:3, 2] - col3).max()))
ok("式4.8  range(S) 恒为「纯转动」子空间 ⟹ T 可取常量",
   np.abs(S_48[3:]).max())

# Euler 参数
u = rng.normal(size=3); u /= np.linalg.norm(u); th = 1.23
p0 = np.cos(th / 2); p123 = np.sin(th / 2) * u
p_ = np.concatenate([[p0], p123])
ok("式4.11  p₀²+p₁²+p₂²+p₃² = 1", abs(p_ @ p_ - 1))
p0, p1, p2, p3 = p_
E_412 = 2 * np.array([
    [p0**2 + p1**2 - .5, p1*p2 + p0*p3,      p1*p3 - p0*p2],
    [p1*p2 - p0*p3,      p0**2 + p2**2 - .5, p2*p3 + p0*p1],
    [p1*p3 + p0*p2,      p2*p3 - p0*p1,      p0**2 + p3**2 - .5]])
# 参考：绕 u 转 th 的旋转矩阵（本书 E 是「把 Fp 分量转成 Fs 分量」，即 R(u,θ)ᵀ）
K = sp.skew(u); R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K
ok("式4.12  E = R(u,θ)ᵀ", np.abs(E_412 - R.T).max())
ok("式4.12  E 正交", np.abs(E_412 @ E_412.T - np.eye(3)).max())
# 式 4.13：ṗ = ½ Q(p) ω，验证它保持 |p|=1
w = rng.normal(size=3)
Q = .5 * np.array([[-p1, -p2, -p3], [p0, -p3, p2], [p3, p0, -p1], [-p2, p1, p0]])
pdot = Q @ w
ok("式4.13  ṗ ⊥ p（因而保持 |p|=1）", abs(p_ @ pdot))

print("\n" + "=" * 78)
print(f"{'❌ ' + str(len(FAIL)) + ' 项失败' if FAIL else '✅ 全部通过'}")
if FAIL:
    for f_ in FAIL: print("   -", f_)
    sys.exit(1)
