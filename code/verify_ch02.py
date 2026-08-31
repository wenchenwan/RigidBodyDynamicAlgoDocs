"""第 2 章重点结论的数值验证——特别是 PDF 中三处批注的困惑点。

运行: python3 code/verify_ch02.py
配套 docs/ch02-spatial-vector-algebra.md
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, spatial as sp

rng = np.random.default_rng(23)
FAIL = []
def ok(n, e, t=1e-10):
    good = e < t
    if not good: FAIL.append(n)
    print(f"{'✅' if good else '❌'} {n:<52} err={e:.2e}")

print("=" * 74)
print("§2.5  线向量 + 自由向量的分解（PDF p23 批注：这里的理解仍然有困惑）")
print("=" * 74)
s_hat = rng.normal(size=6); s, sO = s_hat[:3], s_hat[3:]

# (a) 要求线向量过指定点 P —— 分解唯一
P = rng.normal(size=3)
line = np.concatenate([s, np.cross(P, s)])
free = np.concatenate([np.zeros(3), sO - np.cross(P, s)])
ok("过指定点 P 的分解求和 = 原向量", np.abs(line + free - s_hat).max())
ok("  线部分满足线向量判据 s·s_O = 0", abs(line[:3] @ line[3:]))
ok("  自由部分满足自由向量判据 s = 0", np.abs(free[:3]).max())

# (b) 要求自由向量平行于 s —— Chasles 螺旋分解，也唯一
h = (s @ sO) / (s @ s)
line2 = np.concatenate([s, sO - h * s])
free2 = np.concatenate([np.zeros(3), h * s])
ok("平行（螺旋）分解求和 = 原向量", np.abs(line2 + free2 - s_hat).max())
ok("  线部分满足 s·s_O = 0", abs(line2[:3] @ line2[3:]))
ok("  自由部分平行于 s", float(np.linalg.norm(np.cross(free2[3:], s))))
axis_pt = np.cross(s, sO) / (s @ s)
ok("  螺旋轴过点 (s×s_O)/(s·s)", np.abs(np.cross(axis_pt, s) - (sO - h * s)).max())

print("\n" + "=" * 74)
print("§2.13 I = Σ g_i g_i· 中 g_i 的含义（PDF p40 批注：g_i 的实际含义是什么？？？）")
print("=" * 74)
m, c, Ib = 2.4, rng.normal(size=3), np.diag([.4, .6, .9])
I = sp.rbi(m, c, Ib)

L = np.linalg.cholesky(I)
ok("Cholesky 的 6 个列向量就是一组 g_i", np.abs(L @ L.T - I).max())
w, V = np.linalg.eigh(I); G2 = V * np.sqrt(w)
ok("特征分解给出另一组 g_i，同样成立", np.abs(G2 @ G2.T - I).max())
print(f"     ⚠️ 两组 g₁ 完全不同：{L[:,0][:3].round(3)} vs {G2[:,0][:3].round(3)}")
print("        ⟹ g_i 不唯一（任意正交 Q，G 与 GQ 给出同一个 I）")

# 物理含义：质点系。点质量 m_k 在 p_k 处，J_k 把空间速度映到该点线速度
pts, ms = rng.normal(size=(8, 3)), rng.uniform(.1, 1, 8)
I_pts, dyads = np.zeros((6, 6)), []
for pk, mk in zip(pts, ms):
    Jk = np.hstack([-sp.skew(pk), np.eye(3)])          # v_P = J_k v̂
    I_pts += mk * Jk.T @ Jk
    for a in range(3):
        ea = np.eye(3)[a]
        g = np.sqrt(mk) * np.concatenate([np.cross(pk, ea), ea])  # 过 p_k 沿 e_a 的单位力
        dyads.append(np.outer(g, g))
ok("质点系 I = Σ_k m_k J_kᵀJ_k = Σ g g ᵀ", np.abs(I_pts - sum(dyads)).max())
M = ms.sum(); C = (ms[:, None] * pts).sum(0) / M
Ibar = sum(mk * ((pk - C) @ (pk - C) * np.eye(3) - np.outer(pk - C, pk - C))
           for pk, mk in zip(pts, ms))
ok("  且它确实是合法刚体惯性 rbi(M, C, Ī)", np.abs(I_pts - sp.rbi(M, C, Ibar)).max())
print("     ⟹ 每个 g = √mₖ ×「过质点 pₖ、沿 eₐ 的单位力」的 Plücker 坐标")

v = rng.normal(size=6)
lhs = sp.crf(v) @ I - I @ sp.crm(v)
rhs = sum(np.outer(sp.crf(v) @ L[:, i], L[:, i])
          - np.outer(L[:, i], sp.crm(v).T @ L[:, i]) for i in range(6))
ok("用 dyad 表示推出 İ = v×*I − I v×（式 2.65）", np.abs(lhs - rhs).max())

print("\n" + "=" * 74)
print("§2.13 为什么刚体惯性只需 10 个参数（PDF p42 批注：why???）")
print("=" * 74)
ok("右下块 = m·1（6 个对称元 → 1 个自由，省 5）", np.abs(I[3:, 3:] - m * np.eye(3)).max())
ok("右上块 = m c× 反对称（9 元 → 3 自由，省 6）", np.abs(I[:3, 3:] + I[:3, 3:].T).max())
ok("整体对称 I[3:,:3] = I[:3,3:]ᵀ", np.abs(I[3:, :3] - I[:3, 3:].T).max())
print("     1 (质量) + 3 (质心) + 6 (绕质心惯量) = 10")
print("     一般对称 6×6 有 21 个自由参数，刚体多出 11 个约束：21 − 11 = 10 ✓")
A = rng.normal(size=(6, 6)); IA = A @ A.T
print(f"     对比铰接体惯性（第 7 章），右下块不是 m·1：\n     {np.round(IA[3:,3:],2).tolist()}")
print("     ⟹ 铰接体惯性需要完整的 21 个参数")

print("\n" + "=" * 74)
print("§2.11 空间加速度 vs 经典加速度")
print("=" * 74)
w_ = 2.1
v0 = np.array([0, 0, w_, 0., 0, 0]); a0 = np.zeros(6)   # 绕定轴匀速转动
P = np.array([0.5, 0., 0.])
XP = sp.xlt(P); vP, aP = XP @ v0, XP @ a0
ok("匀速转动：空间加速度 ≡ 0", np.abs(aP).max())
ok("  但 P 点经典加速度 = −ω²r（向心）",
   np.abs(sp.classical_accel(vP, aP) - (-w_**2 * P)).max())
print(f"     空间 a 的线分量 = {aP[3:]}，经典加速度 = {sp.classical_accel(vP,aP)}")

print("\n" + "=" * 74)
print("§2.14 例 2.6：在质心系展开 f = Ia + v×*Iv，还原牛顿+欧拉方程")
print("=" * 74)
m, Ibar = 2.4, np.diag([.31, .52, .44]); I = sp.rbi(m, [0, 0, 0], Ibar)
w_, vC = rng.normal(size=3), rng.normal(size=3)
wd, alin = rng.normal(size=3), rng.normal(size=3)
v = np.concatenate([w_, vC]); a = np.concatenate([wd, alin])
f = I @ a + sp.crf(v) @ I @ v
ok("力分量 = m·c̈（牛顿，式 2.70）", np.abs(f[3:] - m * (alin + np.cross(w_, vC))).max())
ok("矩分量 = Ī ω̇ + ω×Īω（欧拉，式 2.71）",
   np.abs(f[:3] - (Ibar @ wd + np.cross(w_, Ibar @ w_))).max())

print("\n" + "=" * 74)
print("§2.15 逆惯性 Φ = I⁻¹")
print("=" * 74)
I = sp.rbi(2.4, rng.normal(size=3), np.diag([.4, .6, .9]))
Phi = np.linalg.inv(I)
ok("Φ 对称", np.abs(Phi - Phi.T).max())
ok("Φ I = 1", np.abs(Phi @ I - np.eye(6)).max(), 1e-9)
print("     无约束刚体：rank(Φ) = 6 = 运动自由度数")
print("     受约束刚体：rank(Φ) < 6，此时 Φ 奇异、I 不存在（式 2.72 仍可用）")

print("\n" + "=" * 74)
print(f"{'❌ ' + str(len(FAIL)) + ' 项失败' if FAIL else '✅ 全部通过'}")
if FAIL:
    for f_ in FAIL: print("   -", f_)
    sys.exit(1)
