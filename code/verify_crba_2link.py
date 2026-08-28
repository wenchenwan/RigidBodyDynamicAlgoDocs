"""
2-连杆平面机械臂：CRBA 与 H_ij 公式的数值验证。

配套 docs/ch06-forward-dynamics-crba.md 的「$H$ 的元素：完整解析」一节。
纯 numpy，无外部依赖。运行：python3 code/verify_crba_2link.py

验证四件事：
  1. CRBA 算出的 H == 教材上 2R 机械臂的解析 H
  2. H_ij = S_i' Ic_i  iX_j  S_j 的下标条件是「j 是 i 的祖先，Ic 取后代的」
     （反过来写会得到一个对称、量纲正确、但错误的数）
  3. H == sum_k J_k' I_k J_k                （定义 C：Jacobian 视角）
  4. H == d^2T/(dq̇_i dq̇_j)                 （定义 A：动能 Hessian，数值二阶差分）

第 4 项完全不碰 Ic / X* / 下标条件，是唯一能可靠检出「下标写反」的对拍方式。

约定与书一致：空间向量分量排列为 (角, 线)，即 v = [ω; v_O]、f = [n_O; f]。
"""
import numpy as np
np.set_printoptions(precision=10, suppress=True)

def skew(v):
    x,y,z = v
    return np.array([[0,-z,y],[z,0,-x],[-y,x,0]])

def Xrotz(t):
    c,s = np.cos(t), np.sin(t)
    E = np.array([[c,s,0],[-s,c,0],[0,0,1]])       # rotates parent coords into child coords
    X = np.zeros((6,6)); X[:3,:3]=E; X[3:,3:]=E
    return X

def Xtrans(r):
    X = np.eye(6); X[3:,:3] = -skew(r)
    return X

def rbi(m, c, Ibar):
    C = skew(c)
    I = np.zeros((6,6))
    I[:3,:3] = Ibar + m*C@C.T
    I[:3,3:] = m*C
    I[3:,:3] = m*C.T
    I[3:,3:] = m*np.eye(3)
    return I

# ---- 2-link planar arm ----
m1,m2 = 2.7, 1.3
l1     = 0.45
r1,r2  = 0.21, 0.17
Iz1,Iz2 = 0.09, 0.04

def build(q1,q2):
    S = [None, np.array([0,0,1,0,0,0.0]), np.array([0,0,1,0,0,0.0])]
    parent = [None, 0, 1]
    Xtree = [None, np.eye(6), Xtrans(np.array([l1,0,0.0]))]
    XJ    = [None, Xrotz(q1), Xrotz(q2)]
    Xup   = [None] + [XJ[i] @ Xtree[i] for i in (1,2)]
    Ibar1 = np.diag([0.0,0.0,Iz1]); Ibar2 = np.diag([0.0,0.0,Iz2])
    I = [None, rbi(m1, np.array([r1,0,0.0]), Ibar1),
               rbi(m2, np.array([r2,0,0.0]), Ibar2)]
    return S, parent, Xup, I

def crba(q1,q2):
    S,parent,Xup,I = build(q1,q2)
    NB = 2
    Ic = [None] + [I[i].copy() for i in (1,2)]
    H = np.zeros((NB,NB))
    for i in range(NB,0,-1):
        if parent[i] != 0:
            Ic[parent[i]] += Xup[i].T @ Ic[i] @ Xup[i]
        F = Ic[i] @ S[i]
        H[i-1,i-1] = S[i] @ F
        j = i
        while parent[j] != 0:
            F = Xup[j].T @ F
            j = parent[j]
            H[i-1,j-1] = F @ S[j]
            H[j-1,i-1] = H[i-1,j-1]
    return H, Ic, Xup, S

def H_analytic(q2):
    H11 = Iz1 + m1*r1**2 + Iz2 + m2*(l1**2 + r2**2 + 2*l1*r2*np.cos(q2))
    H12 = Iz2 + m2*(r2**2 + l1*r2*np.cos(q2))
    H22 = Iz2 + m2*r2**2
    return np.array([[H11,H12],[H12,H22]])

for q2 in (0.0, 0.7, -1.3, 2.4):
    H,_,_,_ = crba(0.3, q2)
    Ha = H_analytic(q2)
    print(f"q2={q2:+.2f}  max|CRBA - 解析| = {np.abs(H-Ha).max():.3e}")

print("\n--- 检验 H_ij = S_i^T Ic_i  iX_j  S_j  的下标条件 ---")
q1,q2 = 0.3, 0.7
H,Ic,Xup,S = crba(q1,q2)
X21 = Xup[2]                       # ^2X_1
# 候选 A: i=2(后代), j=1(祖先)  ->  S_2^T Ic_2 ^2X_1 S_1
candA = S[2] @ Ic[2] @ X21 @ S[1]
# 候选 B: i=1(祖先), j=2(后代)  ->  S_1^T Ic_1 ^1X_2 S_2
X12 = np.linalg.inv(X21)
candB = S[1] @ Ic[1] @ X12 @ S[2]
print(f"H[1,2] (真值)                       = {H[0,1]:.10f}")
print(f"用 Ic_后代, j=祖先  (S2' Ic2 X21 S1) = {candA:.10f}   <-- 正确")
print(f"用 Ic_祖先, j=后代  (S1' Ic1 X12 S2) = {candB:.10f}   <-- 错误")

print("\n--- 检验 ^2X_1 S_1 的解析形式 [0,0,1, l1 sin q2, l1 cos q2, 0] ---")
print("数值 :", X21 @ S[1])
print("解析 :", np.array([0,0,1, l1*np.sin(q2), l1*np.cos(q2), 0]))

print("\n--- 检验 H = sum_k J_k^T I_k J_k ---")
# J_k 的第 l 列 = ^kX_l S_l  若 l 在 k 的支撑集上
_,parent,Xup,I = build(q1,q2)
X = {1: Xup[1], 2: Xup[2]@Xup[1]}          # ^kX_0
Hj = np.zeros((2,2))
for k in (1,2):
    J = np.zeros((6,2))
    for l in (1,2):
        if l<=k:                            # 串联链: 支撑集 = {1..k}
            Xkl = X[k] @ np.linalg.inv(X[l])
            J[:,l-1] = Xkl @ S[l]
    Hj += J.T @ I[k] @ J
print("max|H - sum J'IJ| =", np.abs(H-Hj).max())

print("\n--- 检验 H_ij = d^2T/dq̇_i dq̇_j (数值二阶差分) ---")
def kinetic(qd):
    _,parent,Xup,I = build(q1,q2)
    v = [None]*3; v[0]=np.zeros(6); T=0.0
    for i in (1,2):
        v[i] = Xup[i] @ v[parent[i]] + S[i]*qd[i-1]
        T += 0.5 * v[i] @ I[i] @ v[i]
    return T
h=1e-5; Hn=np.zeros((2,2))
for a in range(2):
    for b in range(2):
        ea=np.zeros(2); ea[a]=h; eb=np.zeros(2); eb[b]=h
        Hn[a,b]=(kinetic(ea+eb)-kinetic(ea-eb)-kinetic(-ea+eb)+kinetic(-ea-eb))/(4*h*h)
print("max|H - d2T/dqd2| =", np.abs(H-Hn).max())
