"""
空间向量代数核心运算 (Featherstone 约定)。

约定
----
分量排列为 (角, 线):
    运动向量 v = [ω; v_O] ∈ M⁶
    力向量   f = [n_O; f] ∈ F⁶
配对 f·v = fᵀv 给出功率。

坐标变换
--------
X  = ᴮX_A  作用于运动向量: ᴮv = X ᴬv
X* = ᴮX*_A 作用于力向量:   ᴮf = X* ᴬf,  且 X* = X⁻ᵀ

配套 docs/ch02-spatial-vector-algebra.md。纯 numpy。
"""
import numpy as np


# ────────────────────────── 3D 辅助 ──────────────────────────

def skew(v):
    """3D 向量的反对称矩阵，满足 skew(a) @ b == cross(a, b)。"""
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotx(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def roty(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])


def rotz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])


# ─────────────────── 6D Plücker 变换 (运动向量) ───────────────────

def plux(E, r):
    """由 3×3 旋转 E 和位置 r 构造运动向量变换 ᴮX_A。

    E = ᴮR_A 把 A 系分量转为 B 系分量；r 是 B 系原点在 A 系中的位置。
    """
    X = np.zeros((6, 6))
    X[:3, :3] = E
    X[3:, :3] = -E @ skew(r)
    X[3:, 3:] = E
    return X


def xlt(r):
    """纯平移的运动向量变换。"""
    return plux(np.eye(3), np.asarray(r, dtype=float))


def Xrotx(t):
    return plux(rotx(t), np.zeros(3))


def Xroty(t):
    return plux(roty(t), np.zeros(3))


def Xrotz(t):
    return plux(rotz(t), np.zeros(3))


def Xstar(X):
    """由运动向量变换得到对应的力向量变换: X* = X⁻ᵀ。"""
    return np.linalg.inv(X).T


# ────────────────────────── 空间叉乘 ──────────────────────────

def crm(v):
    """v× : 作用于运动向量 (M⁶ → M⁶)。"""
    w, vo = v[:3], v[3:]
    M = np.zeros((6, 6))
    M[:3, :3] = skew(w)
    M[3:, :3] = skew(vo)
    M[3:, 3:] = skew(w)
    return M


def crf(v):
    """v×* : 作用于力向量 (F⁶ → F⁶)，等于 -(v×)ᵀ。"""
    return -crm(v).T


# ────────────────────────── 空间惯性 ──────────────────────────

def rbi(m, c, Ibar):
    """刚体空间惯性。m 质量, c 质心位置(相对参考点), Ibar 绕质心的 3×3 转动惯量。"""
    c = np.asarray(c, dtype=float)
    C = skew(c)
    I = np.zeros((6, 6))
    I[:3, :3] = Ibar + m * C @ C.T
    I[:3, 3:] = m * C
    I[3:, :3] = m * C.T
    I[3:, 3:] = m * np.eye(3)
    return I


def transform_inertia(X, I):
    """把在 A 系中表示的惯性 I 变换到 B 系，X = ᴮX_A。

    ᴮI = ᴮX*_A ᴬI ᴬX_B = X⁻ᵀ I X⁻¹
    """
    Xi = np.linalg.inv(X)
    return Xi.T @ I @ Xi


def parent_from_child(X, I):
    """把子坐标系中的惯性搬到父坐标系。X = ⁱX_λ(i)，I 在 i 系中。

    等价于 transform_inertia(inv(X), I)，但避免求逆：λI = Xᵀ I X。
    """
    return X.T @ I @ X


# ────────────────────────── 关节模型 ──────────────────────────

_JOINT_S = {
    'Rx': np.array([1.0, 0, 0, 0, 0, 0]),
    'Ry': np.array([0, 1.0, 0, 0, 0, 0]),
    'Rz': np.array([0, 0, 1.0, 0, 0, 0]),
    'Px': np.array([0, 0, 0, 1.0, 0, 0]),
    'Py': np.array([0, 0, 0, 0, 1.0, 0]),
    'Pz': np.array([0, 0, 0, 0, 0, 1.0]),
}


def jcalc(jtype, q):
    """返回 (X_J, S)。常见 1-DoF 关节在关节坐标系中 S 为常量、c_J = 0。"""
    if jtype == 'Rx':
        return Xrotx(q), _JOINT_S['Rx']
    if jtype == 'Ry':
        return Xroty(q), _JOINT_S['Ry']
    if jtype == 'Rz':
        return Xrotz(q), _JOINT_S['Rz']
    if jtype == 'Px':
        return xlt([q, 0, 0]), _JOINT_S['Px']
    if jtype == 'Py':
        return xlt([0, q, 0]), _JOINT_S['Py']
    if jtype == 'Pz':
        return xlt([0, 0, q]), _JOINT_S['Pz']
    raise ValueError(f'unknown joint type: {jtype}')


# ─────────────────────── 经典 / 空间加速度 ───────────────────────

def classical_accel(v, a):
    """由空间加速度 a 得到「与原点重合的物质点」的经典加速度。

    a_classical = a_linear + ω × v_O
    """
    w, vo = v[:3], v[3:]
    return a[3:] + np.cross(w, vo)
