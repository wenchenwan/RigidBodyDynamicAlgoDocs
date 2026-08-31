"""示例机构模型。

model 字典的字段与书中的 model 结构对应（1-based 下标，索引 0 保留给基座）:
    NB       刚体数
    parent   λ(i)，满足正则编号 λ(i) < i
    jtype    关节类型
    Xtree    X_T(i)，常量树变换
    I        body i 的空间惯性（在 body i 自己的坐标系中，因而是常量）
    gravity  重力加速度 3D 向量
"""
import numpy as np
import spatial as sp


def _check_regular(parent):
    for i in range(1, len(parent)):
        assert parent[i] < i, f'违反正则编号 λ({i})={parent[i]} 不小于 {i}'


def arm2r(m1=2.7, m2=1.3, l1=0.45, r1=0.21, r2=0.17, Iz1=0.09, Iz2=0.04):
    """平面 2R 机械臂。两关节轴均沿 z，坐标系 i 在关节 i 处、x 轴沿连杆。"""
    parent = [None, 0, 1]
    _check_regular([0] + parent[1:])
    return {
        'NB': 2,
        'parent': parent,
        'jtype': [None, 'Rz', 'Rz'],
        'Xtree': [None, np.eye(6), sp.xlt([l1, 0, 0])],
        'I': [None,
              sp.rbi(m1, [r1, 0, 0], np.diag([0.0, 0.0, Iz1])),
              sp.rbi(m2, [r2, 0, 0], np.diag([0.0, 0.0, Iz2]))],
        'gravity': np.array([0.0, -9.81, 0.0]),
        'params': dict(m1=m1, m2=m2, l1=l1, r1=r1, r2=r2, Iz1=Iz1, Iz2=Iz2),
    }


def arm3r_spatial(seed=1):
    """3R 空间机械臂：转轴方向不共面，用于检验一般情形。"""
    rng = np.random.default_rng(seed)
    NB = 3
    I = [None]
    for _ in range(NB):
        A = rng.normal(size=(3, 3))
        Ibar = A @ A.T + 3 * np.eye(3)          # 保证正定
        I.append(sp.rbi(abs(rng.normal()) + 0.5, rng.normal(size=3) * 0.2, Ibar))
    return {
        'NB': NB,
        'parent': [None, 0, 1, 2],
        'jtype': [None, 'Rz', 'Ry', 'Rx'],
        'Xtree': [None, np.eye(6),
                  sp.plux(sp.rotx(0.3), [0.4, 0.0, 0.1]),
                  sp.plux(sp.roty(-0.5), [0.35, 0.05, 0.0])],
        'I': I,
        'gravity': np.array([0.0, 0.0, -9.81]),
    }


def branched3(seed=2):
    """三刚体分支树：body1 为「肩」，body2 与 body3 是两条独立支链。

            body0
              │ joint1
            body1
            ╱     ╲
       joint2      joint3
       body2       body3

    用于展示 H 的结构性零：ν(2)∩ν(3)=∅ ⟹ H[2,3] ≡ 0。
    """
    rng = np.random.default_rng(seed)
    I = [None]
    for _ in range(3):
        A = rng.normal(size=(3, 3))
        I.append(sp.rbi(abs(rng.normal()) + 0.5, rng.normal(size=3) * 0.2,
                        A @ A.T + 3 * np.eye(3)))
    return {
        'NB': 3,
        'parent': [None, 0, 1, 1],
        'jtype': [None, 'Rz', 'Ry', 'Rx'],
        'Xtree': [None, np.eye(6),
                  sp.plux(sp.rotz(0.2), [0.3, 0.1, 0.0]),
                  sp.plux(sp.rotz(-0.6), [0.3, -0.1, 0.0])],
        'I': I,
        'gravity': np.array([0.0, 0.0, -9.81]),
    }


def chain(n, seed=3):
    """n 连杆串联链，用于复杂度/稀疏性对比（树深度 d = n）。"""
    rng = np.random.default_rng(seed)
    I = [None]
    for _ in range(n):
        A = rng.normal(size=(3, 3))
        I.append(sp.rbi(abs(rng.normal()) + 0.5, rng.normal(size=3) * 0.2,
                        A @ A.T + 3 * np.eye(3)))
    return {
        'NB': n,
        'parent': [None] + list(range(n)),
        'jtype': [None] + ['Rz' if k % 2 else 'Ry' for k in range(n)],
        'Xtree': [None, np.eye(6)] + [sp.xlt([0.3, 0, 0]) for _ in range(n - 1)],
        'I': I,
        'gravity': np.array([0.0, 0.0, -9.81]),
    }


def nu(model, i):
    """子树 ν(i)：以 i 为根（含 i）的所有 body。"""
    out, stack = set(), [i]
    while stack:
        k = stack.pop(); out.add(k)
        stack += [j for j in range(1, model['NB'] + 1) if model['parent'][j] == k]
    return out


def kappa(model, k):
    """支撑集 κ(k)：从 k 到根路径上的所有关节。"""
    out, j = set(), k
    while j != 0:
        out.add(j); j = model['parent'][j]
    return out
