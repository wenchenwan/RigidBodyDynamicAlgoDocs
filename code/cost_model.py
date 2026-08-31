"""第 10 章表 10.1 的运算次数公式，以及 ABA / O(nd²) 路线的交叉点分析。

运行: python3 code/cost_model.py
配套 docs/ch10-accuracy-efficiency.md

基准系统集合是原书的 GU(n)：无分支、全转动关节、一般惯性与几何参数的链。
正弦/余弦计算按原书惯例不计入（对需要它们的算法都一样）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ── 表 10.1 中每个算法的最佳已发表版本 ──────────────────────────
def rnea(n):
    """Balafoutis & Patel (1991), Algorithm 5.7"""
    return (93 * n - 108, 81 * n - 100)


def crba(n):
    """Featherstone (2005) 式 20 / 原书式 10.3"""
    return (10 * n * n + 22 * n - 32, 6 * n * n + 37 * n - 43)


def factor_solve(n):
    """F&S：分解并求解 Hq̈ = τ − C（原书表 6.6 配式 6.28）"""
    return (n**3 / 6 + 1.5 * n * n - 2 * n / 3, n**3 / 6 + n * n - 7 * n / 6)


def aba(n):
    """McMillan et al. (1995), 表 II（固定基）"""
    return (224 * n - 259, 205 * n - 248)


def crba_general_tree(D0_dh, D0_g, D1_dh, D1_g):
    """原书式 10.12：一般运动学树（全转动关节）上的 CRBA 代价。

    代价原语来自附录 A.3/A.4：
        ra      = 10a                 （惯性相加）
        rx_dh   = 32m + 33a           （DH 变换下的 λX*I^c X）
        rx_g    = 47m + 48a           （一般变换）
        vx_dh   = 20m + 12a           （DH 变换下的 λX* F）
        vx_g    = 24m + 18a
    """
    m = 32 * D0_dh + 47 * D0_g + 20 * D1_dh + 24 * D1_g
    a = (10 + 33) * D0_dh + (10 + 48) * D0_g + 12 * D1_dh + 18 * D1_g
    return (m, a)


def D0_D1(model):
    """原书式 10.4、10.5：D0 = n − |µ(0)|，D1 = Σ(|κ(i)| − 1)。"""
    import model as md
    n = model['NB']
    mu0 = sum(1 for i in range(1, n + 1) if model['parent'][i] == 0)
    D0 = n - mu0
    D1 = sum(len(md.kappa(model, i)) - 1 for i in range(1, n + 1))
    return D0, D1


def _tot(f, n):
    m, a = f(n)
    return m + a


def crossover(lo=2, hi=40):
    """找出 O(nd²) 路线（无分支时即 O(n³)）与 ABA 的代价交叉点。"""
    prev = None
    for n in range(lo, hi + 1):
        o3 = _tot(rnea, n) + _tot(crba, n) + _tot(factor_solve, n)
        ab = _tot(aba, n)
        cur = 'O(n³)' if o3 < ab else 'ABA'
        if prev and cur != prev:
            return n
        prev = cur
    return None


if __name__ == '__main__':
    print("=" * 72)
    print("表 10.1：GU(n) 上的运算次数（乘法 m + 加法 a，正余弦不计）")
    print("=" * 72)
    print(f"{'n':>3} {'RNEA':>8} {'CRBA':>9} {'F&S':>9} {'O(n³)合计':>11} {'ABA':>9} {'比值':>7}  谁快")
    for n in [2, 4, 6, 7, 8, 9, 10, 12, 15, 18, 20, 30]:
        r, c, f, ab = _tot(rnea, n), _tot(crba, n), _tot(factor_solve, n), _tot(aba, n)
        o3 = r + c + f
        print(f"{n:>3} {r:>8.0f} {c:>9.0f} {f:>9.0f} {o3:>11.0f} {ab:>9.0f} "
              f"{o3/ab:>7.3f}  {'O(n³)' if o3 < ab else 'ABA'}")

    x = crossover()
    print(f"\n➜ 交叉点：n ≤ {x-1} 时 O(n³) 路线更快，n ≥ {x} 时 ABA 更快")
    print(f"   （原书：'the O(n³) algorithm is slightly faster than ABA for n ≤ 8,")
    print(f"     and has risen to about 1.6 times the cost of ABA by n = 18'）")
    print(f"   n=18 实测比值 = {(_tot(rnea,18)+_tot(crba,18)+_tot(factor_solve,18))/_tot(aba,18):.3f} ✓")

    print("\n" + "=" * 72)
    print("O(n³) 路线内部：哪一项主导？")
    print("=" * 72)
    for n in [5, 8, 10, 20, 30, 45, 50]:
        r, c, f = _tot(rnea, n), _tot(crba, n), _tot(factor_solve, n)
        dom = max([('RNEA', r), ('CRBA', c), ('F&S', f)], key=lambda t: t[1])[0]
        print(f"  n={n:>3}: RNEA {r/(r+c+f)*100:>5.1f}%  CRBA {c/(r+c+f)*100:>5.1f}%  "
              f"F&S {f/(r+c+f)*100:>5.1f}%   主导={dom}")
    print("   原书：CRBA 从 n=8 起主导，直到远超 n=30；F&S 曲线约在 n=45 处与 CRBA 相交 ✓")

    print("\n" + "=" * 72)
    print("分支的影响（式 10.4-10.12）：本仓库四个示例模型")
    print("=" * 72)
    import model as md
    for name, M in [('2R 平面臂', md.arm2r()), ('3R 空间臂', md.arm3r_spatial()),
                    ('三体分支树', md.branched3()), ('6 连杆串联链', md.chain(6)),
                    ('12 连杆串联链', md.chain(12))]:
        n = M['NB']; D0, D1 = D0_D1(M)
        d = max(len(md.kappa(M, i)) for i in range(1, n + 1))
        # 全部按 DH 变换估计（乐观上界）
        m, a = crba_general_tree(D0, 0, D1, 0)
        m_ub, a_ub = crba(n)
        print(f"  {name:<14} n={n:>2} d={d:>2} D0={D0:>2} D1={D1:>3}  "
              f"CRBA≈{m+a:>6.0f}  (无分支上界 {m_ub+a_ub:>6.0f})")
    print("\n  D0、D1 的界（式 10.7-10.9）：0 ≤ D0 ≤ n−1，0 ≤ D1 ≤ min(n(n−1)/2, n(d−1))")
    print("  下界在「每个刚体都直接连到基座」时取到——此时 H 是对角且常量，算它的代价为零")
