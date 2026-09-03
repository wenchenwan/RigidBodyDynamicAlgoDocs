# 验证代码

配套各章笔记的**可运行验证**。纯 numpy，无其他依赖。

```bash
pip install numpy
python3 code/verify_all.py        # 跨章节一致性套件（第 2、3/4、5、6、7、10 章）
python3 code/verify_ch02.py       # 第 2 章：§2.5 分解、§2.13 g_i 与 10 参数、§2.11 加速度、例 2.6、§2.15
python3 code/verify_ch03.py       # 第 3 章：式 3.1/3.2、约束的两种描述与两种施加方式、例 3.1、§3.6 三种解法、式 3.74
python3 code/verify_ch04.py       # 第 4 章：kappa/nu、体雅可比、几何分解、DH、关节极性、球面运动
python3 code/verify_crba_2link.py # 第 6 章：2R 机械臂 H 逐元素核对 + 下标条件
python3 code/cost_model.py        # 第 10 章：表 10.1 的代价公式、交叉点、分支的影响
python3 code/lint_docs.py         # Markdown 渲染规则自检（无输出即通过）
```

| 文件 | 内容 |
|---|---|
| `spatial.py` | 空间向量代数核心：`plux` `xlt` `Xrotx/y/z` `Xstar` `crm` `crf` `rbi` `transform_inertia` `jcalc` `classical_accel` |
| `model.py` | 示例机构（表 4.3 的字段，1-based 父数组）：`arm2r` 2R 平面臂、`arm3r_spatial` 3R 空间臂、`branched3` 三体分支树、`chain(n)` n 连杆串联链；拓扑集合 `nu()` `kappa()` |
| `algorithms.py` | `rnea` `crba` `aba`（与表 5.1、6.2、7.1 逐行对应）；`H_via_rnea` `H_via_jacobians` `H_via_energy` 三种独立的 $H$ 构造；`expanded_parent` `ltl_factor` `ltdl_factor` `sparse_solve_L/LT` `fd_crba_sparse`（表 6.3–6.5） |
| `cost_model.py` | 表 10.1 的运算次数公式；`crossover()` 找 $O(n^3)$ 路线与 ABA 的交叉点；`crba_general_tree()` / `D0_D1()` 算分支树的 $D_0,D_1$ |
| `verify_all.py` | 跨章节套件：第 2 章代数恒等式、$H$ 的性质、重力技巧 vs 逐体外力、CRBA/RNEA/ABA 互验、LTL/LTDL、数值性质、能量守恒 |
| `verify_ch02.py` `verify_ch03.py` `verify_ch04.py` | 按章逐式验证，每条都标注原书式号；`verify_ch02.py` 覆盖 PDF 上三处批注的困惑点 |
| `verify_crba_2link.py` | 第 6 章算例：2R 机械臂 $H$ 与解析解逐元素比对，并演示下标条件写反的后果 |
| `lint_docs.py` | 检查 Markdown 的 GitHub 渲染规则（`CONVENTIONS.md` §2、§4）：一个 `$$` 块只能有一个 `\tag`、行内 `$...$` 不能有 `\tag`、空/未闭合的公式块、代码块内的非 ASCII 上下标与箭头、表格行内公式里的裸 `\|`、失效的相对链接 |

**实现风格**：与书中伪代码逐行对应，**不做任何性能优化**——目的是可读与可验证，
不是快。真要用请看 Pinocchio / RBDL / `spatial_v2`。

## 笔记中所有数值断言的来源

各章笔记里出现的具体数字（如"正确值 0.15363，写反得 0.10755"、"$n=18$ 时比值 1.638"）
都出自这些脚本。改了笔记里的公式，**先跑一遍这里**；改了 Markdown，先跑 `lint_docs.py`。

## 最有价值的三条对拍

| 对拍 | 检出什么 |
|---|---|
| `rnea(q, qd, aba(q, qd, tau)) == tau` | ABA 或 RNEA 的任何一处错误（两算法完全独立） |
| `crba(q) == H_via_energy(q)` | $H$ 的下标条件写反（这条完全不碰 $I^c$ / $X^*$） |
| 能量守恒 | 积分器与整体自洽性 |
