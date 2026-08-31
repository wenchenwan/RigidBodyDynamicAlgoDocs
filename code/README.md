# 验证代码

配套各章笔记的**可运行验证**。纯 numpy，无其他依赖。

```bash
pip install numpy
python3 code/verify_all.py        # 跨章节一致性套件
python3 code/verify_crba_2link.py # 第 6 章 2R 机械臂的 H 逐元素核对
```

| 文件 | 内容 |
|---|---|
| `spatial.py` | 空间向量代数核心：`plux` `xlt` `Xrotz` `Xstar` `crm` `crf` `rbi` `transform_inertia` `jcalc` |
| `model.py` | 示例机构：2R 平面臂、3R 空间臂、三体分支树、n 连杆串联链；以及 `nu()` `kappa()` |
| `algorithms.py` | `rnea` `crba` `aba`，外加 `H_via_rnea` `H_via_jacobians` `H_via_energy` 三种独立的 $H$ 构造 |
| `verify_all.py` | 跨章节套件：第 2 章的代数恒等式、$H$ 的性质、RNEA/CRBA/ABA 的互验、数值性质 |
| `verify_crba_2link.py` | 第 6 章算例：2R 机械臂 $H$ 与解析解逐元素比对 |

**实现风格**：与书中伪代码逐行对应，**不做任何性能优化**——目的是可读与可验证，
不是快。真要用请看 Pinocchio / RBDL / `spatial_v2`。

## 笔记中所有数值断言的来源

各章笔记里出现的具体数字（如"正确值 0.15363，写反得 0.10755"）都出自这两个脚本。
改了笔记里的公式，**先跑一遍这里**。

## 最有价值的三条对拍

| 对拍 | 检出什么 |
|---|---|
| `rnea(q, q̇, aba(q, q̇, τ)) == τ` | ABA 或 RNEA 的任何一处错误（两算法完全独立） |
| `crba(q) == H_via_energy(q)` | $H$ 的下标条件写反（这条完全不碰 $I^c$ / $X^*$） |
| 能量守恒 | 积分器与整体自洽性 |
