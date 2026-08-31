# Rigid Body Dynamics Algorithms — 逐章精读笔记

> 对 Roy Featherstone《**Rigid Body Dynamics Algorithms**》(Springer, 2008) 的
> **逐节精读笔记**，已按仓库内的原书 PDF 全文核对。
> **中文为主，专业术语保留英文原词**，便于与原书对照。

---

## 一、原书信息

| 项目 | 内容 |
|---|---|
| 书名 | Rigid Body Dynamics Algorithms |
| 作者 | Roy Featherstone（The Australian National University） |
| 出版 | Springer, 2008（ISBN 978-0-387-74314-1） |
| 篇幅 | **正文 11 章 + 附录 A**，共 279 页 |
| 前作 | *Robot Dynamics Algorithms* (Kluwer, 1987)，原书明言"已被本书取代" |
| 配套代码 | 作者主页的 `spatial_v2` MATLAB 工具箱（本书不附带源码） |

本仓库含原书 PDF：`[2008 Featherstone] Rigid Body Dynamics Algorithms.pdf`

---

## 二、仓库结构

```
.
├── README.md                     # 本文件：导航 + 阅读进度
├── CONVENTIONS.md                # 笔记撰写与更新约定
├── docs/                         # 逐章笔记（主体）
│   ├── 00-roadmap.md                          # 阅读路线图 + 章节依赖 + 自测题
│   ├── ch01-introduction.md                   # 第 1 章 绪论            pp.1-6
│   ├── ch02-spatial-vector-algebra.md         # 第 2 章 空间向量代数 ★  pp.7-38
│   ├── ch03-rigid-body-system-dynamics.md     # 第 3 章 刚体系统动力学  pp.39-64
│   ├── ch04-modelling.md                      # 第 4 章 刚体系统建模    pp.65-88
│   ├── ch05-inverse-dynamics-rnea.md          # 第 5 章 逆动力学 RNEA ★ pp.89-100
│   ├── ch06-forward-dynamics-crba.md          # 第 6 章 CRBA ★          pp.101-118
│   ├── ch07-forward-dynamics-aba.md           # 第 7 章 ABA ★           pp.119-140
│   ├── ch08-closed-loop-systems.md            # 第 8 章 闭环系统        pp.141-170
│   ├── ch09-hybrid-dynamics.md                # 第 9 章 混合动力学等    pp.171-194
│   ├── ch10-accuracy-efficiency.md            # 第 10 章 精度与效率     pp.195-212
│   ├── ch11-contact-and-impact.md             # 第 11 章 接触与碰撞     pp.213-240
│   └── appendixA-spatial-vector-arithmetic.md # 附录 A 空间向量运算     pp.241-256
├── reference/                    # 速查资料（跨章节）
│   ├── notation.md               # 符号表
│   ├── glossary.md               # 术语中英对照
│   ├── formula-cheatsheet.md     # 空间向量代数公式速查
│   └── algorithm-cards.md        # 全书算法伪代码卡片
├── code/                         # 可运行的 Python 验证实现
│   ├── README.md
│   ├── spatial.py                # 空间向量代数核心
│   ├── model.py                  # 示例机构
│   ├── algorithms.py             # RNEA / CRBA / ABA / LTL / LTDL
│   ├── cost_model.py             # 第 10 章表 10.1 的代价公式与交叉点分析
│   ├── verify_all.py             # 跨章节一致性套件
│   ├── verify_ch02.py            # 第 2 章（含 PDF 中三处批注的困惑点）
│   ├── verify_ch03.py            # 第 3 章
│   ├── verify_ch04.py            # 第 4 章
│   └── verify_crba_2link.py      # 第 6 章 2R 机械臂算例
├── notes/
│   ├── questions.md              # 疑问与待办清单
│   └── derivations.md            # 推导补充
└── templates/
    └── chapter-note.md           # 章节笔记模板
```

---

## 三、阅读进度

图例：`[ ]` 未读 · `[~]` 在读 · `[x]` 已读并整理 · `[✓]` 已读 + 已做习题/复现代码

| 章 | 标题 | 页 | 节数 | 状态 | 笔记 |
|---|---|---|---|---|---|
| 1 | Introduction | 1–6 | 5 | `[ ]` | [ch01](docs/ch01-introduction.md) |
| 2 | **Spatial Vector Algebra** ★ | 7–38 | 17 | `[~]` | [ch02](docs/ch02-spatial-vector-algebra.md) |
| 3 | Dynamics of Rigid Body Systems | 39–64 | 7 | `[ ]` | [ch03](docs/ch03-rigid-body-system-dynamics.md) |
| 4 | Modelling Rigid Body Systems | 65–88 | 6 | `[ ]` | [ch04](docs/ch04-modelling.md) |
| 5 | **Inverse Dynamics** ★ | 89–100 | 5 | `[ ]` | [ch05](docs/ch05-inverse-dynamics-rnea.md) |
| 6 | **FD — Inertia Matrix Methods** ★ | 101–118 | 6 | `[ ]` | [ch06](docs/ch06-forward-dynamics-crba.md) |
| 7 | **FD — Propagation Methods** ★ | 119–140 | 5 | `[ ]` | [ch07](docs/ch07-forward-dynamics-aba.md) |
| 8 | Closed Loop Systems | 141–170 | 13 | `[ ]` | [ch08](docs/ch08-closed-loop-systems.md) |
| 9 | Hybrid Dynamics and Other Topics | 171–194 | 7 | `[ ]` | [ch09](docs/ch09-hybrid-dynamics.md) |
| 10 | Accuracy and Efficiency | 195–212 | 4 | `[ ]` | [ch10](docs/ch10-accuracy-efficiency.md) |
| 11 | Contact and Impact | 213–240 | 9 | `[ ]` | [ch11](docs/ch11-contact-and-impact.md) |
| A | Spatial Vector Arithmetic | 241–256 | 5 | `[ ]` | [appA](docs/appendixA-spatial-vector-arithmetic.md) |

**当前进度：读到第 2 章 §2.13 附近**（据 PDF 中的批注位置推断）

---

## 四、📌 你在 PDF 上标注的困惑点

PDF 中共有 5 处批注，全部集中在第 2 章。三处明确的困惑已写成专门小节：

| PDF 页 | 你的批注 | 对应笔记 | 结论摘要 |
|---|---|---|---|
| p.23 §2.5 | 「这里的理解仍然有困惑」 | [ch02 §2.5-C](docs/ch02-spatial-vector-algebra.md) | 原书两条分解回答**不同的问题**：指定点则唯一；要求平移平行于轴则得 Chasles 螺旋分解，后者是前者取 $P=(s\times s_O)/(s\cdot s)$ 的特例 |
| p.25 §2.8 | 「对运动向量和力向量分别进行坐标变换的关系」 | [ch02 §2.1](docs/ch02-spatial-vector-algebra.md) | $X^{*}=X^{-\mathsf T}$ 是**纯代数结论**（对偶坐标），物理上对应功率守恒 |
| p.32 表 2.3 | 「参考图 2.5a」 | [ch02 §2.9](docs/ch02-spatial-vector-algebra.md) | 叉乘表的每个元素 = 顶行基向量在左列基向量速度下的时间导数 |
| p.40 §2.13 | 「$g_i$ 的实际含义是什么？？？」 | [ch02 §2.13-B](docs/ch02-spatial-vector-algebra.md) | $g_i\in F^6$ 且**不唯一**；一组物理取法是 $\sqrt{m_k}\times$（过质点 $p_k$ 沿 $e_a$ 的单位力）；引入它是为了推出 $\dot I=v\times^{*}I-Iv\times$ |
| p.42 §2.13 | 「why???」（为何 10 个参数） | [ch02 §2.13-C](docs/ch02-spatial-vector-algebra.md) | $10=1+3+6$；也可由 $21-11$ 得到（右下块必为 $m\mathbf 1$ 省 5、右上块必反对称省 6） |

以上结论全部有可运行验证：`python3 code/verify_ch02.py`

---

## 五、可运行的验证代码

```bash
pip install numpy
python3 code/verify_all.py        # 跨章节一致性套件（全部通过）
python3 code/verify_ch02.py       # 第 2 章，含三处困惑点
python3 code/verify_ch03.py       # 第 3 章
python3 code/verify_ch04.py       # 第 4 章
python3 code/verify_crba_2link.py # 第 6 章 2R 算例
python3 code/cost_model.py        # 第 10 章代价与交叉点分析
```

**笔记中所有数值断言都出自这些脚本。** 详见 [`code/README.md`](code/README.md)。

**最有价值的三条对拍**：

| 对拍 | 检出什么 |
|---|---|
| `rnea(q, q̇, aba(q, q̇, τ)) == τ` | ABA 或 RNEA 的任何错误（两算法完全独立） |
| `crba(q) == H_via_energy(q)` | $H$ 下标条件写反（这条完全不碰 $I^c$ / $X^*$） |
| 能量守恒 | 积分器与整体自洽性 |

---

## 六、怎么用这个仓库

1. 读之前看 [`docs/00-roadmap.md`](docs/00-roadmap.md)，确认这一章的位置和前置依赖。
2. 读的时候把 [`reference/notation.md`](reference/notation.md) 和
   [`reference/formula-cheatsheet.md`](reference/formula-cheatsheet.md) 放在手边。
3. 每章笔记开头都有**原书节次对照表**（节号 + 原书标题 + 页码），
   可以直接对照 PDF 定位。
4. 读完一章，补全笔记末尾的「✍️ 我的理解 / ❓ 疑问与待办 / 📌 与原文的出入」三节，
   并更新上面的进度表。
5. 卡住的点丢进 [`notes/questions.md`](notes/questions.md)，不要停下来死磕。

---

## 七、笔记的质量说明

- **所有章节已按仓库内的原书 PDF 全文逐节核对**：节号、标题、页码、
  公式编号、伪代码（表 5.1、6.1、6.2、6.3、7.1、8.2、8.3、9.1、9.2、9.3、11.1 等）
  均取自原书。
- **公式与算法有数值验证**：见第五节。
- **引用原书原文的地方用引用块标出**，其余为笔记作者的整理、解释与补充。
- 标有 💡 的是笔记补充的解释或实践建议，标有 🔑 的是原书中特别关键的论断。

**已修正的重要错误**（早期凭既有知识搭建框架时产生）：

| 错误 | 修正 |
|---|---|
| 以为全书 10 章 | **实为 11 章 + 附录**；第 11 章 Contact and Impact 是独立一章 |
| $H_{ij}=S_i^{\mathsf T}I^c_i{}^iX_jS_j$ 的下标条件写反 | 应为 **$i\in\nu(j)$（$j$ 是祖先，$I^c$ 取后代）**，已由原书式 6.14 确认 |
| 称本书"比 DH 更通用"而未提 DH | 原书 **§4.3 整节讲 DH 参数** |
| 遗漏关节极性 | 补上 §4.1.3 与表 4.2 的极性反转 |
| 齿轮"在对角元加 $\rho^2I_{rotor}$" | 原书把转子建模成**独立刚体**，用 §8.11 的闭环函数技术 |
| 遗漏 §5.1、§5.2、§9.7、§10.2 等整节 | 已补 |
| 第 10 章代价表标〔待填〕 | 已按**表 10.1** 填入，交叉点精确到 $n\le8$ / $n\ge9$ |

---

## 八、许可

个人学习笔记。原书版权归 Springer 与作者所有。
