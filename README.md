# Rigid Body Dynamics Algorithms — 阅读笔记

> 对 Roy Featherstone《**Rigid Body Dynamics Algorithms**》(Springer, 2008) 的
> 逐章梳理笔记。**中文为主，专业术语保留英文原词**，便于与原书对照查阅。

---

## 一、原书信息

| 项目 | 内容 |
|---|---|
| 书名 | Rigid Body Dynamics Algorithms |
| 作者 | Roy Featherstone |
| 出版 | Springer, 2008（ISBN 978-0-387-74314-1） |
| 篇幅 | 正文 10 章 + 附录 |
| 配套代码 | `spatial_v2` MATLAB 工具箱（作者主页提供） |
| 前作 | *Robot Dynamics Algorithms* (Kluwer, 1987) |

这本书是刚体动力学算法领域的标准参考，核心贡献是用**空间向量代数
(spatial vector algebra)** 把 3D 的力/力矩、角速度/线速度统一成 6D 向量，
从而让 RNEA / CRBA / ABA 这些递推算法的推导和实现都大幅简化。

---

## 二、仓库结构

```
.
├── README.md                     # 本文件：导航 + 阅读进度
├── CONVENTIONS.md                # 笔记撰写与更新约定
├── docs/                         # 逐章笔记（主体）
│   ├── 00-roadmap.md             # 全书阅读路线图 + 章节依赖关系
│   ├── ch01-introduction.md      # 第 1 章 绪论
│   ├── ch02-spatial-vector-algebra.md   # 第 2 章 空间向量代数 ★核心
│   ├── ch03-rigid-body-system-dynamics.md # 第 3 章 刚体系统动力学
│   ├── ch04-modelling.md         # 第 4 章 刚体系统建模
│   ├── ch05-inverse-dynamics-rnea.md    # 第 5 章 逆动力学 RNEA ★核心
│   ├── ch06-forward-dynamics-crba.md    # 第 6 章 正动力学·惯性矩阵法 CRBA ★核心
│   ├── ch07-forward-dynamics-aba.md     # 第 7 章 正动力学·传播法 ABA ★核心
│   ├── ch08-closed-loop-systems.md      # 第 8 章 闭环系统
│   ├── ch09-hybrid-dynamics.md          # 第 9 章 混合动力学与其他专题
│   └── ch10-accuracy-efficiency.md      # 第 10 章 精度与效率
├── reference/                    # 速查资料（跨章节）
│   ├── notation.md               # 符号表
│   ├── glossary.md               # 术语中英对照
│   ├── formula-cheatsheet.md     # 空间向量代数公式速查
│   └── algorithm-cards.md        # RNEA / CRBA / ABA 伪代码卡片
├── notes/                        # 个人笔记
│   ├── questions.md              # 疑问与待办清单
│   └── derivations.md            # 自己补的推导
├── code/                         # 可运行的验证脚本
│   └── verify_crba_2link.py      # 2R 机械臂：CRBA / H_ij 公式的三种数值对拍
└── templates/
    └── chapter-note.md           # 章节笔记模板
```

---

## 三、阅读进度

图例：`[ ]` 未读 · `[~]` 在读 · `[x]` 已读并整理 · `[✓]` 已读 + 已做习题/复现代码

| 章 | 标题 | 状态 | 笔记 | 备注 |
|---|---|---|---|---|
| 1 | Introduction | `[ ]` | [ch01](docs/ch01-introduction.md) | 快速过，建立全局观 |
| 2 | Spatial Vector Algebra | `[ ]` | [ch02](docs/ch02-spatial-vector-algebra.md) | ★ 全书地基，必须吃透 |
| 3 | Dynamics of Rigid Body Systems | `[ ]` | [ch03](docs/ch03-rigid-body-system-dynamics.md) | 约束与运动方程的抽象框架 |
| 4 | Modelling Rigid Body Systems | `[ ]` | [ch04](docs/ch04-modelling.md) | 数据结构与关节模型 |
| 5 | Inverse Dynamics | `[ ]` | [ch05](docs/ch05-inverse-dynamics-rnea.md) | ★ RNEA |
| 6 | Forward Dynamics — Inertia Matrix Methods | `[ ]` | [ch06](docs/ch06-forward-dynamics-crba.md) | ★ CRBA |
| 7 | Forward Dynamics — Propagation Methods | `[ ]` | [ch07](docs/ch07-forward-dynamics-aba.md) | ★ ABA |
| 8 | Closed-Loop Systems | `[ ]` | [ch08](docs/ch08-closed-loop-systems.md) | 闭链/并联机构 |
| 9 | Hybrid Dynamics and Other Topics | `[ ]` | [ch09](docs/ch09-hybrid-dynamics.md) | 浮动基、碰撞、齿轮 |
| 10 | Accuracy and Efficiency | `[ ]` | [ch10](docs/ch10-accuracy-efficiency.md) | 代价对比与数值问题 |

**当前进度：0 / 10**

---

## 四、怎么用这个仓库

1. 读之前先看 [`docs/00-roadmap.md`](docs/00-roadmap.md)，确认这一章在全书中的位置和前置依赖。
2. 读的时候把 [`reference/notation.md`](reference/notation.md) 和
   [`reference/formula-cheatsheet.md`](reference/formula-cheatsheet.md) 放在手边——
   本书符号密度很高，符号表能省掉大量翻书时间。
3. 读完一章，回到对应的 `docs/chXX-*.md` 补全「我的理解 / 疑问 / 与原文的出入」三节，
   并在上面的进度表里更新状态。
4. 卡住的点丢进 [`notes/questions.md`](notes/questions.md)，不要停下来死磕，往后读常常自解。

---

## 五、关于本仓库内容的说明

> ⚠️ **重要**：当前各章笔记是在开始精读之前搭建的**知识框架**，内容基于对本书
> 的既有了解整理而成，**尚未逐句与纸质/PDF 原书核对**。因此：
>
> - **章级标题结构**可信度高；
> - **小节编号（如 §2.13）** 请以你手上的版本为准，文中标了 `〔待核对〕` 的地方尤其需要确认；
> - **公式与伪代码**按标准形式书写，与原书可能存在符号或排列上的细微差异。
>
> 阅读过程中发现不一致，直接在对应章节文件里改，并在 commit message 里注明
> `fix(chXX): ...`。这正是本仓库设计成"持续更新"的原因。

---

## 六、许可

个人学习笔记。原书版权归 Springer 与作者所有，本仓库不包含原书任何原文扫描或大段引用。
