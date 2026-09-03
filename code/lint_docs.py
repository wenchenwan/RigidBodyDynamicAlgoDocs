#!/usr/bin/env python3
"""检查仓库里的 Markdown 在 GitHub 上的渲染问题（规则见 CONVENTIONS.md §2、§4）。

用法：
    python3 code/lint_docs.py            # 检查整个仓库
    python3 code/lint_docs.py docs/ch03-rigid-body-system-dynamics.md

无输出且退出码为 0 即通过。检查项：
  1. 一个 $$ 公式块里出现多个 \\tag（KaTeX/MathJax 会把整块渲染成红字）
  2. 行内 $...$ 或单行 $$...$$ 里用了 \\tag
  3. 空的 / 首尾带空行的 / 未闭合的 $$ 块（用户自动修复脚本留下的残片）
  4. 围栏代码块里的非 ASCII 上下标、组合字符、制表符、箭头、希腊字母、数学符号
     （等宽字体没有这些字形，会把并排两栏推歪）；语言标为 text / mermaid 的围栏视为示意图，豁免
  5. 表格行的行内公式里有裸 |（会被当成单元格分隔符）
  6. 行内 $ 数量为奇数（公式未闭合）
  7. 失效的相对链接
"""
import pathlib
import re
import sys
import unicodedata

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXEMPT_FENCES = {"text", "mermaid"}
CODE_SPAN = re.compile(r"`[^`\n]*`")
LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")
DISPLAY_INLINE = re.compile(r"\$\$(.+?)\$\$")
SINGLE_DOLLAR = re.compile(r"(?<!\$)\$(?!\$)")
INLINE_MATH = re.compile(r"(?<!\$)\$(?!\$)([^$\n]+?)\$(?!\$)")
QUOTE_PREFIX = re.compile(r"^(\s*>\s?)+")


def classify(ch):
    """返回该字符在代码块里的问题类别；ASCII、CJK、全角标点等返回 None。"""
    o = ord(ch)
    if o < 128:
        return None
    cat = unicodedata.category(ch)
    if cat.startswith("M"):
        return "组合字符（如 q̈ 的两点）"
    if 0x2070 <= o <= 0x209F or ch in "¹²³":
        return "上下标"
    if 0x2500 <= o <= 0x257F:
        return "制表符"
    if 0x2190 <= o <= 0x21FF or 0x27F0 <= o <= 0x27FF or 0x2900 <= o <= 0x297F:
        return "箭头"
    if 0x2200 <= o <= 0x22FF:
        return "数学符号"
    if 0x0370 <= o <= 0x03FF:
        return "希腊字母"
    if 0x2100 <= o <= 0x214F or 0x1D400 <= o <= 0x1D7FF:
        return "数学字母"
    if ch in "′″‴":
        return "撇号"
    if cat == "So":
        return "符号/emoji"
    return None


def lint_file(path):
    problems = []
    lines = path.read_text(encoding="utf-8").split("\n")
    in_fence = fence_exempt = False
    fence_start = 0
    in_math = False
    math_start = 0
    math_body = []
    for i, raw in enumerate(lines, 1):
        s = raw.strip()
        # ---- 围栏代码块 ----
        if s.startswith("```"):
            if not in_fence:
                in_fence, fence_start = True, i
                info = s[3:].strip().split(" ")[0] if s[3:].strip() else ""
                fence_exempt = info in EXEMPT_FENCES
            else:
                in_fence = False
            continue
        if in_fence:
            if not fence_exempt:
                seen = set()
                for ch in raw:
                    kind = classify(ch)
                    if kind and ch not in seen:
                        seen.add(ch)
                        problems.append((i, f"代码块内的 {kind} {ch!r} 会破坏等宽对齐；示意图请把围栏标成 ```text"))
            continue
        # ---- $$ 公式块（允许在 > 引用块里）----
        body = QUOTE_PREFIX.sub("", raw)
        if body.strip() == "$$":
            if not in_math:
                in_math, math_start, math_body = True, i, []
            else:
                in_math = False
                text = "\n".join(math_body)
                if not text.strip():
                    problems.append((math_start, "空的 $$ 公式块"))
                else:
                    if not math_body[0].strip():
                        problems.append((math_start, "$$ 块开头有空行"))
                    if not math_body[-1].strip():
                        problems.append((math_start, "$$ 块结尾有空行"))
                    if any(not b.strip() for b in math_body[1:-1]):
                        problems.append((math_start, "$$ 块内部有空行（GitHub 会把块拆成两段）"))
                    ntag = text.count("\\tag")
                    if ntag > 1:
                        problems.append((math_start, f"一个 $$ 块里有 {ntag} 个 \\tag，只能有一个：合并写成 \\tag{{a, b}}"))
            continue
        if in_math:
            math_body.append(body)
            continue
        # ---- 普通文本：先去掉行内代码 ----
        clean = CODE_SPAN.sub("", raw)
        for m in DISPLAY_INLINE.finditer(clean):
            if "\\tag" in m.group(1):
                problems.append((i, "单行 $$...$$ 里用了 \\tag，请改成独立的多行公式块"))
        rest = DISPLAY_INLINE.sub("", clean)
        if "$$" in rest:
            problems.append((i, "行内出现未配对的 $$"))
            continue
        is_table_row = rest.lstrip().startswith("|")
        for m in INLINE_MATH.finditer(rest):
            if "\\tag" in m.group(1):
                problems.append((i, "行内公式 $...$ 里用了 \\tag"))
            if is_table_row and "|" in m.group(1):
                problems.append((i, "表格行的行内公式里有裸 |，请用 \\lvert \\rvert 或 \\vert"))
        if len(SINGLE_DOLLAR.findall(rest)) % 2:
            problems.append((i, "行内 $ 数量为奇数，公式可能未闭合"))
        for m in LINK.finditer(clean):
            target = m.group(1)
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            rel = target.split("#")[0]
            if rel and not (path.parent / rel).exists():
                problems.append((i, f"相对链接失效：{target}"))
    if in_math:
        problems.append((math_start, "$$ 块未闭合"))
    if in_fence:
        problems.append((fence_start, "代码围栏未闭合"))
    return problems


def main(argv):
    if argv:
        files = [pathlib.Path(a) for a in argv]
    else:
        files = sorted(p for p in ROOT.rglob("*.md") if ".git" not in p.parts)
    total = 0
    for f in files:
        for line, msg in lint_file(f):
            print(f"{f.relative_to(ROOT) if f.is_absolute() else f}:{line}: {msg}")
            total += 1
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
