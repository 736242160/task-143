#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
层序仲裁工具（Stratigraphy Arbiter）

纯 Python 标准库实现。输入地层定义、叠压声明、发掘操作流，输出地层状态、
可追溯的发掘历史与错误清单。

输入文件格式（每行一条指令，空行与 '#' 开头的行被忽略，按空白分词）：

    地层 <层名> <年代> <探方>
        定义一个地层。年代为数值，数值越大表示年代越晚（越年轻）。

    叠压 <上覆层> <下伏层> [探方]
        声明“上覆层压在下伏层之上”。探方可省略，缺省取上覆层所属探方。

    发掘 <探方> <层名>
        记录一次发掘操作。状态在操作之间延续。

校验规则：
    1. 叠压关系：上覆层年代不得早于下伏层（上覆层年代 >= 下伏层年代）。
    2. 发掘顺序：必须先上后下；同探方内仍有未发掘的上覆层（含传递叠压）
       就发掘下伏层，报告探方与缺失层名。
    3. 跨探方矛盾：两个探方对同一对地层声明了相反的叠压方向。
    4. 叠压成环：报告环上全部地层（Tarjan 强连通分量）。
    5. 发掘引用了不存在的探方或地层、地层不属于该探方。
    6. 重复发掘已发掘地层。

用法：
    python3 stratigraphy.py 输入文件
    python3 stratigraphy.py < 输入文件
    python3 stratigraphy.py --help

有错误时进程退出码为 1，否则为 0。
"""

import sys
from collections import defaultdict


class Layer:
    __slots__ = ("name", "date", "trench", "excavated_at")

    def __init__(self, name, date, trench):
        self.name = name
        self.date = date
        self.trench = trench
        self.excavated_at = None  # 发掘操作序号；None 表示未发掘

    @property
    def excavated(self):
        return self.excavated_at is not None


class Arbiter:
    def __init__(self):
        self.layers = {}                       # 层名 -> Layer
        self.trenches = set()                  # 已登记探方
        self.above = defaultdict(set)          # 下伏层 -> {上覆层}
        self.below = defaultdict(set)          # 上覆层 -> {下伏层}
        self.declarations = []                 # (上覆层, 下伏层, 探方)
        self.errors = []                       # (行号, 错误码, 说明)
        self.history = []                      # (序号, 探方, 层名, 结果)
        self.seq = 0

    def error(self, line_no, code, message):
        self.errors.append((line_no, code, message))

    # ---------- 指令处理 ----------

    def define_layer(self, line_no, name, date_s, trench):
        if name in self.layers:
            self.error(line_no, "LAYER_DUP", "地层重复定义：%s" % name)
            return
        try:
            date = float(date_s)
        except ValueError:
            self.error(line_no, "DATE_BAD", "年代不是数值：%s（地层 %s）" % (date_s, name))
            return
        self.layers[name] = Layer(name, date, trench)
        self.trenches.add(trench)

    def declare(self, line_no, upper, lower, trench=None):
        missing = [n for n in (upper, lower) if n not in self.layers]
        if missing:
            self.error(line_no, "LAYER_UNKNOWN",
                       "叠压声明引用了未定义地层：%s" % "、".join(missing))
            return
        if upper == lower:
            self.error(line_no, "STRAT_SELF", "地层不能叠压自身：%s" % upper)
            return
        up, low = self.layers[upper], self.layers[lower]
        if trench is None:
            trench = up.trench
        self.trenches.add(trench)

        # 规则 1：上覆层年代不早于下伏层
        if up.date < low.date:
            self.error(line_no, "STRAT_DATE",
                       "叠压关系与年代矛盾：%s(年代%s) 叠压 %s(年代%s)，上覆层年代更早"
                       % (upper, fmt_date(up.date), lower, fmt_date(low.date)))

        # 规则 3：跨探方相反声明
        for (u2, l2, t2) in self.declarations:
            if u2 == lower and l2 == upper and t2 != trench:
                self.error(line_no, "STRAT_CONFLICT",
                           "跨探方叠压矛盾：地层 %s、%s；探方 %s 声明 %s 叠压 %s，"
                           "探方 %s 声明相反（%s 叠压 %s）"
                           % (upper, lower, t2, lower, upper, trench, upper, lower))

        if lower in self.below[upper]:
            return  # 完全重复的声明，幂等忽略
        self.below[upper].add(lower)
        self.above[lower].add(upper)
        self.declarations.append((upper, lower, trench))

    def excavate(self, line_no, trench, name):
        self.seq += 1
        seq = self.seq

        # 规则 5：探方 / 地层存在性
        if trench not in self.trenches:
            self.error(line_no, "TRENCH_UNKNOWN",
                       "发掘引用了不存在的探方：%s（地层 %s）" % (trench, name))
            self.history.append((seq, trench, name, "失败：探方不存在"))
            return
        if name not in self.layers:
            self.error(line_no, "LAYER_UNKNOWN",
                       "发掘引用了不存在的地层：%s（探方 %s）" % (name, trench))
            self.history.append((seq, trench, name, "失败：地层不存在"))
            return

        layer = self.layers[name]
        if layer.trench != trench:
            self.error(line_no, "TRENCH_MISMATCH",
                       "地层 %s 属于探方 %s，不属于探方 %s"
                       % (name, layer.trench, trench))
            self.history.append((seq, trench, name, "失败：地层不属于该探方"))
            return

        # 规则 6：重复发掘
        if layer.excavated:
            self.error(line_no, "EXC_DUP",
                       "重复发掘：探方 %s 地层 %s 已在操作 #%d 发掘"
                       % (trench, name, layer.excavated_at))
            self.history.append((seq, trench, name, "失败：重复发掘"))
            return

        # 规则 2：先上后下（同探方、传递叠压、仍未发掘）
        missing = sorted(
            n for n in self._ancestors(name)
            if n != name  # 成环时地层会成为自身的“祖先”，需排除
            and self.layers[n].trench == trench
            and not self.layers[n].excavated
        )
        if missing:
            self.error(line_no, "EXC_ORDER",
                       "发掘顺序错误：探方 %s 跳过未发掘的上覆层 [%s] 直接发掘 %s"
                       % (trench, "、".join(missing), name))
            self.history.append(
                (seq, trench, name, "已发掘，但顺序错误（缺失上覆层：%s）" % "、".join(missing)))
        else:
            self.history.append((seq, trench, name, "已发掘"))

        layer.excavated_at = seq

    # ---------- 图分析 ----------

    def _ancestors(self, name):
        """所有（传递地）叠压在 name 之上的地层。"""
        seen, stack = set(), [name]
        while stack:
            cur = stack.pop()
            for up in self.above.get(cur, ()):
                if up not in seen:
                    seen.add(up)
                    stack.append(up)
        return seen

    def find_cycles(self):
        """规则 4：Tarjan 强连通分量，返回每个环上的地层名列表。"""
        sys.setrecursionlimit(max(10000, len(self.layers) * 4 + 1000))
        index_of, low, on_stack = {}, {}, set()
        stack, groups = [], []
        counter = [0]

        def strongconnect(v):
            index_of[v] = low[v] = counter[0]
            counter[0] += 1
            stack.append(v)
            on_stack.add(v)
            for w in self.below.get(v, ()):
                if w not in index_of:
                    strongconnect(w)
                    low[v] = min(low[v], low[w])
                elif w in on_stack:
                    low[v] = min(low[v], index_of[w])
            if low[v] == index_of[v]:
                comp = []
                while True:
                    w = stack.pop()
                    on_stack.discard(w)
                    comp.append(w)
                    if w == v:
                        break
                groups.append(comp)

        nodes = set(self.below) | set(self.above)
        for v in sorted(nodes):
            if v not in index_of:
                strongconnect(v)

        cycles = []
        for comp in groups:
            if len(comp) > 1 or (len(comp) == 1 and comp[0] in self.below.get(comp[0], ())):
                cycles.append(sorted(comp))
        return sorted(cycles, key=lambda c: c[0])

    # ---------- 输出 ----------

    def render(self):
        # 成环属于声明阶段的全局问题，在汇总时检测
        for cyc in self.find_cycles():
            self.errors.append(("-", "STRAT_CYCLE",
                                "叠压关系成环，环上地层：%s" % "、".join(cyc)))

        lines = []
        lines.append("===== 地层状态 =====")
        lines.append("层名    年代    探方    状态")
        name_w = max([len("层名")] + [display_width(n) for n in self.layers])
        trench_w = max([2] + [display_width(t) for t in self.trenches])
        for name in sorted(self.layers):
            ly = self.layers[name]
            status = "已发掘(操作#%d)" % ly.excavated_at if ly.excavated else "未发掘"
            lines.append("%s%s  %s  %s%s  %s"
                         % (name, " " * (name_w - display_width(name)),
                            fmt_date(ly.date),
                            ly.trench, " " * (trench_w - display_width(ly.trench)),
                            status))
        if not self.layers:
            lines.append("（无地层定义）")

        lines.append("")
        lines.append("===== 发掘历史 =====")
        if self.history:
            for seq, trench, name, result in self.history:
                lines.append("#%-3d 发掘 探方 %s 地层 %s -> %s" % (seq, trench, name, result))
        else:
            lines.append("（无发掘操作）")

        lines.append("")
        lines.append("===== 错误清单 =====")
        if not self.errors:
            lines.append("无错误。")
        else:
            lines.append("共 %d 条错误：" % len(self.errors))
            for line_no, code, msg in self.errors:
                where = "行 %s" % line_no if line_no != "-" else "全局"
                lines.append("[%s] %s: %s" % (where, code, msg))

        return "\n".join(lines)


def fmt_date(value):
    if value == int(value):
        return str(int(value))
    return str(value)


def display_width(text):
    """CJK 字符按 2 个显示宽度计，用于简易对齐。"""
    width = 0
    for ch in text:
        width += 2 if ord(ch) > 0x2E7F else 1
    return width


def run(stream):
    arbiter = Arbiter()
    for line_no, raw in enumerate(stream, 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        cmd, args = parts[0], parts[1:]
        if cmd == "地层":
            if len(args) != 3:
                arbiter.error(line_no, "PARSE",
                              "地层定义格式应为：地层 <层名> <年代> <探方>")
            else:
                arbiter.define_layer(line_no, *args)
        elif cmd == "叠压":
            if len(args) not in (2, 3):
                arbiter.error(line_no, "PARSE",
                              "叠压声明格式应为：叠压 <上覆层> <下伏层> [探方]")
            else:
                arbiter.declare(line_no, args[0], args[1],
                                args[2] if len(args) == 3 else None)
        elif cmd == "发掘":
            if len(args) != 2:
                arbiter.error(line_no, "PARSE",
                              "发掘操作格式应为：发掘 <探方> <层名>")
            else:
                arbiter.excavate(line_no, args[0], args[1])
        else:
            arbiter.error(line_no, "UNKNOWN_CMD",
                          "无法识别的指令：%s（可用：地层 / 叠压 / 发掘）" % cmd)
    return arbiter


def main(argv):
    if argv and argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    if argv:
        try:
            with open(argv[0], "r", encoding="utf-8") as f:
                arbiter = run(f)
        except OSError as exc:
            print("无法打开输入文件：%s" % exc, file=sys.stderr)
            return 2
    else:
        arbiter = run(sys.stdin)

    print(arbiter.render())
    return 1 if arbiter.errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
