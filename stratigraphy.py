#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
层序仲裁工具（stratigraphy.py）—— 考古地层叠压关系与发掘顺序校验

仅依赖 Python 标准库，单文件，开箱即用。

输入格式（文本文件或标准输入，每行一条指令，# 之后为注释，空行忽略）：

    地层 <层名> <年代> <探方>        定义地层；年代为距今年数(BP)，数值越大越古老
    叠压 <探方> <上覆层> <下伏层>    声明某探方内观察到的叠压关系
    发掘 <探方> <层名>              一次发掘操作，按出现顺序执行，状态跨操作延续

指令亦可用英文别名：layer / over / dig。

用法：
    python3 stratigraphy.py 输入文件
    python3 stratigraphy.py            # 从标准输入读取
    python3 stratigraphy.py --demo     # 运行内置示例

校验规则：
  1. 叠压关系校验：上覆层年代不得早于下伏层（上覆层年代值须 <= 下伏层年代值）。
  2. 发掘顺序校验：同探方内须先上后下；跳过未发掘的上覆层直接挖下层时，
     报告探方与被跳过的缺失层（含传递意义上的所有上覆层）。
  3. 跨探方叠压矛盾：两个探方对同一对地层声明了相反方向的叠压关系时，
     报告该地层对与两个探方。
  4. 叠压成环：报告每个环（强连通分量）上的全部地层。
  5. 发掘引用不存在的探方或地层时报告。
  6. 重复发掘已发掘地层时报告。
  7. 发掘历史全程可追溯，输出地层状态、发掘历史与错误清单。
"""

import sys
from collections import defaultdict


# ---------------------------------------------------------------- 数据模型

class Layer:
    """一个地层，身份为 (探方, 层名)。"""
    __slots__ = ("name", "age", "trench", "excavated", "seq")

    def __init__(self, name, age, trench):
        self.name = name
        self.age = age          # 距今年数 BP，越大越古老
        self.trench = trench
        self.excavated = False  # 是否已发掘
        self.seq = None         # 发掘次序（1 起）

    def label(self):
        return "%s:%s" % (self.trench, self.name)


class Site:
    """整个遗址的状态：地层表、叠压关系、发掘历史、错误清单。"""

    def __init__(self):
        self.layers = {}                 # (trench, name) -> Layer
        self.trenches = set()            # 已知探方
        # 同探方叠压邻接表：over_key -> [under_key, ...]（上 -> 下）
        self.adj = defaultdict(list)
        # 按无序层名对记录各探方的方向声明：frozenset({a,b}) -> [(trench, over, under)]
        self.decls_by_pair = defaultdict(list)
        self.errors = []                 # (类别, 行号或None, 消息)
        self.history = []                # (序号, 探方, 层名, 结果描述)
        self._reported_conflicts = set() # 已报告的矛盾，避免重复
        self._dig_seq = 0

    # ------------------------------------------------------------ 错误记录
    def err(self, category, lineno, message):
        self.errors.append((category, lineno, message))

    # ------------------------------------------------------------ 指令：地层
    def define_layer(self, name, age, trench, lineno):
        key = (trench, name)
        if key in self.layers:
            self.err("定义错误", lineno,
                     "地层重复定义：探方 %s 的层 %s" % (trench, name))
            return
        self.layers[key] = Layer(name, age, trench)
        self.trenches.add(trench)

    # ------------------------------------------------------------ 指令：叠压
    def declare_over(self, trench, over, under, lineno):
        if trench not in self.trenches:
            self.err("叠压错误", lineno,
                     "叠压声明引用了不存在的探方：%s" % trench)
            return
        over_key, under_key = (trench, over), (trench, under)
        missing = [n for k, n in ((over_key, over), (under_key, under))
                   if k not in self.layers]
        if missing:
            self.err("叠压错误", lineno,
                     "叠压声明引用了不存在的地层：探方 %s 的 %s"
                     % (trench, "、".join(missing)))
            return
        if over == under:
            self.err("叠压成环", lineno,
                     "地层 %s 自我叠压，构成环" % self.layers[over_key].label())

        # 规则 1：年代校验 —— 上覆层年代不得早于下伏层
        o_age, u_age = self.layers[over_key].age, self.layers[under_key].age
        if o_age > u_age:
            self.err("叠压校验", lineno,
                     "年代倒置：%s（年代 %g）叠压于 %s（年代 %g）之上，"
                     "上覆层年代早于下伏层"
                     % (self.layers[over_key].label(), o_age,
                        self.layers[under_key].label(), u_age))

        # 规则 3：跨探方叠压矛盾 —— 同一层名对的方向冲突
        pair = frozenset((over, under))
        for t2, o2, u2 in self.decls_by_pair[pair]:
            same_dir = (o2 == over and u2 == under)
            if t2 == trench and same_dir:
                self.err("叠压错误", lineno,
                         "重复声明：探方 %s 已声明 %s 叠压 %s"
                         % (trench, over, under))
            elif not same_dir:
                sig = (pair, frozenset((trench, t2)))
                if sig not in self._reported_conflicts:
                    self._reported_conflicts.add(sig)
                    kind = "同探方自相矛盾" if t2 == trench else "跨探方叠压矛盾"
                    self.err(kind, lineno,
                             "%s：地层对 {%s, %s}，探方 %s 声明 %s 叠压 %s，"
                             "探方 %s 声明 %s 叠压 %s"
                             % (kind, over, under,
                                t2, o2, u2, trench, over, under))
        self.decls_by_pair[pair].append((trench, over, under))
        self.adj[over_key].append(under_key)

    # ------------------------------------------------------------ 指令：发掘
    def excavate(self, trench, name, lineno):
        self._dig_seq += 1
        seq = self._dig_seq
        key = (trench, name)

        # 规则 5：引用不存在的探方或地层
        if trench not in self.trenches:
            self.err("发掘错误", lineno, "发掘引用了不存在的探方：%s" % trench)
            self.history.append((seq, trench, name, "失败：探方不存在"))
            return
        if key not in self.layers:
            self.err("发掘错误", lineno,
                     "发掘引用了不存在的地层：探方 %s 的 %s" % (trench, name))
            self.history.append((seq, trench, name, "失败：地层不存在"))
            return

        layer = self.layers[key]

        # 规则 6：重复发掘
        if layer.excavated:
            self.err("发掘错误", lineno,
                     "重复发掘：%s 已在第 %d 次操作中发掘"
                     % (layer.label(), layer.seq))
            self.history.append((seq, trench, name,
                                 "失败：重复发掘（首次于第 %d 次）" % layer.seq))
            return

        # 规则 2：顺序校验 —— 同探方内所有（传递）上覆层须已发掘
        missing = [self.layers[k] for k in self._overlying_keys(key)
                   if k != key and not self.layers[k].excavated]
        if missing:
            names = "、".join(sorted(l.name for l in missing))
            self.err("顺序错误", lineno,
                     "探方 %s：发掘 %s 时跳过未发掘的上覆层：%s"
                     % (trench, name, names))
            self.history.append((seq, trench, name,
                                 "完成（违规：跳过 %s）" % names))
        else:
            self.history.append((seq, trench, name, "完成"))

        # 状态跨操作延续：即使违规也记录为已发掘
        layer.excavated = True
        layer.seq = seq

    def _overlying_keys(self, key):
        """同探方内所有位于 key 之上的地层（沿 上->下 边反向可达，传递闭包）。"""
        # 反向邻接：under -> [over]
        rev = defaultdict(list)
        trench = key[0]
        for o, unders in self.adj.items():
            if o[0] != trench:
                continue
            for u in unders:
                if u[0] == trench:
                    rev[u].append(o)
        seen, stack = set(), [key]
        while stack:
            cur = stack.pop()
            for up in rev.get(cur, ()):
                if up not in seen:
                    seen.add(up)
                    stack.append(up)
        return seen

    # ------------------------------------------------------------ 规则 4：成环检测
    def check_cycles(self):
        """Tarjan 强连通分量；大小>1 的分量或自环即叠压环。"""
        index_of, lowlink, on_stack = {}, {}, set()
        stack, scc_stack, counter = [], [], [0]
        cycles = []

        for root in list(self.adj.keys()) + list(self.layers.keys()):
            if root in index_of:
                continue
            work = [(root, iter(self.adj.get(root, ())))]
            index_of[root] = lowlink[root] = counter[0]
            counter[0] += 1
            stack.append(root)
            on_stack.add(root)
            while work:
                node, it = work[-1]
                advanced = False
                for nxt in it:
                    if nxt not in self.layers:
                        continue
                    if nxt not in index_of:
                        index_of[nxt] = lowlink[nxt] = counter[0]
                        counter[0] += 1
                        stack.append(nxt)
                        on_stack.add(nxt)
                        work.append((nxt, iter(self.adj.get(nxt, ()))))
                        advanced = True
                        break
                    elif nxt in on_stack:
                        lowlink[node] = min(lowlink[node], index_of[nxt])
                if advanced:
                    continue
                work.pop()
                if work:
                    parent = work[-1][0]
                    lowlink[parent] = min(lowlink[parent], lowlink[node])
                if lowlink[node] == index_of[node]:
                    scc = []
                    while True:
                        w = stack.pop()
                        on_stack.discard(w)
                        scc.append(w)
                        if w == node:
                            break
                    if len(scc) > 1:
                        cycles.append(scc)

        for scc in cycles:
            labels = " -> ".join(sorted(self.layers[k].label() for k in scc))
            self.err("叠压成环", None, "环上地层：%s" % labels)


# ---------------------------------------------------------------- 解析

CMD_ALIASES = {
    "地层": "layer", "layer": "layer",
    "叠压": "over", "over": "over",
    "发掘": "dig", "dig": "dig",
}


def parse_and_run(text, site):
    dig_lines = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        cmd = CMD_ALIASES.get(parts[0])
        if cmd is None:
            site.err("解析错误", lineno, "无法识别的指令：%s" % parts[0])
            continue
        if cmd == "layer":
            if len(parts) != 4:
                site.err("解析错误", lineno,
                         "地层指令格式应为：地层 <层名> <年代> <探方>")
                continue
            _, name, age_s, trench = parts
            try:
                age = float(age_s)
            except ValueError:
                site.err("解析错误", lineno, "年代不是数值：%s" % age_s)
                continue
            site.define_layer(name, age, trench, lineno)
        elif cmd == "over":
            if len(parts) != 4:
                site.err("解析错误", lineno,
                         "叠压指令格式应为：叠压 <探方> <上覆层> <下伏层>")
                continue
            _, trench, over, under = parts
            site.declare_over(trench, over, under, lineno)
        else:  # dig —— 先收集，待全部定义与声明解析完、环检测后再按序执行
            if len(parts) != 3:
                site.err("解析错误", lineno,
                         "发掘指令格式应为：发掘 <探方> <层名>")
                continue
            dig_lines.append((lineno, parts[1], parts[2]))

    site.check_cycles()
    for lineno, trench, name in dig_lines:
        site.excavate(trench, name, lineno)


# ---------------------------------------------------------------- 输出

def render(site):
    out = []
    out.append("=" * 56)
    out.append("地层状态")
    out.append("=" * 56)
    if not site.layers:
        out.append("（无地层定义）")
    by_trench = defaultdict(list)
    for layer in site.layers.values():
        by_trench[layer.trench].append(layer)
    for trench in sorted(by_trench):
        out.append("探方 %s：" % trench)
        # 年代小（年轻）在上，自上而下排列
        for layer in sorted(by_trench[trench], key=lambda l: l.age):
            if layer.excavated:
                state = "已发掘（第 %d 次操作）" % layer.seq
            else:
                state = "未发掘"
            out.append("  %-12s 年代 %8g BP   %s" % (layer.name, layer.age, state))

    out.append("")
    out.append("=" * 56)
    out.append("发掘历史")
    out.append("=" * 56)
    if not site.history:
        out.append("（无发掘操作）")
    for seq, trench, name, result in site.history:
        out.append("  第 %-3d 次  发掘 %s:%s  ——  %s" % (seq, trench, name, result))

    out.append("")
    out.append("=" * 56)
    out.append("错误清单（共 %d 条）" % len(site.errors))
    out.append("=" * 56)
    if not site.errors:
        out.append("无错误。")
    for i, (category, lineno, message) in enumerate(site.errors, 1):
        where = "第 %d 行" % lineno if lineno is not None else "全局"
        out.append("  %2d. [%s]（%s）%s" % (i, category, where, message))
    return "\n".join(out)


# ---------------------------------------------------------------- 内置示例

DEMO = """\
# 地层定义：层名 年代(BP) 探方
地层 L1 300  T1
地层 L2 900  T1
地层 L3 1500 T1
地层 L1 320  T2
地层 L2 950  T2
地层 C1 100  T3
地层 C2 200  T3
地层 C3 300  T3
# 叠压声明：探方 上覆层 下伏层
叠压 T1 L1 L2
叠压 T1 L2 L3
叠压 T2 L2 L1        # 与 T1 方向相反 -> 跨探方矛盾
叠压 T3 C1 C2
叠压 T3 C2 C3
叠压 T3 C3 C1        # C1->C2->C3->C1 成环（且 C3 叠压 C1 年代倒置）
叠压 T1 L3 L9        # 引用不存在的地层
# 发掘操作流
发掘 T1 L2           # 跳过未发掘的 L1 -> 顺序错误
发掘 T1 L2           # 重复发掘
发掘 T1 L1
发掘 T9 L1           # 不存在的探方
发掘 T2 L2           # 按 T2 声明先上后下
发掘 T2 L1
"""


# ---------------------------------------------------------------- 入口

def main(argv):
    if "--demo" in argv:
        text = DEMO
    elif len(argv) >= 2 and argv[1] != "-":
        with open(argv[1], "r", encoding="utf-8") as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    site = Site()
    parse_and_run(text, site)
    print(render(site))
    return 1 if site.errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
