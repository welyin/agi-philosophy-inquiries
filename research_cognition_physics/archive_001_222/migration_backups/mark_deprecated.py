# -*- coding: utf-8 -*-
"""给已废止的旧版论文 docx 在文首插入醒目的废止标记。"""
from docx import Document
from docx.shared import Pt, RGBColor

NOTICE_COMMON = (
    "【已废止】本文档为历史旧稿，不代表当前研究结论。其中部分声称（如“已推导出量子概率”"
    "“已完成引力与量子的共同母方程”）已在后续研究中被撤回、证伪或降级为条件性结果。"
    "请以最终稿《8. 认知论视角下的量子信息与引力涌现（假说稿）》及 research_cognition_physics/、"
    "research_information_geometry/、research_physics_construction/ 中的各轮研究笔记为准。"
)

NOTICE_V3_EXTRA = (
    "特别警示（v3 特有的已知错误）：(1) §3.2 性质(3)“折叠操作不增加认知熵”在数学上不成立，"
    "最终稿第 4.2 节已给出显式反例；(2) §4.5.4 的引力符号是以“只有正号才能给出吸引引力”为由事后翻转的，"
    "最终稿已明确废弃该做法；(3) 预言 2 的质量差公式存在量纲错误（多出一个 1/c 因子），最终稿 §6.4 已修正。"
)

TARGETS = {
    "8. 从认知论本体出发：引力量子统一的数学框架.docx": [NOTICE_COMMON],
    "8. 从认知论本体出发：引力量子统一的数学框架（修正版）.docx": [NOTICE_COMMON],
    "8. 从认知论本体出发：引力量子统一的数学框架（v2）.docx": [NOTICE_COMMON],
    "8. 从认知论本体出发：引力量子统一的数学框架（v3）.docx": [NOTICE_V3_EXTRA, NOTICE_COMMON],
}

for path, notices in TARGETS.items():
    doc = Document(path)
    first = doc.paragraphs[0]
    # 逆序插入，使列表首条最终位于最前
    for text in reversed(notices):
        para = first.insert_paragraph_before("")
        run = para.add_run(text)
        run.bold = True
        run.font.size = Pt(12)
        run.font.color.rgb = RGBColor(0xC0, 0x00, 0x00)
    doc.save(path)
    print(f"marked: {path} ({len(notices)} notice(s))")
