"""Apply the user's priority correction to living navigation only."""
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
PATHS = [RESEARCH / n for n in ("README.md", "research_direction.md", "RESEARCH_STATE.md")]
PATHS += [STAGE / n for n in ("README.md", "文件索引.md", "阶段成果总览.md",
    "跨阶段主题索引.md", "_shared/notes/unified_physics_condition_ledger_current.md")]
TITLE = "## 953工作期执行更正：先判共同接合，暂缓动态细化"

def main():
    original = {p: p.read_bytes() for p in PATHS}
    updates = {}
    for p, raw in original.items():
        s = raw.decode("utf-8-sig").replace("\r\n", "\n")
        assert TITLE not in s, f"Already updated: {p}"
        pre = "archive_764_/" if p.parent == RESEARCH else "../../" if p == PATHS[-1] else ""
        audit = pre + "953/drafts/priority_reaudit.md"
        old = pre + "953/drafts/STATUS.md"
        # Only the prospective 953 clauses change; old scientific summaries remain.
        s, count = re.subn(
            r"接\[953\]\([^)]*953/drafts/STATUS\.md\)[^。；\n]*",
            lambda m: f"按[953取舍复核]({audit})先判共同接合，暂缓原定动态细化",
            s,
        )
        assert count >= 1, p
        s = s.replace("接953核静态准备与有限时间过程",
            f"按[953复核]({audit})先判共同接合，暂缓动态细化")
        if p == PATHS[4]:
            s = s.replace("最新取舍见[952范围决定](952/drafts/shared_source_decision.md)",
                f"最新取舍见[953复核]({audit})")
        block = (
            TITLE + "\n\n"
            + f"[本次取舍记录]({audit})保留948—952成果，先检验有限操作材料与较广物理部门能否共同实现；"
            + "不把候选全部自洽性或无截断完成设为纲领门槛。"
            + f"[原953入口]({old})保留，原定静态／动态细化暂缓。"
            + "本次不新增科学轮次，正式952／累计3737；应用目标与定时设置不变。\n\n"
        )
        first, rest = s.split("\n", 1)
        s = first + "\n\n" + block + rest.lstrip("\n")
        if b"\r\n" in raw:
            s = s.replace("\n", "\r\n")
        updates[p] = (b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"") + s.encode("utf-8")
    assert all(p.read_bytes() == raw for p, raw in original.items()), "Concurrent navigation edit"
    for p, data in updates.items():
        p.write_bytes(data)
    print("Updated eight living navigation files; scientific round remains 952.")

if __name__ == "__main__":
    main()
