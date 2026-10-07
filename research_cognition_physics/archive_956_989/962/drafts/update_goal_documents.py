"""Apply the requested research-goal revision to living docs, not app state."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
goal = (HERE / "revised_goal_20261007.txt").read_text("utf-8")
direction = RESEARCH / "research_direction.md"
raw = direction.read_bytes()
text = raw.decode("utf-8-sig").replace("\r\n", "\n")
start = text.index("## 当前目标正文：从认知操作压缩物理的独立假设")
end = text.index("\n## 952：", start)
history = HERE / "previous_goal_section.md"
assert not history.exists()
history.write_text(text[start:end] + "\n", encoding="utf-8")
section = (
    "## 当前目标正文：从认知操作压缩物理的独立假设\n\n"
    "2026-10-07按用户“修改一下目标，然后开始”的授权修订。"
    "以下是本任务现行研究目标，优先于历史轮次的下一步。"
    "完整可粘贴版本见[目标正文](archive_764_/962/drafts/revised_goal_20261007.txt)，"
    "[旧目标段落](archive_764_/962/drafts/previous_goal_section.md)保留。"
    "本轮目标接口不能改写应用正文或恢复执行；读取时应用目标为paused，"
    "故仅仓库目标已更新，应用状态不冒称已同步。当前研究由用户本轮授权继续。\n\n"
    + goal.replace("# 目标：", "### 目标：", 1)
)
text = text[:start] + section + text[end:]
if b"\r\n" in raw:
    text = text.replace("\n", "\r\n")
assert direction.read_bytes() == raw
direction.write_bytes((b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"") + text.encode("utf-8"))
state = RESEARCH / "RESEARCH_STATE.md"
raw = state.read_bytes()
text = raw.decode("utf-8-sig").replace("\r\n", "\n")
heading = "## 目标修订：整体机制与决定性检验（2026-10-07）"
assert heading not in text
block = (
    heading + "\n\n用户已授权修改目标并开始研究。"
    "[现行目标](research_direction.md#当前目标正文从认知操作压缩物理的独立假设)"
    "采用机制草图→冲突检查→最小检验→修订→系统验证；"
    "候选内部技术债仅在影响所采信的有限共同预测时成为验收要求。"
    "维持共同模型、全部物理部门和有效描述的阶段目标。"
    "[可粘贴正文](archive_764_/962/drafts/revised_goal_20261007.txt)已保存；"
    "应用目标读取为paused，现有工具不能改正文或恢复，未声称应用已更新。"
    "当前由本轮用户指令启动962，不新建任务或定时任务。\n\n"
)
text = text.replace("\n\n", "\n\n" + block, 1)
if b"\r\n" in raw:
    text = text.replace("\n", "\r\n")
assert state.read_bytes() == raw
state.write_bytes((b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"") + text.encode("utf-8"))
print("Research goal updated in direction/state; app goal not modified.")
