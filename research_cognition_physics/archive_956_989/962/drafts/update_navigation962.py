"""Publish a direction-selecting result and acknowledge the user's app edit."""
from pathlib import Path
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
paths = [RESEARCH/n for n in ("README.md", "research_direction.md", "RESEARCH_STATE.md")] + [
    STAGE/n for n in ("README.md", "文件索引.md", "阶段成果总览.md", "跨阶段主题索引.md",
                     "_shared/notes/unified_physics_condition_ledger_current.md")]
ack = json.loads((HERE/"app_goal_confirmation.json").read_text("utf-8"))
goal = (HERE/"revised_goal_20261007.txt").read_text("utf-8")
assert ack["status"] == "active" and ack["objective"].strip() == goal.strip()
original = {p:p.read_bytes() for p in paths}
updates = {}
for p, raw in original.items():
    s = raw.decode("utf-8-sig").replace("\r\n", "\n")
    pre = "archive_764_/" if p.parent == RESEARCH else "../../" if p == paths[-1] else ""
    heading = "## 962：钟保护、关系反馈与稳定性的量词"
    old = "## 961：整体机制图与跨部门筛选优先"
    assert heading not in s and s.count(old) == 1
    block = (heading + "\n\n"
        + f"[962报告]({pre}research_note_962.md)直接复用430／438／449／459，"
        + "以同一原始材料区分裸钟精确隔离与有来源的关系反馈；"
        + "采用M2*澄清：保持声明能力与实际读数，不要求任意未知裸钟绝对隔离。"
        + "精确交织、非零响应解析证书及内部能源流已核，不是几何或引力推导。"
        + f"[结果]({pre}962/clock_relation_screen_results.json) · "
        + f"[核验]({pre}962/research_round_962_checks.json)。"
        + "正式962／累计3747，整体目标未完成。\n\n"
        + f"用户已更新并启用[新目标]({pre}962/drafts/revised_goal_20261007.txt)，"
        + f"[应用读取确认]({pre}962/drafts/app_goal_confirmation.json)。"
        + "停止本轮钟与装置优化，"
        + f"接[963]({pre}963/drafts/STATUS.md)先审共同几何与普通介质／材料变化的判别，"
        + "沿整体机制→冲突→最小检验推进。\n\n")
    s = s.replace(old, block + old, 1)
    s = s.replace("001—961轮共961份", "001—962轮共962份")
    if p == paths[1]:
        before = ("本轮目标接口不能改写应用正文或恢复执行；读取时应用目标为paused，"
                  "故仅仓库目标已更新，应用状态不冒称已同步。当前研究由用户本轮授权继续。")
        after = ("初次读取应用目标为paused；随后用户已写入新目标。"
                 "本轮读取确认应用正文与下文一致、状态active，"
                 "[确认记录](archive_764_/962/drafts/app_goal_confirmation.json)已保存。")
        assert before in s
        s = s.replace(before, after, 1)
    if p == paths[2]:
        before = ("应用目标读取为paused，现有工具不能改正文或恢复，未声称应用已更新。"
                  "当前由本轮用户指令启动962，不新建任务或定时任务。")
        after = ("随后用户已将同一正文写入应用，读取确认active，"
                 "[确认记录](archive_764_/962/drafts/app_goal_confirmation.json)已保存。"
                 "962已按新目标完成，不新建任务或定时任务。")
        assert before in s
        s = s.replace(before, after, 1)
    if p == paths[4]:
        s = s.replace("当前正式961／累计3746，961已结项", "当前正式962／累计3747，962已结项")
    if p == paths[5]:
        s = s.replace("# 231—961轮阶段成果总览", "# 231—962轮阶段成果总览", 1)
        s = s.replace("231—961的731份", "231—962的732份", 1)
    if p == paths[-1]:
        s = s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至961）",
                      "## 六条共同协议：全局缺口对应与检验优先级（截至962）", 1)
        rows = {
            "C09": "|C09 钟尺、参考态与可访问性|1、3、6；522—530、941、962|962同一三qubit材料中的钟与关系联合交织；反馈不能与任意未知裸钟精确隔离并要|采用M2*按任务与访问验收；轴、频率及物理读出仍为输入，不是普适钟或度规|",
            "C21": "|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—962|旧重组与条件处理兼容直接复用；962区分同材料功能分工与裸边缘保护|整体酉不保证局部任务；停止本候选钟优化，下一项按共同几何判别选择|",
            "C22": "|C22 应力、守恒与反作用|1、3、4、6、H2；935—962|962同一钟—关系项同时给相位、可辨反馈和内部能源流|有限模型参数源不是度规应力；共同几何及普适耦合仍开放|",
            "C27": "|C27 共同数据与可区分预测|1、3、5、6；945、958—962|962同一模型的钟可见度与关系分支区别共同核验；修正机制稳定性的过强量词|不是新宇宙预测；963先复用旧锥／钟结果，选择可区分共同几何与介质的有限任务|"}
        for key, row in rows.items():
            s, n = re.subn(r"^\|"+key+r" [^\n]*$", lambda _:row, s, count=1, flags=re.M)
            assert n == 1
    if b"\r\n" in raw:
        s = s.replace("\n", "\r\n")
    updates[p] = (b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"") + s.encode("utf-8")
# This is a newly copied, not previously frozen, old section. Adjust only
# links after moving it four levels below the research root.
history = HERE/"previous_goal_section.md"
oldhistory = history.read_bytes()
history_text = oldhistory.decode("utf-8").replace("](archive_764_/", "](../../")
assert all(p.read_bytes() == raw for p, raw in original.items())
assert history.read_bytes() == oldhistory
for p, data in updates.items():
    p.write_bytes(data)
history.write_text(history_text, encoding="utf-8")
print("Updated eight living navigation files; user app goal active; frozen science unchanged.")
