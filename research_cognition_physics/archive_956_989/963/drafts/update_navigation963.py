"""Publish finite scale consistency and the conditional EPS bridge."""
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
RESEARCH = STAGE.parent
paths = [RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")] + [
    STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                     "_shared/notes/unified_physics_condition_ledger_current.md")]
original = {p:p.read_bytes() for p in paths}
updates = {}
for p,raw in original.items():
    s = raw.decode("utf-8-sig").replace("\r\n","\n")
    pre = "archive_764_/" if p.parent == RESEARCH else "../../" if p == paths[-1] else ""
    heading = "## 963：共同钟标定的有限证书与EPS接口"
    old = "## 962：钟保护、关系反馈与稳定性的量词"
    assert heading not in s and s.count(old) == 1
    block = (heading+"\n\n"
        +f"[963报告]({pre}research_note_963.md)把多类钟的共同尺度写成有限差分约束，"
        +"给精确赋值或负环证书；物种差与路径差的两个合成菜单分别有最小log容差1/100、1/75。"
        +"对接EPS时明示光／自由落体／钟输运条件，保留实际物质耦合的区别；"
        +"不同介质声锥、累计时间差和内部相位不自动否定共同几何。"
        +f"[结果]({pre}963/common_scale_certificate_results.json) · "
        +f"[核验]({pre}963/research_round_963_checks.json)。正式963／累计3748，整体目标未完成。\n\n"
        +"停止标定工具优化，EPS只作可选桥梁；"
        +f"接[964]({pre}964/drafts/STATUS.md)回到同一候选的物理恢复与来源，"
        +"不把这份条件证书升级为所有有效模型的门槛。应用目标保持active，任务与定时设置保持。\n\n")
    s = s.replace(old,block+old,1)
    s = s.replace("001—962轮共962份","001—963轮共963份")
    if p == paths[4]:
        s = s.replace("当前正式962／累计3747，962已结项","当前正式963／累计3748，963已结项")
    if p == paths[5]:
        s = s.replace("# 231—962轮阶段成果总览","# 231—963轮阶段成果总览",1)
        s = s.replace("231—962的732份","231—963的733份",1)
    if p == paths[-1]:
        s = s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至962）",
                      "## 六条共同协议：全局缺口对应与检验优先级（截至963）",1)
        rows = {
          "C09":"|C09 钟尺、参考态与可访问性|1、3、6；522—530、941、962—963|963给多钟有限共同尺度的赋值／负环证书；材料差与路径差分别可判|合成合同不是实际钟实现；周期—联络识别额外列明，保持M2*与访问边界|",
          "C10":"|C10 共同传播几何与普适耦合|1、4、5；324、334、352、924、930、941、963|EPS提供光、自由落体与钟的条件桥；963明确共同尺度的独立验收|不同介质声锥不自动违反共同度规；钟定义不可混用，未推出普适耦合或Einstein方程|",
          "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—963|963有限区间数据可给兼容尺度或明确冲突，不要求精确零误差|有限菜单不外推全尺度；实际候选来源与误差仍须同账，不将EPS／UV升为普遍门槛|",
          "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—963|963合成菜单区分物种尺度差和共同路径差，给达到的兼容误差边界|不是宇宙观测或几何唯一性；停止工具扩展，964回同一物理候选的共同实现|"}
        for key,row in rows.items():
            s,n = re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert n == 1
    if b"\r\n" in raw:
        s = s.replace("\n","\r\n")
    updates[p] = (b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes() == raw for p,raw in original.items())
for p,data in updates.items():
    p.write_bytes(data)
print("Updated eight living entries; all prior scientific files preserved.")
