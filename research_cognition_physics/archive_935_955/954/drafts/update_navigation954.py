"""Living navigation for the shared portal recovery window."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
    STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
    "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 954：接合信号与物质恢复的共同窗口"
    old="## 953：有效物体的协变作用与共同能源交换"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[954报告]({pre}research_note_954.md)给同一门户核的接合—物质响应精确关系，"
        +"α族在明确容差下有非空共同系数窗；原949整体弱耦合不减小混合。"
        +f"[结果]({pre}954/portal_recovery_window_results.json) · [核验]({pre}954/research_round_954_checks.json)。"
        +"正式954／累计3739。只签收树级空间类响应，未认证完整协议／SM／引力。"
        +f"停止门户扫描，接[955]({pre}955/drafts/STATUS.md)判断系数窗与实际有限协议能否共同采用。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—953轮共953份","001—954轮共954份")
    s=s.replace("953共同作用接口已核，接[954]("+pre+"954/drafts/STATUS.md)审门户物理恢复",
        "954共同系数窗已核，接[955]("+pre+"955/drafts/STATUS.md)审同一有限任务的采用条件")
    s=s.replace("正式953／累计3738，整体目标未完成","正式954／累计3739，整体目标未完成")
    if p==paths[4]:
        s=s.replace("当前正式953／累计3738，953已结项","当前正式954／累计3739，954已结项")
    if p==paths[5]:
        s=s.replace("# 231—953轮阶段成果总览","# 231—954轮阶段成果总览",1)
        s=s.replace("231—953的723份","231—954的724份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至953）",
                    "六条共同协议：全局缺口对应与检验优先级（截至954）",1)
        rows={
            "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—954|954同一源核约束接合与物质自身偏移，α族有非空树级空间类共同窗|不是完整规范散射／共振／实验拟合；群和场表仍输入|",
            "C17":"|C17 质量、耦合、稳定物质结构|2、3、6；941—954|954实际门户极点／残数与物质恢复共同登记；可选α但不能独立调响应|质量机制及参数是输入；未由认知唯一选择，未认证现实全部物质|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、854、898、944—954|954解析共同系数窗同时保非零接合和物质响应容差；953有限尺寸字典边界保持|核系数窗不代替947全交互仪器界，下一项先判同一有限任务能否采用|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 952[^\n]*$",lambda m:
            "4. 953协变材料来源、954共同系数窗已核；停止世界线和门户细化。"
            +"接[955](../../955/drafts/STATUS.md)判断同一参数／准备下的有限任务窗，不能把静态系数当概率或直接领取旧误差。",s,count=1,flags=re.M)
        assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 954.")
