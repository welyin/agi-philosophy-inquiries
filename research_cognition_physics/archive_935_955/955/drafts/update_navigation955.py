"""Living navigation for the shared finite protocol, with an adoption gate."""
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
    heading="## 955：物质恢复点与有限协议共同采用"
    old="## 954：接合信号与物质恢复的共同窗口"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[955报告]({pre}research_note_955.md)在同一ζ=.1与原资源下，"
        +"以指定效果误差界接通954物质恢复点与完整有限协议；记录、径向物质及弱引力任务均有正下界。"
        +f"[结果]({pre}955/joint_effective_window_results.json) · [核验]({pre}955/research_round_955_checks.json)。"
        +"正式955／累计3740。有限共存见证成立，完整统一目标未完成。"
        +f"按[范围与投入决定]({pre}955/drafts/common_model_scope_v0_5.md)，停止信号及精度优化；"
        +f"接[956]({pre}956/drafts/STATUS.md)先比较共同模型接口的价值，再决定是否计算，"
        +"不将候选边界自动变成全纲领门槛。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—954轮共954份","001—955轮共955份")
    s=s.replace("954共同系数窗已核，接[955]("+pre+"955/drafts/STATUS.md)审同一有限任务的采用条件",
        "955共同有限任务已核，接[956]("+pre+"956/drafts/STATUS.md)先比较共同模型接口")
    s=s.replace("正式954／累计3739，整体目标未完成","正式955／累计3740，整体目标未完成")
    if p==paths[4]:
        s=s.replace("当前正式954／累计3739，954已结项","当前正式955／累计3740，955已结项")
    if p==paths[5]:
        s=s.replace("# 231—954轮阶段成果总览","# 231—955轮阶段成果总览",1)
        s=s.replace("231—954的724份","231—955的725份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至954）",
                    "六条共同协议：全局缺口对应与检验优先级（截至955）",1)
        rows={
            "C01":"|C01 状态、概率、仪器、复合|1、3、4；001—230、923、947、951、955|955在新K重新运输947完整cq仪器、未知输入与旧参考，界≤.000621169|保有限合同；不制造全部仪器，不认证完整父物理匹配|",
            "C19":"|C19 真空、热态、非平衡准备|1、3、6、H2；既有准备、955|955明示采用新K无源真空，有限协议和物质窗使用同一准备规则|新旧真空不等同；制造装置与真空完整几何来源未认证|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、854、898、944—955|955在同一恢复点保完整协议及三项正信号；指定效果误差不替代全仪器界|有限见证有价值，父模型共同匹配仍开放；先比较接口，停止Gaussian精度修补|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—955|953单极Ward及955同H能量／动量控制各在原范围成立，真实运动进入效果误差|协变单极与Gaussian正过程尚非已认证等价；完整应力不是自动续修清单|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 953[^\n]*$",lambda m:
            "4. 955已接通共同有限任务；停止世界线、门户和信号优化。"
            +"接[956](../../956/drafts/STATUS.md)先比较共同模型接口及真实冲突，只为能改变采用决定的问题安排检验。",s,count=1,flags=re.M)
        assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 955.")
