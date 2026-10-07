"""Register the exact constraint interface and full recovery scope."""
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
    heading="## 959：Gauss父表示与共同恢复表"
    old="## 958：材料间的电荷关联写入与同一来源"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[959报告]({pre}research_note_959.md)给958实际材料的精确静电Gauss父表示，"
        +"保未知准备、完整仪器和来源；旧记录界保持。重复计入原静电能会使原概率改变约.68143，独立零通量准备也不等于物理准备。"
        +f"[共同恢复表v0.6]({pre}959/drafts/common_recovery_v0_6.md)区分已核静电部门与尚未实现的完整父物理；"
        +f"[结果]({pre}959/gauss_material_parent_results.json) · [核验]({pre}959/research_round_959_checks.json)。"
        +"正式959／累计3744，整体目标未完成。停止转子／器件扩建，"
        +f"接[960]({pre}960/drafts/STATUS.md)核全部物理部门的共同有效域。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—958轮共958份","001—959轮共959份")
    s=s.replace("958自主材料写入已核，接[959]("+pre+"959/drafts/STATUS.md)整合共同父描述与恢复范围",
        "959静电父接口已核，接[960]("+pre+"960/drafts/STATUS.md)核完整共同父描述的有效域")
    if p==paths[4]:
        s=s.replace("当前正式958／累计3743，958已结项","当前正式959／累计3744，959已结项")
    if p==paths[5]:
        s=s.replace("# 231—958轮阶段成果总览","# 231—959轮阶段成果总览",1)
        s=s.replace("231—958的728份","231—959的729份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至958）",
                    "六条共同协议：全局缺口对应与检验优先级（截至959）",1)
        rows={
            "C01":"|C01 状态、概率、仪器、复合|1、3、4；001—230、947、955、958—959|959固定Gauss等距保958完整未知输入、仪器和旧记录界|只核静电约束父表示；完整父物理和六协议仍未共同实现|",
            "C02":"|C02 子系统、约束、边界拼接|1、3、4、5；360—365、705—706、958—959|959未知电荷与电通量共同准备，局部仪器带本链路；裸部分迹不是原材料字典|本有限中性边界及Gauss部门；不是所有规范区域分解定理|",
            "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—959|959在固定U(1)链路中精确实现实际材料过程的约束与密度来源|群及电场核为输入；两链路不含传播光子，更不等于完整SM|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—959|959静电父表示L1与958过程L2精确保源、保读数，旧有限误差不增|完整L0到L1／L2匹配未证；不要求所有候选先去除UV截断|",
            "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；434—459、929、939、951、956—959|已有自主材料写入可置于同一Gauss物理空间，不需要另配独立空白场|准备、检测、寿命和六协议仍条件；停止转子和器件扩建|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—959|959同一电场能与材料源精确交织；双计静电项改变原有限概率约.68143|参数源不等于全应力，完整父表示需传播场／离子／几何的共同有效字典|",
            "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—959|959在原准备与时刻量化双计静电能的记录差，Gauss准备负对照另列|内部模型对照不是新宇宙预测；完整共同恢复按v0.6逐部门验收|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 957[^\n]*$",lambda m:
            "4. 957假说v0.2保持；959给958的精确静电父表示。"
            +"[共同恢复表v0.6](../../959/drafts/common_recovery_v0_6.md)保全部部门；"
            +"接[960](../../960/drafts/STATUS.md)审完整父物理的同一有效域，停止两转子扩建。",
            s,count=1,flags=re.M);assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 959.")
