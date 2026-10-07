"""Register the material-write interface without expanding the device branch."""
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
    heading="## 958：材料间的电荷关联写入与同一来源"
    old="## 957：任务域预算与共同假说v0.2"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[958报告]({pre}research_note_958.md)在固定电子密度有效模型中，"
        +"把原关系内容自主写入另一份材料；一阶平均电荷为零而真实电荷关联及来源非零。"
        +"全未知输入误差≤.017060，固定读数窗内容差≥.945881；等待时间约126339，现实寿命及准备／末读仍为输入。"
        +f"[结果]({pre}958/capacitive_material_write_results.json) · [核验]({pre}958/research_round_958_checks.json)。"
        +"正式958／累计3743，完整共同模型未完成。停止电容与器件优化；"
        +f"接[959]({pre}959/drafts/STATUS.md)回到同一父描述及全部物理恢复。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—957轮共957份","001—958轮共958份")
    s=s.replace("957任务域合同已核，接[958]("+pre+"958/drafts/STATUS.md)核同一材料的控制与来源",
        "958自主材料写入已核，接[959]("+pre+"959/drafts/STATUS.md)整合共同父描述与恢复范围")
    if p==paths[4]:
        s=s.replace("当前正式957／累计3742，957已结项","当前正式958／累计3743，958已结项")
    if p==paths[5]:
        s=s.replace("# 231—957轮阶段成果总览","# 231—958轮阶段成果总览",1)
        s=s.replace("231—957的727份","231—958的728份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至957）",
                    "六条共同协议：全局缺口对应与检验优先级（截至958）",1)
        s=s.replace("[存在性草案](../../935/drafts/unified_operation_hypotheses_v0_1.md)",
                    "[现行存在性草案v0.2](../../957/drafts/unified_operation_hypotheses_v0_2.md)",1)
        s=s.replace("全部目标部门、A1和应用目标保持",
                    "全部目标部门和应用目标保持；A1_D沿957，旧强版保留")
        rows={
            "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、923、929、947、956—958|958真实密度作用将关系内容写入第二材料，末读效果不需绝对自旋轴|实际制备／仪器与寿命仍输入；不与955或E_rec证书直接拼接|",
            "C04":"|C04 同一内部动力学与控制|1、6、H2；929、935、947、955—958|958固定H不需脉冲产生条件写入，对全部未知输入及参考有全时间误差界|不是完整控制器或全部六协议；不再展开器件编译|",
            "C17":"|C17 质量、耦合、稳定物质结构|2、3、6；434、941—958|958相同电子材料以普通密度耦合产生非加性条件能源，旧交换关系码直接复用|轨道／势阱／参数／屏蔽为输入；现实寿命和父物理匹配未核|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—958|958用精确本征等距核真实乘积准备的统一误差，保留虚电荷与关联来源|有限材料证书不代替父场匹配；来源与概率分别控制|",
            "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；434—459、929、939、951、956—958|原生材料不只存在单体效果，958同密度作用确能自主传关系内容|有效准备、隔离、旁观自旋和末读列输入；停止电容／读口优化，回到共同父描述|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—958|958同H的三项能源流抵消，电荷关联来源非零，独立平均电荷会漏同阶内容|静态参数源不是全应力／动态几何；离子、控制及父场字典仍需共同登记|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 957[^\n]*$",lambda m:
            "4. 957假说v0.2保持；958核原生材料的自主电荷关联写入与同源能源。"
            +"接[959](../../959/drafts/STATUS.md)整合共同父描述及全部物理恢复，停止电容、读口和寿命优化。",
            s,count=1,flags=re.M);assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 958.")
