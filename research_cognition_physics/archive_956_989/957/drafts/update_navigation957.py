"""Publish a declared working-hypothesis revision; preserve frozen evidence."""
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
    heading="## 957：任务域预算与共同假说v0.2"
    old="## 956：原生电子材料的交换、记录与来源"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[957报告]({pre}research_note_957.md)明确区分全空间范数与实际任务预算；"
        +"955原自主运行的全部已列H项有有限矩预算，956完整六维材料保原强界。"
        +f"[现行草案v0.2]({pre}957/drafts/unified_operation_hypotheses_v0_2.md)采用明示A1_D工作修订，"
        +"原A1_∞及依赖原前提的定理保留；不自动领取传播、空间、无界来源或读后继承。"
        +f"[结果]({pre}957/task_domain_budget_results.json) · [核验]({pre}957/research_round_957_checks.json)。"
        +"正式957／累计3742，整体目标未完成。"
        +f"接[958]({pre}958/drafts/STATUS.md)核共同材料的控制／读取与来源合同，停止预算及器件优化。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—956轮共956份","001—957轮共957份")
    s=s.replace("956原生材料接口已核，接[957]("+pre+"957/drafts/STATUS.md)落实共同对象与权限合同",
        "957任务域合同已核，接[958]("+pre+"958/drafts/STATUS.md)核同一材料的控制与来源")
    s=s.replace("正式956／累计3741，整体目标未完成","正式957／累计3742，整体目标未完成")
    # The living execution decision points to the new draft; old reports remain intact.
    s=s.replace("[草案v0.1与构造规格]("+pre+"935/drafts/unified_operation_hypotheses_v0_1.md)已保存；它尚非已验证模型。",
        "[现行草案v0.2与构造规格]("+pre+"957/drafts/unified_operation_hypotheses_v0_2.md)已保存；"
        +"A1任务域修订明示，旧版保留，尚非已验证完整模型。")
    s=s.replace("全部共同菜单、A1、应用目标及定时设置保持。",
        "全部共同菜单、应用目标及定时设置保持；A1按957明示任务域版本试行，旧强版保留。")
    if p==paths[4]:
        s=s.replace("当前正式956／累计3741，956已结项","当前正式957／累计3742，957已结项")
    if p==paths[5]:
        s=s.replace("# 231—956轮阶段成果总览","# 231—957轮阶段成果总览",1)
        s=s.replace("231—956的726份","231—957的727份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至956）",
                    "六条共同协议：全局缺口对应与检验优先级（截至957）",1)
        rows={
            "C04":"|C04 同一内部动力学与控制|1、6、H2；929、935、947、955—957|957明示A1_D预算修订，955完整自主运行已列项有限；956六维材料保旧强界|准备／终端及内部控制仍须各给矩继承和来源，不能自动领取全部权限|",
            "C05":"|C05 因果序、局域性、传播界|1、4；465—466、933、957|旧模型传播界保持原条件；957任务预算控制指定量子变化|A1_D不是一般空间局域性／光锥证明，不能代替旧支撑及传播前提|",
            "C19":"|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957|955新K真空及实际未知准备给957统一运行态族预算|读后若继续须保真实矩；未证明任意场仪器、复位和真空制造|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—957|957按同一任务域控制有界任务的动力学差；强预算仅为更强充分路径|源导数、无界应力和不同候选匹配须独立控制；不假定精确无截断父对象|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 956[^\n]*$",lambda m:
            "4. 957保存共同假说v0.2并明示A1_D修订；原生材料优先接合，B_eff保留。"
            +"接[958](../../958/drafts/STATUS.md)核控制／读取和来源的真实权限，不续修预算、二聚体或全部仪器。",s,count=1,flags=re.M)
        assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 957.")
