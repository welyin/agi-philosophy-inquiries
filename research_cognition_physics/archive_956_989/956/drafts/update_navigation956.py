"""Publish the native material adoption comparison without merging candidates."""
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
    heading="## 956：原生电子材料的交换、记录与来源"
    old="## 955：物质恢复点与有限协议共同采用"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[956报告]({pre}research_note_956.md)先比较候选，再以同一Hubbard材料"
        +"接通旧434关系码、交换、电荷效果及静态来源；此接口无需新增基本门户。"
        +f"[结果]({pre}956/native_material_interface_results.json) · [核验]({pre}956/research_round_956_checks.json)。"
        +"正式956／累计3741。测量后电荷态须保留；完整阵列、六协议和父SM／GR匹配未认证。"
        +f"按[采用比较]({pre}956/drafts/adoption_comparison.md)，保留B_eff见证，优先接合现成物质；"
        +f"接[957]({pre}957/drafts/STATUS.md)形成共同对象与权限合同，停止二聚体及器件优化。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—955轮共955份","001—956轮共956份")
    s=s.replace("955共同有限任务已核，接[956]("+pre+"956/drafts/STATUS.md)先比较共同模型接口",
        "956原生材料接口已核，接[957]("+pre+"957/drafts/STATUS.md)落实共同对象与权限合同")
    s=s.replace("正式955／累计3740，整体目标未完成","正式956／累计3741，整体目标未完成")
    if p==paths[4]:
        s=s.replace("当前正式955／累计3740，955已结项","当前正式956／累计3741，956已结项")
        s=s.replace("最新取舍见[953复核](953/drafts/priority_reaudit.md)",
                    "最新取舍见[956比较](956/drafts/adoption_comparison.md)")
    if p==paths[5]:
        s=s.replace("# 231—955轮阶段成果总览","# 231—956轮阶段成果总览",1)
        s=s.replace("231—955的725份","231—956的726份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至955）",
                    "六条共同协议：全局缺口对应与检验优先级（截至956）",1)
        rows={
            "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、923、929、947、956|956同一电子材料的电荷效果读取434关系内容，无需绝对自旋方向；保全部电荷后态|不等于尖锐原读口或自主制造；不与955不同候选的仪器证书合并|",
            "C17":"|C17 质量、耦合、稳定物质结构|2、3、6；434、941—956|956原生电子低能交换、双占据和参数响应同源，提供替代材料接合依据|轨道、束缚、参数和有效域为物理输入；并非完整现实材料或质量起源|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—956|955共同有限窗保留；956精确衣着身份保自旋演化与静态来源，测量后态边界已核|两候选不合并；956整网长时间、动态来源及父场误差未认证|",
            "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；434—459、929、939、951、956|956使旧关系交换接到现成电子／电荷材料，此接口不必新添基本味／门户|原生材料优先接合；控制、准备、读取及完整六协议尚需共同合同，停止器件优化|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—956|953协变单极、955同H守恒保原范围；956材料力与静态源含标量偏置／虚跃迁|956参数源不是全应力；不同候选不能无字典拼成量子几何反馈|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 955[^\n]*$",lambda m:
            "4. 956比较后优先接合现成电子材料，B_eff保留为有限见证。"
            +"接[957](../../957/drafts/STATUS.md)落实共同对象与权限；停止二聚体、门户和仪器工程，不把范围边界自动生成待办。",s,count=1,flags=re.M)
        assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 956.")
