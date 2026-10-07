"""Install the whole-mechanism priority; preserve all frozen science."""
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
    heading="## 961：整体机制图与跨部门筛选优先"
    old="## 960：原生材料的电磁色散接口与有限有效域"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[961整体机制图v0.1]({pre}research_note_961.md)覆盖量子、时空、物质、相互作用、引力、热与时间箭头、宇宙演化，"
        +"区分认知机制解释物理与给定物理实现协议。机制空白、不相容和反直觉分开登记。"
        +"首项筛选排除全体速率重标／只数节点就认定膨胀；相同往返仍可能有不同局部机制。"
        +f"[结果]({pre}961/operational_scale_screen_results.json) · [核验]({pre}961/research_round_961_checks.json)。"
        +"正式961／累计3746，整体目标未完成。\n\n"
        +"**现行执行顺序：整体机制草图 → 概念冲突检查 → 最小数学检验 → 修订假说 → 系统验证。** "
        +"暂缓自动接续局部器件、控制、来源和误差优化；已有结果按原范围保存。"
        +f"[961现行入口]({pre}961/drafts/mechanism_priority_entry.md)接替冻结原入口；"
        +f"接[962]({pre}962/drafts/STATUS.md)优先筛选同一内部关系转换、局部材料保持及来源的联合机制。"
        +"应用目标、任务与定时设置保持。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—960轮共960份","001—961轮共961份")
    s=s.replace("960领先材料／Maxwell接口已核，接[961]("+pre+"961/drafts/STATUS.md)选择完整共同父描述的有效域",
        "961整体机制图已登记，接[962]("+pre+"962/drafts/STATUS.md)按跨部门机制收益选择检验")
    s=s.replace("## 当前阶段执行裁定：先构造一个共同实现（2026-10-07）",
        "## 阶段交付要求：共同实现，按961机制优先顺序推进（2026-10-07）",1)
    if p==paths[4]:
        s=s.replace("当前正式960／累计3745，960已结项","当前正式961／累计3746，961已结项")
    if p==paths[5]:
        s=s.replace("# 231—960轮阶段成果总览","# 231—961轮阶段成果总览",1)
        s=s.replace("231—960的730份","231—961的731份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至960）",
                    "六条共同协议：全局缺口对应与检验优先级（截至961）",1)
        rows={
            "C04":"|C04 同一内部动力学与控制|1、6、H2；929、935、947、955—958、961|已有自主局部过程保持；961要求先说明整套关系更新、材料和资源机制|共同自然规则未定；暂缓器件／控制优化，不把局部成功代替全局说明|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—961|既有误差与字典保留；961先按整体机制选择真正需验证的跨尺度关系|允许受控有效描述；不把任一候选的全部连续／UV问题升级为纲领门槛|",
            "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；434—459、929、939、951、956—961|961把D/R/C/L/Q/B作为同一系统的角色，M1—M4列待筛选机制要求|同一关系更新下的局部能力保持与实际传递仍需联合检验，不免费冻结材料|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—961|已有来源身份按原范围保存；961把关系、材料及交互放在同账机制中|守恒不选共同度规；原961原子／维里来源局部检验暂缓，优先整体机制|",
            "C23":"|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961|961区分过程顺序、钟读数与热箭头；内部遗忘、环境转存和有限空白分别记账|热化、非平衡边界及共同箭头未生成；不由有限可重复直接推出|",
            "C24":"|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、325、331—340、961|961给关系转换→传播／钟尺→尺度的候选链，明确局部材料与频移源要求|完整膨胀机制、分支和暗部门仍空白；参数时间变慢或节点增多均不充分|",
            "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—961|961五快照以往返钟周期、局部谱比和单链路读数区分错误机制叙述|不是自主膨胀模型或新宇宙预测；先整体机制与冲突，再选择最小物理检验|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 957[^\n]*$",lambda m:
            "4. 957假说v0.2保持；[961整体机制图](../../research_note_961.md)为现行解释入口，"
            +"与[959物理恢复表](../../959/drafts/common_recovery_v0_6.md)分开使用。"
            +"按草图→冲突→最小检验→修订→系统验证推进；接[962](../../962/drafts/STATUS.md)，暂缓局部优化。",
            s,count=1,flags=re.M);assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living files: mechanism map is now the current priority.")
