"""Publish 968 once; preserve frozen inputs and register decisive relay checks."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"967/verify_round967.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"collective_material_results.json")',1)[0]
prefix=prefix.replace("for n in range(776,967)","for n in range(776,968)")
prefix=prefix.replace("Delivery checks for 967","Delivery checks for 968")
prefix=prefix.replace("research_round_967_checks.json","research_round_968_checks.json")
checks=r'''    result=read(HERE/"internal_relay_results.json")
    assert result["round"]==968 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"]
    assert (p["physical_total_dimension"],p["exact_invariant_total_dimension"])==(13824,64)
    assert p["eta"]==.05 and p["kappa_AC"]==0.
    assert max(result["exact_restriction_checks"].values())<1e-14
    assert math.isclose(result["relay_initial_participation_acceleration"],3*.05**2/8)
    assert result["endpoint_record_projectors_conserved"]
    dark,active=result["fixed_menu_rows"]
    assert dark["contrast"]<1e-8 and not dark["accepted_nonzero_window"]
    assert active["window_contrast_lower"]>.338 and active["accepted_nonzero_window"]
    assert active["relay_participation_current_norm"]>.02
    assert abs(active["details"][0]["relay_spin_singlet_probability"]-
               active["details"][1]["relay_spin_singlet_probability"])>.068
    for row in (dark,active):
        assert row["certificate"]["total_operator_error_bound"]<1e-4
        assert row["energy_current_balance"]<1e-12
        for d in row["details"]:
            assert abs(sum(d["before_energies"])-sum(d["after_energies"]))<1e-11
    assert all(c["exact_contrast"]==0. for c in result["edge_cut_controls"])
    assert result["scope"]["signal_window_adopted"]
    for key in ("external_time_dependent_coupling_used","new_graph_register_added",
                "geometry_change_derived","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core968",HERE/"internal_relay.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_968.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是低能截断","固定接触布局仍然是输入",
                 "标准舍入模型","固定功能关系","时间箭头未恢复"):
        assert term in prose,term
    newfiles=[note,HERE/"internal_relay.py",HERE/"internal_relay_results.json",
        Path(__file__),HERE/"drafts/internal_relation_decision.md",
        HERE/"drafts/internal_organization_increment.md",HERE/"drafts/publish968.py",
        STAGE/"969/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
tail=tail.replace("截至967","截至968").replace("n<=967","n<=968")
tail=tail.replace("list(range(1,968))","list(range(1,969))")
tail=tail.replace("001—967轮共967份","001—968轮共968份")
tail=tail.replace("231—967的737份","231—968的738份")
tail=tail.replace("967：原生材料的集体作用与共同来源","968：原生材料的内部组织更新与自主传递")
tail=tail.replace('STAGE/"966/research_round_966_checks.json"','STAGE/"967/research_round_967_checks.json"')
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=968,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=968,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_967=prev["cumulative_numbered_test_groups_from_966"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_invariant_material_sector=True,
        autonomous_internal_organization_updates=True,
        fixed_effect_nonzero_window_certified=True,
        numerical_certificate_scope="standard rounding model; not formal library verification",
        no_external_time_dependent_contact=True,
        same_H_energy_and_parameter_sources=True,
        geometry_change_derived=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"968/verify_round968.py";assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 968：原生材料的内部组织更新与自主传递"
    oldheading="## 967：原生材料的集体作用与共同来源"
    assert heading not in s and s.count(oldheading)==1
    block=(heading+"\n\n"
        +f"[968报告]({pre}research_note_968.md)复用原材料的既有自旋：固定内部交换让中继自行参与电荷传递，"
        +"同一远端效果的概率差在有限窗口内>.338；切断任一边解析上无信号。"
        +"组织回作用、能量流与来源在同一H中；不是空间几何变化。"
        +f"[结果]({pre}968/internal_relay_results.json) · "
        +f"[核验]({pre}968/research_round_968_checks.json)。正式968／累计3753，整体目标未完成。\n\n"
        +f"[机制增量]({pre}968/drafts/internal_organization_increment.md)采用内部组织改变功能关系，"
        +"保留固定接触、交换规则、准备与长时间寿命的输入边界。停止中继和读取优化，"
        +f"接[969]({pre}969/drafts/STATUS.md)回整体机制，优先比较与已有传播、钟尺及物质来源的共同桥接。"
        +"应用目标及任务设置保持。\n\n")
    s=s.replace(oldheading,block+oldheading,1)
    s=s.replace("001—967轮共967份","001—968轮共968份")
    if p==paths[4]:
        s=s.replace("当前正式967／累计3752，967已结项","当前正式968／累计3753，968已结项")
    if p==paths[5]:
        s=s.replace("# 231—967轮阶段成果总览","# 231—968轮阶段成果总览",1)
        s=s.replace("231—967的737份","231—968的738份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至967）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至968）",1)
        rows={
        "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—968|968同一材料的端点记录部门保持，中继组织更新，固定远端效果可区分源|不是未知裸态隔离或复制；初态、效果与参考权限输入|",
        "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—968|967低部门全时间界保留；968直接用64维精确不变部门及长时间残差界|两种部门不可混用；未把965场或964背景自动并入|",
        "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—968|968原有自旋的固定交换让同材料中继参与真实远端通信，不增图寄存器|布局、交换和准备输入，真实寿命未认证；停止中继工程|",
        "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—968|968同一H给组织反馈、内部能流及交换／电荷来源；总能量守恒|固定接触未移动；来源不是共同时空或Einstein动力学|",
        "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—968|968固定窗口概率差>.338，切任一边解析无信号；保留967集体谱判别|采用功能关系机制；尚非几何、经验确认或完整统一|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        before="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[966表示筛选](../../966/drafts/mechanism_map_v0_2.md)及[967共同机制增量](../../967/drafts/collective_mechanism_increment.md)与[965物理恢复表](../../965/drafts/common_recovery_v0_8.md)分开使用。967已接通同材料的集体比较与来源；接[968](../../968/drafts/STATUS.md)先审实际内部关系变量，不扩建三体或器件。"
        after="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[966表示筛选](../../966/drafts/mechanism_map_v0_2.md)、[967集体机制](../../967/drafts/collective_mechanism_increment.md)及[968组织增量](../../968/drafts/internal_organization_increment.md)与[965物理恢复表](../../965/drafts/common_recovery_v0_8.md)分开使用。968采用同材料组织改变功能关系；接[969](../../969/drafts/STATUS.md)比较共同传播、钟尺与来源桥接，停止局部中继优化。"
        assert before in s;s=s.replace(before,after,1)
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 968 preservation verifier.")
