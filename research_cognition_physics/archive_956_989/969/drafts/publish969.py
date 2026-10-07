"""Publish 969 once, preserving frozen science and the existing objective."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"968/verify_round968.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"internal_relay_results.json")',1)[0]
prefix=prefix.replace("for n in range(776,968)","for n in range(776,969)")
prefix=prefix.replace("Delivery checks for 968","Delivery checks for 969")
prefix=prefix.replace("research_round_968_checks.json","research_round_969_checks.json")
checks=r'''    result=read(HERE/"local_field_cosmology_results.json")
    assert result["round"]==969 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    prior=read(STAGE/"964/material_cosmology_results.json")
    assert result["shared_background_parameters"]["mu"]==prior["cosmology"]["mu"]
    assert result["shared_background_parameters"]["V0"]==prior["cosmology"]["coordinate_volume"]
    assert .25<result["local_sources"]["mass_increment"]<.25001
    assert result["local_sources"]["minimum_mass_lower"]>99
    assert result["local_sources"]["energy_change_error"]<1e-12
    bg=result["background"]
    assert 2.0006<bg["scale_at_T"]<2.0007
    assert bg["integration_error"]<2e-9 and bg["Raychaudhuri_constraint_max"]<2e-10
    joint=result["joint_task_transport"]
    assert joint["exact_factorization"] and joint["infinite_occupation_primary"]
    assert joint["inherited_receiver_contrast_lower"]>.9987
    assert joint["inherited_polarization_effect_contrast_lower"]>.0498
    assert joint["global_entropy_initial"]==joint["global_entropy_final"]
    bad=result["wrong_dictionary"]
    assert bad["nonzero_t_squared_coefficient"]<0
    assert bad["source_energy_offset_at_T"]<-.12
    for row in bad["initial_field_derivatives"]:
        assert abs(row["first"])<1e-15 and row["second"]>3.2e-7
    assert result["scope"]["one_nonselected_mean_background"]
    assert not any(result["scope"][key] for key in (
        "microscopic_cavity_matching","microscopic_Einstein_derivation",
        "external_radiation_absorption_recovered","irreversible_arrow_derived","full_goal_completed"))
    import importlib.util
    spec=importlib.util.spec_from_file_location("core969",HERE/"local_field_cosmology.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_969.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不重调原引力系数","不是完整QED",
                 "精确分解","有限尺度","只反证上述混用","未完成"):
        assert term in prose,term
    newfiles=[note,HERE/"local_field_cosmology.py",HERE/"local_field_cosmology_results.json",
        Path(__file__),HERE/"drafts/common_background_decision.md",
        HERE/"drafts/common_recovery_v0_9.md",HERE/"drafts/publish969.py",
        STAGE/"970/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for before,after in [
    ("截至968","截至969"),("n<=968","n<=969"),
    ("list(range(1,969))","list(range(1,970))"),
    ("001—968轮共968份","001—969轮共969份"),
    ("231—968的738份","231—969的739份"),
    ("968：原生材料的内部组织更新与自主传递","969：局部量子交互与宇宙传播的共同来源"),
    ('STAGE/"967/research_round_967_checks.json"','STAGE/"968/research_round_968_checks.json"')]:
    tail=tail.replace(before,after)
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=969,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=969,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_968=prev["cumulative_numbered_test_groups_from_967"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        same_965_preparation_interaction_and_effects=True,
        same_964_gravity_coefficients=True,
        complete_local_energy_in_mean_source=True,
        wrong_bound_free_source_dictionary_refuted=True,
        explicit_stable_support_and_monopole_inputs=True,
        real_cavity_matching=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"969/verify_round969.py";assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 969：局部量子交互与宇宙传播的共同来源"
    previous="## 968：原生材料的内部组织更新与自主传递"
    assert heading not in s and s.count(previous)==1
    block=(heading+"\n\n"
        +f"[969报告]({pre}research_note_969.md)在明确稳定支撑／单极平均输入下，"
        +"将965完整局部交互、准备与读口接到964同一引力参数及自由热光子。"
        +"完整局部能源使a(T)从2变为2.000613，原记录和场效果界保留；"
        +"只红移内部光子来源而保留原动力学会违反守恒。"
        +f"[结果]({pre}969/local_field_cosmology_results.json) · "
        +f"[核验]({pre}969/research_round_969_checks.json)。正式969／累计3754，整体目标未完成。\n\n"
        +f"[共同恢复v0.9]({pre}969/drafts/common_recovery_v0_9.md)保留局域／自由模式区分、"
        +"平均背景及微观匹配边界；外来光吸收与不可逆箭头未恢复。"
        +f"停止背景和腔体优化，接[970]({pre}970/drafts/STATUS.md)先审局部材料与外界传播的共同交互。"
        +"应用目标及任务设置保持。\n\n")
    s=s.replace(previous,block+previous,1)
    s=s.replace("001—968轮共968份","001—969轮共969份")
    if p==paths[4]:
        s=s.replace("当前正式968／累计3753，968已结项","当前正式969／累计3754，969已结项")
    if p==paths[5]:
        s=s.replace("# 231—968轮阶段成果总览","# 231—969轮阶段成果总览",1)
        s=s.replace("231—968的738份","231—969的739份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至968）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至969）",1)
        rows={
        "C09":"|C09 钟尺、参考态与可访问性|1、3、6；522—530、941、962—969|969同一解内局部完整谱与自由光频率不同缩放，原读取概率继承|稳定支撑与单极输入；外来发射／吸收比较仪器未实现|",
        "C19":"|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957、964—969|969同一准备含局部相干模式及不同的自由Gibbs模式，温度降低而总熵保持|模式菜单须无重叠；未热化，真空与支撑匹配为输入|",
        "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—969|969精确因子分解运输965原效果和无限占据界，能源另由初态解析守恒运输|本单极模型不自动给真实腔体曲率误差；968交换未并入|",
        "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—969|969完整局部内能为同一平均引力源；错误只红移内部光子能源破坏连续性|给定Einstein均匀作用；不是微观应力、非均匀或量子引力完成|",
        "C23":"|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—969|969实际局部交互与膨胀共存，总熵仍保持；964关联回落证据保留|非平衡容量与稳定记录机制仍需交代；没有宏观不可逆箭头|",
        "C24":"|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、325、331—340、961、964—969|969沿964固定引力参数给965完整交互／读口与自由光子的同一膨胀解|稳定单极／3+1／初态输入；外界吸收、暗部门及生成机制仍缺|",
        "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—969|969原交互能源使尺度终点变2.000613而旧读口界保留，错误来源混用可排除|有效模型内条件结果；不是经验确认或完整统一|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        pattern=r"^4\. 957假说v0\.2保持；[^\n]*$"
        replacement="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[966表示筛选](../../966/drafts/mechanism_map_v0_2.md)、[967集体机制](../../967/drafts/collective_mechanism_increment.md)与[968组织增量](../../968/drafts/internal_organization_increment.md)保留。[969共同恢复v0.9](../../969/drafts/common_recovery_v0_9.md)接通965局部交互和964均匀背景；接[970](../../970/drafts/STATUS.md)审外界传播与局部材料的共同交互，不扩建腔体或背景精度。"
        s,count=re.subn(pattern,lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 969 preservation verifier.")
