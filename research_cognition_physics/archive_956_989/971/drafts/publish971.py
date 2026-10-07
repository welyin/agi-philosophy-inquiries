"""Publish 971 once: parent-interface evidence is not a complete common model."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"970/verify_round970.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"native_photon_scattering_results.json")',1)[0]
for a,b in [("range(776,970)","range(776,971)"),("Delivery checks for 970","Delivery checks for 971"),
            ("research_round_970_checks.json","research_round_971_checks.json")]:
    prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"native_covariant_source_results.json")
    assert result["round"]==971 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["fixed_native_inputs"];r=result["recoil"];c=result["covariance"]
    assert math.isclose(p["electric_dipole_matrix_element"]**2,p["d"],rel_tol=1e-13)
    assert math.isclose(p["electronic_mass_gap"],p["excited_mass"]-p["ground_mass"],abs_tol=1e-13)
    assert math.isclose(r["no_recoil_mass_shell_residual"],p["electronic_mass_gap"]**2,abs_tol=5e-12)
    assert r["recoil_shift_over_benchmark_width"]>120
    assert r["bare_mass_recoil_difference"]<.05*r["benchmark_970_half_width"]
    assert r["benchmark_is_not_a_waveguide_to_vacuum_identification"]
    assert c["maximum_lorentz_vertex_error"]<1e-13
    assert c["maximum_gauge_shift_error"]<1e-13
    assert c["maximum_ward_residual"]<1e-13
    assert c["maximum_mass_shell_residual"]<2e-11
    assert math.isclose(c["rows"][1]["electric_only_relative_error"],.5,abs_tol=1e-12)
    m=result["matching"]
    assert abs(m["actual_CAR_energy_curvature"]-m["resolved_polarizability"])<1e-7
    assert abs(m["wrong_explicit_plus_full_contact_curvature"]-2*m["resolved_polarizability"])<1e-7
    assert not m["self_polarization_of_965_removed"]
    for key in ("actual_transition_rates_certified","real_rotational_selection_rules_matched",
                "full_SM_to_material_matching","full_quantum_positive_parent_constructed",
                "round970_disproved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core971",HERE/"native_covariant_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_971.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","970没有被推翻","不删除965的自极化项",
                 "运动学必要条件","不是已经认证的完整共同模型"):
        assert term in prose,term
    newfiles=[note,HERE/"native_covariant_source.py",HERE/"native_covariant_source_results.json",
        Path(__file__),HERE/"drafts/native_parent_decision.md",
        HERE/"drafts/native_parent_map_v0_1.md",HERE/"drafts/publish971.py",
        STAGE/"972/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ("截至970","截至971"),("n<=970","n<=971"),
    ("list(range(1,971))","list(range(1,972))"),
    ("001—970轮共970份","001—971轮共971份"),
    ("231—970的740份","231—971的741份"),
    ("970：原生电荷与传播光子的共同记录接口","971：原生材料的共同父接口与重复计数检验"),
    ('STAGE/"969/research_round_969_checks.json"','STAGE/"970/research_round_970_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=971,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=971,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_970=prev["cumulative_numbered_test_groups_from_969"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        same_native_mass_gap_and_electric_vertex=True,
        covariant_vertex_and_kinematic_recoil_verified=True,
        finite_matching_double_counting_counterexample=True,
        real_transition_rates_verified=False,
        full_quantum_parent_verified=False,full_SM_matching_verified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"971/verify_round971.py";assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 971：原生材料的共同父接口与重复计数检验"
    previous="## 970：原生电荷与传播光子的共同记录接口"
    assert heading not in s and s.count(previous)==1
    block=(heading+"\n\n"
        +f"[971报告]({pre}research_note_971.md)选择原生材料的共同父有效接口："
        +"同一h和Q约束惯性、运动电偶极及真空反冲；同一CAR响应再次作为完整接触项加入会重复计数。"
        +"这是匹配与运动学身份，不是实际跃迁率或全部SM匹配。"
        +f"[结果]({pre}971/native_covariant_source_results.json) · "
        +f"[核验]({pre}971/research_round_971_checks.json)。正式971／累计3756，整体目标未完成。\n\n"
        +f"[共同父接口表]({pre}971/drafts/native_parent_map_v0_1.md)区分E_rec、P948与原生材料层级，"
        +"不直接叠加各小模型；970波导没有被自由靶反冲检验推翻。"
        +f"停止光学／反冲优化，接[972]({pre}972/drafts/STATUS.md)先回查共同材料的热与有限记录机制。"
        +"应用目标和957假说保持。\n\n")
    s=s.replace(previous,block+previous,1)
    s=s.replace("001—970轮共970份","001—971轮共971份")
    if p==paths[4]:
        s=s.replace("当前正式970／累计3755，970已结项","当前正式971／累计3756，971已结项")
    if p==paths[5]:
        s=s.replace("# 231—970轮阶段成果总览","# 231—971轮阶段成果总览",1)
        s=s.replace("231—970的740份","231—971的741份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至970）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至971）",1)
        rows={
        "C10":"|C10 共同传播几何与普适耦合|1、4、5；324、334、352、924、930、941、963、971|EPS与963条件桥保留；971在共同度规输入下核对同一材料的运动顶点和反冲|共同度规与Lorentz结构是输入；未推出普适耦合或Einstein方程|",
        "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—971|971明示低能替代与残余响应匹配，CAR例证实重复加响应会翻倍|尚无全SM到材料匹配；970时窗未进入969共同历史|",
        "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—971|原生父接口以同一h和Q连接内部态、电荷、质量与传播顶点|实际体取向／支撑及热记录共同域待核，不增万能装置|",
        "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—971|971同一协变偶极给极化电流与壳上来源交换；核对运动顶点和质量壳|形式Noether身份不等于完整量子应力；固定支撑尚非闭系统|",
        "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—971|971错误运动顶点在指定boost偏差50%；自由靶反冲与重复极化提供可复算检验|仅声明有限有效情形；不推翻970波导，不是经验确认|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        pattern=r"^4\. 957假说v0\.2保持；[^\n]*$"
        replacement="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[967集体机制](../../967/drafts/collective_mechanism_increment.md)、[968组织增量](../../968/drafts/internal_organization_increment.md)和[969共同恢复v0.9](../../969/drafts/common_recovery_v0_9.md)保留。[971父接口](../../971/drafts/native_parent_map_v0_1.md)确定材料、运动与来源的共同匹配规则；不直接累加局部模型。停止光学优化，接[972](../../972/drafts/STATUS.md)先回查同一材料的热与记录共同域。"
        s,count=re.subn(pattern,lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 971 preservation verifier.")
