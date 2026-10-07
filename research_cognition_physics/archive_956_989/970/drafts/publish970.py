"""Publish 970 once. Preserve science; conditional scattering is not unification."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"969/verify_round969.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"local_field_cosmology_results.json")',1)[0]
for a,b in [("range(776,969)","range(776,970)"),("Delivery checks for 969","Delivery checks for 970"),
            ("research_round_969_checks.json","research_round_970_checks.json")]:
    prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"native_photon_scattering_results.json")
    assert result["round"]==970 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["analytic_bounds"]
    assert p["band_min"]>0 and p["lambda_charge"]==.05
    assert math.isclose(p["G"]**2,.05**2*p["d"],rel_tol=1e-14)
    assert max(result[k] for k in ("real_space_matching_error","two_port_unitarity_error",
                                   "record_probability_identity_error"))<1e-12
    assert b["reflection_and_asymptotic_record_lower"]==100/101
    assert b["endpoint_state_error_bound"]<b["reported_endpoint_error"]==.001
    assert b["finite_record_contrast_lower"]>.986
    for row in result["packet_rows"]:
        assert abs(row["normalization"]-1)<1e-12
        assert abs(row["mean_photon_energy"]-p["Delta"])<1e-12
        assert row["record_contrast"]>.999
    r=result["resources"]
    assert r["incoming_parities_same_free_photon_energy_distribution"]
    assert not r["finite_bare_preparations_full_spectral_equality_claimed"]
    assert not r["mechanical_recoil_dynamics_included"]
    assert r["total_energy_current_balance"]<1e-14
    assert result["overlap_with_969"]["cook_duration_over_969_T"]>9000
    for key in ("spatial_dimension_derived","same_full_965_Hamiltonian",
                "Pauli_Fierz_matching_error_proved","same_969_cosmology_window_proved",
                "irreversible_arrow_proved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core970",HERE/"native_photon_scattering.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_970.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==18
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,10)]
    for term in ("整体目标未完成","机械反作用没有完成","不是965完整",
                 "最短时间","自由光子","旋波本身是新增有效模型输入"):
        assert term in prose,term
    newfiles=[note,HERE/"native_photon_scattering.py",HERE/"native_photon_scattering_results.json",
        Path(__file__),HERE/"drafts/scattering_adoption_decision.md",
        HERE/"drafts/scattering_mechanism_increment.md",HERE/"drafts/publish970.py",
        STAGE/"971/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ("截至969","截至970"),("n<=969","n<=970"),
    ("list(range(1,970))","list(range(1,971))"),
    ("001—969轮共969份","001—970轮共970份"),
    ("231—969的739份","231—970的740份"),
    ("969：局部量子交互与宇宙传播的共同来源","970：原生电荷与传播光子的共同记录接口"),
    ('STAGE/"968/research_round_968_checks.json"','STAGE/"969/research_round_969_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=970,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=970,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_969=prev["cumulative_numbered_test_groups_from_968"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        native_charge_vertex_and_old_code=True,
        positive_band_exact_single_excitation_scattering=True,
        finite_wavepacket_Cook_bound=True,
        phase_inputs_same_free_energy_distribution=True,
        common_969_history_verified=False,
        rotating_wave_parent_error_verified=False,mechanical_recoil_verified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"970/verify_round970.py";assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 970：原生电荷与传播光子的共同记录接口"
    previous="## 969：局部量子交互与宇宙传播的共同来源"
    assert heading not in s and s.count(previous)==1
    block=(heading+"\n\n"
        +f"[970报告]({pre}research_note_970.md)以原生电荷矩阵元和旧关系码接成熟单光子散射："
        +"正频旋波模型内，等自由能谱的偶／奇波包留下>0.986的有限时间记录差下界，"
        +"同一S保概率和能源。固定支撑、父作用误差及969共同历史尚未核。"
        +f"[结果]({pre}970/native_photon_scattering_results.json) · "
        +f"[核验]({pre}970/research_round_970_checks.json)。正式970／累计3755，整体目标未完成。\n\n"
        +f"[机制增量]({pre}970/drafts/scattering_mechanism_increment.md)保留传播相位写入材料相干记录的接口，"
        +"不将长的充分时间界当最短时间或不可能定理；机械反冲未完成。"
        +f"停止波导和读取优化，接[971]({pre}971/drafts/STATUS.md)回父有效描述与共同任务域的整体比较。"
        +"应用目标及任务设置保持。\n\n")
    s=s.replace(previous,block+previous,1)
    s=s.replace("001—969轮共969份","001—970轮共970份")
    if p==paths[4]:
        s=s.replace("当前正式969／累计3754，969已结项","当前正式970／累计3755，970已结项")
    if p==paths[5]:
        s=s.replace("# 231—969轮阶段成果总览","# 231—970轮阶段成果总览",1)
        s=s.replace("231—969的739份","231—970的740份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至969）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至970）",1)
        rows={
        "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—970|970外来偶／奇波包写入旧关系码相干性，人口保持且同一效果可读|共同相位及末读输入；不是未知态复制或不可逆记忆|",
        "C05":"|C05 因果序、局域性、传播界|1、4；465—466、933、957、970|970声明正频一维通道，给实际单光子散射及有限波包尾界|通道和旋波是输入；不等于3+1光锥或新的空间局域性证明|",
        "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—970|969原任务运输保持；970以Cook尾界控制有限开始／结束的传播记录|970时窗不能直接输运到969背景；父旋波误差未核，不把充分界当下限|",
        "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—970|970电荷顶点由原CAR谱固定，传播相位可写入旧关系码，不增记录粒子|为956单体接口，不是965完整H；波导／准备／支撑输入|",
        "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—970|969完整平均源保持；970同一光学H保能源并给λ来源|固定支撑尚无机械反冲；光学能流不是全应力或量子引力证书|",
        "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—970|970同自由能谱两相位输入的有限记录差>.986，原生顶点决定散射|非完整总制备成本比较；969共同域未签收，非经验确认|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        pattern=r"^4\. 957假说v0\.2保持；[^\n]*$"
        replacement="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[967集体机制](../../967/drafts/collective_mechanism_increment.md)、[968组织增量](../../968/drafts/internal_organization_increment.md)与[969共同恢复v0.9](../../969/drafts/common_recovery_v0_9.md)保留。[970传播接口](../../970/drafts/scattering_mechanism_increment.md)有条件价值但未进入同一宇宙历史；接[971](../../971/drafts/STATUS.md)回父有效描述和共同域，不扩建波导。"
        s,count=re.subn(pattern,lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 970 preservation verifier.")
