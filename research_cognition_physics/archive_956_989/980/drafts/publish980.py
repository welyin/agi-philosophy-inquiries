"""Publish the finite native thermal-record interface; preserve frozen evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'979/verify_round979.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"spectral_mass_bridge_results.json")',1)[0]
for a,b in [('range(776,979)','range(776,980)'),('Delivery checks for 979','Delivery checks for 980'),
            ('research_round_979_checks.json','research_round_980_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"finite_thermal_records_results.json")
    assert result["round"]==980 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];c=result["channel"];r=result["finite_resources"]
    assert p["eta"]==.05 and p["theta"]==1 and p["Q_star_center"]==2
    assert c["contacts"]==600000 and c["full_charge_reset_bound"]<.001029
    assert c["comparison_thermal_reset_bound"]==13/14000
    assert c["arbitrary_unknown_input_and_passive_reference"]
    assert c["one_contact_dilation_error"]<1e-12 and c["stationary_error"]<1e-12
    assert c["final_choi_reset_bound"]<c["comparison_thermal_reset_bound"]
    cert=result["eigensystem_certificate"]["coarse_analytic_certificate"]
    assert cert["analytic_spectral_gap_lower"]>1/62500
    assert cert["nonresonant_gap_rational_lower"]>.03
    assert cert["minimum_population_diagonal_lower"]>.5
    assert cert["exp_9_6_rational_lower"]>14000
    assert r["prepared_thermal_copies"]==600000 and r["no_bath_reused_or_reset"]
    assert r["total_contact_time"]>5e19 and r["initial_bath_energy_above_ground"]>1e5
    assert r["thermal_state_to_best_pure_blank_distance"]>.58
    assert r["switching_work_absolute_bound"]<6.3e-9
    assert r["controller_preparation_not_constructed"]
    therm=result["thermodynamic_check"]
    assert therm["energy_balance_error"]<1e-12 and therm["landauer_identity_error"]<1e-12
    assert therm["endpoint_entropy_production"]>0
    for row in result["averaging_diagnostics"]:assert row["full_unitary_error"]<row["averaging_bound"]
    for key in ("infinite_bath_required","new_primitive_reset_operator_inserted",
        "constant_autonomous_controller_constructed","pure_blank_reset_proved",
        "universe_permanent_arrow_proved","actual_geometric_feedback_computed",
        "whole_lifecycle_with_976_verified","full_SM_matching_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core980",HERE/"finite_thermal_records.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_980.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是宇宙永久时间箭头，也不是纯空白复位",
        "小开关平均功不是小控制成本证明","不保证实际每一步严格单调",
        "自主控制器未构造"):
        assert term in prose,term
    newfiles=[note,HERE/"finite_thermal_records.py",HERE/"finite_thermal_records_results.json",
        Path(__file__),HERE/"drafts/finite_thermal_decision.md",
        HERE/"drafts/thermal_arrow_adoption.md",HERE/"drafts/publish980.py",
        STAGE/"981/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至979','截至980'),('n<=979','n<=980'),('list(range(1,980))','list(range(1,981))'),
    ('001—979轮共979份','001—980轮共980份'),('231—979的749份','231—980的750份'),
    ('979：完整协议与实际质量的谱模拟桥接','980：有限热材料与有效遗忘'),
    ('STAGE/"978/research_round_978_checks.json"','STAGE/"979/research_round_979_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=980,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=980,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_979=prev["cumulative_numbered_test_groups_from_978"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        finite_native_thermal_reset_with_reference_verified=True,
        analytic_contraction_and_finite_averaging_verified=True,
        full_bath_and_mean_switching_work_accounts_verified=True,
        autonomous_controller_or_pure_blank_reset_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'980/verify_round980.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 980：有限热材料与有效遗忘'
    previous='## 979：完整协议与实际质量的谱模拟桥接'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[980报告]({pre}research_note_980.md)以972同种材料及原电荷接触，给有限热复位合同：'
        +'60万份预备热副本下，全部未知输入与参考的误差<.001029；整体资料仍在，热量及开关平均功共同入账。'
        +f'[结果]({pre}980/finite_thermal_records_results.json) · '
        +f'[核验]({pre}980/research_round_980_checks.json)。正式980／累计3765，整体目标未完成。\n\n'
        +f'[采用范围]({pre}980/drafts/thermal_arrow_adoption.md)保留大而有限的时间／热资源、接触及控制输入；'
        +'不是纯空白、完整生命周期或宇宙永久箭头。'
        +f'停止热接触器件优化，接[981]({pre}981/drafts/STATUS.md)回共同父模型及全部门有效域。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—979轮共979份','001—980轮共980份')
    if p==paths[4]:s=s.replace('当前正式979／累计3764，979已结项','当前正式980／累计3765，980已结项')
    if p==paths[5]:
        s=s.replace('# 231—979轮阶段成果总览','# 231—980轮阶段成果总览',1)
        s=s.replace('231—979的749份','231—980的750份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至979）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至980）',1)
        rows={
        'C03':'|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—980|980原电荷有限热接触使局部旧记录接近同一热态并与参考解关联；整体资料转存|初始热副本、接触及控制输入；不是全局信息消失或与976共用完整生命周期|',
        'C19':'|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—980|980全部60万份热副本预先入账，不免费复用；热与开关平均功有精确身份|准备及控制资源仍输入；混合热态距纯空白约.582，不能当免费再初始化|',
        'C23':'|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—980|980有限碰撞比较通道有严格收缩，原电荷全过程热复位误差<.001029；整体幺正保持|有限接触顺序提供方向性，实际逐步单调和宇宙永久箭头未证；停止热化器件优化|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—980|980有限未知量子输入／参考的热复位有解析统一误差，原电荷与共振比较经有限平均界相连|充分资源和时间大，不是现实最优效率；没有把热复位签成纯空白或完整父物理|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[980采用增量](../../980/drafts/thermal_arrow_adoption.md)给有限原生热复位与全资源边界，未领取永久箭头或纯空白。停止热接触、模拟器与器件优化，接[981](../../981/drafts/STATUS.md)回共同父模型和全部门有效范围。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 980 verifier.')
