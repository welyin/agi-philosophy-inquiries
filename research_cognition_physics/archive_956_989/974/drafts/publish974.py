"""Publish 974 once; preserve previous evidence and the full active goal."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'973/verify_round973.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"common_source_distribution_results.json")',1)[0]
for a,b in [('range(776,973)','range(776,974)'),('Delivery checks for 973','Delivery checks for 974'),
            ('research_round_973_checks.json','research_round_974_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"joint_energy_source_results.json")
    assert result["round"]==974 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];cert=result["certificate"];reports=result["unconditional_reports"]
    assert p["g"]==.003 and p["omega"]==.5 and p["m0"]==100
    assert p["computational_cutoff"]==10
    assert abs(p["proper_time_difference"]-math.pi)<1e-14
    assert cert["transported_bound"]==2e-6 and cert["old_uniform_isometry_bound"]<2e-7
    assert cert["complete_square_mass_lower"]>99
    assert abs(cert["material_charge_commutator"]-.4)<1e-12
    assert cert["current_sum_residual"]<1e-12
    assert cert["source_derivative_residual"]<1e-9
    assert cert["arbitrary_passive_reference"] and cert["coherent_path_input_preserved"]
    assert reports["receiver_contrast_lower"]>.9947
    assert reports["physical_mode_contrast_lower"]>.0224
    assert result["rows"][0]["infinite_model_path_difference_lower"]>.14709
    for row in result["rows"]:
        assert .70<row["visibility"]<.71
        assert row["energy_variance"]>.062
        assert row["visibility_endpoint_difference"]<2e-10
        for energy in row["energies"]:assert abs(energy["joint_energy_change"])<1e-10
    for key in ("fixed_path_supports_certified","full_moving_backreaction_certified",
                "quantum_metric_dynamics_constructed","post_Newton_remainder_certified",
                "full_SM_to_material_matching","cosmology_973_repaired",
                "macroscopic_time_arrow_proven","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core974",HERE/"joint_energy_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_974.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","固定路径有效模型","不是本项目新发现的一般定理",
                 "同一准备、同一H和原末读","无须先令实际电荷Q与材料h对易"):
        assert term in prose,term
    newfiles=[note,HERE/"joint_energy_source.py",HERE/"joint_energy_source_results.json",
        Path(__file__),HERE/"drafts/positive_source_decision.md",
        HERE/"drafts/common_positive_process_increment.md",HERE/"drafts/publish974.py",
        STAGE/"975/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ('截至973','截至974'),('n<=973','n<=974'),('list(range(1,974))','list(range(1,975))'),
    ('001—973轮共973份','001—974轮共974份'),('231—973的743份','231—974的744份'),
    ('973：共同实现的准备审计与平均来源边界','974：完整交互能源与原生记录的共同量子来源'),
    ('STAGE/"972/research_round_972_checks.json"','STAGE/"973/research_round_973_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=974,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=974,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_973=prev["cumulative_numbered_test_groups_from_972"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        original_material_field_preparation_and_effects_retained=True,
        coherent_fixed_positive_source_process=True,
        infinite_occupation_certificate_transported=True,
        supports_and_post_Newton_remainder_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'974/verify_round974.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 974：完整交互能源与原生记录的共同量子来源'
    previous='## 973：共同实现的准备审计与平均来源边界'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[974报告]({pre}research_note_974.md)以965完整材料—场能源作质量来源，'
        +'在固定静路径Newton模型中保原准备和读口：记录差>.99479、场效果差>.02241，'
        +'来源分布报告与均值替代差>.14709；不要求原生电荷与材料h对易。'
        +f'[结果]({pre}974/joint_energy_source_results.json) · '
        +f'[核验]({pre}974/research_round_974_checks.json)。正式974／累计3759，整体目标未完成。\n\n'
        +f'[共同实现增量]({pre}974/drafts/common_positive_process_increment.md)'
        +'给同一固定正过程的有限接口，运输无限占据误差；支撑、运动、后Newton和动态度规未认证。'
        +f'停止路径优化，接[975]({pre}975/drafts/STATUS.md)回全目标审稳定记录与非平衡资源。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—973轮共973份','001—974轮共974份')
    if p==paths[4]:s=s.replace('当前正式973／累计3758，973已结项','当前正式974／累计3759，974已结项')
    if p==paths[5]:
        s=s.replace('# 231—973轮阶段成果总览','# 231—974轮阶段成果总览',1)
        s=s.replace('231—973的743份','231—974的744份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至973）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至974）',1)
        rows={
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—974|974在同一原生材料／场准备上运输965无限占据误差到两路径正过程|仅声明输入及静路径域；不借其他材料或支撑的误差签收全理论|',
        'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—974|974完整H_L同时给内时序、量子质量来源与路径相干反作用，内部交换守恒|未认证运动、支撑应力、后Newton或动态度规；973宇宙边界仍保留|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—974|974同一准备保旧记录／场效果，并以固定路径效果检出均值替代差>.14709|有效模型的算符报告；不等于实际引力检测器或全模型证伪|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[973共同实现审计](../../973/drafts/common_realization_audit_v1.md)区分不同准备和联合任务。[974正过程增量](../../974/drafts/common_positive_process_increment.md)接965实际交互、旧记录与完整量子质量来源，范围是静路径Newton单极。停止路径优化；接[975](../../975/drafts/STATUS.md)审稳定记录、非平衡资源与方向性，不重做一般熵身份或修单模热化率。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 974 verifier.')
