"""Publish 973 once: fixed mean-field claims are not full common statistics."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'972/verify_round972.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"native_thermal_record_results.json")',1)[0]
for a,b in [('range(776,972)','range(776,973)'),('Delivery checks for 972','Delivery checks for 973'),
            ('research_round_972_checks.json','research_round_973_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"common_source_distribution_results.json")
    assert result["round"]==973 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    inputs=result["frozen_inputs"];stats=result["source_statistics"]
    assert inputs["modes"]==12 and inputs["mean_R"]==100
    assert .289<stats["coefficient_of_variation"]<.290
    assert .144<stats["initial_H_squared_relative_std"]<.145
    assert abs(stats["exact_variance"]-stats["numerical_variance"])<1e-8
    summ=result["summation"];diag=result["monopole_branch_diagnostic"]
    assert summ["chernoff_tail"]<6e-25
    assert abs(summ["normalization"]-1)<2e-12
    assert diag["report_certified_lower"]>.6725
    assert diag["report_under_mean_source"]==0
    assert diag["endpoint_time_margins"][0]>6 and diag["endpoint_time_margins"][1]>1
    assert diag["constraint_residual_at_mean_initial_momentum"]>1
    assert .0008<diag["mean_redshift_difference"]<.001
    assert result["mixture_witness"]["mixture_report"]==1
    assert result["mixture_witness"]["mean_source_report"]==0
    assert result["mixture_witness"]["initial_geometry_correlated_with_energy"]
    for key in ("full_stress_noise_or_Einstein_Langevin_solved",
                "original_969_mean_field_result_refuted","all_unknown_quantum_inputs_supported",
                "actual_geometric_detector_constructed","full_common_positive_quantum_model","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core973",HERE/"common_source_distribution.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_973.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是完整应力噪声或Einstein–Langevin解",
                 "保留969原结论","共同理论允许不同任务使用不同准备和边界",
                 "能源与初始几何相关的准备"):
        assert term in prose,term
    newfiles=[note,HERE/"common_source_distribution.py",HERE/"common_source_distribution_results.json",
        Path(__file__),HERE/"drafts/common_realization_decision.md",
        HERE/"drafts/common_realization_audit_v1.md",HERE/"drafts/publish973.py",
        STAGE/"974/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ('截至972','截至973'),('n<=972','n<=973'),('list(range(1,973))','list(range(1,974))'),
    ('001—972轮共972份','001—973轮共973份'),('231—972的742份','231—973的743份'),
    ('972：原生记录的热接触、保护扇区与内部交换','973：共同实现的准备审计与平均来源边界'),
    ('STAGE/"971/research_round_971_checks.json"','STAGE/"972/research_round_972_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=973,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=973,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_972=prev["cumulative_numbered_test_groups_from_971"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        actual_frozen_thermal_source_distribution_verified=True,
        bounded_report_in_monopole_diagnostic=True,
        source_geometry_constraint_correlation_retained=True,
        global_common_realization_audited=True,
        original_mean_field_refuted=False,full_stochastic_gravity_solved=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'973/verify_round973.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 973：共同实现的准备审计与平均来源边界'
    previous='## 972：原生记录的热接触、保护扇区与内部交换'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[973报告]({pre}research_note_973.md)核对964／969实际12模热准备：能源相对涨落约28.9%，'
        +'大占据数不使其趋零。另在明示逐能源壳单极诊断中，有界几何报告与均值源替代相差>.6725；'
        +'不是完整应力噪声或对969原均值解的反证。'
        +f'[结果]({pre}973/common_source_distribution_results.json) · '
        +f'[核验]({pre}973/research_round_973_checks.json)。正式973／累计3758，整体目标未完成。\n\n'
        +f'[共同实现审计]({pre}973/drafts/common_realization_audit_v1.md)区分同一理论的不同准备、重叠任务和同一次联合过程，'
        +'不强迫所有历史数值共处一次事件，也不拼接不同候选结项。'
        +f'停止本均值诊断优化，接[974]({pre}974/drafts/STATUS.md)选择保实际源和未知量子资料的固定正过程。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—972轮共972份','001—973轮共973份')
    if p==paths[4]:s=s.replace('当前正式972／累计3757，972已结项','当前正式973／累计3758，973已结项')
    if p==paths[5]:
        s=s.replace('# 231—972轮阶段成果总览','# 231—973轮阶段成果总览',1)
        s=s.replace('231—972的742份','231—973的743份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至972）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至973）',1)
        rows={
        'C19':'|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957、964—973|973原12模Gibbs能源统计精确给28.9%相对涨落，固定模数大占据不消除此项|不等于所有热源都不能经典化；单极诊断未解全应力噪声|',
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—973|973区分同理论不同准备与联合任务；均值源替代须按实际报告验收|969只保原均值合同；全未知量子资料与几何统计尚未共同签收|',
        'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—973|973逐能源壳诊断保来源与初始几何约束关联，均值量不能任意代替二阶约束|不是完整Einstein-Langevin；未把原平均场解判错|',
        'C24':'|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、325、331—340、961、964—973|969平均膨胀共存保留；973实际热准备及有界统计报告明确其扩展边界|分支是均匀单极附加诊断；相干来源、暗部门与初态生成仍开放|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—973|973单极分支与均值源对有界报告差>.6725；精确旧源统计和诊断预测分列|不是实际宇宙噪声预测、实际检测器或全模型证伪|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[971父接口](../../971/drafts/native_parent_map_v0_1.md)与[972热记录连接](../../972/drafts/thermal_record_mechanism.md)保留。[973共同实现审计](../../973/drafts/common_realization_audit_v1.md)区分不同准备和联合任务；969仅保均值范围，不冒领全来源统计。接[974](../../974/drafts/STATUS.md)复用940等固定正过程，保原生源算符与未知资料，不继续修12模宇宙。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 973 verifier.')
