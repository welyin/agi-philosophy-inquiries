"""Publish the P981 common curvature/radiation-order adoption result."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'982/verify_round982.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"native_radiative_channels_results.json")',1)[0]
for a,b in [('range(776,982)','range(776,983)'),('Delivery checks for 982','Delivery checks for 983'),
            ('research_round_982_checks.json','research_round_983_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"sm_curvature_radiation_results.json")
    assert result["round"]==983 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    m=result["inherited_matter"];a=result["analytic"];b=result["benchmark"]
    assert (m["real_scalars"],m["Weyl_components"],m["gauge_vectors"])==(4,45,12)
    assert m["a"]["exact"]=="1991/720" and m["c"]["exact"]=="283/120"
    assert m["thermal_gstar"]["exact"]=="427/4"
    assert m["beta_curvature_times_16pi2"]==dict(bR="0",bW="283/120",bE="-1991/720")
    assert a["wrong_radiation_pressure_conservation_coefficient"]=="-4"
    assert a["missing_pressure_conservation_coefficient"]=="-5"
    assert a["proper_time_error_upper"]["exact"]=="63/25600000000"
    assert 0<b["reduced_vs_low_time_difference"]<a["proper_time_error_upper"]["decimal"]
    assert b["shared_equations_max_residual"]<2e-14
    assert float(b["branches"]["z_times_high"])>.98
    for key in ("full_interacting_SM_thermal_state","all_quantum_source_to_classical_geometry_error",
        "full_physical_truncation_error_certified","high_curvature_root_adopted",
        "cosmic_arrow_generated","full_goal_completed"):assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core983",HERE/"sm_curvature_radiation.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_983.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不重新算成新发现","β_bR=0不表示有限bR=0",
        "不是全部物理遗漏的上界","同一份量子物质修正", "停止反常宇宙分支"):
        assert term in prose,term
    newfiles=[note,HERE/"sm_curvature_radiation.py",HERE/"sm_curvature_radiation_results.json",
        Path(__file__),HERE/"drafts/curvature_adoption_decision.md",
        HERE/"drafts/parent_recovery_scope_v2.md",HERE/"drafts/publish983.py",STAGE/"984/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至982','截至983'),('n<=982','n<=983'),('list(range(1,983))','list(range(1,984))'),
    ('001—982轮共982份','001—983轮共983份'),('231—982的752份','231—983的753份'),
    ('982：同一材料的电磁级联与动态几何通道','983：标准模型曲率来源与宇宙有效阶'),
    ('STAGE/"981/research_round_981_checks.json"','STAGE/"982/research_round_982_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=983,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=983,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_982=prev["cumulative_numbered_test_groups_from_981"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        inherited_SM_curvature_coefficients_mapped=True,
        same_source_density_pressure_trace_and_geometry_verified=True,
        finite_low_branch_time_bound_verified=True,
        full_physical_error_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'983/verify_round983.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')

paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 983：标准模型曲率来源与宇宙有效阶'
    previous='## 982：同一材料的电磁级联与动态几何通道'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[983报告]({pre}research_note_983.md)复用553／601，将P981自由共形SM阶的曲率修正'
        +'共同接到能量、压力、trace与FLRW几何；错配压力产生一阶守恒缺陷。'
        +'指定有限膨胀区间内，降阶与小根的无量纲时间差<2.461×10⁻⁹。'
        +f'[结果]({pre}983/sm_curvature_radiation_results.json) · '
        +f'[核验]({pre}983/research_round_983_checks.json)。正式983／累计3768，整体目标未完成。\n\n'
        +f'[采用范围v2]({pre}983/drafts/parent_recovery_scope_v2.md)明确这是两个保留方程的比较，'
        +'没有认证完整SM热相互作用或量子到经典几何误差。'
        +f'停止特殊宇宙背景优化，接[984]({pre}984/drafts/STATUS.md)回实际共同来源与经典几何任务。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—982轮共982份','001—983轮共983份')
    if p==paths[4]:s=s.replace('当前正式982／累计3767，982已结项','当前正式983／累计3768，983已结项')
    if p==paths[5]:
        s=s.replace('# 231—982轮阶段成果总览','# 231—983轮阶段成果总览',1)
        s=s.replace('231—982的752份','231—983的753份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至982）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至983）',1)
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；553、601、971、974—983|983同一曲率来源的能量、压力、trace及几何保留阶相容；错配压力是一阶缺陷|保留方程的自洽不提供完整量子几何通道，也不覆盖实际物理全部余项|',
        'C21':'|C21 完整物质与相互作用|1、4、5；531—553、939、951、957、971、981、983|P981三代SM表保持；983沿553纯SM的4实标量／45Weyl／12矢量曲率系数|只用自由共形活跃物质阶；质量、真实热相互作用和维五插入未恢复，旧扩展物种不可移植|',
        'C24':'|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、553、601、964—973、983|983同一SM自由共形阶的有限膨胀时间和曲率来源共同采用，低曲率降阶有明确范围|均匀平均、Λ=0和共形热态为输入；源噪声、量子几何、暗部门及初态生成仍开放|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—983|982通道竞争保留；983有限固有时间给解析比较上界，遗漏压力有独立负对照|两轮背景和输入不同，误差不可直接相加；都未领取全部父物理匹配|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[983共同恢复范围v2](../../983/drafts/parent_recovery_scope_v2.md)把982量子通道与983半经典物质来源分清。停止三模与反常宇宙优化，接[984](../../984/drafts/STATUS.md)核实际正量子源到有限经典几何任务的适用范围；不复证一般均场边界。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 983 verifier.')
