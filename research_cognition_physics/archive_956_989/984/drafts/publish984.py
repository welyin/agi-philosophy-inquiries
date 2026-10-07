"""Publish one bounded source-to-geometry adoption result; retain frozen history."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'983/verify_round983.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"sm_curvature_radiation_results.json")',1)[0]
for a,b in [('range(776,983)','range(776,984)'),('Delivery checks for 983','Delivery checks for 984'),
            ('research_round_983_checks.json','research_round_984_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"smeared_source_geometry_results.json")
    assert result["round"]==984 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["inherited_matter"]["c"]["exact"]=="283/120"
    assert result["analytic"]["filter_integral"]["exact"]=="6883/113513400"
    assert result["rational_certificate"]["all_exact_inequalities_passed"]
    assert result["rational_certificate"]["Chebyshev_upper"]["decimal"]<6.118e-6
    control=result["negative_control"]
    assert abs(control["shared_conservation_variance"])<1e-15
    assert control["independent_component_conservation_variance"]>.011
    for key in ("actual_SM_thermal_state","Gaussian_stress_process_claimed",
        "finite_spacetime_detector_certified","full_quantum_geometry_error_certified",
        "full_physical_remainder_certified","UV_completion_assumed",
        "minimum_physical_scale_assumed","full_goal_completed"):assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core984",HERE/"smeared_source_geometry.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_984.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不是完整非线性Einstein方程的算符恒等式",
        "不必假定应力是Gaussian随机过程","有限带宽当作有限时空装置",
        "停止谱窗口、探测器及响应精度优化"):
        assert term in prose,term
    newfiles=[note,HERE/"smeared_source_geometry.py",HERE/"smeared_source_geometry_results.json",
        Path(__file__),HERE/"drafts/adoption_decision.md",HERE/"drafts/geometry_adoption_scope.md",
        HERE/"drafts/publish984.py",STAGE/"985/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至983','截至984'),('n<=983','n<=984'),('list(range(1,984))','list(range(1,985))'),
    ('001—983轮共983份','001—984轮共984份'),('231—983的753份','231—984的754份'),
    ('983：标准模型曲率来源与宇宙有效阶','984：共同量子来源与有限尺度的经典几何报告'),
    ('STAGE/"982/research_round_982_checks.json"','STAGE/"983/research_round_983_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=984,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=984,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_983=prev["cumulative_numbered_test_groups_from_982"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,actual_free_matter_source_kernel_mapped=True,
        finite_band_classical_report_bound_verified=True,
        same_source_response_and_conservation_checked=True,
        full_quantum_geometry_certified=False,full_physical_error_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'984/verify_round984.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')

paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 984：共同量子来源与有限尺度的经典几何报告'
    previous='## 983：标准模型曲率来源与宇宙有效阶'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[984报告]({pre}research_note_984.md)把P981自由共形物质真空的实际应力核'
        +'接到同源低频响应与谱过滤的线性Einstein报告；指定带宽与匹配下，'
        +'绝对容差10⁻⁵的超差概率上界<6.118×10⁻⁶。'
        +f'[结果]({pre}984/smeared_source_geometry_results.json) · '
        +f'[核验]({pre}984/research_round_984_checks.json)。正式984／累计3769，整体目标未完成。\n\n'
        +f'[采用范围]({pre}984/drafts/geometry_adoption_scope.md)区分受限几何报告、完整时空及有限实体装置；'
        +'没有认证真实SM或全量子几何总误差。停止频率窗口／探测器优化，'
        +f'接[985]({pre}985/drafts/STATUS.md)回共同模型的整体采用判断。应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—983轮共983份','001—984轮共984份')
    if p==paths[4]:s=s.replace('当前正式983／累计3768，983已结项','当前正式984／累计3769，984已结项')
    if p==paths[5]:
        s=s.replace('# 231—983轮阶段成果总览','# 231—984轮阶段成果总览',1)
        s=s.replace('231—983的753份','231—984的754份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至983）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至984）',1)
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；553、631、899、940、971、974—984|984实际物质应力核与低频几何响应共用c；交叉关联保守恒，删交叉有有限频率负对照|只在自由共形真空与线性诱导部门，非全SM／引力总余项|',
        'C21':'|C21 完整物质与相互作用|1、4、5；531—553、939、951、957、971、981、983—984|P981三代SM字段保持；983曲率系数与984自由真空谱共用4／45／12计数|质量、真实相互作用与热态核未恢复；不能将自由阶称完整物理误差证书|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—984|984同一源的指定几何报告有有理认证的超差概率界；来源与响应误差不分家|不认证有限时空实体装置、完整经典时空或所有多时分布；不同准备数值不可相加|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[984几何采用范围](../../984/drafts/geometry_adoption_scope.md)填入实际自由物质源与有限带宽报告，保完整时空和实体装置边界。停止本例优化，接[985](../../985/drafts/STATUS.md)先核共同模型的整体采用价值及真实交叠缺口，不自动继续修候选技术债。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 984 verifier.')
