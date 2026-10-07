"""Publish the common-parent matching result, preserving all earlier receipts."""
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'980/verify_round980.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"finite_thermal_records_results.json")',1)[0]
for a,b in [('range(776,980)','range(776,981)'),('Delivery checks for 980','Delivery checks for 981'),
            ('research_round_980_checks.json','research_round_981_checks.json')]:prefix=prefix.replace(a,b)
prefix=prefix.replace('    for path,keys in extra.items():',
    '    extra["981/drafts/common_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")\n'
    '    for path,keys in extra.items():')
checks=r'''    result=read(HERE/"native_response_source_results.json")
    assert result["round"]==981 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];r=result["response"];pulse=result["finite_pulse"]
    assert p["native_dimension"]==6 and p["U"]==1 and p["v"]==.1 and p["eta"]==0
    assert r["spectral_error"]<1e-13 and r["lapse_five_point_error"]<1e-10
    assert abs(r["frozen_proper_frequency_missing_fraction_at_point6"]-9/17)<1e-13
    assert r["derivative_expansion"][1]["alpha_relative_remainder"]==.0001
    assert r["derivative_expansion"][1]["lapse_relative_remainder"]<.0005
    assert pulse["certified_contrast_lower"]>5.1e-5 and pulse["lapse_slope_difference"]<1e-12
    assert result["static_source"]["hellmann_feynman_error"]<1e-11
    for key in ("full_SM_matching_completed","Einstein_background_solved",
        "dynamic_quantum_geometry_constructed","all_unknown_input_pulse_contrast",
        "finite_pulse_is_autonomous_instrument","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core981",HERE/"native_response_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    spec=importlib.util.spec_from_file_location("cert981",HERE/"certify_response_bounds.py")
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    assert cert.run()==read(HERE/"response_bounds_certificate.json")
    note=STAGE/"research_note_981.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","不使用上述低频导数截断",
        "不是总Maxwell—材料能源或完整T00","不是新发现Kubo公式",
        "停止本混合响应的器件和更高频优化"):
        assert term in prose,term
    newfiles=[note,HERE/"native_response_source.py",HERE/"native_response_source_results.json",
        HERE/"certify_response_bounds.py",HERE/"response_bounds_certificate.json",Path(__file__),
        HERE/"drafts/parent_matching_decision.md",HERE/"drafts/common_parent_contract_v1.md",
        HERE/"drafts/publish981.py",STAGE/"982/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至980','截至981'),('n<=980','n<=981'),('list(range(1,981))','list(range(1,982))'),
    ('001—980轮共980份','001—981轮共981份'),('231—980的750份','231—981的751份'),
    ('980：有限热材料与有效遗忘','981：共同父描述与电磁—时间来源匹配'),
    ('STAGE/"979/research_round_979_checks.json"','STAGE/"980/research_round_980_checks.json"')]:tail=tail.replace(a,b)
tail=tail.replace('oldfiles=[STAGE/"980/research_round_980_checks.json",HERE/"drafts/STATUS.md"]',
    'oldfiles=[STAGE/"980/research_round_980_checks.json",HERE/"drafts/STATUS.md",\n'
    '        HERE/"drafts/common_adoption_audit_checks.json"]')
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=981,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=981,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_980=prev["cumulative_numbered_test_groups_from_979"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        common_parent_contract_specified=True,
        native_causal_response_and_lapse_insertion_verified=True,
        low_frequency_relative_remainders_verified=True,
        rational_finite_pulse_contrast_verified=True,
        dynamic_quantum_geometry_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'981/verify_round981.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 981：共同父描述与电磁—时间来源匹配'
    previous='## 980后采用复核：先合并父模型，停止局部扩建'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[981报告]({pre}research_note_981.md)列明共同父描述的物种、作用阶与材料匹配，'
        +'核出同一原生频率响应的时间来源必须保留频率导数；有限读数差有严格下界>5.1×10⁻⁵。'
        +f'[结果]({pre}981/native_response_source_results.json) · '
        +f'[核验]({pre}981/research_round_981_checks.json)。正式981／累计3766，整体目标未完成。\n\n'
        +f'[父合同]({pre}981/drafts/common_parent_contract_v1.md)保留SM／Einstein与匹配输入；'
        +'本轮只签收电磁—lapse—材料响应，没有完成动态量子几何或全部门恢复。'
        +f'停止混合响应和脉冲优化，接[982]({pre}982/drafts/STATUS.md)回共同动态几何采用域。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—980轮共980份','001—981轮共981份')
    if p==paths[4]:s=s.replace('当前正式980／累计3765，980已结项','当前正式981／累计3766，981已结项')
    if p==paths[5]:
        s=s.replace('# 231—980轮阶段成果总览','# 231—981轮阶段成果总览',1)
        s=s.replace('231—980的750份','231—981的751份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至980）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至981）',1)
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；952、953、971、974—981|981同一材料因果核与lapse插入共同匹配；频率导数不可任意漏掉|指定坐标场源的材料响应，不是完整应力或自洽Einstein解；父匹配仍须保场和支撑|',
        'C21':'|C21 完整物质与相互作用|1、4、5；531—544、939、951、957、971、981|P981明列三代SM、必要时维五中微子项、有效系数与复合材料替代规则|字段表不是全部门正过程；旧E_rec单代、p/R及特定参数结果不可直接移植|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—981|981原生谱响应、混合来源与有限脉冲共用h/Q；读数差有理下界>5.1e−5|宽频脉冲不用低频截断；实体源／支撑与父物理误差未认证，不能由局部证书领取全统一|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[981父合同](../../981/drafts/common_parent_contract_v1.md)明列物种与匹配边界，原生响应新增共同时间来源证书。停止极化、脉冲、热接触与模拟器优化，接[982](../../982/drafts/STATUS.md)回动态几何的同过程采用域；全SM与宇宙未完成处保留。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 981 verifier.')
