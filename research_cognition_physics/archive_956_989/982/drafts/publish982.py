"""Publish one bounded common-channel result; preserve prior frozen records."""
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'981/verify_round981.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"native_response_source_results.json")',1)[0]
for a,b in [('range(776,981)','range(776,982)'),('Delivery checks for 981','Delivery checks for 982'),
            ('research_round_981_checks.json','research_round_982_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"native_radiative_channels_results.json")
    assert result["round"]==982 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];a=result["analytic"]
    assert p["native_dimension"]==6 and p["exact_active_dimension"]==3
    assert p["U"]==1 and p["v"]==.1 and p["eta"]==0
    assert p["g_EM"]==p["g_TT"] and p["g_EM"]==p["J"]/40000
    assert max(result["native_intertwiner_errors"].values())<1e-13
    assert a["full_Fock_isometry_bound"]<a["rational_isometry_upper"]==.00372
    assert .38<a["rwa_graviton_probability"]<.4
    assert a["isolated_RWA_error_lower"]==.59628
    assert [r["full_dimension"] for r in result["matrix_witnesses"]]==[192,375]
    for row in result["matrix_witnesses"]:
        assert row["isometry_error"]<a["full_Fock_isometry_bound"]
        assert row["energy_error"]<1e-10 and row["source_equation_residual"]<1e-12
    for key in ("actual_mode_environment_matching_certified","qstar_derived_from_flat_material",
        "full_SM_matching_completed","nonlinear_Einstein_dynamics_completed",
        "arbitrary_quantum_source_to_classical_geometry","real_graviton_detection_claim",
        "full_goal_completed"):assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core982",HERE/"native_radiative_channels.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    spec=importlib.util.spec_from_file_location("cert982",HERE/"certify_radiative_bounds.py")
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    assert cert.run()==read(HERE/"radiative_bounds_certificate.json")
    note=STAGE/"research_note_982.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","不成为新的认知原则","相等耦合是本轮参数输入",
        "不是实际SM复合物体的四极矩构造","停止本辐射见证的效率",
        "不包括实际父理论到三个保留模式的误差δ_match"):
        assert term in prose,term
    newfiles=[note,HERE/"native_radiative_channels.py",HERE/"native_radiative_channels_results.json",
        HERE/"certify_radiative_bounds.py",HERE/"radiative_bounds_certificate.json",Path(__file__),
        HERE/"drafts/radiative_adoption_decision.md",HERE/"drafts/radiative_adoption.md",
        HERE/"drafts/publish982.py",STAGE/"983/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至981','截至982'),('n<=981','n<=982'),('list(range(1,982))','list(range(1,983))'),
    ('001—981轮共981份','001—982轮共982份'),('231—981的751份','231—982的752份'),
    ('981：共同父描述与电磁—时间来源匹配','982：同一材料的电磁级联与动态几何通道'),
    ('STAGE/"980/research_round_980_checks.json"','STAGE/"981/research_round_981_checks.json"')]:tail=tail.replace(a,b)
tail=tail.replace('oldfiles=[STAGE/"981/research_round_981_checks.json",HERE/"drafts/STATUS.md",\n'
        '        HERE/"drafts/common_adoption_audit_checks.json"]',
        'oldfiles=[STAGE/"981/research_round_981_checks.json",HERE/"drafts/STATUS.md"]')
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=982,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=982,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_981=prev["cumulative_numbered_test_groups_from_980"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        same_native_electromagnetic_and_TT_channels_verified=True,
        rational_full_Fock_finite_time_bound_verified=True,
        mode_Hamiltonian_energy_and_backreaction_verified=True,
        actual_parent_to_mode_matching_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'982/verify_round982.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')

paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 982：同一材料的电磁级联与动态几何通道'
    previous='## 981：共同父描述与电磁—时间来源匹配'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[982报告]({pre}research_note_982.md)保留同一原生材料的电磁级联和显式匹配的TT四极通道，'
        +'同一次过程的几何模读数约.386406；删电磁的孤立旋波预测为1。'
        +'指定有限时间的无限Fock等距误差<.00372。'
        +f'[结果]({pre}982/native_radiative_channels_results.json) · '
        +f'[核验]({pre}982/research_round_982_checks.json)。正式982／累计3767，整体目标未完成。\n\n'
        +f'[采用范围]({pre}982/drafts/radiative_adoption.md)保留四极矩、模式环境、支撑和耦合比输入；'
        +'三模内部误差不是实际父理论匹配误差。'
        +f'停止辐射与器件优化，接[983]({pre}983/drafts/STATUS.md)回共同恢复范围和量子／经典几何接口。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—981轮共981份','001—982轮共982份')
    if p==paths[4]:s=s.replace('当前正式981／累计3766，981已结项','当前正式982／累计3767，982已结项')
    if p==paths[5]:
        s=s.replace('# 231—981轮阶段成果总览','# 231—982轮阶段成果总览',1)
        s=s.replace('231—981的751份','231—982的752份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至981）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至982）',1)
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；952、953、971、974—982|981共同电磁／时间核；982同一H保原电磁级联、TT来源和反作用及能源|三模内部证书不替代父物理匹配；高阶来源和真实环境未领取|',
        'C22':'|C22 引力来源与物质反作用|1、3、4、H2；940—947、953、971、974—982|978—979质量与Newton接口保持；982加入明确四极匹配下的动态TT联合过程|平直h/Q不确定质量四极矩；q*、支撑、模式选择及非线性GR范围仍输入／未核|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—982|982同一次三模过程读数约.386406，删电磁孤立预测1；有限时间等距误差<.00372|只签收所列Hamiltonian内界，实际模式环境匹配及全物理共同恢复未完成|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[981父合同](../../981/drafts/common_parent_contract_v1.md)与[982采用记录](../../982/drafts/radiative_adoption.md)区分共同通道与实际物理匹配。停止三模、极化、热接触与模拟器优化，接[983](../../983/drafts/STATUS.md)回共同恢复和量子／经典几何适用范围；全SM与宇宙缺口保持。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 982 verifier.')
