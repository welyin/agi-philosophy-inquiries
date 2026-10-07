"""Publish the two-active-source subinterface; preserve the preceding audit."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'984/verify_round984.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"smeared_source_geometry_results.json")',1)[0]
for a,b in [('range(776,984)','range(776,985)'),('Delivery checks for 984','Delivery checks for 985'),
            ('research_round_984_checks.json','research_round_985_checks.json')]:prefix=prefix.replace(a,b)
prefix=prefix.replace('    for path,keys in extra.items():',
    '    extra["985/drafts/common_model_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")\n'
    '    for path,keys in extra.items():')
checks=r'''    result=read(HERE/"two_active_sources_results.json")
    assert result["round"]==985 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("core985",HERE/"two_active_sources.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    budget=result["uniform_budget"]
    assert Fraction(budget["exact_sum_upper"]["exact"])<Fraction(1721,10**6)
    assert Fraction(budget["record_contrast_lower"]["exact"])>Fraction(482506,10**6)
    assert result["initial_moments"]["unknown_joint_two_sender_input_and_passive_reference"]
    assert result["initial_moments"]["all_Fock_occupations_retained_in_theorem"]
    assert result["stability"]["relative_negative_form_bound"]<.00125
    for key in ("entire_U1_completed","radiative_TT_matched","full_SM_matching",
                "GR_total_physical_error_certified","physical_minimum_scale_assumed",
                "microscopic_continuity_assumed","UV_completion_claimed",
                "unbounded_source_error_deduced_from_probability","full_goal_completed"):
        assert result["scope"][key] is False
    note=STAGE/"research_note_985.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","停止双体扩建和高能延拓优化",
        "能量无下界不等于量子态不正或不存在酉演化","本轮误差不是距完整GR的误差"):
        assert term in prose,term
    newfiles=[note,HERE/"two_active_sources.py",HERE/"two_active_sources_results.json",
        Path(__file__),HERE/"drafts/two_active_sources_decision.md",HERE/"drafts/two_source_adoption.md",
        HERE/"drafts/publish985.py",STAGE/"986/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至984','截至985'),('n<=984','n<=985'),('list(range(1,985))','list(range(1,986))'),
            ('001—984轮共984份','001—985轮共985份'),('231—984的754份','231—985的755份'),
            ('984：共同量子来源与有限尺度的经典几何报告','985：两个活动来源与有限域共同过程'),
            ('STAGE/"983/research_round_983_checks.json"','STAGE/"984/research_round_984_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=985,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=985,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_984=prev["cumulative_numbered_test_groups_from_983"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,two_active_native_sources_joined=True,
        finite_domain_extension_and_record_bounds_verified=True,
        entire_U1_completed=False,full_physical_error_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'985/verify_round985.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')

paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 985：两个活动来源与有限域共同过程'
    previous='## 984后整体采用审计：先接共同物理过程'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[985报告]({pre}research_note_985.md)将978的第二常数内部质量提升为同型活动量子模块；'
        +'双方完整来源、运动和原A记录在同一Newton有效过程中共存。'
        +'全联合态误差<.001721，记录对比>.482506；高能延拓影响有实际尾界。'
        +f'[结果]({pre}985/two_active_sources_results.json) · '
        +f'[核验]({pre}985/research_round_985_checks.json)。正式985／累计3770，整体目标未完成。\n\n'
        +f'[采用范围]({pre}985/drafts/two_source_adoption.md)只关闭U1的双方活动来源子接口；'
        +'没有认证完整传播、GR／SM匹配或全部门共同模型。停止双体和延拓优化，'
        +f'接[986]({pre}986/drafts/STATUS.md)回传播／约束／材料的共同保留阶。应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—984轮共984份','001—985轮共985份')
    if p==paths[4]:s=s.replace('当前正式984／累计3769，984已结项','当前正式985／累计3770，985已结项')
    if p==paths[5]:
        s=s.replace('# 231—984轮阶段成果总览','# 231—985轮阶段成果总览',1)
        s=s.replace('231—984的754份','231—985的755份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至984）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至985）',1)
        rows={
          'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；553、631、899、940、971、974—985|985将两份活动内部模块接到同一Newton过程，双方完整源与相反力共用H；保留窗内不受所选高能延拓影响|只关闭U1子接口，完整传播／约束／材料共同匹配仍待；不把内部误差当全GR误差|',
          'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—985|985对未知双发送资料及旁观参考给统一态界，A旧记录对比>.482506；有限域外延拓有实际尾界|巨大基线与原准备仍输入；无界来源误差不能由概率界代替，各准备数值不能拼成全部门证书|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[985采用范围](../../985/drafts/two_source_adoption.md)关闭两个活动内部来源的Newton子接口，整体U1—U3仍未结项。停止双体与高能延拓优化，接[986](../../986/drafts/STATUS.md)回同一保留阶的传播／约束／材料匹配，不以局部扩建替代整体恢复。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 985 verifier.')
