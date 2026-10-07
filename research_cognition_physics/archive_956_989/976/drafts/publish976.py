"""Publish 976 finite-reference result without changing frozen evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'975/verify_round975.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"record_resource_audit_results.json")',1)[0]
for a,b in [('range(776,975)','range(776,976)'),('Delivery checks for 975','Delivery checks for 976'),
            ('research_round_975_checks.json','research_round_976_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"finite_internal_reference_results.json")
    assert result["round"]==976 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    contract=result["fixed_relational_readout"];cert=result["interval_certificate"]
    assert contract["effect"]=="(I+SWAP_BR)/2"
    assert contract["free_common_invariance_residual"]<1e-14
    assert contract["reference_preparation_is_resource"] and contract["isolation_is_additional_input"]
    assert contract["external_time_dependent_effect"] is False
    assert cert["whole_interval_not_only_grid_certified"]
    assert cert["uniform_infinite_contrast_lower"]>.48594
    assert cert["grid_points"]==257 and cert["mesh_gap_bound"]<.000956
    assert cert["binary_information_lower"]>.1232
    assert result["infinite_occupation_certificate"]["total_state_vector_error"]<2.733e-8
    account=result["reference_account"]
    assert account["source_includes_complete_reference_energy"]
    assert account["source_mass_positive_lower"]>199.68
    assert account["old_endpoint_check"]<1e-9
    assert account["old_path_predictions_unchanged"] is False
    assert account["internal_energy_variance"]>.00037
    for row in account["source_rows"]:
        assert row["old_path_visibility"]-row["new_path_visibility"]>.0012
        assert row["mass_mean_with_reference"]>200
    for key in ("phase_record_converted_to_protected_population","actual_joint_apparatus_derived",
        "arbitrary_noise_stability_proven","repeatable_nondestructive_readout_proven",
        "macroscopic_irreversible_memory_proven","full_microscopic_SM_matching_proven",
        "moving_support_or_quantum_metric_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core976",HERE/"finite_internal_reference.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_976.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","隔离及共同支撑是输入","不是新的普遍认知公理",
        "尚未给具体仪器、读后能源账或复位","不得把其器件缺口升为全纲领门槛"):
        assert term in prose,term
    newfiles=[note,HERE/"finite_internal_reference.py",HERE/"finite_internal_reference_results.json",
        Path(__file__),HERE/"drafts/finite_reference_decision.md",
        HERE/"drafts/finite_record_adoption.md",HERE/"drafts/publish976.py",
        STAGE/"977/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至975','截至976'),('n<=975','n<=976'),('list(range(1,976))','list(range(1,977))'),
    ('001—975轮共975份','001—976轮共976份'),('231—975的745份','231—976的746份'),
    ('975：原生记录的资源账与保持接口','976：有限内部参考与关系记录'),
    ('STAGE/"974/research_round_974_checks.json"','STAGE/"975/research_round_975_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=976,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=976,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_975=prev["cumulative_numbered_test_groups_from_974"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        finite_internal_reference_and_whole_window_verified=True,
        reference_mass_and_source_predictions_updated=True,
        macroscopic_irreversible_memory_proven=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'976/verify_round976.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 976：有限内部参考与关系记录'
    previous='## 975：原生记录的资源账与保持接口'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[976报告]({pre}research_note_976.md)在974原局部过程加入同种材料参考，固定关系读口在整个预定有限区间的标签对比>.485948。'
        +'参考完整能源进入质量源，旧材料边缘保持，路径相干预测随之变化。'
        +f'[结果]({pre}976/finite_internal_reference_results.json) · '
        +f'[核验]({pre}976/research_round_976_checks.json)。正式976／累计3761，整体目标未完成。\n\n'
        +f'[采用条件]({pre}976/drafts/finite_record_adoption.md)保留参考准备、隔离、支撑和联合读口输入；'
        +'不等于P人口转换、实体读后能源账或宏观不可逆记忆。'
        +f'停止保持方案优化，接[977]({pre}977/drafts/STATUS.md)回共同模型的量子来源与经典几何有效域。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—975轮共975份','001—976轮共976份')
    if p==paths[4]:s=s.replace('当前正式975／累计3760，975已结项','当前正式976／累计3761，976已结项')
    if p==paths[5]:
        s=s.replace('# 231—975轮阶段成果总览','# 231—976轮阶段成果总览',1)
        s=s.replace('231—975的745份','231—976的746份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至975）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至976）',1)
        rows={
        'C19':'|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—976|976参考的纯态准备、相干及完整能源入账；975比较热参考保持原范围|隔离和对齐准备是输入，非实际热库；973热源统计边界保留|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—976|976同种有限参考及固定关系读口实现整段时间可读，区间对比>.485948|不等于P人口保护；实际仪器／读后能源／宏观保持仍未签收|',
        'C22':'|C22 引力来源与物质反作用|1、3、4、H2；940—947、953、971、974、976|976完整参考能源进入974同一量子质量源，路径相干随新增参考改变|静态支撑、弱场输入；自引力、运动支撑、后Newton、动态度规未认证|',
        'C23':'|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—976|975资源／关联账接976有限关系记录，整体幺正不妨碍声明有限可读时间|无宏观不可逆记忆或完整循环；不把永久寿命、无限重复设成前置|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—976|976全区间记录对比>.485948，已含无限占据残差和采样间隙；参考源同步更新|误差相对明示H；未代替全部父理论遗漏或实际仪器验证|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[973共同实现审计](../../973/drafts/common_realization_audit_v1.md)保全部部门。[976有限记录采用](../../976/drafts/finite_record_adoption.md)以同种参考和完整来源闭合一次有限合同；停止保持优化，接[977](../../977/drafts/STATUS.md)回量子来源、统计／条件化报告与经典几何的共同有效域。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 976 verifier.')
