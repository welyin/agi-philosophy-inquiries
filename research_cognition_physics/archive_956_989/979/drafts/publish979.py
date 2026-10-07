"""Publish the source-aware spectral bridge; preserve frozen research evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'978/verify_round978.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"moving_complete_source_results.json")',1)[0]
for a,b in [('range(776,978)','range(776,979)'),('Delivery checks for 978','Delivery checks for 979'),
            ('research_round_978_checks.json','research_round_979_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"spectral_mass_bridge_results.json")
    assert result["round"]==979 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    c=result["finite_contract"];p=result["protocol_transport"]
    assert c["noncommuting_error"]>1e-6
    assert c["low_mass_intertwining_error"]<1e-11 and c["minimum_mass"]>0
    assert c["actual_full_H_error"]<=c["full_H_error_bound"]
    assert c["translation_symmetry_error"]<1e-11
    assert c["exact_compressed_force_error"]<=c["exact_compressed_force_bound"]
    assert c["force_finite_difference_error"]<1e-10 and c["equal_opposite_force_difference"]<1e-10
    for row in c["rows"]:assert row["all_input_isometry_difference"]<=row["analytic_uniform_bound"]
    assert p["time"]==30 and p["G"]==1e-30
    assert p["gaussian_initial_bound"]<p["initial_H_minus_mu_norm_bound"]
    assert p["moving_simulation_error"]<.000148
    assert p["full_protocol_output_bound"]<.000293
    assert p["unknown_input_and_passive_reference"] and p["uniform_time_interval"]
    assert p["physical_simulator_generated"] is False
    b=result["resource_boundaries"]
    assert b["baseline_omission_probability_gap"]>.999999
    assert b["tail_examples"][-1]["excess_mass_mean"]>9999
    assert b["tail_examples"][-1]["all_time_trace_distance_upper"]<.021
    for key in ("exchange_coupling_list_generated","native_976_charge_model_identified_with_exchange_simulator",
        "all_SM_matching_completed","full_GR_completed","macro_arrow_completed",
        "resolved_gravity_signal_claimed_in_protocol_example","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core979",HERE/"spectral_mass_bridge.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_979.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","没有生成最终交换网络","数学Newton H的存在不自动认证",
        "不是947／976／978旧装置的新误差","高能辅助态的资源矩"):
        # The resource-tail phrase occurs in the next entry; prose carries the same explicit boundary.
        assert term in prose or term in (STAGE/"980/drafts/STATUS.md").read_text("utf-8"),term
    newfiles=[note,HERE/"spectral_mass_bridge.py",HERE/"spectral_mass_bridge_results.json",
        Path(__file__),HERE/"drafts/spectral_mass_decision.md",
        HERE/"drafts/protocol_parent_adoption.md",HERE/"drafts/publish979.py",
        STAGE/"980/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至978','截至979'),('n<=978','n<=979'),('list(range(1,979))','list(range(1,980))'),
    ('001—978轮共978份','001—979轮共979份'),('231—978的748份','231—979的749份'),
    ('978：完整原生能源与双体运动','979：完整协议与实际质量的谱模拟桥接'),
    ('STAGE/"977/research_round_977_checks.json"','STAGE/"978/research_round_978_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=979,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=979,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_978=prev["cumulative_numbered_test_groups_from_977"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        continuum_spectral_mass_transport_proved=True,
        existing_complete_protocol_budget_checked=True,
        baseline_and_high_energy_resource_boundaries_checked=True,
        actual_exchange_material_or_SM_matching_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'979/verify_round979.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 979：完整协议与实际质量的谱模拟桥接'
    previous='## 978：完整原生能源与双体运动'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[979报告]({pre}research_note_979.md)复用436静态模拟与929完整协议，证明同一编码可运输完整质量、逆质量动能及Newton来源。'
        +'明示质量匹配下，有限协议联合误差<.000293；能量基线与高能资源尾不能删掉。'
        +f'[结果]({pre}979/spectral_mass_bridge_results.json) · '
        +f'[核验]({pre}979/research_round_979_checks.json)。正式979／累计3764，整体目标未完成。\n\n'
        +f'[采用范围]({pre}979/drafts/protocol_parent_adoption.md)允许同一理论中的不同复合材料，'
        +'不强迫两轨道器件承担所有协议；未生成实际交换材料或完成SM／GR匹配。'
        +f'停止模拟器与器件优化，接[980]({pre}980/drafts/STATUS.md)回共同父层级和有限记录生命周期。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—978轮共978份','001—979轮共979份')
    if p==paths[4]:s=s.replace('当前正式978／累计3763，978已结项','当前正式979／累计3764，979已结项')
    if p==paths[5]:
        s=s.replace('# 231—978轮阶段成果总览','# 231—979轮阶段成果总览',1)
        s=s.replace('231—978的748份','231—979的749份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至978）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至979）',1)
        rows={
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、625、758—761、854、898、955—979|979把静态谱合同接到连续位置及逆质量动能；929完整协议条件预算<.000293|有界输出误差不自动给裸准备能源矩误差；物理父层及全部局域源仍需各自匹配|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—979|979复用436与929，允许不同复合材料；完整协议／质量有同编码条件桥接|交换权重、准备、稳定复合物及质量匹配仍输入；未等同976实际Q，未完成宏观箭头|',
        'C22':'|C22 引力来源与物质反作用|1、3、4、H2；940—947、953、971、974—979|979同一实际H及能量基线进入惯性、Newton源和相反力；近低能谱可运输运动|领先Newton单极输入；高能尾资源及GR有效范围不能从近编码自动领取|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—979|979全协议误差<.000293来自条件性合同；完整质量遗漏可产生有限读数差|小矩阵是合同核查，非实际交换网络；示例极弱引力未认证可分辨信号|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[979采用增量](../../979/drafts/protocol_parent_adoption.md)接通静态谱模拟与完整质量的条件合同，不强迫单一两轨道材料承担所有功能。停止模拟器优化，接[980](../../980/drafts/STATUS.md)回共同父层级、有限记录生命周期及全部门有效范围。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 979 verifier.')
