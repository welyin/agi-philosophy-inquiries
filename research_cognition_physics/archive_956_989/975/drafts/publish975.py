"""Publish the 975 same-process resource audit, preserving all frozen evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'974/verify_round974.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"joint_energy_source_results.json")',1)[0]
for a,b in [('range(776,974)','range(776,975)'),('Delivery checks for 974','Delivery checks for 975'),
            ('research_round_974_checks.json','research_round_975_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"record_resource_audit_results.json")
    assert result["round"]==975 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["unchanged_physical_hamiltonian_preparation_and_endpoint"]
    ref=result["reference"];bal=result["resource_balance"];rec=result["record"]
    assert abs(ref["beta"]-2*math.log(2))<1e-14
    assert ref["reference_is_actual_joint_equilibrium"] is False
    assert ref["physical_heat_bath_added"] is False
    assert 6.78<bal["initial_readiness"]<6.79
    assert 5.25<bal["final_local_readiness"]<5.26
    assert 1.53<bal["final_total_correlation"]<1.531
    assert abs(bal["identity_residual"])<1e-10
    assert bal["interaction_energy_change"]<-1.6e-5
    assert rec["certified_binary_information_lower"]>.67507
    assert rec["sender_retains_label_exactly"]
    assert rec["arbitrary_unknown_quantum_state_copied"] is False
    assert max(abs(x-.5) for x in rec["receiver_conserved_sector_probabilities"])<1e-11
    assert rec["after_sector_dephasing_infinite_distance_upper"]<.000106
    assert rec["sector_dephasing_is_diagnostic_only"]
    errs=result["error_bounds"]
    assert errs["oscillator_mean_number_bound"]==1
    assert errs["derived_sqrt_mean_number_bound"]<.797
    assert errs["total_correlation"]<.000183 and errs["local_readiness"]<.000184
    assert errs["infinite_oscillator_entropy_controlled"]
    for key in ("global_entropy_production","actual_heat_or_work_cost_identified",
                "stable_macroscopic_record_dynamics_completed","full_cyclic_reset_completed",
                "cosmic_low_entropy_origin_derived","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core975",HERE/"record_resource_audit.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_975.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","乘积τ不是完整相互作用H的平衡态",
                 "没有把此退相位当成已经发生的环境噪声","不改变974的H、准备",
                 "新记录尚未写入现有作用所保护的变量"):
        assert term in prose,term
    newfiles=[note,HERE/"record_resource_audit.py",HERE/"record_resource_audit_results.json",
        Path(__file__),HERE/"drafts/record_resource_decision.md",
        HERE/"drafts/record_arrow_adoption.md",HERE/"drafts/publish975.py",
        STAGE/"976/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ('截至974','截至975'),('n<=974','n<=975'),('list(range(1,975))','list(range(1,976))'),
    ('001—974轮共974份','001—975轮共975份'),('231—974的744份','231—975的745份'),
    ('974：完整交互能源与原生记录的共同量子来源','975：原生记录的资源账与保持接口'),
    ('STAGE/"973/research_round_973_checks.json"','STAGE/"974/research_round_974_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=975,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=975,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_974=prev["cumulative_numbered_test_groups_from_973"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        unchanged_joint_process_resource_balance_verified=True,
        same_record_distinguished_from_protected_sector=True,
        infinite_oscillator_entropy_error_controlled=True,
        global_entropy_increase_or_macroscopic_arrow_proven=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'975/verify_round975.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 975：原生记录的资源账与保持接口'
    previous='## 974：完整交互能源与原生记录的共同量子来源'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[975报告]({pre}research_note_975.md)保持974同一过程，核出局部资源转成约1.53自然信息单位的关联，'
        +'总熵及完整能源保持；原发送标签仍在，接收信息下界>.67507。'
        +'但新记录不在受保护P人口中，扇区退相位后可区分度<.000106。'
        +f'[结果]({pre}975/record_resource_audit_results.json) · '
        +f'[核验]({pre}975/research_round_975_checks.json)。正式975／累计3760，整体目标未完成。\n\n'
        +f'[采用条件]({pre}975/drafts/record_arrow_adoption.md)区分转存、新记录、保持及复位；'
        +'比较热参考不是新热库，关联预算不等于宏观箭头。'
        +f'停止资源／单模式优化，接[976]({pre}976/drafts/STATUS.md)筛选新记录到有限保持载体的真实接口。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—974轮共974份','001—975轮共975份')
    if p==paths[4]:s=s.replace('当前正式974／累计3759，974已结项','当前正式975／累计3760，975已结项')
    if p==paths[5]:
        s=s.replace('# 231—974轮阶段成果总览','# 231—975轮阶段成果总览',1)
        s=s.replace('231—974的744份','231—975的745份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至974）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至975）',1)
        rows={
        'C19':'|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957、964—975|975核974实际准备的局部非平衡指标、关联及交互能；未新增热库|比较Gibbs参考不是实际联合平衡；973热源统计边界保持|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—975|975原发送P标签保持、接收形成第二份经典记录，但新记录主要在跨P相位|P保护不等于新记录稳定；实际转换／有限保持接口待验|',
        'C23':'|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—975|975同过程资源账闭合，总关联约1.53而全局熵不增；不把局部关联当全局耗散|没有稳定宏观记录、循环复位或箭头；有限合同不要求永久、持续散热|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—975|975旧末读有>.67507标签信息，受保护人口却不区分；P退相位诊断对比<.000106|退相位只作诊断，未声称发生实际噪声或测得寿命|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[973共同实现审计](../../973/drafts/common_realization_audit_v1.md)保全部部门。[975记录采用条件](../../975/drafts/record_arrow_adoption.md)在974实际过程中闭合资源账并区分相位记录与P保护；接[976](../../976/drafts/STATUS.md)先筛选记录写入与有限保持的共同权限，停止单模式、容量和寿命优化。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 975 verifier.')
