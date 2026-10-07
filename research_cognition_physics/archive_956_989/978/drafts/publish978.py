"""Publish 978 movement/source bridge and preserve old research evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'977/verify_round977.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"coarse_source_record_results.json")',1)[0]
for a,b in [('range(776,977)','range(776,978)'),('Delivery checks for 977','Delivery checks for 978'),
            ('research_round_977_checks.json','research_round_978_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"moving_complete_source_results.json")
    assert result["round"]==978 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["uniform_moving_error"];a=result["operator_bounds"]
    assert p["G"]==1e-30 and p["reference_radius"]==1e10 and p["position_sigma"]==1e3
    assert p["mu"]>2e35
    assert a["internal_occupation_not_physically_cutoff"]
    assert a["initial_rms_mass_upper"]==204 and a["analytic_triangle_mass_upper"]<204
    assert a["weak_potential_global_upper"]<.003
    assert a["full_H_lower_excluding_mu"]>198
    assert b["total"]<.001721
    assert b["actual_moving_record_contrast_lower"]>.4825
    assert b["moving_output_to_977_cq_bound"]<.010383
    assert b["arbitrary_sender_and_passive_reference"] and b["interval_starts_at_zero"]
    assert b["tail_probability_upper"]>0 and b["minimum_tail_exponent"]>1e11
    c=result["conservation"]
    for key in ("total_momentum_conserved_exactly","total_energy_conserved_exactly",
        "full_internal_mass_conserved_exactly","both_centres_are_dynamical",
        "source_is_complete_internal_energy"):assert c[key]
    assert c["force_identity_difference"]<1e-12 and c["speed_norm_upper"]<1e-5
    for row in result["spatial_checks"]:
        assert row["quadrature_second_moment"]<row["analytic_second_moment_upper"]
        assert row["force_Y_is_exact_opposite"]
    for key in ("actual_wavepacket_evolution_computed","force_signal_above_measurement_error_claimed",
        "full_post_Newtonian_budget","dynamical_Einstein_metric_derived",
        "entire_SM_matching_completed","actual_recombination_device_completed","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core978",HERE/"moving_complete_source.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_978.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==12
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,7)]
    for term in ("整体目标未完成","不直接继承945的旧误差数值",
        "尚未计入实体装置及读后能源","不是完整Einstein动态时空",
        "引力作用是物理输入，不是认知原则生成的定理"):
        assert term in prose,term
    newfiles=[note,HERE/"moving_complete_source.py",HERE/"moving_complete_source_results.json",
        Path(__file__),HERE/"drafts/moving_source_decision.md",
        HERE/"drafts/common_recovery_audit_v2.md",HERE/"drafts/publish978.py",
        STAGE/"979/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至977','截至978'),('n<=977','n<=978'),('list(range(1,978))','list(range(1,979))'),
    ('001—977轮共977份','001—978轮共978份'),('231—977的747份','231—978的748份'),
    ('977：粗来源标签与块内量子记录','978：完整原生能源与双体运动'),
    ('STAGE/"976/research_round_976_checks.json"','STAGE/"977/research_round_977_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=978,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=978,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_977=prev["cumulative_numbered_test_groups_from_976"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        native_full_energy_and_two_moving_positions_verified=True,
        uniform_moment_based_transport_bound_verified=True,
        common_momentum_energy_and_force_identity_verified=True,
        full_GR_or_post_Newtonian_matching_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'978/verify_round978.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 978：完整原生能源与双体运动'
    previous='## 977：粗来源标签与块内量子记录'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[978报告]({pre}research_note_978.md)将976完整原生能源接到两个动态位置的同一Newton有效H，'
        +'利用实际质量二阶矩给全时间联合态误差<.001721；记录对比>.482507，来源、相反力和总动量／能源共用规则。'
        +f'[结果]({pre}978/moving_complete_source_results.json) · '
        +f'[核验]({pre}978/research_round_978_checks.json)。正式978／累计3763，整体目标未完成。\n\n'
        +f'[全部门审计v2]({pre}978/drafts/common_recovery_audit_v2.md)保物质、热／箭头、动态几何及宇宙缺口；'
        +'巨大但有限的第二质量、波包、软核和重合读口仍为输入，不领取全GR误差。'
        +f'停止运动方案优化，接[979]({pre}979/drafts/STATUS.md)回共同父描述及保留阶的整体合并。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—977轮共977份','001—978轮共978份')
    if p==paths[4]:s=s.replace('当前正式977／累计3762，977已结项','当前正式978／累计3763，978已结项')
    if p==paths[5]:
        s=s.replace('# 231—977轮阶段成果总览','# 231—978轮阶段成果总览',1)
        s=s.replace('231—977的747份','231—978的748份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至977）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至978）',1)
        rows={
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、625、758—761、854、898、955—978|978用实际M二阶矩将976运输到移动双体，全态误差<.001721；同菜单到977 CQ误差<.010383|界相对Newton模型，未给全GR保留阶及全部无界应力误差|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—978|978原材料与参考共同运动后，旧区间记录对比>.482507；未改Q／h实际交互|波包、第二有限物体和重合权限为输入；不是实体装置或宏观记录箭头|',
        'C22':'|C22 引力来源与物质反作用|1、3、4、H2；940—947、953、971、974—978|978完整原生M同时进入惯性、软核Newton源和相反力；两个位置动态，总能源／动量守恒|领先Newton物理输入；未领取后Newton、动态Einstein度规及全部父理论匹配|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—978|978全联合态与旧固定路径误差<.001721，记录界共同运输，初始力源身份独立校准|力的非零数值不是具体仪器可分辨证书；全GR长时相位仍无总界|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[978全部门审计v2](../../978/drafts/common_recovery_audit_v2.md)合并近期证据并保全部部门。978已给原生完整能源的移动双体界；停止运动优化，接[979](../../979/drafts/STATUS.md)审共同父描述、保留阶与真正剩余接口。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 978 verifier.')
