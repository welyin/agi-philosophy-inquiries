"""Publish the 977 joint-output source coarse representation."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'976/verify_round976.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"finite_internal_reference_results.json")',1)[0]
for a,b in [('range(776,976)','range(776,977)'),('Delivery checks for 976','Delivery checks for 977'),
            ('research_round_976_checks.json','research_round_977_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"coarse_source_record_results.json")
    assert result["round"]==977 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    p=result["parameters"];b=result["uniform_budget"];a=result["algebra"]
    assert p["bin_width"]==.01 and p["joint_internal_dimension"]==352
    assert p["minimal_finite_spectrum_distance_to_bin_edge"]>2e-4
    assert b["uniform_on_entire_interval"]
    assert b["maximum_joint_output_distance"]<.008663
    assert b["record_contrast_lower"]>.46862
    assert b["infinite_occupation_certificate"]["total_state_vector_error"]<2.711e-8
    assert a["source_derivative_error"]<.0025
    assert a["full_energy_commutes_with_cq_generator_exactly"]
    assert a["finite_spectral_polar_basis_defines_positive_candidate"]
    assert a["minimum_cq_generator"]>199.67
    for row in result["rows"]:
        assert row["same_energy_mean_and_distribution_under_pinching"]
        assert row["pinching_menu_uniform_distance_bound"]<1.335e-5
        for endpoint in row["endpoints"]:
            assert endpoint["actual_total_joint_distance"]<row["full_interval_joint_distance_bound"]
    assert min(result["cq_endpoint_record_contrasts"])>.48
    assert max(result["rank_one_energy_diagonal_endpoint_record_contrasts"])<3.04e-6
    for row in result["characteristic_functions"]:
        assert .704<row["coarse_visibility"]<.705
        assert row["characteristic_error"]<math.pi*.01/2
    scope=result["scope"]
    assert scope["classical_bin_and_quantum_fibre_adopted"]
    assert scope["arbitrary_sender_reference_on_declared_joint_output"]
    for key in ("pinching_claimed_as_physical_irreversibility","all_internal_quantum_information_preserved",
        "exact_infinite_H_spectral_projectors_approximated","repeat_measurement_or_future_control_covered",
        "full_source_replaced_by_single_deterministic_mean","dynamical_classical_metric_derived",
        "Einstein_Langevin_noise_kernel_matched","full_goal_completed"):
        assert scope[key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core977",HERE/"coarse_source_record.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_977.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==14
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,8)]
    for term in ("整体目标未完成","不是声称原宇宙实际执行了谱测量或退相干",
        "并没有声称已逼近无限H的不连续谱投影" if False else "没有声称已逼近无限H的不连续谱投影",
        "本轮没有声称恢复Einstein—Langevin方程","宽度与时间没有在看到结果后调整"):
        assert term in prose,term
    newfiles=[note,HERE/"coarse_source_record.py",HERE/"coarse_source_record_results.json",
        Path(__file__),HERE/"drafts/coarse_source_decision.md",
        HERE/"drafts/coarse_source_adoption.md",HERE/"drafts/publish977.py",
        STAGE/"978/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至976','截至977'),('n<=976','n<=977'),('list(range(1,977))','list(range(1,978))'),
    ('001—976轮共976份','001—977轮共977份'),('231—976的746份','231—977的747份'),
    ('976：有限内部参考与关系记录','977：粗来源标签与块内量子记录'),
    ('STAGE/"975/research_round_975_checks.json"','STAGE/"976/research_round_976_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=977,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=977,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_976=prev["cumulative_numbered_test_groups_from_975"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        common_record_and_path_output_bound_verified=True,
        source_derivative_error_and_energy_conservation_verified=True,
        classical_bin_quantum_fibre_distinguished_from_full_dephasing=True,
        dynamical_classical_metric_derived=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'977/verify_round977.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 977：粗来源标签与块内量子记录'
    previous='## 976：有限内部参考与关系记录'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[977报告]({pre}research_note_977.md)在976同一过程采用经典来源箱与箱内量子状态，'
        +'一次记录及任意路径末读的联合误差在旧全区间<.008663，记录对比>.468622。'
        +'同一来源误差≤.0025且完整能源守恒；删除全部块内相干会丢失本记录。'
        +f'[结果]({pre}977/coarse_source_record_results.json) · '
        +f'[核验]({pre}977/research_round_977_checks.json)。正式977／累计3762，整体目标未完成。\n\n'
        +f'[采用增量]({pre}977/drafts/coarse_source_adoption.md)只领取有限来源／记录字典，'
        +'不等于动态经典度规、应力噪声核或宇宙完成。'
        +f'停止粗化优化，接[978]({pre}978/drafts/STATUS.md)把974—977并回全部门共同验收，优先判断剩余跨部门接口。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—976轮共976份','001—977轮共977份')
    if p==paths[4]:s=s.replace('当前正式976／累计3761，976已结项','当前正式977／累计3762，977已结项')
    if p==paths[5]:
        s=s.replace('# 231—976轮阶段成果总览','# 231—977轮阶段成果总览',1)
        s=s.replace('231—976的746份','231—977的747份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至976）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至977）',1)
        rows={
        'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、625、758—761、854、898、955—977|977从976同一准备到来源箱／量子块，联合记录及路径输出全区间误差<.008663|一次末读菜单；不逼近无限H谱投影，不自动保任意未来内部控制|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—977|977来源摘要保箱内量子记录，区间对比>.468622；完全谱基退相位会丢记录|保持实际材料与参考；不等于人口转换、实体读后账或循环协议完成|',
        'C22':'|C22 引力来源与物质反作用|1、3、4、H2；940—947、953、971、974—977|977同一有效H给来源导数误差≤.0025与完整能源守恒，量子来源分布不换成单均值|差分静态路径报告；动态度规、局域应力响应及噪声核仍未恢复|',
        'C26':'|C26 宏观经典、流体等|1、3、4、6；350、354、758—761、973—977|977实际材料支持经典来源箱＋量子内层的受控表示；记录与路径统计共同保留|经典标签不是完整经典时空，未签收所有态、任意后续任务或流体恢复|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—977|977原生联合输出误差<.008663；纯能量基概率记录差约3.04e-6，而所采CQ表示仍>.468622|一次声明菜单的有限对照，不用一般经典模型或全纲领反证作解释|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[973共同实现审计](../../973/drafts/common_realization_audit_v1.md)保全部部门。[977采用增量](../../977/drafts/coarse_source_adoption.md)把来源摘要和量子记录接成同一有限字典；停止粗化优化，接[978](../../978/drafts/STATUS.md)并回全部门验收，确定真正剩余的共同接口。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 977 verifier.')
