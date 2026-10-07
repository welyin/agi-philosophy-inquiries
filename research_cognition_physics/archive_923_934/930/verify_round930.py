"""930 delivery: conditional cone-response bridge and finite polarization/source comparison."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_930_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,930):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'cone_response_bridge_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['round']==930 and result['all_scientific_checks_passed']
    assert result['common_cone_symbol_samples']==32
    assert (result['null_symbol_rank'],result['off_cone_symbol_rank'])==(1,3)
    assert result['metric_dilaton_axion_coefficient_jacobian_rank']==11
    assert result['fixed_common_cone_remaining_coefficients']==2
    for key in ('symbol_identity_max_error','conformal_response_max_error','jacobian_smallest_singular_value','common_lagrangian_energy_error','common_lagrangian_poynting_error','same_power_pair_difference','same_force_pair_difference'):assert result[key]<1e-12,key
    rows=result['interfaces'];assert len(rows)==5
    for row in rows:
        for key in ('boundary_solve_error','flux_error','full_scattering_unitarity_error'):assert row[key]<1e-12,key
        assert abs(row['reflection_power']+row['transmission_power']-1)<1e-13
        assert abs(row['interface_momentum_per_incoming_energy']-2*row['reflection_power'])<1e-13
        assert abs(row['reconstructed_lambda_right']-row['lambda_right'])<1e-13
        assert abs(row['reconstructed_axion_jump']-row['axion_jump'])<1e-13
        assert abs(sum(row['outgoing_port_polarization_probabilities'])-1)<1e-13
    assert abs(rows[1]['reflection_power']-1/9)<1e-13
    assert abs(rows[2]['reflection_power']-1/9)<1e-13
    assert abs(result['same_power_pair_full_polarization_total_variation']-16/81)<1e-13
    assert result['transparency_budget']==.01
    assert abs(result['lambda_ratio_bounds'][0]-9/11)<1e-13
    assert abs(result['lambda_ratio_bounds'][1]-11/9)<1e-13
    assert abs(result['normalized_axion_jump_bound']-(1/99)**.5)<1e-13
    for key in ('transparent_exactly_iff_response_parameters_match_in_declared_class','full_response_readout_recovers_remaining_relative_parameters','nonbirefringent_classification_is_cited_not_numerically_proved','common_cone_matching_and_transparent_replacement_are_extra_hypotheses','no_independent_surface_layers_assumed'):assert result[key]
    for key in ('physical_volume_or_absolute_EM_normalization_generated','interface_carrier_dynamics_and_closed_backreaction_constructed','all_six_protocols_realized_in_EM_candidate','spacetime_dimension_Maxwell_kinematics_or_GR_generated','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_930.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==22 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,12)]
    for phrase in ('停止本界面模型的技术扩展','成熟分类','本轮没有构造它的完整动力学','16/81','分量界','不是已采用的六条本身','整体目标未完成'):
        assert phrase in prose,phrase
    working=HERE/'drafts/bridge_selection_working.md'
    assert '有限预算与停止' in working.read_text('utf-8')
    newdocs=[note,working,STAGE/'931/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至930）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    assert '当前认知假说补充：共同指认的六条协议' in nav[1].read_text('utf-8-sig')
    assert '已采用六条共同认知协议' in nav[2].read_text('utf-8-sig')
    links=0
    for doc in newdocs+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=930)==list(range(1,931))
    assert '001—930轮共930份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—930的700份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '930：共同光锥约束响应' in s and '931/drafts/STATUS.md' in s
    oldfiles=[STAGE/'929/research_round_929_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'cone_response_bridge.py',HERE/'cone_response_bridge_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'929/research_round_929_checks.json')
    out=dict(round=930,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=930,
      cumulative_numbered_test_groups_from_929=previous['cumulative_numbered_test_groups_from_928']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      mature_classification_cited_with_scope_not_claimed_new=True,
      principal_cone_full_polarization_and_source_compared_in_same_object=True,
      cone_matching_and_active_transparent_replacement_kept_extra=True,
      full_relative_response_readout_and_finite_error_bound_checked=True,
      same_power_and_total_force_do_not_fix_full_response_checked=True,
      conformal_volume_and_absolute_normalization_not_claimed_fixed=True,
      full_carrier_dynamics_and_cognitive_EM_generation_not_claimed=True,
      interface_technical_branch_stopped=True,
      universal_program_failure_not_claimed=True,
      full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
