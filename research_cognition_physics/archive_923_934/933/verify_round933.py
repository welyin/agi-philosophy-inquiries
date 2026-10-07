"""933 delivery: unknown-load link contract and full coherent-loop boundary."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_933_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,933):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'transport_holonomy_contract_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==933 and result['all_scientific_checks_passed']
    for link in result['isolated_edge_checks']:
        for row in link:assert row['equal_rate_formula_error']<2e-14
    witness=result['arrival_witness']
    assert abs(witness[0]['target_probability'])<2e-14
    assert abs(witness[1]['target_probability']-1)<2e-14
    assert [r['initial_shifted_H_moments_1_to_4'][-1] for r in witness]==[4.,8.]
    for r in witness:
        assert r['initial_full_energy_mean']==3 and r['initial_full_energy_variance']==2
        assert r['total_energy_conservation_error']<2e-13
    assert result['smallest_full_H_eigenvalue']>.99
    for key in ('full_H_unitarity_error','blocked_sector_squared_identity_error','open_sector_sin4_formula_error','continuity_error','node_frame_H_covariance_error','loop_conjugacy_error','node_frame_actual_probability_error','central_connection_factorization_error','central_loop_scalar_error','central_connection_arrival_effect_spread'):assert result[key]<2e-13,key
    for r in result['finite_precision_rows']:assert r['operator_distance_to_selected_scalar_phase']<=r['analytic_bound']+1e-14
    assert result['one_quadrature_can_hide_internal_information']==[0.,1.]
    for key in ('all_path_quantifier_and_visible_full_state_family_are_extra_contracts','edge_rate_and_internal_transport_constrained_together','noncentral_loop_not_claimed_nonAbelian_group_derivation','initial_mean_and_variance_equal_but_full_energy_distributions_differ','full_global_unknown_state_preservation_does_not_imply_receiver_only_preservation','candidate_detail_optimization_stopped'):assert result[key]
    for key in ('coordinate_dimension_metric_or_Gauss_dynamics_generated','graph_and_control_permissions_derived_from_six','whole_program_disproved','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_933.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('一般共同空间不要求这项加强','全部可见闭路','完整能量分布不同','手征投影跳跃','成熟','停止本图和控制器细化','整体目标保持，未结项'):assert phrase in prose,phrase
    working=HERE/'drafts/transport_holonomy_contract_working.md'
    assert '明确的额外模型与预算' in working.read_text('utf-8')
    newdocs=[note,working,HERE/'drafts/geometry_source_priority_review.md',STAGE/'934/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至933）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=933)==list(range(1,934))
    assert '001—933轮共933份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—933的703份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '933：单边载荷独立与完整端点过程' in s and '934/drafts/STATUS.md' in s
    oldfiles=[STAGE/'932/research_round_932_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'transport_holonomy_contract.py',HERE/'transport_holonomy_contract_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'932/research_round_932_checks.json')
    out=dict(round=933,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=933,
      cumulative_numbered_test_groups_from_932=previous['cumulative_numbered_test_groups_from_931']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      old_endpoint_control_geometry_and_internal_transport_results_reused=True,
      scalar_times_unitary_link_classification_and_full_loop_quantifier_distinguished=True,
      actual_same_H_zero_one_arrival_and_current_checked=True,
      complete_initial_energy_distribution_not_claimed_equal=True,
      mature_vector_bundle_connection_definitions_attributed=True,
      finite_two_quadrature_error_bound_not_claimed_global_physics_budget=True,
      endpoint_autonomy_not_imposed_as_definition_of_common_space=True,
      graph_and_controller_technical_branch_stopped=True,
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
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
