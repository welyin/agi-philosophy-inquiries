"""932 delivery: covariant protection contract and bounded effective realization."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_932_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,932):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'covariant_protection_contract_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==932 and result['all_scientific_checks_passed']
    assert [c['n'] for c in result['cases']]==[2,3,4,6,4]
    for c in result['cases']:
        for k in ('isometry_error','full_U2_covariance_sample_error','total_generator_intertwining_error','total_occupation_error'):assert c[k]<2e-13,k
        assert abs(c['weighted_leakage_sum']-1)<2e-13
        for e in c['erasures']:
            for k in ('recovery_support_error','trace_preserving_error','choi_formula_error'):assert e[k]<2e-13,k
            assert abs(e['entanglement_infidelity']-e['probability'])<2e-13
            assert abs(e['entangled_witness_half_trace_distance']-e['probability'])<2e-13
            assert e['entanglement_infidelity']+2e-13>=e['general_recovery_infidelity_lower_bound']
    assert max(result['internal_erasure'].values())<1e-14
    for key in ('exact_erasure_plus_nontrivial_transversal_continuity_excluded_conditionally','approximate_accepted_contract_has_finite_dimension_for_each_nonzero_accuracy','generator_span_is_not_consumed_energy','mature_theorem_and_generalized_W_construction_not_claimed_original','active_permissions_and_single_share_erasure_are_extra_inputs','prepared_blanks_access_rules_and_encoding_decoding_permissions_are_inputs','erased_information_kept_inside_overall_system','branch_stops_without_code_optimization'):assert result[key]
    for key in ('natural_autonomous_controller_derived','active_permissions_forced_by_six_protocols','all_six_protocols_physically_implemented','spacetime_dimension_gauge_group_or_GR_derived','exact_strong_contract_failure_is_whole_program_failure','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_932.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for phrase in ('谱宽不是消耗能量','额外主动权限','成熟','不取平方根','停止协变码和控制器细化','整体目标保持，未结项'):assert phrase in prose,phrase
    working=HERE/'drafts/covariant_protection_working.md'
    assert '唯一检验与预期判据' in working.read_text('utf-8')
    newdocs=[note,working,HERE/'drafts/operation_selection_review.md',STAGE/'933/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至932）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=932)==list(range(1,933))
    assert '001—932轮共932份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—932的702份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '932：连续主动权限与量子保护' in s and '933/drafts/STATUS.md' in s
    oldfiles=[STAGE/'931/research_round_931_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'covariant_protection_contract.py',HERE/'covariant_protection_contract_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'931/research_round_931_checks.json')
    out=dict(round=932,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=932,
      cumulative_numbered_test_groups_from_931=previous['cumulative_numbered_test_groups_from_930']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      original_KL_autonomous_preparation_and_commutant_results_reused=True,
      mature_Eastin_Knill_and_generalized_W_results_attributed=True,
      exact_additional_contract_failure_separated_from_six_protocols=True,
      complete_channel_partial_trace_and_entangled_error_cross_checked=True,
      finite_accuracy_allowed_without_microscopic_continuum_or_minimum_scale=True,
      physical_controller_permissions_not_claimed_derived=True,
      actual_generator_bandwidth_and_consumed_energy_distinguished=True,
      error_correction_technical_branch_stopped=True,
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
