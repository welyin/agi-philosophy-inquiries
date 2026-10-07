"""923 delivery: cognitive reporting/memory contract, scoped counterexample and common realization."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_923_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,923):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'direction_memory_contract_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['all_scientific_checks_passed'] and not result['sharp_report_cp']
    assert result['cp_range_for_nonnegative_visibility']==[0.,1/3]
    assert result['sharp_normalized_reference_witness_minimum']==-.25
    assert result['legal_boundary_report_kraus_rank']==7
    assert result['finite_task_algebra_dimensions']==dict(instantaneous_closed_direction_algebra=8,mass_mixed_future_algebra=16)
    assert result['initial_report_equal'] and result['initial_classical_chirality_summary_equal']
    assert abs(result['future_report_trace_distance']-1/3)<1e-14
    assert result['original_future_probability_prediction_minimax_error_lower_bound']==.5
    for k in ('kraus_completeness_error','kraus_channel_identity_error','isometry_error','same_Hamiltonian_intertwining_error','paired_symbol_square_identity_error'):assert result[k]<1e-14
    for k in ('Dirac_dynamics_derived_from_cognition','apparatus_autonomous_implementation_proved','energy_cost_of_preparing_dilation_computed','whole_stage_completed','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_923.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==22 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,12)]
    newdocs=[note,STAGE/'924/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in newdocs+nav:
        text=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=923)==list(range(1,924))
    assert '001—923轮共923份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—923的693份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:assert '923：方向报告不能取代内部状态' in q.read_text('utf-8-sig')
    oldfiles=[STAGE/'922/research_round_922_checks.json',HERE/'drafts/STATUS.md',STAGE/'_shared/notes/cognitive_generation_pivot_20261006.md']
    newfiles=newdocs+[Path(__file__),HERE/'direction_memory_contract.py',HERE/'direction_memory_contract_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=923,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=923,
      cumulative_numbered_test_groups_from_922=3708,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      explicit_cognitive_contract_checked=True,scoped_counterexample_and_same_process_positive_realization=True,
      known_universal_NOT_result_not_claimed_new=True,report_fidelity_not_optimized=True,
      legacy_candidate_not_prerequisite=True,autonomous_interaction_generation_completed=False,
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
