"""922 delivery: background readout variation and preserved direction-change handoff."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_922_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,922):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'background_readout_variation_results.json');validation=read(HERE/'mixed_derivative_validation_results.json')
    for obj in (result,validation):
        for rel,digest in obj['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['N']==17 and result['projection_quadrature_M']==25
    for k in ('force_is_negative_background_kinematic_residual','all104_variables_propagated','same_original_preparation_and_readout','moving_material_path_source_reference_and_weight_evaluation_included','one_integrable_background_field_used','background_derivative_Hessian_readout_included'):assert result[k]
    assert result['force_slots']==['phi','A'] and not result['independent_local_jet_patch_used']
    assert len(result['deformed_readouts'])==4
    assert result['baseline915_reproduction_error']==0
    x=result['centered_background_derivatives'][-1]['weighted_source']['real'][2]
    assert 8.5e-6<x<8.6e-6
    assert -1.23e-4<result['partial_921_tangent_plus_922_background_weighted']['real'][2]<-1.22e-4
    for k in ('background_induced_D2F_propagation_included','background_induced_receiver_and_stress_changes_included','full_background_residual_included','initial_geometric_constraint_error_included','finite_operator_dictionary_error_certified','physical_GPi_inverse_implemented','full_observable_error_certified','full872_response_computed','full_goal_completed'):assert not result[k]
    assert validation['all_diagnostic_checks_passed']
    assert validation['extrapolated_order_difference_max']<2e-17
    assert validation['weighted_half_step_relative_difference']<1e-6
    assert validation['maximum_Ward_array_identity_error']<1e-16
    assert not validation['full_physical_error_certified']
    note=STAGE/'research_note_922.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    newdocs=[note,HERE/'drafts/background_derivative_proof.md',STAGE/'923/drafts/STATUS.md',STAGE/'923/drafts/legacy_coupled_response_STATUS_before_pivot.md',STAGE/'923/drafts/coupled_kinematic_response_working.md',STAGE/'_shared/notes/cognitive_generation_pivot_20261006.md']
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
    assert sorted(n for n in nums if n<=922)==list(range(1,923))
    assert '001—922轮共922份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—922的692份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:assert '当前主线：认知操作对物理结构的实际约束' in q.read_text('utf-8-sig')
    oldfiles=[STAGE/'921/research_round_921_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/joint_error_transport_working.md']
    newfiles=newdocs+[Path(__file__)]+[HERE/n for n in ('background_readout_variation.py','background_readout_variation_results.json','analyze_mixed_derivative.py','mixed_derivative_validation_results.json')]
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=922,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=922,
      cumulative_numbered_test_groups_from_921=3707,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      same_integrable_background_readout_derivative_computed=True,independent_finite_record_derivative_order_checked=True,
      joint_error_not_certified=True,user_direction_change_applied=True,legacy_unfinished_response_not_automatically_resumed=True,
      physical_GPi_inverse_implemented=False,full872_response_computed=False,full_goal_completed=False,
      visual_checks_performed=False,app_goal_changed=False,
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
