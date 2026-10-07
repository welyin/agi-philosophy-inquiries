"""921 delivery: original residual contribution and original weighted-source response."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_921_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,921):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'kinematic_defect_readout_results.json');validation=read(HERE/'kinematic_response_validation_results.json')
    for obj in (result,validation):
        for rel,digest in obj['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['N']==17 and result['projection_quadrature_M']==25
    for k in ('original_background_source_and_receiver_retained','force_is_negative_tangent_kinematic_residual','all104_dynamical_variables_propagated','integrated_field_first_jets_used','same_source_reference_terms_retained'):assert result[k]
    assert result['force_slots']==['phi','A'] and not result['pointwise_jet_replacement_used']
    assert result['collocation_vs_projection_difference']['phi']>3e-5
    assert max(abs(v) for v in result['original_closed_source_pairing']['total'])<2e-13
    assert abs(result['original_weighted_Ward_source']['real'][2])>1e-4
    assert result['weighted_Ward_identity_error']<1e-18
    assert result['original_closed_source_pairing']['jet_pairing_identity_residual']<1e-23
    assert max(abs(v) for v in result['finite_path_curvature_source_difference_not_error_bound'])<3e-16
    for k in ('background_error_Hessian_and_readout_variation_included','initial_geometric_constraint_error_included','remaining_momentum_and_receiver_residuals_included','continuum_projection_error_certified','physical_GPi_inverse_implemented','actual_full_observable_error_certified','full872_response_computed','full_goal_completed'):assert not result[k]
    assert validation['all_checks_passed']
    assert max(row['maximum_error'] for row in validation['independent_finite_amplitude_rows'])<2e-19
    assert validation['integrated_field_final_value_error']<1e-22
    assert validation['constant_weight_original_source_identity_error']<1e-24
    assert validation['weighted_original_source_reproduction_error']<1e-16
    assert validation['large_validation_amplitude_is_numerical_check_not_physical_preparation']
    assert not validation['full_observable_error_certified'] and not validation['physical_GP_inverse_certified']
    note=STAGE/'research_note_921.md';text=note.read_text('utf-8')
    assert text.count('$$')==12 and re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,7)]
    docs=[note,STAGE/'922/drafts/STATUS.md',STAGE/'922/drafts/joint_error_transport_working.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=921)==list(range(1,922))
    oldfiles=[STAGE/'920/research_round_920_checks.json',HERE/'drafts/STATUS.md']+[STAGE/f'research_note_{n}.md' for n in (899,909,914,915,919,920)]
    newfiles=[note,Path(__file__),STAGE/'922/drafts/STATUS.md',STAGE/'922/drafts/joint_error_transport_working.md']+[HERE/n for n in ('kinematic_defect_readout.py','kinematic_defect_readout_results.json','validate_kinematic_response.py','kinematic_response_validation_results.json')]
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=921,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=921,
      cumulative_numbered_test_groups_from_920=3706,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      original_kinematic_residual_transported_through_all104=True,original_closed_and_weighted_sources_used=True,
      independent_finite_amplitude_variational_check_passed=True,
      isolated_contribution_not_total_error=True,physical_observable_error_certified=False,
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
