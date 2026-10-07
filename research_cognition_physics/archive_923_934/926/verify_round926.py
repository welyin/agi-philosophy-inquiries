"""926 delivery: bounded selection audit, not completion of a preferred physical candidate."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_926_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,926):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'direction_update_selection_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['all_scientific_checks_passed']
    for row in result['finite_operator_classification']:
        m=row['positive_role_multiplicity'];n=row['negative_role_multiplicity']
        assert row['rotation_dimension']==(m+n)**2
        assert row['quadratic_separation_dimension']==2*m*n
        if m==n:
            assert row['rotation_and_role_exchange_dimension']==2*m*m
            assert row['with_quadratic_separation_and_role_exchange_dimension']==m*m
    counter=result['rotational_role_exchange_counterexample']
    assert counter['rotation_error']==counter['role_exchange_error']==0
    assert abs(counter['square_separation_defect']-.48)<1e-14
    assert abs(counter['rest_square']-.65)<1e-14
    assert counter['analytic_spectrum_error']<1e-14
    assert max(counter['mean_reference_cross_terms'])<1e-14
    assert result['minimal_balanced_role_case_square_error']<1e-14
    cases=result['mixing_gap_and_stabilizer_cases']
    assert [r['mass_nullity'] for r in cases]==[0,0,2,4]
    assert [r['internal_stabilizer_lie_dimension'] for r in cases]==[2,4,3,3]
    for row in cases:
        assert row['square_identity_error']<1e-14 and row['comparison_covariance_error']<1e-14
    var=result['distinguishability_additivity_counterexample']
    assert var['anticommutator_norm']==0 and abs(var['covariance']+.5)<1e-14
    assert abs(var['variance_A']-.5)<1e-14 and abs(var['variance_B']-.5)<1e-14
    assert abs(var['variance_sum'])<1e-14 and var['uncentred_second_moment_additivity_error']<1e-14
    assert result['all_pure_variance_additivity_solution_dimension']==1
    for key in ('square_separation_is_additional_input','variance_cost_not_identified_with_physical_energy_cost','classical_dirac_classification_not_claimed_new','full_quantum_unknown_state_contract_kept'):assert result[key]
    for key in ('all_matter_and_gauge_representations_selected','physical_mass_values_selected','Lorentz_covariance_derived_from_cognition','GR_generated','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_926.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==18 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,10)]
    audit=STAGE/'_shared/notes/candidate_priority_audit_20261006.md'
    assert '结束本支线的分类投入' in prose
    assert '值得研究' in audit.read_text('utf-8') and '整个纲领不成立' in audit.read_text('utf-8')
    newdocs=[note,audit]
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=926)==list(range(1,927))
    assert '001—926轮共926份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—926的696份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig')
        assert '926：先判定要求的选择力' in s and 'candidate_priority_audit_20261006.md' in s
    oldfiles=[STAGE/'925/research_round_925_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'direction_update_selection.py',HERE/'direction_update_selection_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'925/research_round_925_checks.json')
    out=dict(round=926,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=926,
      cumulative_numbered_test_groups_from_925=previous['cumulative_numbered_test_groups_from_924']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      rotation_and_role_exchange_not_mistaken_for_general_dispersion=True,
      quadratic_separation_marked_as_input=True,all_pure_state_variance_additivity_proved_analytically=True,
      numerical_rank_not_used_as_universal_proof=True,prior_Dirac_and_space_results_reused=True,
      classification_branch_stopped=True,candidate_value_review_before_technical_extension=True,
      candidate_failure_not_program_failure=True,uncut_continuum_not_required=True,
      physical_gauge_group_or_GR_generated=False,full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
