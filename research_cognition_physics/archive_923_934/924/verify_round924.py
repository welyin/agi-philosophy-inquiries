"""924 delivery: bounded candidate-value test; encoding, transport and composition."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_924_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,924):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'role_transport_selection_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['all_scientific_checks_passed']
    for row in result['classification']:
        r=row['r'];assert row['label_permutation_commutant']==2
        assert row['cyclic_only_commutant']==r and row['two_encoding_commutant']==1
        assert row['transport_real_parameters_before']==4*r*r and row['transport_real_parameters_after']==4
    counter=result['permutation_counterexample']
    assert max(abs(a-4/3) for a in counter['all_classical_label_velocities'])<1e-14
    assert abs(counter['coherent_velocities'][0]-counter['coherent_velocities'][1]-1)<1e-14
    assert counter['permutation_commutator_error']==0 and counter['passive_covariance_error']<1e-14
    pair=result['composite_counterexample']
    assert pair['collective_commutant_dimension']==2 and pair['independent_commutant_dimension']==1
    assert pair['collective_commutator_error']==pair['whole_swap_commutator_error']==0
    assert pair['independent_replacement_violation']>.9
    assert result['identical_symbol_different_full_symmetry_commutant_dimensions']==[4,2]
    finite=result['finite_scope']
    assert max(finite['unitary_distance'],finite['reference_entangled_trace_distance'])<=finite['all_effect_probability_bound']
    assert abs(finite['all_effect_probability_bound']-.0126)<1e-14
    for row in result['weyl_twirl']:
        assert row['twirl_identity_error']<2e-14 and row['actual_distance']<=row['analytic_distance_bound']+2e-14
    for k in ('physical_gauge_group_selected','gauge_dynamics_generated','mass_parameters_selected','common_geometry_for_all_matter_proved','whole_stage_completed','full_goal_completed'):assert not result[k]
    assert not finite['source_or_backreaction_budget_certified']
    note=STAGE/'research_note_924.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==22 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,12)]
    newdocs=[note,STAGE/'925/drafts/STATUS.md']
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
    assert sorted(n for n in nums if n<=924)==list(range(1,925))
    assert '001—924轮共924份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—924的694份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:assert '924：重编码选择共同主部' in q.read_text('utf-8-sig')
    oldfiles=[STAGE/'923/research_round_923_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/role_transport_classification_working.md']
    newfiles=newdocs+[Path(__file__),HERE/'role_transport_selection.py',HERE/'role_transport_selection_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=924,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=924,
      cumulative_numbered_test_groups_from_923=3709,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      candidate_value_before_technical_completion=True,conditional_encoding_transport_bridge=True,
      overstrong_joint_geometry_gauge_selection_disproved_in_declared_scope=True,
      mature_commutant_theorem_not_claimed_original=True,finite_scope_error_bound_not_backreaction_certificate=True,
      legacy_candidate_not_prerequisite=True,autonomous_interaction_generation_completed=False,
      full_goal_completed=False,visual_checks_performed=False,app_goal_changed_by_this_round=False,
      user_updated_app_goal_verified_active=True,
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
