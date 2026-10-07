"""915 delivery audit: immutable evidence, actual error transport, local links."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_915_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def complex_values(d):return np.array(d['real'])+1j*np.array(d['imag'])
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,915):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    budget=read(HERE/'weighted_source_budget_results.json');result=read(HERE/'analytic_weighted_source_results.json');jet=read(HERE/'analytic_jet_checks.json')
    for d in (budget,result,jet):
        for rel,digest in d['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
        assert not d['full_goal_completed']
    for row in budget['rows']:
        assert all(v is None for v in row['unknown_actual_physical_input_errors'].values())
        assert row['independent_formula_error']<1e-20
        assert row['diagnostic_random_error']<=row['diagnostic_bound']
        assert abs(row['diagnostic_aligned_complex_error']-row['diagnostic_bound'])<1e-18
    assert budget['constant_weight_extra_coefficients_max']==0
    assert result['finite_difference_step_in_new_source'] is None
    assert result['old_FD_reference_reproduction_error']<1e-14
    change=complex_values(result['actual_source_change'])
    pred=complex_values(result['exact_linear_error_transport_prediction'])
    assert np.max(abs(change-pred))<1e-17
    assert np.all(abs(change)<=np.array(result['finite_array_componentwise_difference_envelope']))
    assert np.max(abs(change))>6e-11
    assert not result['floating_roundoff_interval_certified'] and not result['physical_A_error_certified']
    assert not result['actual_finite_observable_error_certified'] and not result['full872_response_computed']
    assert jet['all_checks_passed'] and max(jet['original_reference_relative_errors'].values())<1e-10
    note=STAGE/'research_note_915.md';audit=HERE/'drafts/effective_scope_decision_after_915.md'
    text=note.read_text('utf-8');assert text.count('$$')==6
    assert re.findall(r'\\tag\{(\d+)\}',text)==['1','2','3']
    docs=[note,audit,STAGE/'916/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=915)==list(range(1,916))
    oldfiles=[STAGE/'914/research_round_914_checks.json',HERE/'drafts/working_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/research_note_915_working.md']
    newfiles=[note,audit,Path(__file__),STAGE/'916/drafts/STATUS.md']+[HERE/n for n in ('weighted_source_budget.py','weighted_source_budget_results.json','analytic_material_jets.py','analytic_weighted_source.py','analytic_weighted_source_results.json','check_analytic_jets.py','analytic_jet_checks.json')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    d=dict(round=915,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_914=3700,formal_reports=915,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
        independent_derivative_checks=jet['original_reference_relative_errors'],
        original_source_error_transport_residual=result['error_transport_identity_residual'],
        finite_algorithmic_error_separated=True,full_physical_error_certified=False,
        controlled_effective_description_accepted_as_stage_criterion=True,
        microscopic_continuity_assumed=False,physical_minimum_scale_assumed=False,
        full872_response_computed=False,full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        old=read(TARGET)
        for key in ('frozen_inputs','new_scientific_and_entry_files'):assert old[key]==d[key],key
    return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
