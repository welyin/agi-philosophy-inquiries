"""919 delivery: original canonical/continuous constraint dictionary and scope."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_919_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,919):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(HERE/'drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    weak=read(HERE/'weak_gauss_transport_results.json');rec=read(HERE/'canonical_gauss_commutator_results.json')
    for result in (weak,rec):
        for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
        assert not result['full872_response_computed'] and not result['full_goal_completed']
    assert weak['same_background_response_and_receiver'] and not weak['new_physical_preparation']
    assert weak['physical_interpolant_N']==17 and weak['quadrature_sizes']==[17,25]
    assert weak['finite_constraint_tests_not_complete_observable_basis']
    assert max(weak['adapter_direct_Fourier_checks'].values())<2e-10
    assert not weak['actual_finite_observable_error_certified'] and not weak['all_constraints_certified']
    fine=weak['transport'][1]
    assert max(abs(x) for x in fine['weak_balance_defect'])<2e-14
    assert abs(fine['drop_initial_defect'][0])>6e-9
    assert abs(fine['drop_spatial_terms_defect'][0])>5e-11
    assert abs(weak['quadrature_differences_not_error_bounds']['actual_initial'][0])>8e-9
    assert rec['same912_preparation_and_flow'] and rec['physical_initial_surface_time']==0
    first=rec['rows'][0]['derivative_commutator_decomposition']
    assert first['canonical_Gauss']['max']==0. and first['reconstructed_Gauss']['max']>5e-5
    for row in rec['rows']:
        d=row['derivative_commutator_decomposition']
        assert d['canonical_Gauss']['max']<5e-11
        assert d['reconstructed_Gauss']['max']>5e-5
        assert d['decomposition_error']['max']<2e-20
        assert d['algebraic_momentum_mismatch']['max']<3e-17
        assert row['nodal_dictionary_errors']['Euler_density_dictionary']<1e-17
        assert row['independent_canonical_constraint_error']<2e-19
        assert not row['physical_error_certified']
    note=STAGE/'research_note_919.md';text=note.read_text('utf-8')
    assert text.count('$$')==12 and re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,7)]
    docs=[note,STAGE/'920/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=919)==list(range(1,920))
    oldfiles=[STAGE/'918/research_round_918_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/effective_scope_after_918_checks.json',HERE/'drafts/effective_scope_decision_after_918.md',HERE/'verify_effective_scope_after_918.py']+[STAGE/f'research_note_{n}.md' for n in (901,912,916,917,918)]
    newfiles=[note,Path(__file__),STAGE/'920/drafts/STATUS.md']+[HERE/n for n in ('weak_gauss_transport.py','weak_gauss_transport_results.json','canonical_gauss_commutator.py','canonical_gauss_commutator_results.json')]
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    out=dict(round=919,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=919,
      cumulative_numbered_test_groups_from_918=3704,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      same_preparation_initial_linear_canonical_Gauss_zero=True,nonlinear_reconstruction_derivative_commutator_identified=True,
      weak_transport_retains_actual_boundary_values=True,weak_tests_not_complete_observable_basis=True,
      physical_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
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
