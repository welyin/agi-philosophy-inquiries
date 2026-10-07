"""916 delivery verification; historical content preserved, no visual checks."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_916_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,916):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    paths=('covariant_joint_residual_results.json','covariant_residual_validation_results.json','internal_covariance_results.json')
    actual,check,internal=[read(HERE/n) for n in paths]
    for row in (actual,check,internal):
        for rel,digest in row['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
        assert not row['full_goal_completed']
    assert actual['samples']==512 and actual['original_background_response_receiver_retained']
    assert actual['sample_residual_is_not_uniform_bound'] and not actual['physical_solution_error_certified']
    noether=actual['offshell_noether_checks']
    assert noether['diffeomorphism_noether_error']<2e-11 and noether['internal_noether_error']<1e-11
    assert noether['uncancelled_matter_exchange']>1e-5 and noether['uncancelled_scalar_charge']>1e-4
    assert check['all_checks_passed'] and internal['all_checks_passed']
    for key in ('manufactured_equation_checks','canonical_dictionary_checks','independent_affine_diffeomorphism_covariance_errors'):
        assert max(check[key].values())<1e-11
    assert max(internal['constant_internal_gauge_Euler_covariance_errors'].values())<1e-11
    assert check['actual_extractor_background_residual_correction_max']['YM']>.03
    note=STAGE/'research_note_916.md';text=note.read_text('utf-8')
    assert text.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,7)]
    docs=[note,STAGE/'917/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=916)==list(range(1,917))
    oldfiles=[STAGE/'915/research_round_915_checks.json',HERE/'drafts/STATUS.md']+[STAGE/f'research_note_{n}.md' for n in (773,785,897,902,904,905,906,912,913,914,915)]
    newfiles=[note,Path(__file__),STAGE/'917/drafts/STATUS.md']+[HERE/n for n in ('covariant_joint_residual.py','covariant_joint_residual_results.json','validate_covariant_residual.py','covariant_residual_validation_results.json','check_internal_covariance.py','internal_covariance_results.json')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    d=dict(round=916,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_915=3701,formal_reports=916,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
        actual_covariant_joint_residual_computed=True,independent_canonical_stress_and_gauge_covariance_checked=True,
        background_offshell_projection_term_retained=True,noether_identity_errors={k:v for k,v in noether.items() if k.endswith('_error')},
        physical_solution_error_certified=False,full872_response_computed=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==d[k],k
    return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
