"""918 delivery: same-action density contacts and immutable archive verification."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_918_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,918):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    r=read(HERE/'euler_density_contacts_results.json')
    for rel,digest in r['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert r['same_action_background_response_receiver'] and r['physical_error_not_removed_by_representation']
    assert max(r['direct_variation_density_conversion_error'].values())<1e-13
    assert r['volume_log_derivative_identity_error']<1e-14 and r['internal_volume_cancellation_error']<1e-18
    assert r['internal_density_reduction_error']<1e-14
    assert r['raw_tensor_internal_driver_times_volume_max']>1e-7
    assert r['corrected_density_internal_driver_max']<2e-12
    assert r['same_Gauss_variation_sample_max']>5e-5 and r['same_Gauss_variation_identity_error']<1e-17
    assert not r['small_internal_driver_resolved_by_coordinate_difference']
    assert not r['actual_finite_observable_error_certified'] and not r['full872_response_computed'] and not r['full_goal_completed']
    for row in r['independent_coordinate_checks']:
        assert row['diffeomorphism_density_identity_error']<2e-10
        assert row['internal_density_identity_error']<2e-10
    note=STAGE/'research_note_918.md';text=note.read_text('utf-8')
    assert text.count('$$')==14 and re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,8)]
    docs=[note,STAGE/'919/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=918)==list(range(1,919))
    oldfiles=[STAGE/'917/research_round_917_checks.json',HERE/'drafts/STATUS.md']+[STAGE/f'research_note_{n}.md' for n in (765,785,905,916,917)]
    newfiles=[note,Path(__file__),STAGE/'919/drafts/STATUS.md',HERE/'euler_density_contacts.py',HERE/'euler_density_contacts_results.json']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    d=dict(round=918,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=918,
      cumulative_numbered_test_groups_from_917=3703,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      action_density_conversion_verified=True,internal_volume_contact_absorbed=True,same_nonzero_Gauss_residual_preserved=True,
      physical_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
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
