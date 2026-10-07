"""917 delivery: frozen history, original joint drivers, scope and links."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_917_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,917):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    stress=read(HERE/'receiver_stress_derivative_results.json');joint=read(HERE/'joint_constraint_driving_results.json')
    for row in (stress,joint):
        for rel,digest in row['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
        assert not row['full_goal_completed'] and not row['actual_finite_observable_error_certified']
    assert stress['original_receiver_and_full_offshell_stress_retained'] and stress['independent_original_stress_error']<1e-12
    assert joint['same_original_background_response_receiver'] and joint['original_receiver_divergence_reproduced']
    assert joint['full_background_contacts_including_volume_retained'] and not joint['full_support_constraint_driver_certified']
    for row in joint['independent_differentiated_residual_checks']:
        assert row['diffeomorphism_identity_error']<2e-10 and row['internal_identity_error']<2e-10
        assert row['wrong_drop_background_contact_error']>9e-7
        assert row['wrong_drop_receiver_divergence_error']>8e-6
        assert row['wrong_drop_internal_volume_contact_error']>8e-8
    note=STAGE/'research_note_917.md';text=note.read_text('utf-8')
    assert text.count('$$')==10 and re.findall(r'\\tag\{(\d+)\}',text)==[str(i) for i in range(1,6)]
    docs=[note,STAGE/'918/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=917)==list(range(1,918))
    oldfiles=[STAGE/'916/research_round_916_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/working_checks.json',HERE/'drafts/receiver_stress_working.md']+[STAGE/f'research_note_{n}.md' for n in (765,785,905,916)]
    newfiles=[note,Path(__file__),STAGE/'918/drafts/STATUS.md']+[HERE/n for n in ('receiver_stress_derivative.py','receiver_stress_derivative_results.json','joint_constraint_driving.py','joint_constraint_driving_results.json')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    d=dict(round=917,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=917,
      cumulative_numbered_test_groups_from_916=3702,historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      joint_source_background_contacts_computed=True,internal_volume_contact_retained=True,
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
