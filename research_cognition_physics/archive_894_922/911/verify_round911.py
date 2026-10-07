"""911: verify full-source flavour reduction and preserve prior evidence."""
from pathlib import Path
import sys,json,hashlib,ast,re,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import receiver_flavour_source_reduction as science
TARGET=HERE/'research_round_911_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text('utf-8'))
def run(writing=False):
    frozen={}
    def incorporate(table):
        for name,digest in table.items():
            assert name not in frozen or frozen[name]==digest,name
            frozen[name]=digest
    for n in range(776,911):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):incorporate(old[k])
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json'):
        incorporate(read(STAGE/name)['files'])
    audit=read(STAGE/'909/drafts/scope_reaudit_checks.json')
    incorporate(audit['preserved_files']);incorporate(audit['evidence_and_audit_hashes'])
    # The911 checkpoint is a historical pre-completion snapshot; do not execute
    # its then-current assertion that formal911 did not yet exist.
    incorporate(read(HERE/'drafts/finite_scope_checkpoint_checks.json')['evidence_and_audit_hashes'])
    for name,digest in frozen.items():assert sha(ROOT/name)==digest,name
    layout=Layout().verify()
    saved=read(science.TARGET);current=science.run();assert current==saved
    assert saved['cumulative_numbered_groups']==read(STAGE/'910/physical_direction_validation_results.json')['cumulative_numbered_groups']+1==3696
    assert len(saved['exact_rational_cases'])==9
    assert all(r['full_recursion_equals_independent_Euler_solution'] and r['reduced_matches_full'] for r in saved['exact_rational_cases'])
    assert saved['negative_changed_old_flavour_to_identity']['reduced_matches_full'] is False
    note=STAGE/'research_note_911.md';body=note.read_text('utf-8')
    assert body.count('$$')==12 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,7)]
    docs=[note,STAGE/'912/drafts/STATUS.md']+list((HERE/'drafts').glob('*.md'))
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=911)==list(range(1,912))
    oldfiles=[STAGE/'910/research_round_910_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/finite_scope_checkpoint_checks.json']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (785,845,846,853,860,869,870,871,872,899,906,910)]
    oldfiles += [STAGE/n for n in ('846/constrained_six_leg_recursion.py','870/receiver_four_point_backreaction.py','906/receiver_density_pair_identity.py')]
    newfiles=[note,Path(__file__),science.TARGET,HERE/'receiver_flavour_source_reduction.py',STAGE/'912/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    result=dict(round=911,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_910=3696,formal_reports=911,display_equations=6,
        exact_rational_cases_reproduced=9,independent_full_Euler_comparison_passed=True,
        full_four_mode_CAR_identity_reproduced=True,changed_flavour_negative_control_passed=True,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
        argument_scope=saved['actual_theorem_scope'],
        original_continuous_total_source_numerically_evaluated=False,
        actual_finite_observable_error_certified=False,full872_numerical_feedback_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
    if not writing:
        prior=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert result[k]==prior[k],k
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run(a.write)
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
