"""910 freeze the actual original constrained direction and bounded-scope dual bridge."""
from pathlib import Path
import sys,importlib.util,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(STAGE/'909'))
spec=importlib.util.spec_from_file_location('prior909_verifier',STAGE/'909/verify_round909.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
import validate_physical_direction as experiment
TARGET=HERE/'research_round_910_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False,reproduce=False):
    previous=json.loads((STAGE/'909/research_round_909_checks.json').read_text('utf-8'))
    frozen={}
    def incorporate(table):
        for name,digest in table.items():
            assert name not in frozen or frozen[name]==digest,name
            frozen[name]=digest
    for n in range(776,910):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):incorporate(old[key])
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json'):
        incorporate(json.loads((STAGE/name).read_text('utf-8'))['files'])
    audit=json.loads((STAGE/'909/drafts/scope_reaudit_checks.json').read_text('utf-8'))
    incorporate(audit['preserved_files']);incorporate(audit['evidence_and_audit_hashes'])
    for name,digest in frozen.items():assert sha(ROOT/name)==digest,name
    history=prior.Layout().verify()
    saved=json.loads(experiment.TARGET.read_text('utf-8'));current=experiment.run(reproduce)
    for k in current:
        if k!='one_actual_tangent_source_pairing_reproduced':assert current[k]==saved[k],k
    assert saved['one_actual_tangent_source_pairing_reproduced'] is True
    assert saved['cumulative_numbered_groups']==previous['cumulative_numbered_test_groups_from_908']+1==3695
    note=STAGE/'research_note_910.md';body=note.read_text('utf-8')
    assert body.count('$$')==8 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,5)]
    docs=[note,STAGE/'911/drafts/STATUS.md',HERE/'drafts/dual_read_working.md']
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for p in docs:
        text=re.sub(r'\$\$.*?\$\$','',p.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(p.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(p,link)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=910)==list(range(1,911))
    oldfiles=[STAGE/'909/research_round_909_checks.json',HERE/'drafts/STATUS.md']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (765,785,863,869,870,872,899,904,905,908,909)]
    oldfiles += [STAGE/n for n in ('904/coupled_boson_tangent.py','908/relational_loop_source.py','908/material_background.py','908/physical_family_source_checks_results.json')]
    newfiles=[note,STAGE/'911/drafts/STATUS.md',Path(__file__),experiment.TARGET]
    newfiles += [HERE/n for n in ('physical_direction_dual_bridge.py','physical_direction_dual_bridge_results.json','validate_physical_direction.py','drafts/dual_read_working.md')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    return dict(round=910,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_909=3695,one_actual_case_executed_and_reproduced_in_round=True,
        actual_case_rerun_by_this_verifier=reproduce,all_N16_N24_cases_previously_executed=True,
        all_frozen_history_and_saved_prior_evidence_verified=True,
        historical_unique_files_verified=len(frozen),
        historical_manifest_evidence=history,
        local_links_checked=links,formal_reports=910,display_equations=4,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
        argument_scope=saved['argument_scope'],continuous_actual_response_error_certified=False,
        future_test_kernel_array_computed=False,complete872_quantum_feedback=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--reproduce',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    v=verify(a.write,a.reproduce)
    if a.write:TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(TARGET.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert v[k]==old[k]
    print(json.dumps({k:v for k,v in v.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
