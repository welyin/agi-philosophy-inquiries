"""909: freeze actual-source implementation and its explicitly failed observer acceptance."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import validate_weak_propagation as experiment
RECEIPT=HERE/'research_round_909_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False,reproduce=False):
    science=experiment.run(reproduce=reproduce)
    saved=json.loads(experiment.TARGET.read_text('utf-8'))
    for k in science:
        if k!='one_actual_response_reproduced':experiment.close(science[k],saved[k])
    assert saved['one_actual_response_reproduced'] is True
    assert saved['cumulative_numbered_groups']==3694
    frozen={}
    def incorporate(table):
        for path,digest in table.items():
            assert path not in frozen or frozen[path]==digest,path
            frozen[path]=digest
    for n in range(776,909):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files'):incorporate(old[k])
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json'):
        incorporate(json.loads((STAGE/name).read_text('utf-8'))['files'])
    audit=json.loads((HERE/'drafts/scope_reaudit_checks.json').read_text('utf-8'))
    incorporate(audit['preserved_files']);incorporate(audit['evidence_and_audit_hashes'])
    for name,digest in frozen.items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_909.md';body=note.read_text('utf-8')
    assert body.count('$$')==12 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,7)]
    docs=[note,STAGE/'910/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            p=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert p.exists() or (writing and p==RECEIPT.resolve()),(doc,target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=909)==list(range(1,910))
    oldfiles=[STAGE/'908/research_round_908_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/working_checks.json',HERE/'drafts/scope_reaudit_checks.json']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (785,860,863,869,870,872,899,902,904,905,906,907,908)]
    oldfiles += [STAGE/p for p in ('904/coupled_boson_tangent.py','905/canonical_source_bridge.py','907/material_reference_jets.py','908/material_background.py','908/relational_loop_source.py','908/magnetic_reference_source.py')]
    newfiles=[note,STAGE/'910/drafts/STATUS.md',Path(__file__),experiment.TARGET]
    newfiles += [HERE/p for p in ('weak_source_propagation.py','weak_propagation_checks.py','weak_propagation_checks_results.json','observer_chart_check.py','observer_chart_check_results.json','observer_quadrature_check.py','observer_quadrature_check_results.json','validate_weak_propagation.py','drafts/weak_propagation_initial_results.json')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    return dict(round=909,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_908=3694,actual_source_response_executed_and_reproduced_in_round=True,
        actual_case_rerun_by_this_verifier=reproduce,all_six_diagnostics_previously_executed=True,
        historical_manifest_evidence=history,historical_unique_files_verified=len(frozen),historical_science_rerun=False,
        formal_reports=909,local_links_checked=links,display_equations=6,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
        argument_scope=saved['argument_scope'],physical_observer_accuracy_accepted=False,
        physical_retarded_Green_certified=False,original872_quantum_feedback_computed=False,
        continuous_support_and_time_error_proved=False,mathematical_peer_review=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--reproduce',action='store_true');a=p.parse_args()
    if a.write:assert not RECEIPT.exists()
    v=verify(a.write,a.reproduce)
    if a.write:RECEIPT.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert v[k]==old[k]
    print(json.dumps({k:v for k,v in v.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
