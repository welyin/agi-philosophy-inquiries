"""Verify904 full mixed boson tangent and independent constrained-family trajectories; preserve historical evidence."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import validate_boson_tangent as experiment
RECEIPT=HERE/'research_round_904_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    prev=json.loads((STAGE/'903/curved_fermion_validation_results.json').read_text('utf-8'))
    assert science['cumulative_numbered_groups']==prev['cumulative_numbered_groups']+1==3689
    frozen={}
    for n in range(776,904):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():
                assert name not in frozen or frozen[name]==digest
                frozen[name]=digest
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json'):
        old=json.loads((STAGE/name).read_text('utf-8'))
        for path,digest in old['files'].items():
            assert path not in frozen or frozen[path]==digest
            frozen[path]=digest
    for name,digest in frozen.items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_904.md';body=note.read_text('utf-8')
    assert body.count('$$')==16 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,9)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula)
    docs=[note,STAGE/'905/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=904)==list(range(1,905))
    oldfiles=[STAGE/'903/research_round_903_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/preparation_and_relative_source_audit.md',HERE/'drafts/working_checks.json']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (765,785,859,860,869,872,899,901,902,903)]
    oldfiles += [STAGE/'902/common_background_evolution.py',STAGE/'902/joint_background_validation_results.json',STAGE/'903/curved_fermion_validation_results.json']
    oldfiles += [STAGE.parent/'archive_554_584/572/joint_gauss_einstein_initial_data.py']
    newfiles=[note,STAGE/'905/drafts/STATUS.md',HERE/'coupled_boson_tangent.py',HERE/'boson_tangent_checks.py',HERE/'boson_tangent_checks_results.json',HERE/'validate_boson_tangent.py',experiment.TARGET,Path(__file__)]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=904,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_903=3689,saved_validation_result_reproduced=True,
        full_three_grid_time_integration_rerun=False,independent_small_grid_family_trajectories_rerun=True,
        full_time_integration_previously_executed=True,
        historical_manifest_evidence=history,previous_776_through_903_frozen_hashes_verified=True,
        historical_unique_files_verified=len(frozen),prior_working_files_preserved=True,
        historical_science_rerun=False,formal_reports=904,local_links_checked=links,display_equations=8,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
        argument_scope=science['argument_scope'],complete_classical_mixed_tangent_implemented=True,
        physical_sourced_retarded_inverse_constructed=False,receiver_modes_evolved=False,
        quantum_source_budget_evaluated=False,finite_time_error_certified=False,mathematical_peer_review=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not RECEIPT.exists()
    result=verify(a.write)
    if a.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_903'):assert old[key]==result[key],key
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
