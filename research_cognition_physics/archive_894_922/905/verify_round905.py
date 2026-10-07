"""Verify905 common Euler source, full constraints and independent action checks."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import validate_source_bridge as experiment
RECEIPT=HERE/'research_round_905_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    prev=json.loads((STAGE/'904/boson_tangent_validation_results.json').read_text('utf-8'))
    assert science['cumulative_numbered_groups']==prev['cumulative_numbered_groups']+1==3690
    frozen={}
    for n in range(776,905):
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
    note=STAGE/'research_note_905.md';body=note.read_text('utf-8')
    assert body.count('$$')==16 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,9)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula)
    docs=[note,STAGE/'906/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=905)==list(range(1,906))
    oldfiles=[STAGE/'904/research_round_904_checks.json',HERE/'drafts/STATUS.md']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (765,785,860,869,872,899,902,904)]
    oldfiles += [STAGE.parent/'archive_702_741/research_note_732.md']
    oldfiles += [STAGE/'902/common_background_evolution.py',STAGE/'904/coupled_boson_tangent.py',STAGE/'904/boson_tangent_validation_results.json']
    newfiles=[note,STAGE/'906/drafts/STATUS.md',HERE/'canonical_source_bridge.py',HERE/'source_bridge_checks.py',HERE/'source_bridge_checks_results.json',HERE/'validate_source_bridge.py',experiment.TARGET,Path(__file__)]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=905,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_904=3690,saved_validation_result_reproduced=True,
        full_Ward_grid_time_integration_rerun=False,independent_action_and_constraint_checks_rerun=True,
        full_Ward_grid_checks_previously_executed=True,
        historical_manifest_evidence=history,previous_776_through_904_frozen_hashes_verified=True,
        historical_unique_files_verified=len(frozen),prior_working_files_preserved=True,
        historical_science_rerun=False,formal_reports=905,local_links_checked=links,display_equations=8,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
        argument_scope=science['argument_scope'],common_Euler_force_and_normal_constraints_implemented=True,
        physical_retarded_mapping_conditionally_established=True,
        physical_sourced_retarded_solution_numerically_computed=False,
        original872_quantum_source_computed=False,finite_time_error_certified=False,mathematical_peer_review=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not RECEIPT.exists()
    result=verify(a.write)
    if a.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_904'):assert old[key]==result[key],key
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
