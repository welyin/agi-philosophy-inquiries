"""Reproduce875; verify original local graph scope and preserve all evidence."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import local_graph_process_bridge as experiment
RECEIPT=HERE/'research_round_875_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    previous=json.loads((STAGE/'874/clock_history_source_match_results.json').read_text('utf-8'))
    assert science['cumulative_numbered_groups']==previous['cumulative_numbered_groups']+1==3660
    assert science['fresh_numbered_groups']==1
    for n in range(776,875):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    work=json.loads((HERE/'drafts/clock_transport_working_checks.json').read_text('utf-8'))
    for name,digest in work['files'].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_875.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'876/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=875)==list(range(1,876))
    frozen=[STAGE/'874/research_round_874_checks.json',HERE/'drafts/STATUS.md',
        HERE/'drafts/clock_transport_working_checks.json']
    frozen += [STAGE.parent/f'archive_{a}/research_note_{n}.md' for a,n in
       [('585_628',623),('585_628',625),('629_652',637),('702_741',704),('702_741',705),('702_741',706)]]
    frozen += [STAGE/f'research_note_{n}.md' for n in (854,855,873,874)]
    frozen += [STAGE.parent/'archive_702_741/706/joint_local_history_source_limit.py']
    fresh=[note,STAGE/'876/drafts/STATUS.md',HERE/'local_graph_process_bridge.py',experiment.TARGET,Path(__file__)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=875,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_874=3660,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_874_frozen_hashes_verified=True,
        prior875_working_files_preserved=True,historical_science_rerun=False,
        formal_reports=875,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='On the original fixed finite graph and positive compact geometry parameter set, source-independent local product spectral cutoffs preserve original interaction supports, Gauss and chosen bounded invariant even CAR records. Correct form-core density gives the original dynamics, own Gibbs states and first scalar preparation/dynamical source responses. This is not spatial refinement, a strictly one-sided finite channel, or an identification with the continuous E model.',
        original_interaction_and_source_family_retained=True,
        product_preserves_support_and_Gauss=True,
        own_thermal_state_and_first_scalar_real_history_response_proved=True,
        second_real_time_source_limit_for_product_cutoff_proved=False,
        original_Q_continuous_E_bridge_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        previous_goal_turn_classification='progress')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists()
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_874'):
            assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
