"""Verify882 science, immutable evidence, report numbering and live links."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import mass_reference_source_jets as experiment
RECEIPT=HERE/'research_round_882_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    previous=json.loads((STAGE/'881/local_mode_source_capacity_results.json').read_text('utf-8'))
    assert science['cumulative_numbered_groups']==previous['cumulative_numbered_groups']+1==3667
    assert science['fresh_numbered_groups']==1
    for n in range(776,882):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json'):
        old=json.loads((STAGE/name).read_text('utf-8'))
        for path,digest in old['files'].items():assert sha(ROOT/path)==digest,path
    history=Layout().verify()
    note=STAGE/'research_note_882.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'883/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=882)==list(range(1,883))
    frozen=[STAGE/'881/research_round_881_checks.json',HERE/'drafts/STATUS.md',
            HERE/'drafts/common_generator_reuse_audit.md',
            STAGE/'875/drafts/clock_transport_working_checks.json',
            STAGE/'878/drafts/working_checks.json']
    frozen += [STAGE/f'research_note_{n}.md' for n in (803,805,842,843,844,854,876,880,881)]
    fresh=[note,STAGE/'883/drafts/STATUS.md',HERE/'mass_reference_source_jets.py',experiment.TARGET,Path(__file__)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=882,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_881=3667,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_881_frozen_hashes_verified=True,
        prior_working_files_preserved=True,historical_science_rerun=False,
        formal_reports=882,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope=science['argument_scope'],
        original_full_mass_family_preserved=True,
        global_normal_family_under_identity_CAR_ruled_out=True,
        finite_CAR_preparation_and_actual_history_C2_proved=True,
        smooth_background_family_843_invalidated=False,
        all_reference_transports_excluded=False,
        full_mixed_renormalized_source_family_matched=False,
        full_interacting_Q_E_identified=False,
        curved_future_PDE_numerically_solved=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        previous_goal_turn_classification='progress')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists()
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_881'):
            assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
