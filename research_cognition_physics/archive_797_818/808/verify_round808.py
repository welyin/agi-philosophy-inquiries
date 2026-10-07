"""Reproduce 808 and preserve all earlier frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import cutoff_response_contract as experiment
RECEIPT=HERE/'research_round_808_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,808):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_808.md';body=note.read_text('utf-8')
    assert body.count('$$')==14 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,8)]
    docs=[note,STAGE/'809/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text('utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=808)==list(range(1,809))
    frozen=[STAGE/'807/research_round_807_checks.json',STAGE/'research_note_807.md',
            STAGE/'research_note_795.md',STAGE/'research_note_801.md',
            STAGE/'research_note_803.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'809/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=808,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_807=3586,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_807_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=808,local_links_checked=links,display_equations=7,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Within the existing real compact-cutoff process family, a collar variation preserves the old action plateau while, for a suitable fixed physical boson/sterile output menu and fixed original free state, the full leading joint response has nonzero affine slope. Thus all but at most one parameter on that line give a nonzero leading coefficient. This changes the scattering process; common transport of a single process instead preserves its statistics.',
        continuous_nonzero_response_exists_in_allowed_cutoff_process_family=True,
        previously_fixed_cutoff_nonzero_proven=False,
        original_free_state_changed=False,scattering_process_changes_along_cutoff_family=True,
        all_old_interacting_statistics_preserved_claim=False,
        finite_calibration_is_original_PDE=False,
        original_continuous_response_value_computed=False,
        autonomous_preparation_or_finite_coupling_proven=False,
        graph_to_continuum_map_proven=False,spacetime_group_action_remain_inputs=True,
        original_all_physics_unification_complete=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    r=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==r[k],k
    print(json.dumps({k:v for k,v in r.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
