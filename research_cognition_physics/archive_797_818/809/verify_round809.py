"""Reproduce 809, inspect links and retain all previously frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import weyl_wick_extension as experiment
import bounded_record_probe as early
RECEIPT=HERE/'research_round_809_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (experiment,early):assert ex.run()==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,809):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_809.md';body=note.read_text('utf-8')
    assert body.count('$$')==20 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,11)]
    docs=[note,STAGE/'810/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=809)==list(range(1,810))
    frozen=[STAGE/'808/research_round_808_checks.json',STAGE/'research_note_808.md',
            STAGE/'research_note_796.md',STAGE/'research_note_798.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'810/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT),
           HERE/'drafts/research_note_809_working.md',HERE/'drafts/working_checks.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=809,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_808=3587,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_808_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=809,local_links_checked=links,display_equations=10,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='A finite smooth Weyl-multiplier extension of the original physical Wick representation preserves its Gaussian positive state, original polynomial sources and formal scattering action. The nonzero leading field/fermion correlation of one allowed 808 process becomes a nonzero four-outcome formal correlation with a bounded free bosonic readout. No finite-coupling bounded interacting operator or autonomous instrument is proved.',
        original_physical_Wick_Weyl_extension_proven=True,
        same_original_state_action_cutoff_and_sources=True,
        four_outcome_formal_nonzero_correlation_proven=True,
        finite_coupling_bounded_interacting_observables_proven=False,
        all_offshell_time_ordered_exponential_insertions_proven=False,
        original_continuous_signal_numerically_computed=False,
        autonomous_preparation_or_readout_proven=False,graph_to_continuum_map_proven=False,
        spacetime_group_action_remain_inputs=True,original_all_physics_unification_complete=False,
        visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    r=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==r[k],k
    print(json.dumps({k:v for k,v in r.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
