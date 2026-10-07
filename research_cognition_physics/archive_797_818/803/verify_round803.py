"""Reproduce 803 and protect all frozen science and the next research entry."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import cauchy_vertex_transport as experiment
RECEIPT=HERE/'research_round_803_checks.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    result=experiment.run()
    assert result==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,803):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_803.md';body=note.read_text('utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'804/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=803)==list(range(1,804))
    frozen=[STAGE/'802/research_round_802_checks.json',STAGE/'research_note_802.md',
            STAGE/'research_note_801.md',STAGE/'research_note_796.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'804/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT),
           HERE/'drafts/research_note_803_working.md',HERE/'drafts/working_checks.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=803,date='2026-10-05',all_checks_passed=True,fresh_test_groups=2,
        cumulative_numbered_test_groups_from_802=3579,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_802_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=803,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='The full original BFF coefficient is mapped to the background derivative of the original Dirac evolution on a common Cauchy space. The varying positive current density, spin frame, original masses and connections, compact cutoff and endpoint terms are retained. The original leading mixed output reduces to an ordered pairing of two original Cauchy solutions with that complete derivative. The selected centered state has zero first-order outgoing record marginal. No actual continuous response value or nonzero bound is supplied by the finite diagnostics.',
        original_complete_Cauchy_mapping=True,leading_record_marginal_zero_in_selected_state=True,
        fixed_original_state_not_instantaneous_spectrum=True,
        continuum_record_response_evaluated=False,finite_probes_are_continuum_signal=False,
        autonomous_preparation_or_terminal_proven=False,finite_coupling_probability_proven=False,
        original_graph_to_continuum_map_proven=False,spacetime_group_and_action_remain_inputs=True,
        original_all_physics_unification_complete=False,independent_agent_review=False,
        visual_checks_performed=False,app_goal_changed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    report=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
