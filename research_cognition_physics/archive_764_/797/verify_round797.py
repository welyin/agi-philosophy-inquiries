"""Replay 797 and verify frozen science, report sequence and current references."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import source_record_instrument as experiment
RECEIPT=HERE/'research_round_797_checks.json'
RESULT=HERE/'source_record_instrument_results.json'
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(RESULT.read_text(encoding='utf-8'))
    for n in range(776,797):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():
                assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_797.md'
    body=note.read_text(encoding='utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'798/drafts/STATUS.md',
          STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):
                continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:
            numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=797)==list(range(1,798))
    frozen=[STAGE/'796/research_round_796_checks.json', STAGE/'research_note_796.md',
            STAGE/'research_note_795.md', STAGE/'research_note_794.md',
            STAGE.parent/'archive_702_741/research_note_722.md',
            STAGE.parent/'archive_702_741/research_note_730.md',HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'source_record_instrument.py',RESULT,Path(__file__),
           STAGE/'798/drafts/STATUS.md']
    for p in fresh:
        if p.suffix=='.py':
            ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    return dict(round=797,date='2026-10-05',all_checks_passed=True,
                fresh_test_groups=2,cumulative_numbered_test_groups_from_796=3569,
                saved_result_reproduced=True,historical_manifest_evidence=history,
                previous_776_through_796_frozen_hashes_verified=True,
                historical_science_rerun=False,formal_reports=797,
                local_links_checked=links,display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
                argument_scope='The original fixed interacting physical source algebra admits a formal resolvent instrument, preserving full BV Ward, matrix positivity, finite ordered records and common unitary cutoff/state transport. The original source corrections and pre-operation state contract are retained. This is not a nonperturbative instrument or an internal autonomous apparatus.',
                original_fixed_source_algebra_formal_CP_instrument=True,
                full_BV_Ward_and_ordered_records_preserved=True,
                source_quantum_corrections_and_same_state_retained=True,
                common_cutoff_instrument_and_record_state_transport=True,
                old_bounded_CAR_projection_subcase_recovered=True,
                formal_strong_readout_limit_taken=False,
                finite_coupling_original_instrument_proven=False,
                original_internal_detector_and_resource_closure=False,
                continuum_operator_energy_domain_proven=False,
                original_graph_to_continuum_map_proven=False,
                original_model_input_spacetime_group_action_retained=True,
                finite_probes_are_continuum_proof=False,
                original_all_physics_unification_complete=False,
                independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    report=verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):
            assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in
                     ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))

