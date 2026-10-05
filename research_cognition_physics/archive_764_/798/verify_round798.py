"""Replay 798; audit frozen science, report sequence and current links."""
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
import physical_record_matrix_units as experiment
import internal_register_probe as probe
RECEIPT=HERE/'research_round_798_checks.json'
RESULT=HERE/'physical_record_matrix_units_results.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(RESULT.read_text(encoding='utf-8'))
    assert probe.run()==json.loads((HERE/'internal_register_probe_results.json').read_text(encoding='utf-8'))
    for n in range(776,798):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():
                assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_798.md'
    body=note.read_text(encoding='utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'799/drafts/STATUS.md',
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
    assert sorted(n for n in numbers if n<=798)==list(range(1,799))
    frozen=[STAGE/'797/research_round_797_checks.json', STAGE/'research_note_797.md',
            STAGE/'research_note_796.md',STAGE/'research_note_795.md',
            STAGE/'research_note_793.md',STAGE/'research_note_785.md',
            STAGE/'research_note_783.md',
            STAGE.parent/'archive_742_763/research_note_743.md',
            STAGE.parent/'archive_742_763/research_note_746.md',
            HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'physical_record_matrix_units.py',RESULT,
           HERE/'internal_register_probe.py',HERE/'internal_register_probe_results.json',
           HERE/'drafts/research_note_798_working.md',Path(__file__),
           STAGE/'799/drafts/STATUS.md']
    for p in fresh:
        if p.suffix=='.py':
            ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    return dict(round=798,date='2026-10-05',all_checks_passed=True,
                fresh_test_groups=2,cumulative_numbered_test_groups_from_797=3571,
                saved_results_reproduced=True,historical_manifest_evidence=history,
                previous_776_through_797_frozen_hashes_verified=True,
                historical_science_rerun=False,formal_reports=798,
                local_links_checked=links,display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
                argument_scope='Two original physical fermion modes, after a compatible finite source-menu extension and normalized epsilon completion, generate a parity-even unital formal M2 record algebra. Projection and partial-isometry corrections preserve full Ward, conjugation, the old action/menu/state and common cutoff transport. This is not an assertion about uncompleted hbar polynomial projections, a finite-coupling detector, blank preparation, or native Hamiltonian implementation.',
                finite_source_extension_preserves_old_menu=True,
                normalized_epsilon_completion_explicit=True,
                original_uncompleted_polynomial_algebra_contains_projectors_claimed=False,
                parity_even_unital_formal_M2_record_algebra=True,
                original_state_and_source_free_action_retained=True,
                original_record_source_commutation_for_arbitrary_source=False,
                same_state_blank_preparation_proven=False,
                original_native_Hamiltonian_implements_797_instrument=False,
                finite_coupling_record_operator_proven=False,
                original_graph_to_continuum_map_proven=False,
                finite_probes_are_continuum_proof=False,
                original_model_input_spacetime_group_action_retained=True,
                original_all_physics_unification_complete=False,
                independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:
        assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    report=verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):
            assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in
                     ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))

