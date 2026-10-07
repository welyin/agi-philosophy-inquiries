"""Replay 799 and audit historical science, report sequence and live links."""
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
import constrained_cauchy_records as experiment
RECEIPT=HERE/'research_round_799_checks.json'
RESULT=HERE/'constrained_cauchy_records_results.json'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(RESULT.read_text(encoding='utf-8'))
    for n in range(776,799):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():
                assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_799.md'
    body=note.read_text(encoding='utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'800/drafts/STATUS.md',
          STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'): continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m: numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=799)==list(range(1,800))
    frozen=[STAGE/'798/research_round_798_checks.json', STAGE/'research_note_798.md',
            STAGE/'research_note_796.md',STAGE/'research_note_795.md',
            STAGE/'research_note_794.md',STAGE/'research_note_785.md',
            STAGE/'research_note_783.md',HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'constrained_cauchy_records.py',RESULT,
           HERE/'drafts/research_note_799_working.md',Path(__file__),
           STAGE/'800/drafts/STATUS.md']
    for p in fresh:
        if p.suffix=='.py': ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    return dict(round=799,date='2026-10-05',all_checks_passed=True,
                fresh_test_groups=2,cumulative_numbered_test_groups_from_798=3573,
                saved_result_reproduced=True,historical_manifest_evidence=history,
                previous_776_through_798_frozen_hashes_verified=True,
                historical_science_rerun=False,formal_reports=799,
                local_links_checked=links,display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
                argument_scope='The original constrained mixed free physical Wick algebra admits a Cauchy-neighborhood representation preserving the original kernel and mean. Applying it coefficientwise to the existing interacting physical representatives transports the original finite sources and corrected records, with all interaction coefficients retained, to one incoming physical Cauchy representation. This is not a theorem of faithful time-slice isomorphism for arbitrary local interacting BV cohomology, nor an autonomous apparatus construction.',
                original_differential_slice_and_source_dual_preserved=True,
                mixed_Wick_wavefront_and_graded_products_in_scope=True,
                same_W_and_compatible_mean_retained=True,
                original_interacting_coefficients_retained_in_incoming_representation=True,
                original_interacting_early_local_algebra_time_slice_proven=False,
                local_BV_cohomology_faithfulness_proven=False,
                finite_fixed_record_menu_dynamically_complete=False,
                native_continuous_Hamiltonian_simulated=False,
                autonomous_preparation_or_readout_proven=False,
                finite_coupling_probability_proven=False,
                original_graph_to_continuum_map_proven=False,
                finite_probes_are_continuum_proof=False,
                original_model_input_spacetime_group_action_retained=True,
                original_all_physics_unification_complete=False,
                independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write: assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    report=verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):
            assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in
                     ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))

