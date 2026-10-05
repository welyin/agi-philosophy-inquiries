"""Replay the 800 joint diagnostic and audit frozen evidence and live links."""
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
import physical_time_slice_transport as experiment
RECEIPT=HERE/'research_round_800_checks.json'
RESULT=HERE/'physical_time_slice_transport_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    actual=experiment.run()
    assert actual==json.loads(RESULT.read_text(encoding='utf-8'))
    assert actual['finite_causal_BV_probe']==json.loads(
        (HERE/'causal_cancellation_probe_results.json').read_text(encoding='utf-8'))
    for n in range(776,800):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_800.md'
    body=note.read_text(encoding='utf-8')
    assert body.count('$$')==24
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    assert not any(ord(c)<32 and c not in '\r\n' for c in body)
    docs=[note,STAGE/'801/drafts/STATUS.md',
          STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=800)==list(range(1,801))
    frozen=[STAGE/'799/research_round_799_checks.json',STAGE/'research_note_799.md',
            STAGE/'research_note_798.md',STAGE/'research_note_796.md',
            STAGE/'research_note_795.md',STAGE/'research_note_794.md',
            STAGE/'research_note_793.md',STAGE/'research_note_783.md',HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'physical_time_slice_transport.py',RESULT,
           HERE/'causal_cancellation_probe.py',HERE/'causal_cancellation_probe_results.json',
           HERE/'drafts/research_note_800_working.md',Path(__file__),STAGE/'801/drafts/STATUS.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    return dict(round=800,date='2026-10-05',all_checks_passed=True,
        fresh_test_groups=1,cumulative_numbered_test_groups_from_799=3574,
        saved_result_reproduced=True,historical_manifest_evidence=history,
        previous_776_through_799_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=800,
        local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='In the existing physical representation of the original real compact QME branch, a real two-sided fixed zero-strip completion gives a closed common past unitary. Causal cancellation and the constrained mixed free time-slice map cover the declared original source/record family in the explicitly enlarged algebra of admissible supported relative actions and Wick-distribution closure. The state is transported with the same map. This is represented target coverage, not faithful abstract local BV cohomology or a physically implemented cancellation apparatus.',
        real_relative_zero_strip_QME_completion_proven_in_existing_local_class=True,
        supported_relative_action_and_Wick_closure_explicit=True,
        common_closed_unitary_and_same_state_transport=True,
        original_sources_records_and_ordered_products_retained=True,
        arbitrary_abstract_local_BV_faithfulness_proven=False,
        original_finite_fixed_readout_menu_complete=False,
        cancellation_action_physically_implemented=False,
        finite_coupling_probability_proven=False,
        autonomous_preparation_or_readout_proven=False,
        original_graph_to_continuum_map_proven=False,
        finite_probes_are_continuum_proof=False,
        original_model_input_spacetime_group_action_retained=True,
        original_all_physics_unification_complete=False,
        independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    report=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):
            assert saved[key]==report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in
                     ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
