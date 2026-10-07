"""Replay 777 and check the frozen predecessors without historical ZIPs."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import Layout
import joint_auxiliary_bv_reduction as experiment

RECEIPT = HERE/'research_round_777_checks.json'
RESULT = HERE/'joint_auxiliary_bv_reduction_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    result = experiment.run()
    assert result == json.loads(RESULT.read_text(encoding='utf-8'))
    previous_path = STAGE/'776/research_round_776_checks.json'
    previous = json.loads(previous_path.read_text(encoding='utf-8'))
    for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
        for name, digest in previous[section].items():
            assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_777.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 24
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 13)]
    documents = [note, STAGE/'778/drafts/STATUS.md',
                 STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    documents += [STAGE/name for name in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    documents += [STAGE.parent/name for name in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    checked = 0
    for document in documents:
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            path = (document.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path == RECEIPT.resolve()), (str(document),target)
            checked += 1
    frozen = [previous_path, STAGE/'research_note_776.md', HERE/'drafts/STATUS.md',
              STAGE/'776/joint_auxiliary_source_tower.py']
    new = [note, HERE/'joint_auxiliary_bv_reduction.py', RESULT,
           Path(__file__), STAGE/'778/drafts/STATUS.md']
    return dict(round=777,date='2026-10-05',fresh_test_groups=3,
                cumulative_numbered_test_groups_from_776=3516,
                all_checks_passed=True,saved_result_reproduced=True,
                historical_manifest_evidence=history,all_previous_776_receipt_hashes_verified=True,
                historical_science_rerun=False,local_links_checked=checked,display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in new},
                continuum_quantum_pushforward_proven=False,full_N1_N2_realization_proven=False,
                original_full_loop_anomaly_computed=False,actual_instrument_proven=False,
                independent_agent_review=False,visual_checks_performed=False,app_goal_changed=False,
                scope='Exact original auxiliary BV splitting, local contraction and same-state free algebra reduction; full time-ordered quantum pushforward remains open.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    report = verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files','scope'):
            assert saved[key] == report[key],key
    print(json.dumps({k:v for k,v in report.items() if k not in
                      ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
