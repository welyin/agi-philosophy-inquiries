"""Replay 783 and verify frozen evidence, current links and report numbering."""
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
import all_order_local_repair as experiment

RECEIPT = HERE/'research_round_783_checks.json'
RESULT = HERE/'all_order_local_repair_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    assert experiment.run() == json.loads(RESULT.read_text(encoding='utf-8'))
    for number in range(776, 783):
        previous = json.loads((STAGE/f'{number}/research_round_{number}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
            for name, digest in previous[section].items():
                assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_783.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 18
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 10)]
    docs = [note, STAGE/'784/drafts/STATUS.md',
            STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/name for name in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md')]
    docs += [STAGE.parent/name for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]
    links = 0
    for document in docs:
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            path = (document.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path == RECEIPT.resolve()), (str(document), target)
            links += 1
    numbers = []
    for path in STAGE.parent.rglob('research_note_*.md'):
        if path.parent.name.startswith('archive_'):
            match = re.fullmatch(r'research_note_(\d+).md', path.name)
            if match:
                numbers.append(int(match[1]))
    assert sorted(n for n in numbers if n <= 783) == list(range(1, 784))
    frozen = [STAGE/'782/research_round_782_checks.json', STAGE/'research_note_782.md',
              STAGE/'research_note_773.md', STAGE/'research_note_772.md',
              STAGE/'777/joint_auxiliary_bv_reduction.py', HERE/'drafts/STATUS.md']
    fresh = [note, HERE/'all_order_local_repair.py', RESULT, Path(__file__),
             STAGE/'784/drafts/STATUS.md']
    return dict(round=783, date='2026-10-05', all_checks_passed=True,
                fresh_test_groups=2, cumulative_numbered_test_groups_from_782=3534,
                saved_result_reproduced=True, historical_manifest_evidence=history,
                previous_776_through_782_frozen_hashes_verified=True,
                historical_science_rerun=False, formal_reports=783,
                local_links_checked=links, display_equations=9,
                frozen_inputs={str(p.relative_to(ROOT)): sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): sha(p) for p in fresh},
                argument_scope='Local formal all-loop recursion in the regular reference patch, constant-switch region and enlarged smooth inverse-reference coefficient class; no global switching or state theorem.',
                local_formal_all_loop_master_identity=True,
                original_one_loop_counterterm_preserved=True,
                original_anomaly_coefficients_computed=False,
                interacting_charge_constructed=False, interacting_positive_state_proven=False,
                finite_coupling_convergence_proven=False, global_or_UV_result=False,
                independent_agent_review=False, visual_checks_performed=False, app_goal_changed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    report = verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        saved = json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs', 'new_scientific_and_entry_files', 'argument_scope'):
            assert saved[key] == report[key], key
    print(json.dumps({k: v for k, v in report.items() if k not in
                     ('frozen_inputs', 'new_scientific_and_entry_files')}, ensure_ascii=False, indent=2))
