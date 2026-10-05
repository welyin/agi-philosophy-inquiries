"""Replay 776, verify archived evidence in place, and check current local links."""
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
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import Layout
import joint_auxiliary_source_tower as experiment

RECEIPT = HERE / 'research_round_776_checks.json'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(writing=False):
    actual = experiment.run()
    assert actual == json.loads(experiment.RESULT.read_text(encoding='utf-8'))
    history = Layout().verify()
    note = STAGE / 'research_note_776.md'
    text = note.read_text(encoding='utf-8')
    assert text.count('$$') == 22
    assert re.findall(r'\\tag\{(\d+)\}', text) == [str(i) for i in range(1, 12)]
    # New report and live navigation only; old publication receipts remain frozen.
    documents = [note, STAGE / '777/drafts/STATUS.md', STAGE / 'README.md',
                 STAGE / '阶段成果总览.md', STAGE / '跨阶段主题索引.md',
                 STAGE / '文件索引.md'] + [STAGE.parent / name for name in
                  ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]
    links = 0
    for document in documents:
        body = document.read_text(encoding='utf-8-sig')
        for target in re.findall(r'\]\(([^)]+)\)', body):
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            target = target.split('#')[0].strip('<>')
            path = (document.parent / target).resolve()
            if writing and path == RECEIPT.resolve():
                pass
            else:
                assert path.exists(), (str(document), target)
            links += 1
    frozen = [STAGE / f'research_note_{n}.md' for n in (768, 770, 771, 773, 774, 775)]
    frozen += [HERE / 'drafts/STATUS.md',
               STAGE / '775/joint_wick_source_anomaly_results.json',
               STAGE / '775/research_round_775_checks.json']
    new_files = [note, HERE / 'joint_auxiliary_source_tower.py', experiment.RESULT,
                 Path(__file__), STAGE / '777/drafts/STATUS.md']
    return dict(round=776, date='2026-10-05', fresh_test_groups=3,
                diagnostic_subchecks=1386+84+252+63,
                cumulative_numbered_test_groups_from_775_receipt=3513,
                all_fresh_checks_passed=True, saved_result_reproduced=True,
                historical_manifest_evidence=history,
                historical_science_rerun=False,
                display_equations=11, local_links_checked=links,
                frozen_inputs={str(p.relative_to(ROOT)): digest(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): digest(p) for p in new_files},
                full_N1_N2_realization_proven=False,
                original_interacting_QME_unconditionally_proven=False,
                original_continuum_loop_coefficients_computed=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
                visual_checks_performed=False, app_goal_changed=False,
                scope='Conditional all-order basic auxiliary-source lift, original single-insertion Wick transport, and regular-patch distribution-coefficient interface. Full common BV normalization remains open.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    report = verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        saved = json.loads(RECEIPT.read_text(encoding='utf-8'))
        # Later rounds legitimately change live navigation and link counts.
        for key in ('frozen_inputs', 'new_scientific_and_entry_files', 'scope'):
            assert saved[key] == report[key], key
    print(json.dumps({k: v for k, v in report.items() if k not in
                      ('frozen_inputs', 'new_scientific_and_entry_files')}, ensure_ascii=False, indent=2))
