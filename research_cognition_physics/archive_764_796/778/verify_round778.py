"""Replay 778; verify frozen evidence and current navigation without ZIPs."""
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
import joint_propagating_normalization as experiment

RECEIPT = HERE/'research_round_778_checks.json'
RESULT = HERE/'joint_propagating_normalization_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    assert experiment.run() == json.loads(RESULT.read_text(encoding='utf-8'))
    for previous_round in (776, 777):
        previous = json.loads((STAGE/f'{previous_round}/research_round_{previous_round}_checks.json')
                              .read_text(encoding='utf-8'))
        for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
            for name, digest in previous[section].items():
                assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_778.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 20
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 11)]
    docs = [note, STAGE/'779/drafts/STATUS.md',
            STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/name for name in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md')]
    docs += [STAGE.parent/name for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]
    link_count = 0
    for document in docs:
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            path = (document.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path == RECEIPT.resolve()), (str(document), target)
            link_count += 1
    numbers = []
    for path in STAGE.parent.rglob('research_note_*.md'):
        if path.parent.name.startswith('archive_'):
            match = re.fullmatch(r'research_note_(\d+).md', path.name)
            if match:
                numbers.append(int(match[1]))
    assert sorted(numbers) == list(range(1, 779)), sorted(numbers)[-5:]
    frozen = [STAGE/'777/research_round_777_checks.json', STAGE/'research_note_777.md',
              STAGE/'777/joint_auxiliary_bv_reduction.py', HERE/'drafts/STATUS.md']
    fresh = [note, HERE/'joint_propagating_normalization.py', RESULT, Path(__file__),
             STAGE/'779/drafts/STATUS.md']
    return dict(round=778, date='2026-10-05', all_checks_passed=True,
                fresh_test_groups=3, cumulative_numbered_test_groups_from_777=3519,
                saved_result_reproduced=True, historical_manifest_evidence=history,
                previous_776_777_frozen_hashes_verified=True, historical_science_rerun=False,
                formal_reports=778, local_links_checked=link_count, display_equations=10,
                frozen_inputs={str(p.relative_to(ROOT)): sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): sha(p) for p in fresh},
                finite_tests_prove_continuum_normalization=False,
                propagating_normalization_argument='Analytic object mapping in research_note_778.md; diagnostics only test signs, insertion identities and the contact boundary.',
                original_full_N1_N2_proven=False, frame_quantum_transfer_proven=False,
                original_loop_anomaly_computed=False, interacting_positive_state_proven=False,
                independent_agent_review=False, visual_checks_performed=False,
                app_goal_changed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    report = verify(args.write)
    if args.write:
        RECEIPT.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        saved = json.loads(RECEIPT.read_text(encoding='utf-8'))
        for key in ('frozen_inputs', 'new_scientific_and_entry_files', 'propagating_normalization_argument'):
            assert saved[key] == report[key], key
    print(json.dumps({k: v for k, v in report.items() if k not in
                     ('frozen_inputs', 'new_scientific_and_entry_files')}, ensure_ascii=False, indent=2))
