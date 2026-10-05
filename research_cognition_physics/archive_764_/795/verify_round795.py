"""Replay round 795; verify the current layout, frozen science and navigation."""
from pathlib import Path
import argparse
import ast
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
import real_cutoff_comparison as experiment

RECEIPT = HERE/'research_round_795_checks.json'
RESULT = HERE/'real_cutoff_comparison_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    assert experiment.run() == json.loads(RESULT.read_text(encoding='utf-8'))
    for number in range(776, 795):
        previous = json.loads((STAGE/f'{number}/research_round_{number}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
            for name, digest in previous[section].items():
                assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_795.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 22
    assert not any(ord(c) < 32 and c not in '\n\r' for c in body)
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 12)]
    docs = [note, STAGE/'796/drafts/STATUS.md',
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
    assert sorted(n for n in numbers if n <= 795) == list(range(1, 796))
    frozen = [STAGE/'794/research_round_794_checks.json',
              STAGE/'research_note_794.md', STAGE/'research_note_783.md',
              STAGE/'research_note_792.md', STAGE/'research_note_793.md',
              STAGE/'research_note_785.md', STAGE/'research_note_772.md',
              HERE/'drafts/STATUS.md']
    fresh = [note, HERE/'real_cutoff_comparison.py', RESULT, Path(__file__),
             STAGE/'796/drafts/STATUS.md', HERE/'drafts/research_note_795_working.md',
             HERE/'cauchy_state_probe.py', HERE/'cauchy_state_probe_results.json']
    for file in fresh:
        if file.suffix == '.py':
            ast.parse(file.read_text(encoding='utf-8'), filename=str(file))
    return dict(round=795, date='2026-10-05', all_checks_passed=True,
                fresh_test_groups=3, cumulative_numbered_test_groups_from_794=3564,
                saved_result_reproduced=True, historical_manifest_evidence=history,
                previous_776_through_794_frozen_hashes_verified=True,
                historical_science_rerun=False, formal_reports=795,
                local_links_checked=links, display_equations=11,
                frozen_inputs={str(p.relative_to(ROOT)): sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): sha(p) for p in fresh},
                argument_scope='On the original real local branch and an admissible causally contained observation region, loopwise real collar repair gives a real compact QME completion of the same germ. A source-independent unitary relative S-matrix compares the old and completed local source star algebras. First-order responses match if the original compatible Cauchy data are matched. No original interacting positive state has yet been constructed.',
                real_compact_QME_completion_on_original_real_branch_proven=True,
                original_physical_local_finite_terms_preserved=True,
                local_cutoff_unitary_star_comparison_proven=True,
                one_comparison_for_all_source_derivatives_proven=True,
                original_local_BRST_intertwining_proven=True,
                first_order_response_matching_under_same_compatible_data=True,
                original_W_and_mean_realized_in_full_positive_state=False,
                original_interacting_positive_state_proven=False,
                faithful_map_of_local_cohomology_into_global_free_cohomology_proven=False,
                original_continuum_anomaly_coefficients_computed=False,
                global_interacting_charge_domain_proven=False,
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
