"""Replay 794 and verify frozen evidence, current links and report numbering."""
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
import compact_ward_homotopy as experiment

RECEIPT = HERE/'research_round_794_checks.json'
RESULT = HERE/'compact_ward_homotopy_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    assert experiment.run() == json.loads(RESULT.read_text(encoding='utf-8'))
    for number in range(776, 794):
        previous = json.loads((STAGE/f'{number}/research_round_{number}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
            for name, digest in previous[section].items():
                assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_794.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 24
    assert not any(ord(c) < 32 and c not in '\n\r' for c in body)
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 13)]
    docs = [note, STAGE/'795/drafts/STATUS.md',
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
    assert sorted(n for n in numbers if n <= 794) == list(range(1, 795))
    frozen = [STAGE/'793/research_round_793_checks.json', STAGE/'research_note_788.md',
              STAGE/'research_note_793.md', STAGE/'research_note_775.md', STAGE/'research_note_783.md',
              STAGE/'research_note_769.md', STAGE/'research_note_771.md',
              STAGE/'785/local_physical_slice.py',
              STAGE/'research_note_773.md', STAGE/'research_note_772.md',
              STAGE/'research_note_765.md', STAGE/'research_note_767.md',
              STAGE/'777/joint_auxiliary_bv_reduction.py',
              HERE/'drafts/STATUS.md']
    fresh = [note, HERE/'compact_ward_homotopy.py', RESULT, Path(__file__),
             STAGE/'795/drafts/STATUS.md', HERE/'drafts/research_note_794_working.md',
             HERE/'source_transport_probe.py', HERE/'source_transport_probe_results.json']
    return dict(round=794, date='2026-10-05', all_checks_passed=True,
                fresh_test_groups=3, cumulative_numbered_test_groups_from_793=3561,
                saved_result_reproduced=True, historical_manifest_evidence=history,
                previous_776_through_793_frozen_hashes_verified=True,
                historical_science_rerun=False, formal_reports=794,
                local_links_checked=links, display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)): sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): sha(p) for p in fresh},
                argument_scope='A compact QME completion has the original local quantum action/source germ on the observation neighborhood. In this constructed completion, the source-relative S-matrices differ along the quantum flow by exact elements of the same free BV star algebra; ordered products and adjoints share the same cohomology classes. The old arbitrary cutoff, its involution comparison and original actual state remain to be matched.',
                compact_QME_completion_of_original_local_germ_proven=True,
                original_physical_finite_terms_and_source_germ_preserved=True,
                constructed_completion_relative_S_cohomology_identity_proven=True,
                constructed_completion_ordered_products_and_adjoint_classes_preserved=True,
                classical_gauge_fixing_interaction_removed_in_constructed_cohomology=True,
                original_free_quadratic_gauge_unchanged=True,
                original_arbitrary_cutoff_star_dictionary_proven=False,
                original_actual_mean_and_state_matched_across_cutoffs=False,
                original_interacting_positive_state_proven=False,
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
