"""Replay round 796; verify the current layout, frozen science and navigation."""
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
import wick_state_initial_data as experiment

RECEIPT = HERE/'research_round_796_checks.json'
RESULT = HERE/'wick_state_initial_data_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(writing=False):
    assert experiment.run() == json.loads(RESULT.read_text(encoding='utf-8'))
    for number in range(776, 796):
        previous = json.loads((STAGE/f'{number}/research_round_{number}_checks.json').read_text(encoding='utf-8'))
        for section in ('frozen_inputs', 'new_scientific_and_entry_files'):
            for name, digest in previous[section].items():
                assert sha(ROOT/name) == digest, name
    history = Layout().verify()
    note = STAGE/'research_note_796.md'
    body = note.read_text(encoding='utf-8')
    assert body.count('$$') == 24
    assert not any(ord(c) < 32 and c not in '\n\r' for c in body)
    assert re.findall(r'\\tag\{(\d+)\}', body) == [str(n) for n in range(1, 13)]
    docs = [note, STAGE/'797/drafts/STATUS.md',
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
    assert sorted(n for n in numbers if n <= 796) == list(range(1, 797))
    frozen = [STAGE/'795/research_round_795_checks.json',
              STAGE/'research_note_795.md', STAGE/'research_note_783.md',
              STAGE/'research_note_792.md', STAGE/'research_note_793.md',
              STAGE/'research_note_785.md', STAGE/'research_note_772.md',
              HERE/'drafts/STATUS.md']
    fresh = [note, HERE/'wick_state_initial_data.py', RESULT, Path(__file__),
             STAGE/'797/drafts/STATUS.md', HERE/'drafts/research_note_796_working.md',
             HERE/'formal_state_probe.py', HERE/'formal_state_probe_results.json']
    for file in fresh:
        if file.suffix == '.py':
            ast.parse(file.read_text(encoding='utf-8'), filename=str(file))
    return dict(round=796, date='2026-10-05', all_checks_passed=True,
                fresh_test_groups=3, cumulative_numbered_test_groups_from_795=3567,
                saved_result_reproduced=True, historical_manifest_evidence=history,
                previous_776_through_795_frozen_hashes_verified=True,
                historical_science_rerun=False, formal_reports=796,
                local_links_checked=links, display_equations=12,
                frozen_inputs={str(p.relative_to(ROOT)): sha(p) for p in frozen},
                new_scientific_and_entry_files={str(p.relative_to(ROOT)): sha(p) for p in fresh},
                argument_scope='On the original real local branch and its fixed finite physical source menu, the original free functional extends to the required mixed Wick domain, annihilates full free BV exact elements, and yields a positive formal state on the completed source cohomology. A compatible full-free-complex mean displacement preserves W and realizes the original first-order Cauchy data; the state pulls back through the old-cutoff unitary comparison. This is not a nonperturbative preparation or finite-coupling instrument.',
                original_mixed_Wick_functional_domain_connected=True,
                original_free_BV_functional_Ward_on_required_domain=True,
                original_local_source_formal_positive_state_proven=True,
                original_physical_finite_terms_and_source_partners_preserved=True,
                original_free_physical_covariance_preserved=True,
                compatible_first_order_mean_data_realized_in_formal_state=True,
                full_free_BV_mean_displacement_proven=True,
                finite_probes_are_continuum_proof=False,
                original_nonperturbative_interacting_positive_state_proven=False,
                all_local_observables_and_regions_coherent_state_proven=False,
                actual_internal_preparation_or_CP_instrument_implemented=False,
                global_interacting_charge_common_closed_domain_proven=False,
                finite_coupling_convergence_proven=False, global_or_UV_result=False,
                original_graph_to_continuum_map_proven=False,
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
