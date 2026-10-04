"""Integrate two independent rounds without rerunning older experiments."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round420_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round421_422_integration_checks.json'
CONFIG = {421: ('boundary_interaction_energy_audit', 8, 10),
          422: ('translation_clock_interface_audit', 9, 18)}


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    links = 0
    for number, (stem, tests, formulas) in CONFIG.items():
        checked = core.read(HERE / f'research_round_{number}_checks.json')
        assert checked['all_reported_checks_passed']
        assert checked['scientific_base_through_round'] == 420
        assert checked['fresh_tests'] == dict(run=tests, failures=0, errors=0)
        assert len(checked['new_file_hashes']) == 3
        assert len(checked['preserved_draft_hashes']) == (3 if number == 422 else 0)
        for key in ('new_file_hashes', 'preserved_draft_hashes'):
            for name, sha in checked[key].items():
                assert core.digest(HERE / name) == sha, name
        saved = core.read(HERE / (stem + '_results.json'))
        assert (saved['tests_run'], saved['failures'], saved['errors']) == (tests, 0, 0)
        for key in ('three_dimensional_space_derived', 'full_GR_goal_completed',
                    'full_cognitive_countermodel_completed', 'phase_closure_triggered'):
            assert not saved['scope'][key], key
        note = HERE / f'research_note_{number}.md'
        assert core.text_checks(note) == checked['text_checks']
        assert checked['text_checks']['display_formulas'] == formulas
        ast.parse((HERE / (stem + '.py')).read_text(encoding='utf-8'))
        for link in core.link_parser()(note.read_text(encoding='utf-8')):
            dest = (HERE / link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (number, link)
            links += 1
    for name in ('verify_boundary_interface_rounds.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_420') == 571
    assert result.pop('unchanged_prior_science_hashes_231_419') == 568
    assert result['stage_saved_tests'] == 1934
    assert result['total_protected_evidence_hashes'] == 601
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 422 and int(latest[2]) >= 1951
    result.update(date='2026-09-24', rounds=[421, 422],
        execution_mode='two independent complete rounds with separate derivations and implementations',
        scientific_base_through_round_by_round={421: 420, 422: 420},
        fresh_tests_by_round={421: 8, 422: 9}, fresh_tests=17,
        stage_saved_tests=1951, science_hashes_verified_231_422=577,
        unchanged_prior_science_hashes_231_420=571, total_protected_evidence_hashes=610,
        unchanged_prior_evidence_hashes=601, newly_preserved_draft_files=3,
        new_scientific_display_formulas_checked=sum(c[2] for c in CONFIG.values()),
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 2,
        interaction_energy_boundary_audited=True, normal_clock_outgoing_interface_audited=True,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
