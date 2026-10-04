"""Read-only science and scope audit for round 505."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
import verify_round504_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_505_checks.json'
NAMES = ['collective_port_closure_audit.py', 'collective_port_closure_audit_results.json',
         'research_note_505.md']
DRAFT = 'round505_drafts/research_note_505.txt'
BEFORE = 'round505_drafts/research_note_505_before_typesetting_review.txt'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_504'] == 823
    assert frozen['total_protected_evidence_hashes'] == 903
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFT).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('collective_port_closure_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (505, 504)
    for key in ['common_graph_independent_subspaces_classified', 'all_graphs_retained_in_exact_contract',
                'I2_exception_retained', 'two_role_fixed_time_uniform_size_gap',
                'uniform_population_source_can_have_graph_reference_correlations']:
        assert saved['scope'][key]
    for key in ['collective_source_preparation_derived', 'all_approximate_codes_excluded',
                'dimension_three_generated', 'full_GR_goal_completed', 'phase_closure_triggered']:
        assert not saved['scope'][key]
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 14
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-27', round=505, scientific_base_through_round=504,
        frozen_integration_base_through_round=504,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=823,
        previous_protected_evidence_hashes_verified=903, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in [DRAFT, BEFORE]},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent analytic, scope, code and final draft review with read-only reproduction.',
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(round=result['round'], fresh_tests=result['fresh_tests'],
        all_reported_checks_passed=result['all_reported_checks_passed']), ensure_ascii=False))
