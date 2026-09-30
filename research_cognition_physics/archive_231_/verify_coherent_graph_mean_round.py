"""Read-only reproduction and scope audit for round 504."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
import verify_round503_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_504_checks.json'
NAMES = ['coherent_graph_mean_obstruction.py', 'coherent_graph_mean_obstruction_results.json',
         'research_note_504.md']
DRAFT = 'round504_drafts/research_note_504.txt'
BEFORE = 'round504_drafts/research_note_504_before_quantifier_review.txt'
SCOPE_CHECKS = 'internal_implementation_scope_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_503'] == 820
    assert frozen['total_protected_evidence_hashes'] == 892
    scope = core.read(HERE/SCOPE_CHECKS)
    extra = dict(scope['additional_unnumbered_evidence_sha256'])
    for name, sha in extra.items():
        assert core.digest(HERE/name) == sha, name
    extra[SCOPE_CHECKS] = core.digest(HERE/SCOPE_CHECKS)
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes', 'preserved_additional_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFT).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('coherent_graph_mean_obstruction_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (504, 503)
    assert saved['scope']['fixed_time_uniform_size_probability_gap']
    assert saved['scope']['arbitrary_actual_graph_reference']
    for key in ['symmetric_source_preparation_derived', 'other_effective_geometries_excluded',
                'dimension_three_generated', 'full_GR_goal_completed', 'phase_closure_triggered']:
        assert not saved['scope'][key]
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 12
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-27', round=504, scientific_base_through_round=503,
        frozen_integration_base_through_round=503,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=820,
        previous_protected_evidence_hashes_verified=892, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in [DRAFT, BEFORE]},
        preserved_additional_hashes=extra, unnumbered_diagnostic_tests_not_added=2,
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent analytic, scope, final-code and draft review; read-only reproduction.',
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
