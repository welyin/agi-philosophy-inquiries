"""Read-only science, certificate and scope audit for round 512."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
import verify_round511_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_512_checks.json'
NAMES = ['target_identity_delivery.py', 'target_identity_delivery_results.json',
         'research_note_512.md']
DRAFTS = ['round512_drafts/research_note_512.txt',
          'round512_drafts/next_problem_checkpoint.md',
          'round512_drafts/target_delivery_scout.py',
          'round512_drafts/target_delivery_scout_results.json',
          'round512_drafts/target_delivery_aligned_scout_results.json',
          'round512_drafts/target_delivery_scout_before_role_alignment.txt']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_511'] == 844
    assert frozen['total_protected_evidence_hashes'] == 938
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFTS[0]).read_bytes()
    assert core.digest(HERE/NAMES[0]) == '037d6bd6c681e9c033c84ea6caaa727d75b6137ea23cad57bd3d6ef11e13a3ad'
    assert core.digest(HERE/NAMES[2]) == 'ad85e82813785812710abb0814a21a9c8649adcb30ae7204a512cf1169efb6bb'
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('target_identity_delivery_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (7, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (512, 511)
    for name, sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name) == sha, name
    true_keys = ['existing_identity_hamiltonian_reused',
        'exact_query_to_target_bright_space_certificate',
        'exact_integer_finite_arrival_budget',
        'arbitrary_unknown_graph_and_reference_covered',
        'payload_identity_factor_and_erasure_channel_bound']
    false_keys = ['passive_graph_free_evolution_preserved',
        'all_conditional_payload_marginals_unchanged', 'all_network_sizes_proved',
        'autonomous_detector_clock_capture_derived',
        'complete_511_parallel_record_routing_implemented',
        'dimension_three_generated', 'full_GR_goal_completed', 'phase_closure_triggered']
    assert set(saved['scope']) == set(true_keys+false_keys)
    assert all(saved['scope'][k] is True for k in true_keys)
    assert all(saved['scope'][k] is False for k in false_keys)
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 18
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-28', round=512, scientific_base_through_round=511,
        frozen_integration_base_through_round=511,
        fresh_tests=dict(run=7, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=844,
        previous_protected_evidence_hashes_verified=938, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent exact-certificate, channel-scope, code and final 18-equation draft review; read-only reproduction passed.',
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
