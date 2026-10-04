"""Read-only science and scope audit for round 511."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
# The preserved audit chain now spans over 250 rounds of nested imports.
sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
import verify_round510_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_511_checks.json'
NAMES = ['coherent_member_relay.py', 'coherent_member_relay_results.json',
         'research_note_511.md']
DRAFTS = ['round511_drafts/research_note_511.txt',
          'round511_drafts/research_note_511_before_recorded_readout.txt',
          'round511_drafts/research_note_511_before_delivery_wait_review.txt']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_510'] == 841
    assert frozen['total_protected_evidence_hashes'] == 932
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFTS[0]).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('coherent_member_relay_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (8, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (511, 510)
    for key in ['coherent_unknown_membership_supported',
                'full_transcript_independent_of_membership_for_ideal_cat',
                'complete_old_information_preserved_in_joint_isometry',
                'original_h_continues', 'heralded_source_error_composed',
                'recorded_membership_readout_closed_with_declared_controls']:
        assert saved['scope'][key]
    for key in ['old_membership_reduced_state_unchanged',
                'physical_removal_by_unknown_member_list',
                'autonomous_classical_record_routing_derived',
                'shared_resource_and_anchor_eliminated',
                'independent_fixed_cat_readout_formula_inherited',
                'pure_original_swap_only_implementation',
                'multiple_complete_spatial_endpoints_generated', 'dimension_three_generated',
                'full_GR_goal_completed', 'phase_closure_triggered']:
        assert not saved['scope'][key]
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 16
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-27', round=511, scientific_base_through_round=510,
        frozen_integration_base_through_round=510,
        fresh_tests=dict(run=8, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=841,
        previous_protected_evidence_hashes_verified=932, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
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
