"""Read-only science, certificate and scope audit for round 513."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
import verify_round512_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_513_checks.json'
NAMES = ['continuous_relational_confinement.py', 'continuous_relational_confinement_results.json',
         'research_note_513.md']
DRAFTS = ['round513_drafts/research_note_513.txt',
          'round513_drafts/continuous_relational_confinement_before_source_effect_check.txt',
          'round513_drafts/research_note_513_before_transport_contract_review.txt']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_512'] == 847
    assert frozen['total_protected_evidence_hashes'] == 947
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFTS[0]).read_bytes()
    assert core.digest(HERE/NAMES[0]) == '862ff5845f1baf43ba78690b5104b0bddf679ba17cfb668c38105079a5b2374d'
    assert core.digest(HERE/NAMES[2]) == '61a276b182df9f2e794b9df4c00af6a61a951b7f0f0018ec8b7c3a6b2e629d31'
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('continuous_relational_confinement_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (513, 512)
    for name, sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name) == sha, name
    true_keys = ['same_continuous_510_generator', 'finite_time_size_uniform_relative_sector_bound', 'all_unknown_memory_graph_and_reference_for_protection', 'nontrivial_actual_same_class_reception', 'conserved_label_source_boundary_proved']
    false_keys = ['original_Z_history_membership_permanently_frozen', 'membership_partition_dynamically_created', 'connected_macroscopic_spatial_regions_generated', 'infinite_time_protection_at_fixed_g', 'autonomous_readout_and_controller_generated', 'dimension_three_generated', 'full_GR_goal_completed', 'phase_closure_triggered']
    assert set(saved['scope']) == set(true_keys+false_keys)
    assert all(saved['scope'][k] is True for k in true_keys)
    assert all(saved['scope'][k] is False for k in false_keys)
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 16
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-28', round=513, scientific_base_through_round=512,
        frozen_integration_base_through_round=512,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=847,
        previous_protected_evidence_hashes_verified=947, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent exact-certificate, channel-scope, code and final 16-equation draft review; read-only reproduction passed.',
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
