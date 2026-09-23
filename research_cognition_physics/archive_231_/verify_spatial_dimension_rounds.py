"""Verify independent spatial-dimension rounds from frozen baseline 375."""
import argparse
import json
import verify_clifford_connection_round as previous

CONFIG = {
    376: ('spatial_direction_integrability_audit', 14),
    377: ('isotropic_direction_orbit_audit', 13),
}
BASES = {376: 375, 377: 375}
BATCH = [376, 377]
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def dimension_priority_hashes():
    name = 'three_dimensional_stage_priority_checks.json'
    saved = core.read(core.HERE / name)
    assert saved['all_reported_checks_passed']
    assert saved['scientific_baseline_through_round'] == 375
    assert saved['new_scientific_tests'] == 0
    hashes = dict(saved['new_decision_hashes'])
    for path, sha in hashes.items():
        assert core.digest(core.HERE / path) == sha, path
    hashes[name] = core.digest(core.HERE / name)
    return hashes


def verify(number, pending=False):
    priority = dimension_priority_hashes()
    target = core.HERE / f'research_round_{number}_checks.json'
    if target.exists():
        assert core.read(target)['dimension_stage_priority_hashes'] == priority
    result = previous.verify(number, pending)
    result['date'] = '2026-09-23'
    result['parallel_batch'] = BATCH
    result['execution_mode'] = 'independent complete spatial-dimension rounds in parallel'
    result['scientific_base_through_round'] = BASES[number]
    result['additional_frozen_dependency_rounds'] = []
    result['batch_scientific_dependencies'] = []
    result['dimension_stage_priority_hashes'] = priority
    result['priority_maintenance_counted_as_science'] = False
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        path = core.HERE / f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
