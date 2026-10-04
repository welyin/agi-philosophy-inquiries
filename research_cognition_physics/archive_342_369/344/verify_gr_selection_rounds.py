"""Verify complete independent rounds 342--344 against their frozen baseline."""
import argparse
import json
import verify_parallel_physics_rounds as previous

CONFIG = {342: ('cognitive_dynamics_underdetermination_audit', 15),
          343: ('local_quantum_geometry_selection_audit', 14),
          344: ('covariant_gravity_uniqueness_audit', 13)}
core = previous.core
core.CONFIG.update(CONFIG)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result['date'] = '2026-09-23'
    result['parallel_batch'] = [342, 343, 344]
    result['scientific_base_through_round'] = 341
    result['batch_scientific_dependencies'] = []
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        path = core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
