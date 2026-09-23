"""Verify the finite redirection dimension round on frozen baseline 386."""
import argparse
import json
import verify_asymptotic_displacement_round as previous

CONFIG = {387: ('finite_redirection_dimension_audit', 9)}
BASES = {387: 386}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result['date'] = '2026-09-23'
    result['parallel_batch'] = [387]
    result['execution_mode'] = 'complete finite redirection dimension round on frozen round-386 baseline'
    result['scientific_base_through_round'] = BASES[number]
    result['additional_frozen_dependency_rounds'] = []
    result['batch_scientific_dependencies'] = []
    result.pop('independent_of_parallel_rounds', None)
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
