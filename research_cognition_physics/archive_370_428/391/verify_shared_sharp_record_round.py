"""Verify round 391 against the frozen round-390 baseline; no image checks."""
import argparse
import json
import verify_finite_time_trace_round as previous

CONFIG = {391: ('shared_sharp_record_audit', 9)}
BASES = {391: 390}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result.update(date='2026-09-23', parallel_batch=[391],
                  execution_mode='complete shared sharp record audit on frozen round-390 baseline',
                  scientific_base_through_round=390,
                  additional_frozen_dependency_rounds=[], batch_scientific_dependencies=[])
    result.pop('independent_of_parallel_rounds', None)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        with (core.HERE / f'research_round_{args.round}_checks.json').open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
