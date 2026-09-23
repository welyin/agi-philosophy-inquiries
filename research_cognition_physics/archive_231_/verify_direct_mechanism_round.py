"""Verify round 392 on the frozen round-391 baseline; no image checks."""
import argparse
import json
import verify_shared_sharp_record_round as previous

CONFIG = {392: ('direct_mechanism_process_audit', 10)}
BASES = {392: 391}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result.update(date='2026-09-23', parallel_batch=[392],
                  execution_mode='complete direct-mechanism process audit on frozen round-391 baseline',
                  scientific_base_through_round=391,
                  additional_frozen_dependency_rounds=[], batch_scientific_dependencies=[],
                  final_science_review='primary-agent analytic and executable checks; parallel reviewers stopped at usage limit')
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
