"""Verify round 393 on the frozen round-392 baseline; no rendering."""
import argparse
import json
import verify_direct_mechanism_round as previous

CONFIG = {393: ('program_transport_index_audit', 10)}
BASES = {393: 392}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result.update(date='2026-09-23', parallel_batch=[393],
                  execution_mode='program transport index and compensating circuit on frozen round-392 baseline',
                  scientific_base_through_round=392,
                  additional_frozen_dependency_rounds=[], batch_scientific_dependencies=[],
                  final_science_review='primary-agent source, analytic, and executable review; no independent agent final review')
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
