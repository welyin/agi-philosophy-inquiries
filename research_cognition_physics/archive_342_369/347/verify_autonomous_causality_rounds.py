"""Verify independent autonomous-resource and causal-limit rounds 345--347."""
import argparse
import json
import verify_gr_selection_rounds as previous

CONFIG = {345: ('autonomous_program_resource_audit', 12),
          346: ('autonomous_history_clock_audit', 17),
          347: ('strict_causal_limit_audit', 14)}
core = previous.core
core.CONFIG.update(CONFIG)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result['parallel_batch'] = [345, 346, 347]
    result['scientific_base_through_round'] = 341
    result['additional_frozen_dependency_rounds'] = [343] if number == 347 else []
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
