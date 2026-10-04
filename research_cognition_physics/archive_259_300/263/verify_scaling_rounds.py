"""Read-only science, historical hash and text checks for rounds 261 through 263."""
import argparse
import json
import verify_interaction_rounds as core

CONFIG = {261: ('fisher_hanoi_scaling_audit', 9),
          262: ('hanoi_transport_spectrum_audit', 10),
          263: ('fisher_update_schedule_audit', 8)}
core.CONFIG.update(CONFIG)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = core.verify(args.round, args.write_checks)
    if args.write_checks:
        target = core.HERE/f'research_round_{args.round}_checks.json'
        if target.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
