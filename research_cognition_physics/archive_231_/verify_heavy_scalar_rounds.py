"""Read-only scientific and text checks for rounds 328-330."""
import argparse
import json
import verify_species_closure_rounds as previous

CONFIG = {328: ('heavy_scalar_kernel_audit', 10),
          329: ('heavy_scalar_potential_audit', 9),
          330: ('heavy_scalar_dynamics_audit', 10)}
core = previous.core
core.CONFIG.update(CONFIG)


def verify(number, pending=False):
    return previous.verify(number, pending)


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
