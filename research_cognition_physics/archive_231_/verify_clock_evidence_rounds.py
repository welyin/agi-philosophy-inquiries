"""Read-only verification of rounds 295-296 and prior frozen evidence."""
import argparse
import json
import verify_clock_structure_rounds as previous

CONFIG = {295:('finite_clock_witness_audit',8),296:('clock_window_certificate_audit',9)}
previous.core.CONFIG.update(CONFIG)
core = previous.core


def verify(number,pending=False):
    return previous.verify(number,pending)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.round,args.write_checks)
    if args.write_checks:
        path = core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

