"""Read-only audit of rounds 276 and 277, with protected proposal evidence."""
import argparse
import json
import verify_delivery_rounds as previous

CONFIG = {276: ('growing_stream_audit', 8), 277: ('stable_ingress_audit', 9)}
previous.core.CONFIG.update(CONFIG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = previous.verify(args.round,args.write_checks)
    if args.write_checks:
        path = previous.core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Check report exists; use read-only mode.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
