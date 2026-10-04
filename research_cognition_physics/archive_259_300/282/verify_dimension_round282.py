"""Read-only verification for operational dimension round 282."""
import argparse
import json
import verify_conservation_rounds as previous

previous.core.CONFIG[282] = ('operational_dimension_audit', 8)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = previous.verify(282, args.write_checks)
    if args.write_checks:
        path = previous.core.HERE/'research_round_282_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
