"""Read-only verification of rounds 278-279 and their frozen predecessors."""
import argparse
import json
import verify_delivery_rounds as previous

CONFIG = {278: ('defect_capacity_audit',10),279: ('screened_defect_audit',9)}
previous.core.CONFIG.update(CONFIG)


def verify(number,pending=False):
    result = previous.verify(number,pending)
    references = previous.core.read(previous.core.HERE/'round276_277_integration_checks.json')
    for name,sha in references['new_reference_hashes'].items():
        assert previous.core.digest(previous.core.HERE/name)==sha,name
    result['frozen_literature_files_verified'] = len(references['new_reference_hashes'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.round,args.write_checks)
    if args.write_checks:
        path = previous.core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Check report already exists; use read-only mode.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
