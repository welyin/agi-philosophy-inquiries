"""Verify three independent complete rounds sharing frozen results through 338."""
import argparse
import json
import verify_disformal_stability_rounds as previous

CONFIG={339:('multispecies_characteristics_audit',14),
        340:('variable_disformal_coupling_audit',11),
        341:('disformal_validity_window_audit',12)}
core=previous.core
core.CONFIG.update(CONFIG)


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['parallel_batch']=[339,340,341]
    result['scientific_base_through_round']=338
    result['batch_scientific_dependencies']=[]
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.round,args.write_checks)
    if args.write_checks:
        path=core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
