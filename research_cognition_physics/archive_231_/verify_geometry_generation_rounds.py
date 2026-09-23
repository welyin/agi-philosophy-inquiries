"""Verify independent complete rounds 349--352 on the frozen 348 baseline."""
import argparse
import json
import verify_quantum_walk_round as previous

CONFIG={349:('curved_quantum_walk_audit',14),
        350:('quantum_geometry_backreaction_audit',16),
        351:('hypersurface_constraint_selection_audit',13),
        352:('matter_constraint_closure_audit',16)}
core=previous.core
core.CONFIG.update(CONFIG)


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[349,350,351,352]
    result['execution_mode']='independent complete rounds in parallel'
    result['scientific_base_through_round']=348
    result['additional_frozen_dependency_rounds']=[]
    result['batch_scientific_dependencies']=[]
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.round,args.write_checks)
    if args.write_checks:
        target=core.HERE/f'research_round_{args.round}_checks.json'
        if target.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
