"""Verify independent dynamic-source, spin-two, and gauge-region rounds 357--360."""
import argparse
import json
import verify_classical_geometry_rounds as previous

CONFIG={357:('dynamic_two_sector_limit_audit',13),
        358:('spin2_self_coupling_audit',15),
        359:('emergent_spin2_obstruction_audit',15),
        360:('gauge_region_composition_audit',13)}
core=previous.core
core.CONFIG.update(CONFIG)
previous.BASES.update({number:356 for number in CONFIG})


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[357,358,359,360]
    result['execution_mode']='independent complete rounds in parallel'
    result['scientific_base_through_round']=356
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
