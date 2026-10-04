"""Verify independent complete rounds 353--356 and their frozen baselines."""
import argparse
import json
import verify_geometry_generation_rounds as previous

CONFIG={353:('collective_geometry_limit_audit',17),
        354:('collective_poisson_limit_audit',14),
        355:('acoustic_geometry_dynamics_audit',16),
        356:('two_sector_mean_field_audit',12)}
BASES={353:350,354:352,355:352,356:352}
core=previous.core
core.CONFIG.update(CONFIG)


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[353,354,355,356]
    result['execution_mode']='independent complete rounds in parallel'
    result['scientific_base_through_round']=BASES[number]
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
