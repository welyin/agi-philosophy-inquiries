"""Verify independent rounds 373--374 against the frozen round-372 baseline."""
import argparse
import json
import verify_direction_geometry_rounds as previous

CONFIG={373:('propagation_direction_calibration_audit',15),
        374:('stable_two_band_node_audit',13)}
BASES={373:372,374:372}
core=previous.core
core.CONFIG.update(CONFIG)
loaded=previous
while loaded is not None:
    if hasattr(loaded,'BASES'):
        loaded.BASES.update(BASES)
    loaded=getattr(loaded,'previous',None)


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[373,374]
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
