"""Verify independent complete rounds 368--369 against frozen evidence."""
import argparse
import json
import verify_admission_geometry_rounds as previous

CONFIG={368:('quasi_local_encoding_audit',16),
        369:('modular_inclusion_limit_audit',14)}
BASES={368:367,369:367}
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
    result['parallel_batch']=[368,369]
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
