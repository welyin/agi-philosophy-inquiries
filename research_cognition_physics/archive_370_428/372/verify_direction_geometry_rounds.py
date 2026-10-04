"""Verify independent rounds 370--372; prior results are referenced, not rerun."""
import argparse
import json
import verify_transport_modular_rounds as previous

CONFIG={370:('radar_port_dimension_audit',14),
        371:('direction_bit_geometry_audit',14),
        372:('qubit_lorentz_cone_audit',18)}
BASES={370:367,371:367,372:367}
core=previous.core
core.CONFIG.update(CONFIG)
loaded=previous
while loaded is not None:
    if hasattr(loaded,'BASES'):
        loaded.BASES.update(BASES)
    loaded=getattr(loaded,'previous',None)


def preserved_direction_draft():
    name='round371_drafts/direction_bit_geometry_audit_results_pre_dedup.json'
    sha='e3e51d050f2f269d5feee49d7c0c11390bf765dda03eee3eb46b233af3ed6f2b'
    assert core.digest(core.HERE/name)==sha
    before=core.read(core.HERE/name)
    assert before['checks']=={'run':14,'failures':0,'errors':0}
    assert 'natural_dynamics_not_selected' in before
    return {name:sha}


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[370,371,372]
    result['execution_mode']='independent complete rounds in parallel'
    result['scientific_base_through_round']=BASES[number]
    result['additional_frozen_dependency_rounds']=[]
    result['batch_scientific_dependencies']=[]
    result['preserved_round371_draft_hashes']=preserved_direction_draft()
    result['draft_checks_counted_as_new_science']=False
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
