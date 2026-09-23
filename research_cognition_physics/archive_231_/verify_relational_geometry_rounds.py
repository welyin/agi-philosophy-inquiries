"""Verify independent complete rounds 361--364 against frozen evidence."""
import argparse
import json
import verify_gravity_interface_rounds as previous

CONFIG={361:('relational_collective_limit_audit',15),
        362:('boundary_charge_encoding_audit',18),
        363:('emergent_boundary_dynamics_audit',15),
        364:('local_clock_constraint_audit',12)}
BASES={361:357,362:360,363:360,364:360}
core=previous.core
core.CONFIG.update(CONFIG)
# Extend the loaded older wrapper's metadata map without changing its source.
previous.previous.BASES.update(BASES)


def preserved_draft():
    name='round363_drafts/emergent_boundary_dynamics_audit_results_pre_square_check.json'
    sha='ebb4b74a32d69cf1ef9e62e4d574047df13aacde7eaf1bdfe1152058a135425a'
    assert core.digest(core.HERE/name)==sha
    before=core.read(core.HERE/name)
    current=core.read(core.HERE/'emergent_boundary_dynamics_audit_results.json')
    assert before['checks']=={'run':14,'failures':0,'errors':0}
    assert current['checks']=={'run':15,'failures':0,'errors':0}
    reduced=json.loads(json.dumps(current))
    for row in reduced['compression_counterexamples']:
        assert 'observable_square_compression_defect' in row
        del row['observable_square_compression_defect']
    reduced['checks']['run']=14
    assert reduced==before
    return {name:sha}


def verify(number,pending=False):
    result=previous.verify(number,pending)
    result['date']='2026-09-23'
    result['parallel_batch']=[361,362,363,364]
    result['execution_mode']='independent complete rounds in parallel'
    result['scientific_base_through_round']=BASES[number]
    result['additional_frozen_dependency_rounds']=[]
    result['batch_scientific_dependencies']=[]
    if number>=363:
        result['preserved_round363_draft_hashes']=preserved_draft()
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
