"""581: conditional first-order reduction of the same scalar effective coefficient."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_effective_operator_reduction_results.json'


def load(name):
    spec=importlib.util.spec_from_file_location(name,HERE/'round581_drafts'/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def run():
    algebra=load('operator_reduction_entry_probe').run()
    actual=load('field_redefinition_probe').run()
    assert algebra==json.loads((HERE/'round581_drafts/operator_reduction_entry_results.json').read_text('utf8'))
    assert actual==json.loads((HERE/'round581_drafts/field_redefinition_results.json').read_text('utf8'))
    assert algebra['checks_passed']==3 and actual['checks_passed']==2
    evidence=dict(metric_factorization=algebra['metric_factorization'],scalar_chain=algebra['scalar_chain'],
        RS_block=algebra['RS_block'],
        actual_joint_redefinition=dict(initial_action=actual['initial_action'],first_action_variation=actual['first_action_variation'],
            rows=[{k:v for k,v in r.items() if 'length' not in k} for r in actual['rows']],
            successive_action_remainder_ratios=actual['successive_action_remainder_ratios'],
            potential_gradient_error=actual['potential_gradient_error']),
        length_pullback=dict(untransformed_length=actual['untransformed_length'],
            rows=[{k:v for k,v in r.items() if 'length' in k or k=='parameter'} for r in actual['rows']]))
    deps=('joint_curved_quantum_source.py','joint_scalar_effective_geometry.py','research_note_580.md',
          'research_round_580_checks.json','round581_drafts/operator_reduction_entry_probe.py',
          'round581_drafts/operator_reduction_entry_results.json','round581_drafts/field_redefinition_probe.py',
          'round581_drafts/field_redefinition_results.json')
    return dict(round=581,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='same original K and U, four-dimensional Euclidean Einstein-scalar leading action at kappa=1, scalar-loop local coefficient with zero gauge curvature; explicit first-order metric and scalar redefinition plus total derivative; controlled on bounded smooth profiles and finite perturbative coefficients, not a complete independent operator basis, full quantum-frame equivalence, full nonzero-gauge matter source or graph continuum')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
