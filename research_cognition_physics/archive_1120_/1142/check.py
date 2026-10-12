"""1142 formal replay: source hashes, combined exact budgets and diagnostics.
Only --write-results writes this round's result. No continuum simulation.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/'results.json'

def load(name):
    p=STAGE/'_shared'/name
    spec=importlib.util.spec_from_file_location(name.replace('.','_'),p)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def run():
    d=load('after1141_defect_check.py').run()
    p=load('after1141_BT_power_check.py').run()
    assert d==json.loads((STAGE/'_shared/after1141_defect_results.json').read_text(encoding='utf-8'))
    assert p==json.loads((STAGE/'_shared/after1141_BT_power_results.json').read_text(encoding='utf-8'))
    hashes=d['groups']['frozen_source_identity']['sha256'].copy()
    for name in ('_shared/after1141_defect_check.py','_shared/after1141_defect_results.json',
                 '_shared/round1142_admission.md','1142/proof.md'):
        hashes[name]=hashlib.sha256((STAGE/name).read_bytes()).hexdigest()
    return {'round':1142,'status':'PASS','groups':{
        'frozen_sources_and_scope':{'status':'PASS','sha256':hashes,
            'new_cognitive_axioms':0,'scientific_calibrations':0,
            'original_inputs_parameters_quantum_tails_retained':True,
            'full_local_stress_proved':False,'energy_instrument_built':False},
        'combined_exact_budgets':{'status':'PASS',
            'defect':d['groups']['exact_fraction_budgets'],
            'power':p['groups']['exact_fraction_budgets'],
            'corrected_intermediate_coefficients_in_formal_proof':[66,30,142],
            'qualifications_preserved_as_historical_records':True,
            'analytic_proof_not_replaced_by_arithmetic':True},
        'exact_constrained_matrix_sylvester':d['groups']['exact_constrained_matrix_sylvester'],
        'finite_charge_and_semigroup_diagnostics':d['groups']['finite_scalar_diagnostics']},
        'scope':{'fixed_lab_w_epsilon_and_T_5_over_4':True,
            'actual_source_subaccount_requires_analytic_proof_and_review':True,
            'separate_actual_relative_energy_flux_old_bounds_claimed':False,
            'full_material_local_balance_or_GR_claimed':False,'goal_complete':False}}

def main():
    a=argparse.ArgumentParser(); a.add_argument('--write-results',action='store_true')
    args=a.parse_args(); data=run()
    if args.write_results:
        RESULT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert data==json.loads(RESULT.read_text(encoding='utf-8')),'Frozen results differ'
    b=data['groups']['combined_exact_budgets']['defect']['values']
    print(json.dumps({'round':1142,'status':'PASS','mode':'write-results' if args.write_results else 'read-only',
        'groups':4,'frozen_sources':len(data['groups']['frozen_sources_and_scope']['sha256']),
        'defect_integral_upper':b['integrated']['decimal'],
        'actual_fluctuation_energy_upper':b['e_actual']['decimal'],
        'physical_simulation':False,'goal_complete':False},ensure_ascii=False))
if __name__=='__main__': main()
