"""745: reproduce exact native history coefficients and a strict integral sign."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round745_drafts'))
import exact_history_density as exact
import certify_history_integral as bound
TARGET=HERE/'joint_native_history_certificate_results.json'
def run():
    a=exact.run();assert a==json.loads((HERE/'round745_drafts/exact_history_density_results.json').read_text('utf8'))
    b=bound.run();assert b==json.loads((HERE/'round745_drafts/certify_history_integral_results.json').read_text('utf8'))
    deps=('research_note_623.md','research_note_744.md','joint_singlet_common_mass_rg_results.json',
          'joint_fermion_gauss_completion.py','round745_drafts/exact_history_density.py',
          'round745_drafts/exact_history_density_results.json','round745_drafts/certify_history_integral.py')
    return dict(round=745,tests_run=2,failures=0,errors=0,
        checks=['exact_original_onsite_H_CAR_density','outward_integral_with_analytic_remainders'],
        exact_coefficients=dict(zero_orders=a['exact_zero_density_orders'],
            maximum_floating_calibration_difference=a['maximum_difference_from_floating_density']),
        certificate=b,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='The original one-vertex Gaussian preparation with exact-decimal frozen model parameters has a strictly negative sixth parity-history contrast derivative. Therefore arbitrarily short nonzero wait choices with nonzero contrast exist, and at each such choice sufficiently small positive inter-read gaps preserve contrast. No explicit usable time window, connected-graph result, ideal occupation instrument, autonomous detector or gravity completion is claimed.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=745,tests_run=2,strict_negative=r['certificate']['strict_negative'])))
