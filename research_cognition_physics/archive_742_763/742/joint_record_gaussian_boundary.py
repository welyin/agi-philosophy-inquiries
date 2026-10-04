"""742: certify the original-record obstruction, without assuming a detector.

The frozen entry contains the actual full-matter covariance and independent
Fock calculation. Count this once when promoting the complete analytic result.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round742_drafts'))
import gaussian_record_entry as entry
TARGET=HERE/'joint_record_gaussian_boundary_results.json'

def run():
    measured=entry.run()
    assert measured==json.loads(entry.TARGET.read_text('utf8'))
    p=measured['actual_record_occupation_probability']
    assert 0<p<1
    defect=4*p*(1-p)
    bound=defect/14
    assert abs(measured['fourth_moment_Wick_defect']-defect)<4e-12
    assert abs(measured['rigorous_lower_bound_to_any_even_Gaussian_state']-bound)<2e-15
    names=('research_note_633.md','research_note_634.md','research_note_716.md',
           'research_note_730.md','research_note_735.md','research_note_741.md',
           'round742_drafts/gaussian_record_entry.py','round742_drafts/gaussian_record_entry_results.json')
    return dict(round=742,tests_run=1,failures=0,errors=0,
        checks=['original_pure_reference_partner_and_nonGaussian_record_gap'],
        results=measured,
        analytic_constant_reference=dict(inherited_from_round633=True,p='1/2',
            fourth_moment_defect='1',any_even_Gaussian_trace_distance_lower_bound='1/14',
            same_covariance_Gaussian_trace_distance='1/2',new_continuum_simulation=False),
        numerical_bounds_are_formula_evaluations_not_interval_certificates=True,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope='On the original pure quasifree reference with 0<p<1, the nonselective occupation record produces a nonGaussian fourth moment. Any fixed even quadratic fermion dilation with an independent even Gaussian probe and only probe discard has output trace-distance at least 2p(1-p)/7. This excludes that implementation class, not primitive measurement FLO, quantized boson-fermion interactions, the first-order mean-source construction, or the unified objective.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
