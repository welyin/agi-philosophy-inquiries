"""Certify a response coefficient, keeping its scope separate from the total."""
from pathlib import Path
from fractions import Fraction as Q
import json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'round746_drafts'))
import higgs_readout_certificate as local
import remote_hop_coefficient as jets
TARGET=HERE/'remote_hop_sign_certificate_results.json'
def rational_rows(result):
    return {(v['r'],v['s']):Q(int(v['numerator']),int(v['denominator'])) for v in result['exact_coefficients']}
def run():
    initial=jets.run()
    assert initial==json.loads((HERE/'remote_hop_coefficient_results.json').read_text('utf8'))
    original=jets.HOP_MATRIX
    try:
        jets.HOP_MATRIX=(((40,0),(20,13)),((-17,11),(31,0)))
        jets.act_terms.cache_clear();changed=jets.run()
    finally:
        jets.HOP_MATRIX=original;jets.act_terms.cache_clear()
    ratio=Q(177,1000)
    assert rational_rows(changed)=={k:ratio*v for k,v in rational_rows(initial).items()}
    saved=local.coefficients
    try:
        local.coefficients=lambda:{(r//2,s//2):v for (r,s),v in rational_rows(initial).items()}
        certificate=local.run()
    finally:local.coefficients=saved
    return dict(unit_sterile_spin='identity2',original_complex_spin_HS_squared='177/500',
        exact_complex_spin_scale_relative_to_identity='177/1000',
        exact_spin_covariance_coefficient_check=True,
        unit_hop_squared_sixth_unnormalized_interval=certificate['unnormalized_fourth_contrast_interval'],
        quadrature_error_upper=certificate['quadrature_error_upper'],
        sine_error_upper=certificate['sine_error_upper'],tail_error_upper=certificate['tail_error_upper'],
        normalization='Divide by the same positive local Gaussian Z, and by w_B²; hbar=1.',
        full_original_scalar_edge_contribution_not_evaluated=True,
        cannot_infer_total_remote_signal_from_one_coefficient=True,all_checks_passed=True)
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
