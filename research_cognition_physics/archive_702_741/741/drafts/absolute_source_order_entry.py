"""741 entry: original same-past, same-record source and backreaction order.

Only a full-species finite-symbol calibration, never a renormalized continuum
stress computation. The history is actually propagated at every parameter.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_retarded_reference_response as response
import joint_relative_source_development as prior
import joint_local_source_normalization as local
TARGET=HERE/'absolute_source_order_entry_results.json'

def actual(end,gamma):
    P,dP,_=response.flow(end,gamma)
    P=P+response.family.record_change(P)
    dP=dP+response.family.record_change(dP)
    _,_,G,dG=response.matrices(end,gamma)
    _,phi,_,_,_=response.ref.collar(end,.08)
    weight=response.entry.pulse(end);phi=np.exp(gamma*weight)*phi
    N=[prior.assemble([n,n]) for n in local.basis()]
    J=local.jacobian(phi);dJ=local.jacobian_jet(phi,weight*phi)
    sx=np.array([prior.source(P,n) for n in N])
    dsx=np.array([prior.source(dP,n) for n in N])
    source=np.r_[prior.source(P,G),J.T@sx]
    tangent=np.r_[prior.source(dP,G)+prior.source(P,dG),J.T@dsx+dJ.T@sx]
    return source,tangent

def run():
    rows=[]
    for end in (-.07,0.):
        Q0,dQ=actual(end,0.)
        scales=[]
        for eps in (.02,.01,.005):
            Qp,_=actual(end,eps);Qm,_=actual(end,-eps)
            # Odd-in-background source difference times its loop prefactor:
            # [eps Q(+eps) - eps Q(-eps)]/(2 eps^2) -> DQ[history].
            normalized=(eps*Qp-eps*Qm)/(2*eps*eps)
            derivative_error=float(np.max(abs(normalized-dQ)))
            weighted_remainder=float(np.max(abs(eps*Qp-eps*Q0)))
            scales.append(dict(parameter=eps,normalized_second_order_history_coefficient=normalized.tolist(),
                derivative_error=derivative_error,first_order_frozen_source_residual=weighted_remainder,
                residual_divided_by_parameter_squared=weighted_remainder/eps**2))
        errors=[r['derivative_error'] for r in scales]
        ratios=[errors[i]/errors[i+1] for i in range(2)]
        assert all(3.7<r<4.3 for r in ratios),(end,errors,ratios)
        assert np.max(abs(dQ))>1e-4
        if end==0.:
            assert response.entry.pulse(end)==0
            assert np.max(abs(Q0-actual(end,.02)[0]))>1e-6
        rows.append(dict(time=end,source=Q0.tolist(),complete_history_derivative=dQ.tolist(),
            current_background_pulse=response.entry.pulse(end),rows=scales,
            independent_centered_error_ratios=ratios))
    deps=('research_note_731.md','research_note_732.md','research_note_734.md','research_note_735.md','research_note_740.md',
          'joint_retarded_reference_response.py','joint_local_source_normalization.py')
    return dict(entry_round=741,latest_completed_round=740,formal_test_count_unchanged=3419,
        original_dimension=128,original_reference_record_and_parameters_retained=True,checks=rows,
        source_is_complete_finite_symbol_not_renormalized_continuum=True,
        no_numeric_Einstein_solution_or_actual_lambda_one_accuracy_claim=True,
        new_connection='Background history response enters lambda Q[B0+lambda b] at second order, including nonzero same-record memory after the background returns.',
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
