"""749 entry: analytic positive second-order averaging defect in731 constraints.

Prescribed smooth relative energy profiles, not actual measurement branches.
Only the explicitly stated mean-psi/mean-canonical-fields map is tested.
"""
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_source_constraint_response as old
TARGET=HERE/'positive_constraint_defect_results.json'
def reaction(q,psi,A):
    return -8*old.geo.laplace(psi)+q['C']*psi**5-q['B']*psi-A*psi**-7-2*q['Y']*psi**-3
def run():
    d=old.setup();q=d['q'];p=d['psi'];A=np.sum(d['tensor']**2,axis=(-1,-2))+q['pKp']
    # Exact source choice makes the first response u=psi by731's identity.
    rho=2*q['C']+4*A*p**-12+4*q['Y']*p**-8
    S2=10*q['C']*p**5+68*A*p**-7+52*q['Y']*p**-3
    predicted=S2*p**-5
    assert predicted.min()>10*q['C'].min()>0
    L=lambda v:-8*old.geo.laplace(v)+(5*q['C']*p**4-q['B']+7*A*p**-8+6*q['Y']*p**-4)*v
    ucheck=old.maximum(L(p)-2*rho*p**5);assert ucheck<1e-8
    rows=[]
    for eps in (.004,.002,.001):
        vals=[];errs=[]
        for sign in (-1,1):
            z=dict(q);z['C']=q['C']-2*sign*eps*rho
            assert z['C'].min()>0
            pp,tt,_=old.geo.solve_hamiltonian(z,initial=p)
            vals.append(pp);errs.append(old.maximum(reaction(z,pp,A)*pp**-5))
            assert old.maximum(tt-d['tensor'])<1e-12
        mean=(vals[0]+vals[1])/2
        defect=reaction(q,mean,A)*mean**-5
        error=old.maximum(defect/eps**2-predicted)
        rows.append(dict(epsilon=eps,branch_equation_error=max(errs),
            mean_constraint_defect_min=float(defect.min()),mean_constraint_defect_max=float(defect.max()),
            scaled_defect_error=error,mean_psi_shift_min=float((mean-p).min())))
        assert defect.min()>1000*max(errs)
    assert rows[-1]['scaled_defect_error']<rows[0]['scaled_defect_error']/10
    return dict(entry_round=749,latest_completed_round=748,formal_tests_unchanged=3434,
        inherited_conditions='731 original smooth non-flat background, C>0,Y>0,A>=0 and fixed matter/CMC data.',
        source='rho_star = 2C+4A psi^-12+4Y psi^-8; energy-only relative profiles +/- epsilon rho_star.',
        analytic_identity='D2 coefficient after averaging psi is S2=10C psi^5+68A psi^-7+52Y psi^-3 > 0.',
        physical_constraint_leading_coefficient='psi^-5 S2 >= 10 C > 0 pointwise.',
        analytic_counterexample_scope='Equal-weight mean of psi and canonical fields is not closed under the original nonlinear constraint, even when every branch is an admissible prescribed-source solution.',
        source_from_actual_quantum_record_not_proven=True,all_geometric_averages_not_ruled_out=True,
        conditional_probability_not_physical_collapse=True,
        first_response_identity_error=ucheck,leading_coefficient_min=float(predicted.min()),
        leading_coefficient_max=float(predicted.max()),rows=rows,
        next='Full joint contract must retain source-geometry correlations or alter its coarse map. Establish actual common record realization separately; do not equate labels with instruments.')
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
