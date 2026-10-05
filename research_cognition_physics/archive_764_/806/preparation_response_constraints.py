"""806: positive preparation increments link covariance, record and source.

Original full-mass derivatives calibrate the fermion response functional.
Finite CCR covariance and source matrices check the joint identities and their
limits; they are not a physical-solution basis of the full original PDE.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'804'))
import covariance_contract_probe as old
TARGET=HERE/'preparation_response_constraints_results.json'

def run():
    data=old.run()
    j=np.array([-2*r['correct_ordered_covariance'][1] for r in data['five_original_mass_derivatives']])
    r=np.eye(5)[:,4]
    rng=np.random.default_rng(806)
    maximum=0.;samples=[]
    for i in range(12):
        a=rng.normal(size=(5,3));gamma=a@a.T/7
        signal=float(r@gamma@j);vb=float(r@gamma@r);vj=float(j@gamma@j)
        defect=signal**2-vb*vj
        maximum=max(maximum,defect)
        assert defect<2e-14
        samples.append(dict(response_change=signal,old_variance_increase=vb,current_variance_increase=vj))
    projector=np.eye(5)-np.outer(r,r)
    a=rng.normal(size=(5,5));protected=projector@a@a.T@projector
    assert r@protected@r==0 and abs(r@protected@j)<1e-14
    # Positive increment preserves a chosen indefinite source, but changes signal.
    gamma=np.zeros((5,5));gamma[3,3]=gamma[4,4]=1.;gamma[3,4]=gamma[4,3]=.5
    t=np.diag([0,0,0,1.,-1.])
    positive_source=np.diag([1.,1.,1.,2.,3.])
    assert np.trace(t@gamma)==0 and np.linalg.eigvalsh(gamma).min()>=0
    signal=float(r@gamma@j);cost=float(np.trace(positive_source@gamma))
    bound=cost*np.linalg.norm(r/np.sqrt(np.diag(positive_source)))*np.linalg.norm(j/np.sqrt(np.diag(positive_source)))
    assert abs(signal)<=bound+1e-14 and abs(signal)>.1
    # Arbitrary differences of positive states need NOT themselves be positive.
    delta=np.zeros((5,5));delta[3,4]=delta[4,3]=.2
    symplectic=np.block([[np.zeros((5,5)),np.eye(5)],[-np.eye(5),np.zeros((5,5))]])
    covariance0=np.eye(10)
    increment=np.zeros((10,10));increment[:5,:5]=delta
    uncertainty0=float(np.linalg.eigvalsh(covariance0+.5j*symplectic).min())
    uncertainty1=float(np.linalg.eigvalsh(covariance0+increment+.5j*symplectic).min())
    assert uncertainty0>0 and uncertainty1>0
    assert r@delta@r==0 and np.trace(increment)==0 and abs(r@delta@j)>.05
    return dict(round=806,all_checks_passed=True,random_PSD_samples=12,
        original_mass_derivative_response_functional=j.tolist(),
        Cauchy_Schwarz_max_positive_defect=maximum,sample_budget_checks=samples,
        preserved_variance_positive_increment=dict(response_change=float(r@protected@j),variance_change=float(r@protected@r)),
        indefinite_source_counterexample=dict(source_change=float(np.trace(t@gamma)),
            response_change=signal,positive_resource_change=cost,resource_response_upper_bound=float(bound)),
        two_valid_states_with_indefinite_difference=dict(min_uncertainty_eigenvalues=[uncertainty0,uncertainty1],
            old_variance_change=float(r@delta@r),total_quadratic_energy_change=float(np.trace(increment)/2),
            response_change=float(r@delta@j),difference_eigenvalues=np.linalg.eigvalsh(delta).tolist()),
        finite_modes_are_original_continuous_solutions=False,
        full_original_stress_source_computed=False,original_continuous_nonzero_response_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
