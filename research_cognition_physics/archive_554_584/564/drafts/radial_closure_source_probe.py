"""564 preliminary: equal radial quantum marginals need not predict the same radial future.

Uses the unchanged full H and smooth normal Gauss families from 562.
This is an object-connection probe, not a completed numbered result.
"""
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
if str(BASE) not in sys.path:sys.path.insert(0,str(BASE))
import joint_radial_link_alignment as old
TARGET=HERE/'radial_closure_source_probe_results.json'


def radial_record(n):
    x,w=old.legendre(n);r=1+.4*x
    f=np.exp(-1/(1-x*x));weight=.4*w*r**3*f*f;weight/=weight.sum()
    R=(r[:,None]-r[None,:])**2/2;wt=weight[:,None]*weight[None,:]
    sigma=.5;y0=.4
    G=np.array([.5*(1+math.erf((q-y0)/(sigma*math.sqrt(2)))) for q in R.ravel()]).reshape(R.shape)
    Gprime=np.exp(-(R-y0)**2/(2*sigma*sigma))/(sigma*math.sqrt(2*math.pi))
    return dict(source_R_mean=float(np.sum(wt*R)),
        identical_instantaneous_event_probability=float(np.sum(wt*G)),
        radial_RGprime_mean=float(np.sum(wt*R*Gprime)))


def run():
    lam=1.;plus=old.angular(lam);minus=old.angular(-lam)
    assert abs(plus['m']+minus['m'])<1e-13
    assert abs(plus['sphere_kinetic']-minus['sphere_kinetic'])<1e-13
    radial=radial_record(160);other=radial_record(240)
    assert max(abs(radial[k]-other[k]) for k in radial)<1e-12
    Ep=old.source_budget(lam)['full_source_energy']
    Em=old.source_budget(-lam)['full_source_energy']
    # Parity U -> -U fixes the radial algebra and flips V_odd=-k X1.U.X2.
    # Exact difference of second derivatives of a bounded G(R):
    # -2<[V_odd,[T,G]]>_+ = -8*a*k*<z>_+ * <R G'(R)>.
    gap=-8*old.A*old.K*plus['m']*radial['radial_RGprime_mean']
    assert gap<0
    rho=old.radial()
    energy_gap=-2*old.K*rho['mean']**2*plus['m']
    assert abs(Ep-Em-energy_gap)<1e-11
    rng=np.random.default_rng(564);worst=0.
    for _ in range(40):
        X,Y=rng.normal(size=(2,4));r1=np.linalg.norm(X);r2=np.linalg.norm(Y)
        z=X@Y/(r1*r2);R=(r1-r2)**2/2
        gradR=np.r_[(r1-r2)*X/r1,-(r1-r2)*Y/r2]
        gradVodd=-old.K*np.r_[Y,X]
        worst=max(worst,abs(gradR@gradVodd-2*old.K*R*z))
    assert worst<1e-12
    return dict(status='preliminary_564_not_completed',lambda_=lam,angular_mean_plus=plus['m'],
        exact_same_radial_and_singlet_reduced_density_by_product_factorization=True,
        same_initial_radial_current_zero=True,same_angular_kinetic_energy=True,
        full_source_energy_plus=Ep,full_source_energy_minus=Em,
        common_energy_budget=max(Ep,Em),same_total_energy_not_claimed=True,
        radial_record=radial,exact_second_derivative_difference=gap,
        short_source_time_probability_difference_coefficient=gap/2,
        odd_potential_cross_gradient_residual=worst,
        old_family_reused_not_new_preparation_discovery=True,
        interpretation='Only excludes autonomous prediction from the radial/singlet reduced state for all these sources.',
        dependencies={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                      for name in ('joint_radial_link_alignment.py','joint_radial_link_alignment_results.json')},
        no_spacetime_or_unified_completion=True)


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
