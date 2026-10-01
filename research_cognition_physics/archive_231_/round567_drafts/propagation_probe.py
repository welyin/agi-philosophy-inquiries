"""Unfrozen 567 entry probe: actual h+s finite graph versus common gradients.
Fixed-background classical quadratic sector only; not a quantum or GR limit.
"""
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_finite_time_gauge_probe as inherited

def run():
    L,u,_=inherited.constants()
    D=np.diag(np.sqrt(u));mass=2*D@L@D
    m2,basis=np.linalg.eigh(mass)
    ch2=.8
    slopes=ch2*basis[0,:]**2
    rows=[]
    for q2 in (0.,.01,.1,1.,10.,100.):
        old=np.linalg.eigvalsh(mass+np.diag([ch2*q2,0.]))
        shared=np.linalg.eigvalsh(mass+ch2*q2*np.eye(2))
        assert np.max(abs(shared-(m2+ch2*q2)))<1e-12
        rows.append(dict(momentum_squared=q2,original_omega_squared=old.tolist(),
                         added_equal_s_gradient_omega_squared=shared.tolist()))
    assert np.linalg.det(mass)>0 and mass[0,1]!=0
    assert np.max(abs(slopes-slopes[::-1]))>.01
    # Independent Hessian of the original positive on-site potential.
    centre=np.sqrt(u)
    def W(x):
        d=x*x-u;return d@L@d/4
    h=1e-4
    numerical=np.zeros((2,2))
    for i in range(2):
        e=np.eye(2)[i]*h
        numerical[i,i]=(W(centre+e)-2*W(centre)+W(centre-e))/h**2
        for j in range(i):
            f=np.eye(2)[j]*h
            numerical[i,j]=numerical[j,i]=(W(centre+e+f)-W(centre+e-f)-W(centre-e+f)+W(centre-e-f))/(4*h*h)
    residual=float(np.max(abs(numerical-mass)))
    assert residual<2e-8
    return dict(status='unfrozen_probe_not_completed_round',background=centre.tolist(),
      original_mass_squared_matrix=mass.tolist(),zero_momentum_eigenvalues=m2.tolist(),
      low_momentum_velocity_squared_slopes=slopes.tolist(),rows=rows,
      independent_potential_Hessian_residual=residual,
      original_spatial_principal_rank=1,required_two_scalar_common_nonzero_cone_rank=2,
      scalar_s_edge_term_is_additional_model_input=True,
      no_claim_about_radiative_or_gravity_generated_gradients=True,
      no_new_scientific_round_count=True)

if __name__=='__main__':
    result=run()
    with (HERE/'propagation_probe_results.json').open('x',encoding='utf8',newline='\n') as f:
        f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
