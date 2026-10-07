"""887: same EM conditional Gibbs, source and quantum-coordinate energy.
The normal-state/kinetic connection is an explicit additional bridge contract.
No full dynamical Gauss Gibbs diagonal is assumed to be pointwise free Gibbs.
"""
from pathlib import Path
import argparse,json
import numpy as np
import renormalized_weight_variance_probe as working
old=working.old;cut=working.cut
HERE=Path(__file__).resolve().parent
TARGET=HERE/'conditional_thermal_kinetic_certificate_results.json'

def annihilator(j,n=4):
    c=np.zeros((2**n,2**n),complex)
    for b in range(2**n):
        if (b>>j)&1:
            sign=(-1)**((b & ((1<<j)-1)).bit_count())
            c[b^(1<<j),b]=sign
    return c

def electron_Fock_test(N,w,m,q,m2):
    u,wp=cut.solve_u(N,w);alpha=-(.5-u)/6;k=N//2
    ids=[26,27,28,29]
    h=old.H(m,q,N,[k,0,0],alpha)[np.ix_(ids,ids)]
    g=(m.GAMMA[0]@np.diag(q*cut.source(N,k+alpha*q)))[np.ix_(ids,ids)]
    cs=[annihilator(j) for j in range(4)]
    hF=sum(h[i,j]*(cs[i].conj().T@cs[j]) for i in range(4) for j in range(4))
    gF=sum(g[i,j]*(cs[i].conj().T@cs[j]) for i in range(4) for j in range(4))
    ev,vec=np.linalg.eigh(hF);gp=vec.conj().T@gF@vec
    p=np.exp(-old.BETA*(ev-ev.min()));p/=p.sum()
    mean=float(np.dot(p,np.diag(gp).real))
    dp=-old.BETA*p*(np.diag(gp).real-mean)
    fisher=float(np.sum(dp*dp/p))
    e=np.sqrt(w*w+m2);r=np.exp(-old.BETA*e)
    selected=144*r/(1+r)**2*w*w/(w*w+m2)*wp*wp
    expected=old.BETA**2*selected
    # Frechet derivative of normalized Gibbs matrix and its positive square root.
    dg=np.zeros_like(gp)
    for i in range(16):
        for j in range(16):
            if abs(ev[i]-ev[j])<1e-9:
                dg[i,j]=-old.BETA*p[i]*gp[i,j]
            else:
                dg[i,j]=(p[i]-p[j])/(ev[i]-ev[j])*gp[i,j]
        dg[i,i]+=old.BETA*p[i]*mean
    root_derivative=dg/(np.sqrt(p)[:,None]+np.sqrt(p)[None,:])
    purif=float(np.sum(abs(root_derivative)**2))
    rel=abs(fisher-expected)/expected
    assert rel<2e-9 and purif>=fisher/4*(1-1e-12)
    return dict(N=N,w=w,Fock_dimension=16,
        full_original_electron_Fock_spectral_Fisher=fisher,
        beta_squared_original_source_lower=expected,
        relative_identity_error=rel,
        canonical_square_root_derivative_norm=purif,
        universal_eigenvalue_lower=fisher/4,
        eigenvector_extra_nonnegative=purif-fisher/4)

def hole_profile(N,m2,order=96):
    # A scalar profile calibration of the weak-convergence escape, not a
    # substitution for the original continuum conditional matter state.
    length=.2;eps=float(N)**-4;delta=old.ETA/(6*np.sqrt(np.log(N)))
    u,_=cut.solve_u(N,1.);center=-(.5-u)/6
    x,w=np.polynomial.legendre.leggauss(order);z=(x+1)/2;wg=w/2
    step=old.smooth_step(z);der=cut.step_derivative(z)
    amp=np.sqrt(eps)+(1-np.sqrt(eps))*step
    deficit=1-eps+float(np.dot(wg,1-amp*amp))
    normalizer=length-delta*deficit
    gradient=4*(1-np.sqrt(eps))**2*float(np.dot(wg,der*der))
    kinetic=gradient/(delta*normalizer)
    for ww in (.5,2.):
        uu,_=cut.solve_u(N,ww)
        assert abs(-(.5-uu)/6-center)<delta/2
    assert center-delta>-.1 and center+delta<.1
    strip=working.integral(N,m2)*eps/normalizer
    return dict(N=N,depth=eps,width=delta,center=center,
        normalizer=normalizer,L1_distance_upper_bound=2*delta*deficit/normalizer,
        retained_selected_strip_variance=strip,
        marginal_square_root_gradient_energy=kinetic,
        kinetic_times_width=kinetic*delta,
        strip_fully_inside_low_density_plateau=True)

def run():
    prior=working.run()
    assert prior==json.loads(working.TARGET.read_text('utf-8'))
    m=old.load();_,matter=old.charges(m)
    q=working.last.em_charges(m,matter)
    m2=prior['original_mass_squared']
    fs=[electron_Fock_test(N,w,m,q,m2) for N in (17,65,257) for w in (.5,1.,2.)]
    rows=[hole_profile(N,m2) for N in (17,65,257,4097,1048577)]
    assert all(y['L1_distance_upper_bound']<x['L1_distance_upper_bound'] and
               y['retained_selected_strip_variance']<x['retained_selected_strip_variance'] and
               y['marginal_square_root_gradient_energy']>x['marginal_square_root_gradient_energy']
               for x,y in zip(rows,rows[1:]))
    a=hole_profile(257,m2,64);b=hole_profile(257,m2,128)
    qerr=abs(a['marginal_square_root_gradient_energy']-b['marginal_square_root_gradient_energy'])
    assert qerr<1e-7
    return dict(round=887,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3672,
        argument_scope='For883 hopping on the original EM flat family, exact scalar-energy thermal calibration fails the smooth conditional source test. Weak marginal convergence alone allows narrow holes. Any normal quantum-coordinate realization retaining the pointwise full free conditional Gibbs state, a weak positive target marginal, and a uniformly positive kinetic coefficient must have unbounded coordinate kinetic energy. This is a conditional-thermal ansatz obstruction, not the full dynamic Gauss Gibbs theorem.',
        prior_working_evidence_reproduced=True,
        same_background_source_lower_coefficient=prior['coefficient'],
        full64_original_matrix_checks=prior['matrix_lower_bound_checks'],
        exact_electron_Fock_checks=fs,weak_marginal_hole_rows=rows,
        hole_quadrature64_vs128_kinetic_error=qerr,
        normal_state_matrix_density_bound_is_analytic=True,
        bounded_coordinate_kinetic_plus_weak_marginal_implies_uniform_density=True,
        pointwise_free_conditional_Gibbs_is_explicit_extra_contract=True,
        uniform_coordinate_kinetic_budget_is_explicit_extra_contract=True,
        scalar_counterterm_can_change_conditional_connected_variance=False,
        full_dynamic_Gauss_Gibbs_diagonal_identified=False,
        all_preparation_repairs_excluded=False,
        full_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
