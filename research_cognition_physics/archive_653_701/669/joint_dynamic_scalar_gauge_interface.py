"""669: original thermal scalar moments and the common gauge interface.

Numerics verify actual potential/mass identities and conditional gauge tensors.
Thermal moment and full polynomial RP conclusions are analytic, not a simulated
path integral or a full interacting gauge reflection-positivity theorem.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_operator_domain_completion as domain
import joint_mass_auxiliary_reflection as mass
old=mass.old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_dynamic_scalar_gauge_interface_results.json'


def scalar_moment_coefficients_check():
    constants=domain.constants()
    ell=float(max(abs(domain.LINEAR)));h=float(np.max(abs(domain.HESS)))
    cu=(28/3)*(ell**2+h**2)*constants['coercive_control']**2
    rng=np.random.default_rng(66911);ratios=[];errors=[];mass_errors=[]
    mass_basis=[mass.mass_pairing(1,domain.ball(v))[0] for v in np.eye(5)]
    for radius in (0.,.1,.8,2.,8.,25.):
        for _ in range(4):
            x=rng.normal(size=5);x*=radius/np.linalg.norm(x)
            a=x[:4]@x[:4];b=x[4]**2;t=a+b
            u,grad,_=domain.polynomial(a,b)
            exact=4*a*grad[0]**2+4*b*grad[1]**2+(2/3)*(a*grad[0]+b*grad[1])**2
            dx=2*x*np.r_[np.full(4,grad[0]),grad[1]]
            matrix=float(dx@(np.eye(5)+np.outer(x,x)/6)@dx)
            errors.append(float(abs(matrix-exact)/(1+exact)))
            ratio=float(exact/(cu*(1+u)**2));ratios.append(ratio)
            assert ratio<=1+1e-13
            p,_=mass.mass_pairing(1,domain.ball(x))
            predicted=sum(xi*pi for xi,pi in zip(x,mass_basis))
            mass_errors.append(old.norm(p-predicted)/(1+old.norm(p)))
    assert max(errors)<1e-13 and max(mass_errors)<4e-13
    # Coefficients of the actual eigenfunction moment induction, not eigenvalues.
    energy=2.;poly=[1.]
    for k in range(6):poly.append(poly[-1]*(1+energy+cu*k*k/4))
    return dict(original_target_gradient_constant=cu,
        original_coercive_constants=constants,points=len(ratios),
        maximal_sampled_gradient_bound_ratio=max(ratios),
        target_metric_gradient_identity_error=max(errors),
        mass_linear_in_original_global_coordinate_error=max(mass_errors),
        moment_recurrence_diagnostic_energy=energy,moment_polynomials=poly,
        eigenfunctions_and_thermal_integrals_not_numerically_computed=True,
        all_orders_moment_result_depends_on_note_induction=True)


def local_gauge_dictionary_check():
    nx,nt=2,2;count=nx*nt;r=32*count
    phi=np.array([.45,.31,-.17,.26,.62])
    n0,_,select,_=mass.physical(nx,nt,phi)
    s=select[:r,:r]@np.linalg.inv(n0[r:,:r]);kin=np.linalg.inv(s)
    nphys=old.diag(np.zeros_like(kin),np.zeros_like(kin))
    nphys[:r,r:]=-kin.T;nphys[r:,:r]=kin
    reps=[];phis=[];j=mass.dictionary.dictionary()
    for t in range(nt):
        for site in range(nx):
            if site==0:
                c=mass.dictionary.old.gauge.group_exp(np.array([.1,-.04,.08,.02,.03,.05,-.06,.07]),3)
                w=mass.dictionary.old.gauge.group_exp(np.array([.17,-.13,.09]),2);z=np.exp(.12j)
            else:c=np.eye(3);w=np.eye(2);z=1.
            ir=j@mass.dictionary.left_rep(c,w,z)@j.conj().T
            reps.append(np.kron(np.eye(2),ir))
            h=z**3*w@(phi[:2]+1j*phi[2:4]);phis.append(np.r_[h.real,h.imag,phi[4]])
    rlocal=np.zeros((r,r),complex)
    for i,rep in enumerate(reps):rlocal[32*i:32*i+32,32*i:32*i+32]=rep
    g=old.diag(rlocal,rlocal.conj())
    p,_=mass.mass_pairing(count,phi);pg=np.zeros_like(p)
    for i,ph in enumerate(phis):
        pi,_=mass.mass_pairing(1,ph);sl=slice(32*i,32*i+32);sr=slice(r+32*i,r+32*i+32)
        pg[sl,sl]=pi[:32,:32];pg[sr,sr]=pi[32:,32:]
    mass_covariance=old.err(g.T@pg@g-p)
    assert mass_covariance<3e-12
    defect=old.norm(g.T@nphys@g-nphys)
    assert defect>.1
    # Transporting the kinetic links repairs the exact same transformation.
    transported=g.conj()@nphys@g.conj().T
    recovery=old.err(g.T@transported@g-nphys);assert recovery<3e-12
    lam=.37;pf=old.old.old.pfaffian
    original=pf(nphys+lam*p)/pf(nphys)
    wrong=pf(nphys+lam*pg)/pf(nphys)
    correct=pf(transported+lam*pg)/pf(transported)
    assert abs(original-correct)<3e-10 and abs(original-wrong)>1e-4
    return dict(nx=nx,nt=nt,original_physical_Weyl_Nambu_dimension=len(nphys),
        time_independent_spatial_local_SM_change=True,mass_scale=lam,
        original_mass_covariance_error=mass_covariance,
        frozen_free_kinetic_gauge_defect_Frobenius=defect,
        transported_kinetic_gauge_error=recovery,
        original_mass_weight_ratio=old.cpair(original),
        scalar_only_transformed_ratio=old.cpair(wrong),
        common_mass_and_kinetic_transformed_ratio=old.cpair(correct),
        pure_gauge_transport_not_general_curved_gauge_RP=True)


def run():
    deps=('joint_operator_domain_completion.py','joint_mass_auxiliary_reflection.py',
          'joint_spinor_subgroup_mass.py','research_note_603.md','research_note_623.md',
          'research_note_624.md','research_note_643.md','research_note_668.md',
          'round669_drafts/dynamic_scalar_probe.py','round669_drafts/dynamic_scalar_probe_results.json')
    return dict(date='2026-10-02',round=669,tests_run=2,failures=0,errors=0,
        original_scalar_moment_coefficients=scalar_moment_coefficients_check(),
        joint_gauge_interface=local_gauge_dictionary_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Original bosonic heat process controls all finite-time polynomial scalar insertions; a declared free-Weyl joint dynamic-mass functional is integrable and unnormalized RP. Its real normalization polynomial is positive off finitely many possible zeros, not proved at every chosen coupling. Keeping fermion kinetics free while integrating original bosons does not restore local Gauss; actual mass and link dictionary must transform together. No full original Gibbs, general gauge RP, continuum or quantum GR identity.',
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
