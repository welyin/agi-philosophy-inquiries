"""632: shared fermion potential, local Weyl contacts and background conditions.

Uses the original full one-generation mass content. Finite MS prescription
and optional gauge-invariant on-shell matching conditions are explicit inputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_continuum_source_spectrum as previous
import joint_fermion_gauss_completion as matter
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_background_contact_matching_results.json'
L,U0,_=matter.original.lattice.scalar.parameters()
FI=np.array([-1/6,-1/6])


def jordan_eigen_jets(x):
    X,Y=x;a=abs(matter.Y['nu'])**2;b=abs(matter.Y['s'])**2
    P=b*b*Y*Y+4*a*b*X*Y;R=np.sqrt(P)
    p=np.array([4*a*b*Y,2*b*b*Y+4*a*b*X])
    pp=np.array([[0.,4*a*b],[4*a*b,2*b*b]])
    r=p/(2*R);rr=pp/(2*R)-np.outer(p,p)/(4*R**3)
    vals=[];grads=[];hess=[];weights=[]
    for name,d in (('u',3.),('d',3.),('e',1.)):
        c=abs(matter.Y[name])**2
        vals.append(c*X);grads.append([c,0.]);hess.append(np.zeros((2,2)));weights.append(d)
    for sign in (-1,1):
        vals.append(a*X+b*Y/2+sign*R/2)
        grads.append(np.array([a,b/2])+sign*r/2)
        hess.append(sign*rr/2);weights.append(.5)
    return np.array(vals),np.array(grads),np.array(hess),np.array(weights)


def W_jets(x,mu2=1.):
    # W=F^2 Vf is the Jordan numerator of the fixed-Einstein fermion loop.
    t,g,h,d=jordan_eigen_jets(x)
    f=matter.original.M-np.sum(x)/6;log=np.log(t/(f*mu2));c=-1/(16*np.pi**2)
    W=c*np.sum(d*t*t*(log-1.5))
    grad=c*np.sum(d[:,None]*(2*t[:,None]*(log-1)[:,None]*g-t[:,None]**2*FI/f),axis=0)
    H=np.zeros((2,2))
    for ti,gi,hi,di,li in zip(t,g,h,d,log):
        H+=c*di*(2*li*np.outer(gi,gi)+2*ti*(li-1)*hi
                    -2*ti/f*(np.outer(gi,FI)+np.outer(FI,gi))
                    +ti*ti/f**2*np.outer(FI,FI))
    return float(W),grad,H


def canonical_mass_data(x):
    t,_,_,d=jordan_eigen_jets(x)
    f=matter.original.M-np.sum(x)/6
    return np.sqrt(t/f),d


def fermion_scaling_jets(m,d,mu2=1.):
    a=d*m**4/(16*np.pi**2);logs=np.log(m*m/mu2)
    V=-np.sum(a*(logs-1.5));vp=-np.sum(4*a*(logs-1))
    vpp=-np.sum(a*(12*logs-4));C=float(np.sum(d*m**4)/(8*np.pi**2))
    return float(V),float(vp),float(vpp),C


def Vscale(lam,m,d,mu2=1.):
    a=d*(m*lam)**4/(16*np.pi**2)
    return float(-np.sum(a*(np.log((m*lam)**2/mu2)-1.5)))


def hessian(fun,x,step=1e-4):
    ans=np.zeros((len(x),len(x)));basis=np.eye(len(x))*step
    for i in range(len(x)):
        ans[i,i]=(fun(x+basis[i])-2*fun(x)+fun(x-basis[i]))/step**2
        for j in range(i):
            ans[i,j]=ans[j,i]=(fun(x+basis[i]+basis[j])-fun(x+basis[i]-basis[j])
                            -fun(x-basis[i]+basis[j])+fun(x-basis[i]-basis[j]))/(4*step**2)
    return ans


def weyl_contact_check():
    _,m,d,_,_=previous.original_data();q0=previous.QSTAR
    V,vp,vpp,C=fermion_scaling_jets(m,d)
    b=1/np.sqrt(6)/np.tanh(q0/np.sqrt(6));c=1/6
    assert abs(vp-4*V+C)<1e-15
    assert abs(vpp-3*vp+4*C)<1e-15
    rng=np.random.default_rng(632);identity_error=0.
    for _ in range(24):
        eta,sigma=rng.uniform(-.12,.12,2);lam=1+eta
        cov=np.exp(4*sigma)*Vscale(lam,m,d)
        transformed=Vscale(np.exp(sigma)*lam,m,d)+C*sigma*np.exp(4*sigma)*lam**4
        identity_error=max(identity_error,abs(cov-transformed))
    expected=np.array([[b*b*vpp+c*vp,4*b*vp],[4*b*vp,16*V]])
    naive=vpp*np.outer([b,1],[b,1])
    contacts=np.array([[c*vp,b*(vp+4*C)],[b*(vp+4*C),vp+8*C]])
    assert np.max(abs(expected-naive-contacts))<1e-15
    def covariant(z):
        q,sigma=z;lam=np.sinh(q/np.sqrt(6))/np.sinh(q0/np.sqrt(6))
        return np.exp(4*sigma)*Vscale(lam,m,d)
    fd=hessian(covariant,np.array([q0,0.]),step=3e-5)
    error=float(np.max(abs(fd-expected)))
    assert identity_error<1e-16 and error<2e-9
    assert np.linalg.norm(contacts)>.001
    return dict(Vf=V,mass_first_derivative=vp,mass_second_derivative=vpp,
        local_Weyl_log_coefficient=C,original_radial_b=float(b),
        covariant_hessian_q_sigma=expected.tolist(),
        naive_common_mass_hessian=naive.tolist(),
        required_local_contact_matrix=contacts.tolist(),
        exact_finite_Weyl_identity_error=identity_error,
        full_hessian_finite_difference_error=error,
        no_complete_trace_anomaly_or_general_contact_Ward_claim=True)


def background_check():
    W,grad,H=W_jets(U0);f=matter.original.M-sum(U0)/6
    m,d=canonical_mass_data(U0)
    phi=np.array([0.,np.sqrt(U0[0]),0.,0.,np.sqrt(U0[1])])
    h,delta=matter.mass_matrices(phi)
    bdg=np.block([[h,delta],[-delta.conj(),-h.T]])
    expected_m=np.r_[np.repeat(m[:3],[6,6,2]),m[3:]]
    mass_error=float(np.max(abs(np.linalg.eigvalsh(bdg)[32:]-np.sort(np.repeat(expected_m,2)))))
    assert mass_error<3e-15
    # Independent derivative differences validate all eigenvalue/F terms.
    step=2e-5;gfd=[];hfd=[]
    for i in range(2):
        e=np.eye(2)[i]*step
        gfd.append((W_jets(U0+e)[0]-W_jets(U0-e)[0])/(2*step))
        hfd.append((W_jets(U0+e)[1]-W_jets(U0-e)[1])/(2*step))
    error=max(float(np.max(abs(np.array(gfd)-grad))),float(np.max(abs(np.array(hfd).T-H))))
    assert error<2e-10
    # Gauge-invariant finite matching δV_J= -W0-gradW .(x-u).
    # A declared renormalization condition, not a predicted counterterm.
    c0=-W;linear=-grad
    value=(W+c0)/f**2;force=(grad+linear)/f**2
    assert abs(value)<1e-15 and np.max(abs(force))<1e-15
    D=np.diag(2*np.sqrt(U0));K=matter.original.metric(phi)[np.ix_([1,4],[1,4])]
    k,v=np.linalg.eigh(K);ki=(v*(1/np.sqrt(k)))@v.T
    rows=[]
    for lam in (0.,.1,1.):
        radial=D@(L/2+lam*H)@D/f**2
        eigen=np.linalg.eigvalsh(ki@radial@ki)
        rows.append(dict(loop_multiplier=lam,radial_hessian=radial.tolist(),
                         canonical_radial_eigenvalues=eigen.tolist()))
    assert min(rows[-1]['canonical_radial_eigenvalues'])>0
    # Old unnumbered 554 has already derived the generic first-order shift.
    shift=-2*np.linalg.solve(L,grad)
    return dict(original_L=L.tolist(),original_squared_vacuum=U0.tolist(),
        F=float(f),masses=m.tolist(),Dirac_weights=d.tolist(),
        mass_dictionary_error=mass_error,loop_Jordan_numerator=W,
        loop_numerator_gradient=grad.tolist(),loop_numerator_hessian=H.tolist(),
        derivative_check_error=error,
        declared_finite_constant=c0,declared_finite_linear_in_x=linear.tolist(),
        anchored_value=value,anchored_gradient=force.tolist(),
        radial_hessian_rows=rows,
        old_shift_formula_reused_not_new=shift.tolist(),
        constant_only_residual_gradient=(grad/f**2).tolist(),
        scalar_gauge_orbit_zero_modes_are_not_physical_stability_claims=True,
        finite_matching_conditions_are_additional_input=True)


def contact_and_state_check():
    # At constant backgrounds, sqrt(g)U fixes all zero-derivative metric
    # contacts. Test covariant volume variation without choosing a cutoff.
    W,g,H=W_jets(U0);f=matter.original.M-sum(U0)/6
    V=W/f**2
    grad=g/f**2-2*W*FI/f**3
    mixed=4*grad
    single_ct=-W
    residual_value=(W+single_ct)/f**2
    residual_gradient=g/f**2
    # Only cancelling vacuum value leaves mixed scalar/metric contact.
    assert abs(residual_value)<1e-15 and np.linalg.norm(4*residual_gradient)>1e-5
    e=np.diag([1.,-1.,0])/np.sqrt(2)
    eps=2e-4
    jac=lambda z:np.sqrt(np.linalg.det(np.eye(3)+2*z*e))
    fd=(jac(eps)-2*jac(0.)+jac(-eps))/eps**2
    assert abs(fd+2)<1e-7
    # Gauge-invariant finite numerator matching cancels value and gradient
    # exactly at u; scalar Hessian remains changed, checked independently.
    return dict(fermion_only_value=V,fermion_only_x_gradient=grad.tolist(),
        fermion_only_constant_contacts=dict(trace_trace=16*V,
                                            shear_shear=-2*V,
                                            x_trace=mixed.tolist()),
        constant_only_total_mixed_x_trace=(4*residual_gradient).tolist(),
        on_shell_total_trace_trace=0.,on_shell_total_shear_shear=0.,
        on_shell_total_mixed_x_trace=[0.,0.],
        volume_shear_second_derivative=float(fd),
        nonlocal_spectra_unchanged_by_finite_local_counterterms=True,
        reference_branch='Same chosen continuum vacuum, now evaluated at the original tree vacuum. Not the finite graph Gauss thermal state.',
        inherited_630_diagnostic_background_not_used_as_a_solution=True,
        q_sigma_contact_cancellation_does_not_remove_finite_frequency_noise=True)


def run():
    deps=('research_note_552.md','research_note_553.md','research_note_599.md',
          'research_note_600.md','research_note_601.md','research_note_602.md',
          'research_note_630.md','research_note_631.md',
          'round554_drafts/vacuum_matching_scope_review.md',
          'round554_drafts/vacuum_matching_probe.py',
          'joint_continuum_source_spectrum.py','joint_fermion_gauss_completion.py',
          'joint_curved_quantum_source.py')
    return dict(round=632,tests_run=3,failures=0,errors=0,
        finite_Weyl_and_contacts=weyl_contact_check(),
        joint_background_matching=background_check(),
        common_zero_momentum_metric_conditions=contact_and_state_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original complete one-generation fermion determinant in the fixed 3+1 Einstein-frame continuum vacuum branch, finite MS prescription and declared optional gauge-invariant constant/two-mass matching conditions. Explicit zero-derivative Weyl logarithm and source contacts, local homogeneous background and radial Hessian of the one-fermion-loop truncated potential. Not complete SM loops, full momentum-dependent Ward completion, RG-invariant vacuum prediction, nonperturbative stability, graph state matching or GR generation.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))

