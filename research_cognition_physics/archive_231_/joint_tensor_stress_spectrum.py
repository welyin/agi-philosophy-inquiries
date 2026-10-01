"""631: full rest-frame stress cut of the original continuum fermion branch."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_continuum_source_spectrum as old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_tensor_stress_spectrum_results.json'


def spin2(omega,masses,weights):
    omega=np.asarray(omega,float);out=np.zeros_like(omega)
    for m,d in zip(masses,weights):
        mask=omega>2*m;w=omega[mask];r=m*m/(w*w)
        out[mask]+=d*w**4/(240*np.pi)*(1-4*r)**1.5*(3+8*r)
    return out


def dirac_matrices():
    Z=np.zeros((2,2))
    sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    return np.array([np.block([[Z,s],[s,Z]]) for s in sig]),np.diag([1.,1.,-1.,-1.])


def stress_basis():
    e=[np.eye(3)/np.sqrt(3),np.diag([1.,-1.,0])/np.sqrt(2),np.diag([1.,1.,-2.])/np.sqrt(6)]
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3));a[i,j]=a[j,i]=1/np.sqrt(2);e.append(a)
    e=np.array(e)
    assert np.max(abs(np.einsum('aij,bij->ab',e,e)-np.eye(6)))<1e-15
    return e


def angular_directions(n=6):
    u,wu=np.polynomial.legendre.leggauss(n)
    phi=2*np.pi*np.arange(12)/12
    return [(np.array([np.sqrt(1-v*v)*np.cos(p),np.sqrt(1-v*v)*np.sin(p),v]),weight/24)
            for v,weight in zip(u,wu) for p in phi]


def angular_cut(omega,masses,weights):
    alpha,beta=dirac_matrices();basis=stress_basis()
    out=np.zeros((7,7),complex)
    for m,d in zip(masses,weights):
        if omega<=2*m:continue
        E=omega/2;k=np.sqrt(E*E-m*m)
        shell=d*k*E/(2*np.pi)
        for direction,weight in angular_directions():
            kv=k*direction;H=np.einsum('i,ijk->jk',kv,alpha)+m*beta
            plus=(np.eye(4)+H/E)/2;minus=np.eye(4)-plus
            vertices=[m*beta]
            vertices.extend(np.einsum('i,ijk->jk',e@kv,alpha) for e in basis)
            vertices=np.array(vertices)
            left=np.einsum('ij,ajk,kl->ail',minus,vertices,plus)
            gram=np.einsum('aij,bji->ab',left,vertices)
            out+=shell*weight*gram
    return out


def tensor_check():
    _,masses,weights,_,_=old.original_data()
    rows=[]
    for omega in (1.5*min(masses),2.1*min(masses),4*max(masses)):
        direct=angular_cut(omega,masses,weights)
        rho0=float(old.density(np.array([omega]),masses,weights)[0])
        rho2=float(spin2(np.array([omega]),masses,weights)[0])
        expected=np.diag([rho0,rho0/3]+[rho2]*5)
        expected[0,1]=expected[1,0]=-rho0/np.sqrt(3)
        error=float(np.max(abs(direct-expected)))
        eigen=np.linalg.eigvalsh(direct)
        assert error<2e-14 and eigen[0]>-1e-14
        null=np.r_[1.,np.sqrt(3),np.zeros(5)]
        assert abs(null@direct@null)<1e-14
        rows.append(dict(frequency=float(omega),rho_mass=rho0,rho_spin2=rho2,
                         angular_error=error,spectral_eigenvalues=eigen.tolist(),
                         predicted_rank=6 if rho0>0 else 0))
    # Covariant projectors on arbitrary timelike momenta, signature (-+++).
    g=np.diag([-1.,1.,1.,1.]);rng=np.random.default_rng(631);max_ward=0.
    for _ in range(20):
        momentum=rng.normal(size=3);s=float(rng.uniform(.3,2.))
        q=np.r_[np.sqrt(s+momentum@momentum),momentum];ql=g@q
        pi=g+np.outer(ql,ql)/s
        p0=np.einsum('ij,kl->ijkl',pi,pi)/3
        p2=(np.einsum('ik,jl->ijkl',pi,pi)+np.einsum('il,jk->ijkl',pi,pi))/2-p0
        max_ward=max(max_ward,float(np.max(abs(np.einsum('i,ijkl->jkl',q,p0)))),
                     float(np.max(abs(np.einsum('i,ijkl->jkl',q,p2)))))
        assert np.max(abs(np.einsum('ij,ijkl->kl',g,p2)))<5e-13
    assert max_ward<5e-13
    return dict(original_Weyl_components=16,angular_nodes=72,
                rows=rows,max_covariant_projector_Ward_error=max_ward,
                tensor_cut_structure='rho2 P2 + (rho_mass/3) P0',
                five_TT_channels_in_addition_to_trace=True,
                timelike_off_shell_channels_not_five_on_shell_graviton_polarizations=True,
                contact_terms_not_covered_by_cut_Ward_test=True)


def curvature_match_check():
    _,masses,weights,_,_=old.original_data()
    # Independent heat-kernel spin trace, omitting total derivatives only.
    riem=F(4*2,360)-F(1,24)
    ric=F(-4*2,360)
    r2=F(4*5,360)+F(1,8)-F(1,6)
    # riem=bW+bE, ric=-2bW-4bE, r2=bW/3+bE+bR.
    bE=-(ric+2*riem)/2;bW=riem-bE;bR=r2-bW/3-bE
    assert (bW,bE,bR)==(F(-1,20),F(11,360),0)
    D=float(weights.sum());c=-float(bW)*D
    betaW=c/(16*np.pi**2)
    expected=D/(80*np.pi**2)
    cutoff=500*max(masses)
    observed=float(spin2(np.array([cutoff]),masses,weights)[0]/(np.pi*cutoff**4))
    assert abs(observed/expected-1)<2e-5
    assert abs(4*betaW-expected)<1e-16
    # h_ij=2 epsilon e_ij, tr(e)=0, tr(e^2)=1:
    # int W^2 to second order is 2 int (epsilon'')^2, Hessian 4 Q^4.
    assert F(2)*F(2)==4
    zero_mass_spin2=float(spin2(np.array([1.]),np.zeros(16),weights)[0])
    assert abs(zero_mass_spin2-D/(80*np.pi))<1e-16
    assert old.density(np.array([1.]),np.zeros(16),weights)[0]==0
    return dict(Dirac_heat_kernel_Riem2_Ric2_R2=[str(riem),str(ric),str(r2)],
        Dirac_heat_kernel_Weyl2_Euler_R2=[str(bW),str(bE),str(bR)],
        original_one_generation_c=c,beta_bW_fermion_one_generation=betaW,
        shear_source_Weyl2_Hessian_multiplier=4,
        spin2_z4_log_from_cut=observed,spin2_z4_log_from_heat_kernel=expected,
        zero_mass_one_generation_rho2_at_unit_frequency=zero_mass_spin2,
        zero_mass_trace_cut=0.,
        no_new_massless_graviton_pole=True,no_Newton_constant_prediction=True,
        no_full_553_three_generation_scalar_vector_count_substitution=True)


def tensor_dispersion(z,masses,weights,n=192):
    x,w=old.quadrature(n);value=0j
    for m,d in zip(masses,weights):
        value+=d/(960*np.pi**2*m*m)*np.dot(
            w,x*(1-x*x)**1.5*(3+2*x*x)/(1-(z*x/(2*m))**2))
    return z**6*value


def joint_window_check():
    _,masses,weights,_,_=old.original_data();gap=2*min(masses)
    coeff=float(np.sum(weights/masses**2)/(1344*np.pi**2))
    rows=[]
    for ratio in (.2,.6,.9):
        z=ratio*gap
        value=tensor_dispersion(z,masses,weights)
        fine=tensor_dispersion(z,masses,weights,384)
        bound=z**6*coeff/(1-ratio**2)
        assert abs(value-fine)<2e-13 and 0<value.real<=bound*(1+1e-12)
        rows.append(dict(frequency=float(z),sixth_order_remainder=float(value.real),
                         bound=float(bound),quadrature_change=float(abs(value-fine))))
    omega=float(4*max(masses));rho2=float(spin2(np.array([omega]),masses,weights)[0])
    # Gaussian approximation, at each frequency, can reproduce these second
    # moments with one shared mass/trace variable and five traceless variables.
    # No claim that the fermion stress itself is Gaussian or fundamental.
    b=1/np.sqrt(6)/np.tanh(old.QSTAR/np.sqrt(6))
    rho0=float(old.density(np.array([omega]),masses,weights)[0])
    covariance=np.zeros((7,7))
    covariance[:2,:2]=rho0*np.outer([b,1],[b,1])/2
    covariance[2:,2:]=np.eye(5)*rho2/2
    assert np.linalg.matrix_rank(covariance,tol=1e-12)==6
    return dict(common_pair_threshold=float(gap),spin2_sixth_moment=coeff,
        below_threshold_rows=rows,reference_frequency=omega,
        shear_noise_per_unit_normalized_source=rho2/2,
        mass_metric_and_five_shear_covariance_rank=6,
        three_subtractions_required_for_spin2=True,
        independent_finite_coefficients=['constant','frequency_squared','frequency_fourth'],
        only_second_moments_not_full_Gaussian_quantum_process=True,
        fixed_background_not_self_consistent_Einstein_solution=True)


def run():
    deps=('research_note_553.md','research_note_599.md','research_note_601.md',
          'research_note_620.md','research_note_628.md','research_note_630.md',
          'joint_continuum_source_spectrum.py','joint_fermion_gauss_completion.py')
    return dict(round=631,tests_run=3,failures=0,errors=0,
        tensor_cut=tensor_check(),common_curvature_coefficient=curvature_match_check(),
        tensor_effective_window=joint_window_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Same original one-generation continuum free-fermion vacuum branch as 630. Complete parity-even stress two-point cut on flat space and its timelike Lorentz tensor reconstruction, with the shared mass source and one-loop Weyl-squared logarithm. This is off-shell matter response, not a graviton state, gravitational dynamics, full interacting theory, general curved-background Ward construction or graph-continuum matching.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))

