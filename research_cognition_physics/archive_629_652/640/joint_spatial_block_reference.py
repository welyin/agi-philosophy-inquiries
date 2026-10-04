"""640: actual spatial Haar blocking of the original neutral quadratic sector.

Flat constant vacuum, canonical quadratic quantization, hbar=1. Both original
radial masses are retained; no pointer is added. Not the full interacting Gauss
Gibbs state, not a nonperturbative continuum limit. The discarded modes belong
to the same spatial lattice and are not an externally supplied heat bath.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_pointer_matter_compatibility as matter
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_spatial_block_reference_results.json'
N=8
EPS=1.


def mass_data():
    x=matter.VAC;f=original.M-np.dot(x,x)/6
    K=matter.base_metric(x)
    V=2*np.diag(x)@matter.MAT@np.diag(x)/f**2
    masses,B=matter.generalized_modes(V,K)
    return masses,B,K,V


def spatial_data():
    # x dependence of a genuine N^3 torus, transverse momentum zero.
    # Q is the real cosine wave of nearest-pair averages along x.
    L=2*np.eye(N)-np.roll(np.eye(N),1,axis=0)-np.roll(np.eye(N),-1,axis=0)
    nc=N//2;theta=2*np.pi/nc
    qc=np.sqrt(2/nc)*np.cos(theta*np.arange(nc))
    q=np.repeat(qc,2)/np.sqrt(2)
    u=float(q@L@q)
    r=L@q-u*q;r/=np.linalg.norm(r)
    rows=np.stack((q,r))
    block=rows@L@rows.T
    k=theta/2
    lam=np.array([4*np.sin(k/2)**2,4*np.cos(k/2)**2])
    weights=np.array([np.cos(k/2)**2,np.sin(k/2)**2])
    return L,rows,block,lam,weights


def data(m2,beta,gamma=0.):
    _,_,_,lam,a=spatial_data()
    A=np.exp(-6*gamma)
    d=m2+np.exp(-4*gamma)*lam/EPS**2
    om=np.sqrt(d);ct=1/np.tanh(beta*om/2)
    X=float(A*np.sum(a*ct/om)/2)
    P=float(np.sum(a*ct*om)/(2*A))
    sym=np.sqrt(X*P)
    mf_omega=float(np.log((2*sym+1)/(2*sym-1))/beta)
    mf_A=float(mf_omega*np.sqrt(X/P))
    mf_V=float(mf_omega*np.sqrt(P/X))
    logZ=float(-np.sum(np.log(2*np.sinh(beta*om/2))))
    logZm=float(-np.log(2*np.sinh(beta*mf_omega/2)))
    shift=(logZm-logZ)/beta
    u=float(a@d);v=float(a[::-1]@d)
    c2=float(np.prod(a)*(d[1]-d[0])**2)
    return dict(A=A,d=d,omega=om,X=X,P=P,sym=sym,
        mf_omega=mf_omega,mf_A=mf_A,mf_V=mf_V,logZ=logZ,
        scalar_shift=shift,u=u,v=v,c2=c2)


def spatial_quadratic_check():
    m2,B,K,V=mass_data();L,rows,block,lam,a=spatial_data()
    assert np.max(abs(B.T@K@B-np.eye(2)))<2e-14
    assert np.max(abs(B.T@V@B-np.diag(m2)))<2e-14
    invariant=float(np.max(abs(L@rows.T-rows.T@block)))
    assert invariant<2e-14
    assert np.max(abs(rows@rows.T-np.eye(2)))<2e-14
    # The eliminated partner lies entirely in Haar differences.
    assert np.max(abs(rows[1].reshape(-1,2).sum(axis=1)))<2e-14
    ev,U=np.linalg.eigh(block)
    assert np.max(abs(ev-lam))<2e-14
    assert np.max(abs(U[0]**2-a))<2e-14
    # Independent original geodesic and potential Hessian on the N-node slice.
    # N transverse constant layers multiply the action; normalized zero
    # transverse momentum leaves this same matrix.
    base=np.array([0.,matter.VAC[0],0.,0.,matter.VAC[1]])
    rng=np.random.default_rng(640);z=rng.normal(size=(N,2))*.15
    dx=z@B.T
    delta=np.zeros((N,5));delta[:,1]=dx[:,0];delta[:,4]=dx[:,1]
    def energy(step):
        ph=base+step*delta
        return float(np.sum(original.node_potential(ph))+
            np.sum(original.distance_squared(ph,np.roll(ph,-1,axis=0)))/(2*EPS**2))
    predicted=float(np.sum(z*z*m2)+sum(z[:,j]@L@z[:,j]/EPS**2 for j in range(2)))
    errors=[]
    for h in (2e-3,1e-3,5e-4):
        diff=(energy(h)+energy(-h)-2*energy(0))/h**2
        errors.append(float(abs(diff-predicted)))
    assert errors[-1]<errors[0]/10 and errors[-1]<2e-7
    return dict(original_mass_squared=m2.tolist(),original_radial_basis=B.tolist(),
        original_target_and_potential_Hessian_errors=errors,
        actual_spatial_average_row=rows[0].tolist(),discarded_Haar_partner=rows[1].tolist(),
        spatial_block=block.tolist(),fine_laplacian_eigenvalues=lam.tolist(),
        average_spectral_weights=a.tolist(),invariant_block_error=invariant,
        finite_torus_N=N,transverse_momenta_zero=True,no_pointer_added=True)


def reference_and_dynamics_check():
    masses,_,_,_=mass_data();_,_,_,lam,a=spatial_data()
    records=[];maxres=0.;source_records=[]
    for m2 in masses:
        thermal=[]
        for beta in (.5,2.,8.,30.):
            d=data(m2,beta)
            mfct=1/np.tanh(beta*d['mf_omega']/2)
            xx=d['mf_A']*mfct/(2*d['mf_omega'])
            pp=d['mf_V']*mfct/(2*d['mf_omega'])
            maxres=max(maxres,abs(xx-d['X']),abs(pp-d['P']))
            thermal.append(dict(beta=beta,X=d['X'],P=d['P'],
                symplectic_eigenvalue=d['sym'],
                mean_force_frequency=d['mf_omega'],mean_force_A=d['mf_A'],
                mean_force_V=d['mf_V'],partition_scalar_shift=d['scalar_shift']))
        assert abs(thermal[0]['mean_force_A']-thermal[-1]['mean_force_A'])>.1
        d=data(m2,2.)
        defect=d['c2'] # determinant of first three commutator spectral moments
        assert defect>.4
        direct=[];schur_errors=[]
        for nu in (0.,.2,.7,1.5,3.):
            z=nu*nu
            exact=d['A']*np.sum(a/(z+d['d']))
            S=z+d['u']-d['c2']/(z+d['v'])
            schur_errors.append(float(abs(exact-d['A']/S)))
            direct.append(dict(euclidean_frequency=nu,blocked_resolvent=float(exact),
                mean_force_resolvent=float(d['mf_A']/(z+d['mf_omega']**2))))
        # State-independent commutator and thermal symmetric correlator.
        times=np.array([.3,.8,1.4])
        true=np.sum(a[:,None]*np.sin(d['omega'][:,None]*times)/d['omega'][:,None],axis=0)*d['A']
        fit=d['mf_A']*np.sin(d['mf_omega']*times)/d['mf_omega']
        assert np.max(abs(true-fit))>.005 and max(schur_errors)<1e-13
        vacuum_X=np.sum(a/np.sqrt(m2+lam))/2
        vacuum_P=np.sum(a*np.sqrt(m2+lam))/2
        records.append(dict(mass_squared=float(m2),thermal_states=thermal,
            vacuum_symplectic_eigenvalue=float(np.sqrt(vacuum_X*vacuum_P)),
            two_pole_moment_determinant=float(defect),resolvents=direct,
            Schur_resolvent_error=max(schur_errors),
            times=times.tolist(),actual_commutator=true.tolist(),
            mean_force_commutator=fit.tolist(),maximum_commutator_error=float(np.max(abs(true-fit)))))
    assert maxres<1e-12
    return dict(rows=records,thermal_covariance_reconstruction_error=maxres,
        exact_single_time_Gaussian_state_match_not_real_time_generator_match=True)


def geometry_source_check():
    masses,_,_,_=mass_data();_,_,_,lam,a=spatial_data()
    rows=[];h=2e-5;beta=2.
    for m2 in masses:
        d=data(m2,beta);plus=data(m2,beta,h);minus=data(m2,beta,-h)
        der=lambda k:(plus[k]-minus[k])/(2*h)
        # Full actual pair source in fixed original canonical coordinates.
        ct=1/np.tanh(beta*d['omega']/2)
        Xfine=ct/(2*d['omega']);Pfine=ct*d['omega']/2
        source=float(np.sum((-6*Pfine+(2*lam+6*m2)*Xfine)/2))
        free_energy_derivative=-der('logZ')/beta
        mf_without_shift=(der('mf_A')*d['P']+der('mf_V')*d['X'])/2
        mf_with_shift=mf_without_shift+der('scalar_shift')
        assert abs(source-free_energy_derivative)<2e-8
        assert abs(source-mf_with_shift)<2e-8
        assert abs(source-mf_without_shift)>.1
        # Euclidean kernel determinant carries the same geometric derivative.
        detrows=[]
        for nu in (0.,.3,1.,2.):
            z=nu*nu
            def logs(dd):
                S=z+dd['u']-dd['c2']/(z+dd['v'])
                return np.array([np.sum(np.log(z+dd['d'])),np.log(z+dd['v']),np.log(S)])
            ld=(logs(plus)-logs(minus))/(2*h)
            assert abs(ld[0]-ld[1]-ld[2])<1e-9
            detrows.append(dict(nu=nu,full_frequency_determinant_derivative=float(ld[0]),
                eliminated_factor_derivative=float(ld[1]),retained_Schur_derivative=float(ld[2])))
        rows.append(dict(mass_squared=float(m2),beta=beta,original_geometry_source=source,
            free_energy_derivative=float(free_energy_derivative),
            mean_force_source_without_scalar=float(mf_without_shift),
            scalar_normalization_derivative=float(der('scalar_shift')),
            mean_force_source_with_scalar=float(mf_with_shift),
            total_source_error=float(abs(source-mf_with_shift)),
            frequency_determinant_checks=detrows))
    return dict(rows=rows,geometry_is_original_uniform_psi_equals_exp_gamma=True,
        no_arbitrary_vacuum_subtraction=True,geometric_reference_normalization_required=True)


def run():
    deps=('research_note_567.md','research_note_574.md','research_note_576.md',
          'research_note_617.md','research_note_625.md','research_note_639.md',
          'joint_pointer_matter_compatibility.py','joint_curved_quantum_source.py')
    return dict(round=640,tests_run=3,failures=0,errors=0,
        spatial_quadratic=spatial_quadratic_check(),
        state_and_dynamics=reference_and_dynamics_check(),geometry=geometry_source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_two_neutral_quadratic_masses=True,actual_neighbor_pair_blocking=True,
            quadratic_Gibbs_not_full_interacting_Gauss_Gibbs=True,hbar_one_explicit=True,
            exact_state_requires_mean_force_and_geometry_normalization=True,
            same_mean_force_not_exact_real_time_generator=True,
            no_full_gauge_blocking_continuum_limit_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
