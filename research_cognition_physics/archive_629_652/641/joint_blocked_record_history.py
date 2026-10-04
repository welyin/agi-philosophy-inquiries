"""641: original nonlinear binary instrument on the common spatial Gaussian process.

The state is the explicitly declared original tree-level quadratic branch,
not the full nonlinear Gauss Gibbs state. No Gaussian replacement instrument.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_spatial_block_reference as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_blocked_record_history_results.json'
N=8
BETA=2.
MU=float(prior.matter.VAC[1])


def process(t,gamma=0.,beta=BETA):
    """Same site-pair singlet average on an actual N^3 periodic lattice."""
    masses,B,_,_=prior.mass_data()
    kx=2*np.pi*np.arange(N//2)/N
    ky=2*np.pi*np.arange(N)/N
    kz=ky.copy()
    x,y,z=np.meshgrid(kx,ky,kz,indexing='ij')
    a=np.cos(x/2)**2;b=1-a
    lam0=4*(np.sin(x/2)**2+np.sin(y/2)**2+np.sin(z/2)**2)
    lam1=4*(np.cos(x/2)**2+np.sin(y/2)**2+np.sin(z/2)**2)
    A=np.exp(-6*gamma);coef=np.exp(-4*gamma)
    V=0.;C=0.;D=0.;Cm=0.;Dm=0.;blocks=[]
    for j,m2 in enumerate(masses):
        w0=np.sqrt(m2+coef*lam0);w1=np.sqrt(m2+coef*lam1)
        ct0=1/np.tanh(beta*w0/2);ct1=1/np.tanh(beta*w1/2)
        X=A*(a*ct0/w0+b*ct1/w1)/2
        P=(a*ct0*w0+b*ct1*w1)/(2*A)
        sym=np.sqrt(X*P)
        om=np.log((2*sym+1)/(2*sym-1))/beta
        ae=om*np.sqrt(X/P)
        # s_pair=s_star+sum_j B_sj Q_j/sqrt(2), epsilon=1.
        weight=B[1,j]**2/2
        V+=weight*float(np.mean(X))
        C+=weight*float(np.mean(A*(a*ct0*np.cos(w0*t)/w0+b*ct1*np.cos(w1*t)/w1)/2))
        D+=weight*float(np.mean(A*(a*np.sin(w0*t)/w0+b*np.sin(w1*t)/w1)))
        Cm+=weight*float(np.mean(X*np.cos(om*t)))
        Dm+=weight*float(np.mean(ae*np.sin(om*t)/om))
        blocks.append(dict(mass_squared=float(m2),min_frequency=float(min(w0.min(),w1.min()))))
    return dict(V=V,C=C,D=D,Cm=Cm,Dm=Dm,mu=MU,blocks=blocks)


def coefficients(sign,K=16,grid=2048):
    q=2*np.pi*np.arange(grid)/grid
    L=np.sqrt(.5+sign*.25*np.sin(q))
    fft=np.fft.fft(L)/grid
    modes=np.arange(-K,K+1)
    return modes,fft[modes%grid]


def branch_prob(sign,d):
    return .5+sign*.25*np.exp(-d['V']/2)*np.sin(d['mu'])


def sandwich(sign,d,K=16,commute=False):
    k,c=coefficients(sign,K)
    kk=k[:,None];ll=k[None,:];n=kk+ll
    # <e^(ik S0) e^(i St) e^(il S0)>, [S0,St]=i D.
    phase=np.exp(1j*((n+1)*d['mu']-(kk-ll)*(0. if commute else d['D'])/2))
    gaussian=np.exp(-.5*((n*n+1)*d['V']+2*n*d['C']))
    return complex(np.sum(c[:,None]*c[None,:]*phase*gaussian))


def probabilities(d,K=16,commute=False):
    out=np.empty((2,2))
    for i,r in enumerate((1,-1)):
        p=branch_prob(r,d)
        s=sandwich(r,d,K,commute).imag
        out[i]=[.5*p+.25*s,.5*p-.25*s]
    assert abs(out.sum()-1)<1e-13 and out.min()>0
    return out


def history_check():
    rows=[];tail=0.;marginal=0.;coefficient_error=0.;fine_error=0.
    # Analytic strip |Im z|<=1.2 avoids the square-root zeros sin z=+/-2.
    strip=1.2;bound=np.sqrt(.5+.25*np.cosh(strip))
    delta=2*bound*np.exp(-strip*17)/(1-np.exp(-strip))
    probability_tail=.25*delta*(2*np.sqrt(.75)+delta)
    q=np.linspace(-2*np.pi,2*np.pi,311)
    for r in (1,-1):
        k,c=coefficients(r)
        reconstructed=np.exp(1j*q[:,None]*k)@c
        coefficient_error=max(coefficient_error,float(np.max(abs(reconstructed-np.sqrt(.5+r*.25*np.sin(q))))))
    for t in (0.,.4,1.,2.):
        d=process(t);p=probabilities(d);fine=probabilities(d,K=24)
        # Independent unblocked N^3 Fourier sum for the same pair observable.
        axes=2*np.pi*np.arange(N)/N
        x,y,z=np.meshgrid(axes,axes,axes,indexing='ij')
        lam=4*(np.sin(x/2)**2+np.sin(y/2)**2+np.sin(z/2)**2)
        filt=(1+np.cos(x))/2
        masses,B,_,_=prior.mass_data()
        vv=cc=dd=0.
        for j,m2 in enumerate(masses):
            om=np.sqrt(m2+lam);ct=1/np.tanh(BETA*om/2);w=B[1,j]**2*filt
            vv+=float(np.mean(w*ct/(2*om)))
            cc+=float(np.mean(w*ct*np.cos(om*t)/(2*om)))
            dd+=float(np.mean(w*np.sin(om*t)/om))
        fine_error=max(fine_error,abs(vv-d['V']),abs(cc-d['C']),abs(dd-d['D']))
        tail=max(tail,float(np.max(abs(p-fine))))
        dm=dict(d,C=d['Cm'],D=d['Dm'])
        pm=probabilities(dm);pc=probabilities(d,commute=True)
        first=np.array([branch_prob(r,d) for r in (1,-1)])
        marginal=max(marginal,float(np.max(abs(p.sum(axis=1)-first))),float(np.max(abs(pm.sum(axis=1)-first))))
        if t==0:
            # Independent exact integral for commuting repetitions E_r(S)E_s(S).
            sinmean=np.exp(-d['V']/2)*np.sin(MU)
            sin2=(1-np.exp(-2*d['V'])*np.cos(2*MU))/2
            expected=np.array([[.25+(r+s)*sinmean/8+r*s*sin2/16 for s in (1,-1)] for r in (1,-1)])
            assert np.max(abs(p-expected))<2e-12
        rows.append(dict(t=t,variance=d['V'],symmetric_time_covariance=d['C'],
            commutator=d['D'],exact_original_history=p.tolist(),
            mean_force_history=pm.tolist(),commuting_surrogate=pc.tolist(),
            total_variation_from_mean_force=float(np.sum(abs(p-pm))/2),
            total_variation_when_commutator_omitted=float(np.sum(abs(p-pc))/2)))
    assert coefficient_error<5e-12 and tail<2e-12 and marginal<1e-13 and fine_error<1e-13
    assert max(r['total_variation_from_mean_force'] for r in rows)>.0001
    assert max(r['total_variation_when_commutator_omitted'] for r in rows)>1e-7
    return dict(rows=rows,Fourier_function_error=coefficient_error,
        K16_vs_K24_probability_error=tail,shared_first_record_error=marginal,
        independent_unblocked_spatial_covariance_error=fine_error,
        analytic_K16_function_tail_bound=float(delta),
        analytic_per_event_probability_tail_bound=float(probability_tail),
        original_square_root_sine_instrument_not_Gaussian_replacement=True)


def derivative_prob(d,dd,K=16):
    k,c=coefficients(1,K);kk=k[:,None];ll=k[None,:];n=kk+ll
    factor=-.5*(n*n+1)*dd['V']-n*dd['C']-.5j*(kk-ll)*dd['D']
    base=np.exp(1j*((n+1)*d['mu']-(kk-ll)*d['D']/2)
        -.5*((n*n+1)*d['V']+2*n*d['C']))
    out=np.empty((2,2))
    for i,r in enumerate((1,-1)):
        _,c=coefficients(r,K)
        ds=np.sum(c[:,None]*c[None,:]*base*factor).imag
        dp=-r*.125*np.exp(-d['V']/2)*np.sin(MU)*dd['V']
        out[i]=[.5*dp+.25*ds,.5*dp-.25*ds]
    return out


def geometry_record_check():
    rows=[]
    for t in (.4,1.,2.):
        d=process(t)
        # Different steps: differentiate the underlying common covariance first,
        # then the characteristic identity, versus full history differences.
        h=1e-5
        dp=process(t,h);dm=process(t,-h)
        dd={k:(dp[k]-dm[k])/(2*h) for k in ('V','C','D')}
        analytic=derivative_prob(d,dd)
        errors=[]
        for step in (2e-4,1e-4,5e-5):
            fd=(probabilities(process(t,step))-probabilities(process(t,-step)))/(2*step)
            errors.append(float(np.max(abs(fd-analytic))))
        p=probabilities(d);score=analytic/p
        fisher=float(np.sum(analytic**2/p))
        # State-only variation: retain the true equal-time derivative but erase
        # the derivatives of the correlation and commutator kernels.
        incomplete=derivative_prob(d,dict(V=dd['V'],C=0.,D=0.))
        mismatch=float(np.max(abs(analytic-incomplete)))
        assert errors[-1]<2e-8 and abs(analytic.sum())<1e-12
        assert abs(float(np.sum(p*score)))<1e-12 and fisher>0
        assert mismatch>1e-5
        rows.append(dict(t=t,original_geometry_covariance_derivatives=dd,
            joint_probability_derivative=analytic.tolist(),joint_record_score=score.tolist(),
            record_Fisher_information=fisher,finite_difference_errors=errors,
            missing_kernel_derivative_error=mismatch))
    return dict(rows=rows,score_is_record_response_not_full_stress_tensor=True,
        uniform_geometry_keeps_spatial_average_definition_fixed=True)


def run():
    deps=('research_note_586.md','research_note_624.md','research_note_625.md',
          'research_note_639.md','research_note_640.md','joint_spatial_block_reference.py')
    return dict(round=641,tests_run=2,failures=0,errors=0,
        actual_nonlinear_record_history=history_check(),shared_geometry_record_response=geometry_record_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_instrument_function_on_linearized_singlet_average=True,
            actual_N8_three_dimensional_quadratic_lattice=True,
            full_interacting_Gauss_state_not_computed=True,
            same_covariance_and_commutator_control_record_and_geometry_score=True,
            no_autonomous_implementation_continuum_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
