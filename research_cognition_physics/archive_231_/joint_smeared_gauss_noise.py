"""621: thermal form-source noise in finite time windows.

The full Gauss theorem is analytic. Numerical checks use the original
neutral Fock factor and an explicitly separate infinite-dimensional example.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quantum_response_matching as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_smeared_gauss_noise_results.json'


def size(x):return float(np.max(np.abs(x)))


def triangle_hat(omega,tau):
    return np.sinc(omega*tau/(2*np.pi))**2


def original_factor_check():
    _,pairs=old.matrices(old.S0,[24,25,30,31])
    H,J,_=[old.fock(p) for p in pairs]
    e,v,p,rho=old.thermal(H,2.)
    a=e-e.min()+1
    delta=e[:,None]-e[None,:]
    sources=[v.conj().T@g@v for g in (H,J)]
    moments={k:float(p@a**k) for k in (1,2,5,6)}
    rows=[]
    for tau in (.05,.2,.8):
        smeared=[];bounds=[];Ls=[];integration_errors=[]
        for g,center in zip(sources,(0.,.4)):
            L=float(np.linalg.norm(g/np.sqrt(a[:,None]*a[None,:]),2))
            sf=triangle_hat(delta,tau)*np.exp(1j*center*delta)*g
            quad=np.zeros_like(g)
            x,w=np.polynomial.legendre.leggauss(48)
            # Two halves handle the triangular window's cusp exactly.
            for left,right in ((-tau,0.),(0.,tau)):
                ts=left+(x+1)*(right-left)/2
                ws=w*(right-left)/2
                for t,weight in zip(ts,ws):
                    quad+=weight*(1-abs(t)/tau)/tau*np.exp(1j*(t+center)*delta)*g
            err=size(quad-sf);assert err<2e-14
            integration_errors.append(err)
            second=float(np.sum(p[:,None]*np.abs(sf)**2))
            bound=L*L*(moments[2]+2/tau*moments[1])
            assert second<=bound+1e-13
            smeared.append(sf);bounds.append(bound);Ls.append(L)
        means=np.array([np.sum(p*np.diag(g)).real for g in smeared])
        gram=np.array([[np.trace(np.diag(p)@g@h) for h in smeared] for g in smeared])
        centered=gram-np.outer(means,means)
        assert size(centered-centered.conj().T)<1e-14
        assert np.linalg.eigvalsh(centered).min()>-1e-13
        assert abs(centered[0,1])**2 <= centered[0,0].real*centered[1,1].real+1e-13
        # Spectral truncation at fixed state, without rethermalizing.
        cutrows=[]
        for cutoff in (1.25,1.65,2.1):
            keep=a<=cutoff
            errors=[];tail_bounds=[]
            for g,sf,L in zip(sources,smeared,Ls):
                diff=sf-np.outer(keep,keep)*sf
                hs2=float(np.sum(p[:,None]*abs(diff)**2))
                high=a>cutoff/2
                omega=cutoff/2
                m0=16/(tau**4*omega**4)
                m1=16/(tau**4*omega**3)
                bound=L*L*(np.sum(p[high]*a[high]**2)+2/tau*np.sum(p[high]*a[high])
                            +m0*moments[2]+m1*moments[1])
                assert hs2<=bound+1e-13
                errors.append(hs2);tail_bounds.append(float(bound))
            cutrows.append(dict(A_cutoff=cutoff,weighted_matrix_tail_squared=errors,analytic_tail_bounds=tail_bounds))
        rows.append(dict(tau=tau,form_norms=Ls,second_moment_bounds=bounds,
            noise_matrix=centered.real.tolist(),commutator_imaginary_part=centered.imag.tolist(),
            direct_time_integral_errors=integration_errors,cutoff_rows=cutrows))
    # Detailed balance is tested on the unsmeared original scalar matrix.
    delta_pos=e[None,:]-e[:,None]
    g=sources[1]
    lhs=p[None,:]*abs(g)**2
    rhs=np.exp(-2*delta_pos)*p[:,None]*abs(g)**2
    kms_error=size(lhs-rhs);assert kms_error<1e-14
    return dict(original_Fock_dimension=16,beta=2.,A_moments={str(k):v for k,v in moments.items()},rows=rows,
                detailed_balance_error=kms_error,not_a_full_graph_spectrum_computation=True)


def infinite_form_example():
    beta=.7;q=np.exp(-beta);p0=1-q;constant=6/np.pi**2
    m1=1/(1-q);m2=(1+q)/(1-q)**2
    cutoffs=[32,128,512,2048,8192]
    rows=[]
    for cutoff in cutoffs:
        n=np.arange(1,cutoff+1,dtype=float)
        weights=p0*(1+q**n)*constant*(n+1)/n**2
        raw=float(weights.sum())
        smooth=[]
        for tau in (.2,.6,1.5):
            val=float(np.sum(weights*triangle_hat(n,tau)**2))
            tail=16*p0*constant/(tau**4*cutoff**4)
            bound=m2+2/tau*m1
            assert val<=bound
            smooth.append(dict(tau=tau,variance_lower_sum=val,omitted_positive_tail_bound=float(tail),
                               general_form_bound=float(bound)))
        rows.append(dict(cutoff=cutoff,instant_second_moment_lower_sum=raw,smeared=smooth))
    assert all(rows[i+1]['instant_second_moment_lower_sum']>rows[i]['instant_second_moment_lower_sum']+.3 for i in range(4))
    for j in range(3):
        for i in range(4):
            lo=rows[i]['smeared'][j];hi=rows[-1]['smeared'][j]
            assert 0<=hi['variance_lower_sum']-lo['variance_lower_sum']<=lo['omitted_positive_tail_bound']+2e-15
    return dict(beta=beta,source_form_norm=1.,thermal_A_moments=[float(m1),float(m2)],rows=rows,
        analytic_raw_divergence='(6/pi^2)*(1-exp(-beta))*sum(1/n) is a lower bound',
        abstract_example_not_the_original_Gauss_model=True,
        all_cutoffs_use_the_same_infinite_Gibbs_state=True)


def run():
    a=original_factor_check();b=infinite_form_example()
    deps=('research_note_590.md','research_note_591.md','research_note_603.md','research_note_620.md',
          'joint_quantum_response_matching.py','joint_thermal_gauss_source.py')
    return dict(round=621,tests_run=2,failures=0,errors=0,original_factor=a,infinite_example=b,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_fixed_graph_Gauss_statement_proved_analytically=True,
            stationary_thermal_state_and_form_bound_required=True,
            joint_GNS_source_vectors_not_self_adjoint_measurement_construction=True,
            finite_time_window_not_a_continuum_UV_completion=True,
            instantaneous_original_noise_not_proved_divergent=True,
            no_GR_or_internal_background_generation_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
