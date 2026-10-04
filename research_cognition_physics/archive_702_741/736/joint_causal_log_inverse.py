"""736: exact old spectral kernels, regular remainder and a causal inverse.

The short-time theorem is analytic in note736.  Numerical checks retain all
original630/631 masses and verify independently computable integral identities.
This is not a nonlinear semiclassical spacetime solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
import joint_continuum_source_spectrum as scalar
import joint_tensor_stress_spectrum as tensor

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_causal_log_inverse_results.json'


@lru_cache(maxsize=12)
def quad(n=192):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)/2,w/2


def profile(x,channel):
    return (1-x*x)**1.5*(1+2*x*x/3 if channel=='tensor' else 1)


def data(channel):
    _,m,d,_,_=scalar.original_data()
    c=d*m*m/(4*np.pi**2) if channel=='scalar' else d/(80*np.pi**2)
    beta=4/3 if channel=='scalar' else 6/5
    b=float(np.dot(c,-np.log(m)-beta))
    return m,d,c,b


def old_normalized(s,channel,n=256):
    m,d,_,_=data(channel)
    if channel=='scalar':
        return scalar.finite_subtracted(1j*s,m,d,n)/s**2
    return -tensor.tensor_dispersion(1j*s,m,d,n)/s**4


def remainder(s,channel,n=256):
    x,w=quad(n);m,d,c,b=data(channel)
    # 1-a(x), evaluated without catastrophic cancellation near x=0.
    base=np.exp(1.5*np.log1p(-x*x))
    gap=-np.expm1(1.5*np.log1p(-x*x))
    if channel=='tensor':gap-=2*x*x*base/3
    ans=0j
    for mi,ci in zip(m,c):
        ans+=ci*(.5*np.log(1+(2*mi/s)**2)+
                 np.dot(w,gap*(4*mi*mi)/(x*(4*mi*mi+s*s*x*x))))
    return ans


def jet_identity_check():
    rows=[]
    for channel in ('scalar','tensor'):
        m,d,c,b=data(channel)
        maximum=0.;sign_error=0.
        for s in (1.1+.3j,2.7+1.4j,.4+1.1j,3.2-.7j):
            left=old_normalized(s,channel)
            right=c.sum()*np.log(s)+b+remainder(s,channel)
            maximum=max(maximum,abs(left-right))
            # Independent positive-measure identity for w=s^2. This is the
            # key zero-free test, not a finite grid proof of zero-freeness.
            x,q=quad(256);w=s*s
            positive=0.
            for mi,ci in zip(m,c):
                positive+=ci*np.dot(q,profile(x,channel)*4*mi*mi*x/abs(4*mi*mi+w*x*x)**2)
            observed=left.imag/w.imag
            sign_error=max(sign_error,abs(observed-positive))
            assert positive>0
        # The constants are independently checked by the beta-type integral.
        x,q=quad(384)
        constant=float(np.dot(q,(profile(x,channel)-1)/x))
        beta=4/3 if channel=='scalar' else 6/5
        assert abs(constant-(np.log(2)-beta))<2e-12
        assert maximum<3e-12 and sign_error<3e-13
        rows.append(dict(channel=channel,total_log_coefficient=float(c.sum()),
                         exact_local_constant=b,normalized_constant=b/float(c.sum()),
                         complex_spectral_decomposition_error=float(maximum),
                         positive_measure_identity_error=float(sign_error),
                         independent_constant_integral_error=abs(constant-np.log(2)+beta)))
    return dict(rows=rows,all_original_masses_retained=True,
                no_threshold_or_long_memory_replacement=True)


def log_resolvent_integral(s,n=192):
    # 1/log(s) = 1/(s-1) + int_0^inf dr/[(s+r)(log(r)^2+pi^2)].
    # Subtract the positive-y step in r=exp(y); the discarded tails are <1e-17.
    x,q=quad(n);y=40*x
    pos=1/(1+s*np.exp(-y))-1
    neg=np.exp(-y)/(s+np.exp(-y))
    integral=.5+40*np.dot(q,(pos+neg)/(y*y+np.pi**2))
    return 1/(s-1)+integral


def inverse_norm(T,n=192):
    # Exact L1 norm of the positive inverse-log kernel on [0,T], mu=1.
    ell=np.log(1/T);x,q=quad(n)
    minus=-40*x;plus=7*x
    low=-np.expm1(-np.exp(minus))/((minus+ell)**2+np.pi**2)
    high=-np.exp(-np.exp(plus))/((plus+ell)**2+np.pi**2)
    cut=.5-np.arctan(ell/np.pi)/np.pi
    return float(np.expm1(T)+cut+40*np.dot(q,low)+7*np.dot(q,high))


def inverse_norm_upper(T):
    assert 0<T<1
    ell=np.log(1/T)
    return float(np.expm1(T)+np.sqrt(T)+.5-np.arctan(ell/(2*np.pi))/np.pi)


def remainder_norm_upper(T,channel):
    m,d,c,b=data(channel)
    assert 2*max(m)*T<=1
    return float(np.sum(c*m*m*T*T*(5.5+3*np.log(1/(2*m*T)))))


def inverse_contract_check():
    errors=[]
    for s in (2.,3.+2j,1.5-.6j):
        errors.append(abs(log_resolvent_integral(s)-1/np.log(s)))
    assert max(errors)<2e-13
    T=.01
    norm=inverse_norm(T)
    upper=inverse_norm_upper(T)
    assert 0<norm<upper and abs(norm-inverse_norm(T,384))<2e-13
    rows=[]
    for channel in ('scalar','tensor'):
        m,d,c,b=data(channel)
        reg=remainder_norm_upper(T,channel)
        # A sufficient contraction certificate for the exact I_p u = f
        # inverse; all additional physical local terms remain explicit in
        # the theorem, not silently set to these diagnostic values.
        rate=upper*(abs(b)+reg)/c.sum()
        assert rate<1
        rows.append(dict(channel=channel,exact_remainder_L1_upper=reg,
                         sufficient_contraction_upper=float(rate),
                         normalized_inverse_bound=upper/(1-float(rate)),
                         no_claim_that_full_Einstein_local_terms_vanish=True))
    return dict(inverse_log_Laplace_identity_error=float(max(errors)),
                auxiliary_log_inverse_pole=1.,pole_not_a_full_response_instability=True,
                time_window=T,inverse_log_L1_norm=norm,analytic_L1_upper=upper,
                original_exact_kernel_certificates=rows,
                weak_causal_response_inverse_not_nonlinear_spacetime=True)


def run():
    names=('jet_identity_check','inverse_contract_check')
    results={name:globals()[name]() for name in names}
    deps=('research_note_601.md','research_note_630.md','research_note_631.md','research_note_632.md',
          'research_note_735.md','joint_continuum_source_spectrum.py','joint_tensor_stress_spectrum.py',
          'round736_drafts/response_regularity_entry.py','round736_drafts/response_regularity_entry_results.json')
    return dict(round=736,tests_run=2,failures=0,errors=0,checks=list(names),results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original630/631 complete-generation stationary vacuum branch: exact logarithm-plus-L1-memory decomposition and a short-time causal inverse on the active scalar/shear response sectors. Additional finite local operators obey an explicit contraction criterion. No full constraint elimination, general dynamic-background inverse, nonlinear semiclassical solution, or long-time stability is claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as handle:
            handle.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
