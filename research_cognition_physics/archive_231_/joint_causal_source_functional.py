"""622: same-state quadratic CTP coefficients and causal source matching.

Numerical original sector is the 602 neutral factor. The full Gauss
coefficient-convergence statement is proved in the accompanying note.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quantum_response_matching as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_causal_source_functional_results.json'
T=2.4
BETA=2.
_,PAIRS=old.matrices(old.S0,[24,25,30,31])
H,J,C=[old.fock(x) for x in PAIRS]
E,V,P,RHO=old.thermal(H,BETA)
GS=[V.conj().T@g@V for g in (H,J)]
CONTACT=np.array([[0.,np.trace(RHO@J).real],[np.trace(RHO@J).real,np.trace(RHO@C).real]])
GAP=E[:,None]-E[None,:]
WEIGHT=(P[:,None]-P[None,:])*abs(GS[1])**2
F0=float(old.matter.original.F(np.array([0.,old.HIGGS,0.,0.,old.S0])))
ROOT=np.sqrt(F0)
NUM1=ROOT*J-old.S0/(6*ROOT)*H
NUM0=ROOT*H-old.S0*NUM1


def pulse(t):
    f=np.sin(np.pi*t/T)**2
    return np.array([.2*f,f]),np.array([-.1*f,.4*f*(1+.3*np.sin(2*np.pi*t/T))])


def quadrature(order,left=0.,right=T):
    x,w=np.polynomial.legendre.leggauss(order)
    return left+(x+1)*(right-left)/2,w*(right-left)/2


def covariance(a,b):
    ma=np.sum(P*np.diag(a)).real;mb=np.sum(P*np.diag(b)).real
    return float(np.trace(np.diag(P)@a@b).real-ma*mb)


def coefficients(order):
    ts,ws=quadrature(order)
    integrated=[np.zeros_like(H),np.zeros_like(H)]
    contact=0.;response=0.
    for t,w in zip(ts,ws):
        plus,minus=pulse(t);d=plus-minus;avg=(plus+minus)/2
        phase=np.exp(1j*GAP*t)
        for a in range(2):integrated[a]+=w*d[a]*phase*GS[a]
        contact+=w*float(d@CONTACT@avg)
        us,vs=quadrature(order,0.,t)
        for u,vw in zip(us,vs):
            pp,mm=pulse(u);cu=(pp+mm)/2
            chi=float((1j*np.sum(WEIGHT*np.exp(1j*GAP*(t-u)))).real)
            # All thermal commutators involving the energy H vanish.
            response+=w*vw*d[1]*chi*cu[1]
    obs=sum(integrated)
    mean=float(np.sum(P*np.diag(obs)).real)
    noise=covariance(obs,obs)
    mixed=noise-sum(covariance(x,x) for x in integrated)
    return dict(linear=[0.,-mean],quadratic=[-.5*noise,response-contact],
                noise=noise,mixed_noise=mixed,contact=contact,retarded=response)


def coefficient_check():
    a=coefficients(40);b=coefficients(64)
    err=max(abs(a[k]-b[k]) for k in ('noise','mixed_noise','contact','retarded'))
    assert err<1e-13
    assert b['noise']>0 and abs(b['mixed_noise'])>1e-4
    assert abs(b['contact'])>1e-3 and abs(b['retarded'])>1e-5
    # The primitive source second derivative is independently differenced.
    eps=1e-4
    def mass(n,s):return n*(NUM0+s*NUM1)/np.sqrt(F0-(s*s-old.S0**2)/6)
    cfd=(mass(1.,old.S0+eps)-2*H+mass(1.,old.S0-eps))/eps**2
    second_error=float(np.max(abs(cfd-C)));assert second_error<2e-7
    return dict(original_neutral_Fock_dimension=16,beta=BETA,duration=T,
                coefficients=b,quadrature_difference=err,second_source_derivative_error=second_error,
                original_complex_Yukawa_Majorana_and_lapse_retained=True)


def unitary(epsilon,branch,steps):
    u=np.eye(len(H),dtype=complex);dt=T/steps
    for k in range(steps):
        lam=pulse((k+.5)*dt)[branch]
        n=1+epsilon*lam[0];s=old.S0+epsilon*lam[1]
        h=n*(NUM0+s*NUM1)/np.sqrt(F0-(s*s-old.S0**2)/6)
        e,v=np.linalg.eigh(h)
        u=((v*np.exp(-1j*dt*e))@v.conj().T)@u
    return u


def log_overlap(epsilon,steps,swap=False):
    up=unitary(epsilon,0,steps);um=unitary(epsilon,1,steps)
    if swap:up,um=um,up
    z=np.trace(up@RHO@um.conj().T)
    assert abs(z)<=1+2e-12
    return np.log(z),float(abs(np.trace(up@RHO@up.conj().T)-1))


def actual_process_check():
    theory=coefficients(48);q2=complex(*theory['quadratic']);q1=complex(*theory['linear'])
    rows=[]
    for epsilon,steps in ((.04,768),(.02,768),(.01,768),(.01,1536)):
        lp,normal=log_overlap(epsilon,steps)
        lm,_=log_overlap(-epsilon,steps)
        l0,_=log_overlap(0.,steps)
        got2=(lp+lm-2*l0)/(2*epsilon**2)
        got1=(lp-lm)/(2*epsilon)
        err=abs(got2-q2)
        assert err<2e-4 and normal<2e-12
        rows.append(dict(epsilon=epsilon,time_steps=steps,
            quadratic_log_coefficient=[float(got2.real),float(got2.imag)],
            quadratic_error=float(err),linear_error=float(abs(got1-q1)),
            equal_history_normalization_error=normal))
    lp,_=log_overlap(.02,768);swap,_=log_overlap(.02,768,swap=True)
    conjugation_error=float(abs(swap-lp.conjugate()));assert conjugation_error<2e-13
    assert rows[-1]['quadratic_error']<2e-5
    return dict(rows=rows,branch_exchange_conjugation_error=conjugation_error,
        one_fixed_initial_state_used=True,
        finite_amplitude_original_neutral_factor_verified=True,
        full_Gauss_limit_and_source_differentiation_not_interchanged=True)


def high_transition_check():
    beta=.7;p0=1-np.exp(-beta);end=2*np.pi
    ts,ws=quadrature(256,0.,end)
    pulse1=(1-np.cos(ts))/2
    rows=[]
    for n in (2,4,16,64):
        pn=p0*np.exp(-beta*n)
        hat=np.sum(ws*pulse1*np.exp(1j*n*ts))
        # Independent direct integral of the forced-oscillator solution.
        y=1/(2*n)-n*np.cos(ts)/(2*(n*n-1))+np.cos(n*ts)/(2*n*(n*n-1))
        causal_integral=float(np.sum(ws*pulse1*y))
        exact=np.pi*(3*n*n-2)/(4*n*(n*n-1))
        response=2*(p0-pn)*(n+1)*causal_integral
        predicted=2*(p0-pn)*(n+1)*exact
        noise=(p0+pn)*(n+1)*abs(hat)**2
        assert abs(hat)<2e-13 and noise<1e-22
        assert abs(response-predicted)<2e-13 and response>2
        # G_n=A^(1/2)B_n A^(1/2), with B_n a norm-one swap.
        instantaneous=(p0+pn)*(n+1)
        rows.append(dict(high_level=n,form_norm=1.,pulse_noise=float(noise),
            retarded_response=float(response),formula_error=float(abs(response-predicted)),
            instantaneous_second_moment=float(instantaneous)))
    return dict(beta=beta,pulse='sin(t/2)^2 on [0,2*pi]',rows=rows,
        asymptotic_response=float(3*np.pi*p0/2),
        noise_exactly_zero_by_Fourier_orthogonality_for_integer_levels_ge_2=True,
        abstract_uniform_class_counterexample_not_original_Gauss_source=True,
        fixed_H_and_fixed_thermal_state_but_different_sources=True)


def run():
    a=coefficient_check();b=actual_process_check();c=high_transition_check()
    deps=('research_note_591.md','research_note_602.md','research_note_603.md','research_note_620.md',
          'research_note_621.md','joint_quantum_response_matching.py','joint_smeared_gauss_noise.py')
    return dict(round=622,tests_run=3,failures=0,errors=0,same_source_coefficients=a,
        original_process=b,high_transition=c,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_fixed_Gauss_regulated_second_order_coefficients_converge=True,
            arbitrary_full_process_differentiation_not_proved=True,
            same_state_contact_noise_and_retardation_required=True,
            no_continuum_or_GR_completion=True,
            mature_CTP_formula_applied_not_new_fundamental_theory=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
