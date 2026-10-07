"""813: original Nambu transport, local CAR repair and source-support checks.
The Fourier diagnostic uses the original constant auxiliary-past matrices.
It is NOT the 753 nonlinear future, an autonomous instrument, or a GR derivation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
TARGET=HERE/'carrier_localization_results.json'

def charge(f):
    return np.einsum('ab,...b->...a',original.CHARGE,f.conj())

def inner(f,g):
    return np.mean(np.einsum('...a,...a->...',f.conj(),g))

def fft_derivative(f):
    n=f.shape[0];k=np.fft.fftfreq(n,1/n)
    return np.fft.ifft((1j*k).reshape((n,)+(1,)*(f.ndim-1))*np.fft.fft(f,axis=0),axis=0)

def original_packet(n,t):
    x=2*np.pi*np.arange(n)/n
    modes=np.arange(-7,8); coeff=[]
    for k in modes:
        f=np.zeros(64,complex)
        f[30]=np.exp(-.17*k*k)*(1+.09j*k)
        f[31]=.37*np.exp(-.23*(k-1)**2)*np.exp(.41j*k)
        val,vec=np.linalg.eigh(original.hamiltonian(np.array([float(k),0.,0.])))
        coeff.append(vec@(np.exp(-1j*t*val)*(vec.conj().T@f)))
    f=np.exp(1j*x[:,None]*modes[None,:])@np.array(coeff)
    return x,f/np.sqrt(inner(f,f).real)

def localized_mode(f,chi):
    cf=charge(f)
    e=np.stack([(f+cf)/np.sqrt(2),1j*(f-cf)/np.sqrt(2)],axis=-1)
    cut=chi[:,None,None]*e
    gram=np.einsum('nsa,nsb->ab',cut.conj(),cut)/len(f)
    assert np.max(abs(gram.imag))<1e-12
    val,vec=np.linalg.eigh(gram.real)
    inv=(vec/np.sqrt(val))@vec.T
    v=cut@inv
    g=(v[:,:,0]-1j*v[:,:,1])/np.sqrt(2)
    delta=inner(f*(np.sqrt(1-chi**2))[:,None],f*(np.sqrt(1-chi**2))[:,None]).real
    bound=np.sqrt(delta)+1-np.sqrt(1-2*delta) if delta<.5 else None
    naive=chi[:,None]*f
    naive/=np.sqrt(inner(naive,naive).real)
    error=np.sqrt(inner(f-g,f-g).real)
    row=dict(tail=float(delta),gram_min=float(val.min()),
        gram_bound_error=float(max(0,1-2*delta-val.min())),
        naive_CAR_anomaly=float(abs(inner(charge(naive),naive))),
        repaired_CAR_anomaly=float(abs(inner(charge(g),g))),
        repaired_normalization_error=float(abs(inner(g,g)-1)),
        mode_error=float(error),mode_error_bound=None if bound is None else float(bound),
        projection_error_bound=None if bound is None else float(min(1,2*bound)))
    assert max(row['gram_bound_error'],row['repaired_CAR_anomaly'],row['repaired_normalization_error'])<1e-12
    if bound is not None:assert error<=bound+1e-12
    return row

def packet_checks():
    rows=[]
    for n in (64,128):
        x,f=original_packet(n,.7)
        dx=fft_derivative(f);h=-1j*dx@original.GAMMA[0].T+f@original.MASS.T
        ft=-1j*h
        rho=np.sum(abs(f)**2,axis=1)
        dt=2*np.sum(f.conj()*ft,axis=1).real
        j=np.array([np.einsum('na,ab,nb->n',f.conj(),g,f).real for g in original.GAMMA])
        div=2*np.einsum('na,ab,nb->n',f.conj(),original.GAMMA[0],dx).real
        chi=.96+.04*np.cos(x-.27)
        w=chi**2
        dweight=.03*np.sin(x)
        actual=float(np.mean(dweight*rho+w*dt))
        flux=float(np.mean(dweight*rho+fft_derivative(w).real*j[0]))
        errors=dict(continuity=float(np.max(abs(dt+div))),
            causal_current=float(max(0,np.max(np.linalg.norm(j,axis=0)-rho))),
            moving_weight_flux=float(abs(actual-flux)),
            global_CAR_isotropy=float(abs(inner(charge(f),f))))
        assert max(errors.values())<1e-11
        q=np.diag([1.]*32+[-1.]*32)
        charge_rate=float(2*inner(f@q.T,ft).real)
        local=localized_mode(f,chi)
        rows.append(dict(grid=n,time=.7,errors=errors,localization=local,
            same_Nambu_charge_rate=charge_rate,weighted_norm_rate=actual))
    assert rows[-1]['localization']['naive_CAR_anomaly']>1e-7
    assert abs(rows[-1]['same_Nambu_charge_rate'])>1e-7
    return rows

def annihilation(n,index):
    a=np.zeros((2**n,2**n),complex)
    for word in range(2**n):
        if word>>index&1:
            sign=(-1)**((word&((1<<index)-1)).bit_count())
            a[word^(1<<index),word]=sign
    return a

def source_support_check():
    # Independent exact-size CAR check of the support argument, not original P.
    a=[annihilation(2,i) for i in range(2)];xi=a+[x.conj().T for x in a]
    state=np.zeros(4,complex);state[0]=np.sqrt(3)/2;state[3]=.5
    def cov(v):
        return np.array([[np.vdot(v,xi[i]@xi[j].conj().T@v) for j in range(4)] for i in range(4)])
    P=cov(state);proj=np.eye(4)-a[0].conj().T@a[0]
    prob=float(np.vdot(state,proj@state).real)
    v1=proj@state/np.sqrt(prob);v0=(np.eye(4)-proj)@state/np.sqrt(1-prob)
    P1,P0=cov(v1),cov(v0);bar=prob*P1+(1-prob)*P0
    Q=np.diag([1.,0.,1.,0.]);R=np.eye(4)-2*Q
    out=np.array([1,3]);block=lambda a:a[np.ix_(out,out)]
    reflected=.5*(P+R@P@R)
    errors=dict(reflection_formula=float(np.max(abs(bar-reflected))),
        exterior_nonselective_source=float(np.max(abs(block(bar-P)))),
        exterior_conditional_cancellation=float(np.max(abs(prob*block(P1-P)+(1-prob)*block(P0-P)))))
    assert max(errors.values())<1e-12
    branch=float(np.max(abs(block(P1-P))))
    assert branch>.2
    return dict(errors=errors,empty_probability=prob,
        exterior_conditional_covariance_change=branch,
        supports_all_local_quadratic_source_matrices=True,original_continuum_covariance_evaluated=False)

def run():
    return dict(round=813,all_checks_passed=True,original_Nambu_dimension=64,
        packet_checks=packet_checks(),source_support=source_support_check(),
        scope='Original auxiliary-past Fourier coefficient checks plus independent CAR support calibration.',
        original_future_PDE_solved=False,full_interacting_record_retention_proven=False,
        fixed_width_carrier_proven=False,autonomous_readout_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

