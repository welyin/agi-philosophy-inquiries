"""688 executed entry: actual temporal Wilson matrices and the same spherical kernel."""
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_full_holonomy_positive_measure as prior
old=prior.old;single=prior.entry
TARGET=HERE/'temporal_sphere_probe_results.json'


def rising(x,n):
    out=Fraction(1)
    for k in range(n):out*=x+k
    return out


def eigenvalues():
    mu0=rising(Fraction(9,2),8)/rising(Fraction(9),8)
    out=[]
    for ell in range(9):
        mu=mu0*Fraction(math.factorial(8),math.factorial(8-ell))/rising(Fraction(17),ell)
        dim=math.comb(9+ell,ell)-(math.comb(7+ell,ell-2) if ell>=2 else 0)
        out.append((ell,dim,mu))
    assert sum(dim*mu for ell,dim,mu in out)==1
    assert sum(dim*mu**2 for ell,dim,mu in out)==rising(Fraction(9,2),16)/rising(Fraction(9),16)
    return out


def harmonic_characters(O):
    coeff=[1.]
    power=np.eye(10);traces=[0.]
    for k in range(1,9):
        power=power@O;traces.append(float(np.trace(power)))
        coeff.append(sum(traces[j]*coeff[k-j] for j in range(1,k+1))/k)
    return np.array([c-(coeff[k-2] if k>=2 else 0) for k,c in enumerate(coeff)])


def full_sphere_trace(holonomy,length):
    chars=harmonic_characters(single.rotation(holonomy))
    return float(sum(float(mu**length)*chars[ell] for ell,dim,mu in eigenvalues()))


def temporal_case(length,seed):
    rng=np.random.default_rng(seed)
    links=[single.gauge.rep(*single.gauge.group(seed+i,.21+.07*i)) for i in range(length)]
    E=rng.normal(size=(length,10));E/=np.linalg.norm(E,axis=1)[:,None]
    size=16*length
    S=np.zeros((size,size),complex)
    for i,r in enumerate(links):
        j=(i+1)%length;sgn=-1 if i==length-1 else 1
        S[16*i:16*i+16,16*j:16*j+16]=sgn*r
    R=-S
    eye=np.eye(size);a=(eye-R.conj().T)/2;b=(eye+R.conj().T)/2
    u=1j*(np.kron(a,old.VP)+np.kron(b,old.spin.GAMMA[3]@old.VP))
    v=1j*(np.kron(a,old.VM)+np.kron(b,old.spin.GAMMA[3]@old.VM))
    c=(R+R.conj().T)/2;s=(R-R.conj().T)/(2j)
    X=np.kron(c,np.eye(4))-1j*np.kron(s,old.spin.GAMMA[3])
    H=np.kron(eye,old.spin.G5)@X;D=(np.eye(4*size)+X)/2
    assert np.linalg.norm(H@H-np.eye(4*size),2)<3e-13
    assert np.linalg.norm(u@u.conj().T-(np.eye(4*size)-H)/2,2)<3e-13
    Tint=np.zeros((size,size),complex)
    expected=np.zeros((2*size,2*size),complex)
    pf_product=1.
    hol=np.eye(16,dtype=complex)
    for i,r in enumerate(links):
        hol=hol@r
        ti=sum(x*t for x,t in zip(E[i],old.T))
        Tint[16*i:16*i+16,16*i:16*i+16]=ti
        avg=(E[i]+single.rotation(r)@E[(i+1)%length])/2
        et=-np.kron(sum(x*t for x,t in zip(avg,old.T)),old.EPS)
        expected[32*i:32*i+32,32*i:32*i+32]=et
        pf_product*=float((avg@avg)**8)
    actual=u.T@np.kron(Tint,old.B)@u
    error=float(np.linalg.norm(actual-expected,2))
    assert error<3e-13
    pf=old.pfaffian(actual)
    relative_pf=float(abs(pf-pf_product)/pf_product)
    assert relative_pf<2e-11
    bar=np.kron(eye,old.VP.conj().T)
    phase,logW=np.linalg.slogdet(bar@D@v)
    _,logHol=np.linalg.slogdet(np.eye(16)+hol)
    predicted=2*logHol-32*length*np.log(2)
    assert abs(phase-1)<3e-12 and abs(logW-predicted)<3e-11
    trace=full_sphere_trace(hol,length)
    assert trace>0
    if length==1:
        assert abs(trace-single.moment(single.rotation(hol)))<3e-14
    return dict(time_slices=length,full_spin_internal_dimension=4*size,
        original_pairing_error=error,relative_Pfaffian_product_error=relative_pf,
        physical_logdet_error=float(abs(logW-predicted)),
        retained_log_normalization_per_slice=-32*np.log(2),
        full_S9_transfer_trace=trace,raw_physical_weight_log=float(logW),
        independent_temporal_links=True,all_plaquettes_spatially_trivial=True)


def run():
    modes=eigenvalues()
    rows=[temporal_case(n,68800+100*n) for n in (1,2,3,4)]
    return dict(entry_round=688,latest_formal_round=687,not_formal_round=True,
        sphere_dimension=9,kernel_power=8,
        transfer_modes=[dict(ell=ell,multiplicity=dim,eigenvalue=str(mu)) for ell,dim,mu in modes],
        rank=sum(dim for ell,dim,mu in modes),
        exact_trace_at_length1='1',
        exact_trace_at_identity_length2=str(sum(dim*mu**2 for ell,dim,mu in modes)),
        matrix_rows=rows,
        original_full_space_Hb_process_not_identified=True,
        original_physical_Y_source_time_dictionary_not_proved=True)


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=688,all_checks_passed=True,rank=result['rank'],
        max_pairing_error=max(r['original_pairing_error'] for r in result['matrix_rows']),
        two_slice_sphere_trace=result['exact_trace_at_identity_length2'])))

