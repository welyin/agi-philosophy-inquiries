"""672 entry: actual reflection-compatible time-dependent original gauge links.

Necessary conditional diagnostics; no claim about integrated Gauss positivity.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_static_gauge_reflection as prior
old=prior.old;internal=prior.internal;mass=prior.mass


def comm(a,b):return a@b-b@a


def data(amplitude):
    sites=[(t,x) for t in range(4) for x in range(2)];count=len(sites);n=64*count
    base=prior.prior.rep(*prior.prior.group(67201,.08))
    extra=prior.prior.rep(*prior.prior.group(67202,amplitude))
    links=[base,extra@base,extra@base,base]
    tx=np.zeros((n,n),complex);tt=np.zeros_like(tx)
    for i,(t,x) in enumerate(sites):
        j=sites.index((t,1-x));tx[64*i:64*i+64,64*j:64*j+64]=np.kron(np.eye(4),links[t])
        j=sites.index(((t+1)%4,x));tt[64*i:64*i+64,64*j:64*j+64]=(-1 if t==3 else 1)*np.eye(64)
    g0=np.kron(np.kron(np.eye(count),internal.spin.GAMMA[3]),np.eye(16))
    g1=np.kron(np.kron(np.eye(count),internal.spin.GAMMA[0]),np.eye(16))
    g5=np.kron(np.kron(np.eye(count),internal.spin.G5),np.eye(16))
    b=np.eye(n)-(tx+tx.conj().T)/2;c=g1@(tx-tx.conj().T)/2
    s=(tt+tt.conj().T)/2;a=(tt-tt.conj().T)/2
    x=b+c-s+g0@a;hs=g0@(b+c)
    baseq=np.eye(n)+hs@hs-b@s-s@b
    electric=-comm(s,c)+g0@comm(b+c,a)
    actual=x.conj().T@x
    error=old.err(actual-baseq-electric);assert error<4e-13
    h=g5@x;ev,vec=np.linalg.eigh(h);assert min(abs(ev))>.5
    d=(np.eye(n)+g5@(vec*np.sign(ev))@vec.conj().T)/2
    pos=np.arange(n//2)
    refl=np.concatenate([np.arange(64*sites.index((3-t,k)),64*(sites.index((3-t,k))+1)) for t,k in sites[:count//2]])
    cross=-g0[:n//2,:n//2]@d[np.ix_(refl,pos)]
    herm=old.err(cross-cross.conj().T);assert herm<4e-12
    action_values=np.linalg.eigvalsh((cross+cross.conj().T)/2)
    v=vec[:,ev>0];r=v.shape[1]
    jm=np.kron(np.kron(np.eye(count),internal.VM),np.eye(16));jp=np.kron(np.kron(np.eye(count),internal.VP),np.eye(16))
    kl=jp.conj().T@d@v;w=jm.conj().T@v;l=old.diag(w,np.eye(r));z=np.zeros_like(kl)
    n0=np.block([[z,-kl.T],[kl,z]])
    pair,_=mass.mass_pairing(count,mass.car.PHI[0])
    reflect=np.eye(count)[[sites.index((3-t,k)) for t,k in sites]]
    tr=np.kron(reflect,np.eye(32));theta=np.block([[np.zeros_like(tr),tr],[tr,np.zeros_like(tr)]])
    idx=np.flatnonzero(np.tile(np.repeat([t>=2 for t,k in sites],32),2))
    rows=[]
    for lam in (0.,.37):
        nw=n0+lam*l.T@pair@l;g=-l@np.linalg.solve(nw,l.T)
        gram=(theta@g)[np.ix_(idx,idx)];gh=old.err(gram-gram.conj().T);assert gh<4e-12
        values,vectors=np.linalg.eigh((gram+gram.conj().T)/2)
        witness=vectors[:,0];rayleigh=np.vdot(witness,gram@witness)
        rows.append(dict(mass_scale=lam,minimum_physical_Gram=float(values[0]),
            hermiticity_error=gh,witness_rayleigh=old.cpair(rayleigh),
            witness_real=witness.real.tolist(),witness_imag=witness.imag.tolist(),
            eigen_residual=old.norm(gram@witness-values[0]*witness)))
    return dict(amplitude=amplitude,electric_plaquette_defect=old.norm(links[0]@links[1].conj().T-np.eye(16)),
        Wilson_gap=float(min(abs(ev))),exact_dynamic_pencil_error=error,
        omitted_electric_terms_norm=old.norm(electric),
        B_time_commutator_norm=old.norm(comm(b,s)),
        action_reflection_minimum=float(action_values[0]),action_reflection_hermiticity_error=herm,
        physical_rows=rows)


def run():
    static=data(0.);dynamic=data(.07)
    assert static['omitted_electric_terms_norm']<1e-13
    assert dynamic['omitted_electric_terms_norm']>.1 and dynamic['electric_plaquette_defect']>.01
    return dict(date='2026-10-02',entry_for_round=672,formal_round_complete=False,
        original_SM_channels=16,nx=2,nt=4,time_reflection_compatible_links=True,
        static_control=static,dynamic_background=dynamic,
        all_algebra_checks_passed=True,
        conditional_Gram_signs_are_results_not_assumed_positive=True,
        full_dynamic_gauge_or_Gauss_claim_not_made=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in
            ('joint_static_gauge_reflection.py','joint_nonflat_mass_measure.py','research_note_671.md')})


if __name__=='__main__':
    result=run();target=HERE/'dynamic_commutator_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps(dict(entry=672,algebra_passed=True,
        dynamic_physical_minima=[x['minimum_physical_Gram'] for x in result['dynamic_background']['physical_rows']]),ensure_ascii=False))
