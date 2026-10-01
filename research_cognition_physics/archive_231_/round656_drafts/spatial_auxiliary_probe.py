"""656 entry: restore actual spatial Wilson terms in the old auxiliary Pfaffian.

Absolute Pfaffian is frame independent. A pointwise mismatch is not yet a
theorem about the integrated measure or all positive transfer representations.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_subgroup_measure_source as old


def geometry(nx=3,nt=2):
    sites=[(t,x) for t in range(nt) for x in range(nx)];size=len(sites)
    st=np.zeros((size,size));sx=np.zeros_like(st)
    for i,(t,x) in enumerate(sites):
        st[i,sites.index(((t+1)%nt,x))]=-1 if t==nt-1 else 1
        sx[i,sites.index((t,(x+1)%nx))]=1
    return sites,st,sx


def frame(strength,nx=3,nt=2):
    sites,st,sx=geometry(nx,nt);size=len(sites)
    x=-np.eye(4*size,dtype=complex)
    for shift,gamma,weight in ((st,old.spin.GAMMA[3],1.),(sx,old.spin.GAMMA[0],strength)):
        x+=weight*(np.kron(np.eye(size)-(shift+shift.T)/2,np.eye(4))
                   +np.kron((shift-shift.T)/2,gamma))
    g5=np.kron(np.eye(size),old.spin.G5);h=g5@x
    eig,v=np.linalg.eigh(h);gap=float(min(abs(eig)))
    assert gap>.4 and np.max(abs(h-h.conj().T))<1e-13
    u=v[:,eig<0];assert u.shape==(4*size,2*size)
    return np.kron(u,np.eye(16)),gap


def value(u,scalars):
    n=len(scalars);q=np.zeros((64*n,64*n),complex)
    for j,e in enumerate(scalars):
        te=sum(a*t for a,t in zip(e,old.T))
        q[64*j:64*(j+1),64*j:64*(j+1)]=np.kron(old.B,te)
    a=u.T@q@u
    return float(abs(old.pfaffian(a)))


def copied_time_kernel(e,nx=3,nt=2):
    out=1.;sites,_,_=geometry(nx,nt)
    for j,(t,x) in enumerate(sites):
        neighbor=sites.index(((t+1)%nt,x))
        out*=((1+e[j]@e[neighbor])/2)**8
    return float(out)


def run():
    rng=np.random.default_rng(656);base=np.eye(10)[0]
    e=base+.22*rng.normal(size=(6,10));e/=np.linalg.norm(e,axis=1)[:,None]
    factor=copied_time_kernel(e);rows=[]
    for strength in (0.,.25,.5,1.):
        u,gap=frame(strength);actual=value(u,e)
        relative=abs(actual-factor)/factor
        if strength==0:assert relative<2e-11
        rows.append(dict(spatial_Wilson_strength=strength,Wilson_gap=gap,
            actual_auxiliary_Pfaffian_absolute=actual,copied_time_kernel=factor,
            relative_difference=relative))
    assert rows[-1]['relative_difference']>1e-3
    return dict(date='2026-10-02',status='656 entry; incomplete',rows=rows,
        scalar_configuration=e.tolist(),spatial_sites=3,antiperiodic_time_sites=2,
        other_two_spatial_directions_one_site=True,original_m0=1,
        original_point_strength=1.,intermediate_strengths_are_declared_interpolation=True,
        dependency_hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
            ('joint_subgroup_measure_source.py','joint_chiral_fibre_source.py','research_note_653.md')},
        scope='Actual free auxiliary integrand absolute value only; no phase-fixing or S9 integration. No claim of impossibility of a different spatial transfer kernel.')


if __name__=='__main__':
    r=run()
    with (HERE/'spatial_auxiliary_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=r['status'],rows=len(r['rows']),
        original_strength_relative_difference=r['rows'][-1]['relative_difference'])))
