"""679 executed entry: original Weyl source reflection vs finite soft regulator.
Checks reflection covariance, not positivity. No full Haar integral is sampled.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_rational_physical_limit as old
base=old.base
TARGET=HERE/'physical_reflection_probe_results.json'


def reflected(links,e,phi):
    lt=links[:,[2,3,0,1]].copy()
    lt[1]=np.array([a.conj().T for a in links[1]])
    return lt,e[[2,3,0,1]],phi[[2,3,0,1]]


def data(links,e,phi,mode,lam):
    mat=base.fixed_matrices(e,phi)
    u,v,d,h,gap=base.kernel(links)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    if mode=='exact':
        eps=v@v.conj().T-u@u.conj().T
    elif mode=='soft':
        eps,_,_=old.regulate(h,g5,.23,1)
    else:
        eps,_,_=old.regulate(h,g5,.001953125,4194304)
    q=old.soft(eps,mat,lam)
    q.update(mat=mat,gap=gap,weight=base.pf(q['N']))
    return q


def theta_matrix(mat):
    jm,jp=mat[:2]
    perm=np.eye(4)[[2,3,0,1]]
    gamma0=np.kron(perm,np.kron(base.internal.spin.GAMMA[3],np.eye(16)))
    c=jp.conj().T@gamma0@jm
    r=len(c)
    theta=np.block([[np.zeros((r,r)),c.T],[c,np.zeros((r,r))]])
    assert np.max(abs(theta.conj()@theta-np.eye(2*r)))<1e-13
    return theta


def evaluate(links,e,phi,mode,lam,z):
    q=data(links,e,phi,mode,lam)
    qt=data(*reflected(links,e,phi),mode,lam)
    theta=theta_matrix(q['mat'])
    minimum=float(min(np.linalg.svd(q['N'],compute_uv=False)[-1],
                      np.linalg.svd(qt['N'],compute_uv=False)[-1]))
    assert minimum>1e-8
    sources=[]
    for count in (2,4):
        zz=z[:,:count]
        # Theta is antilinear and reverses Grassmann product order.
        zt=theta.T@zz.conj()[:,::-1]
        c=old.source_coeff(q,zz);ct=old.source_coeff(qt,zt)
        sources.append(dict(count=count,original=base.old.cpair(c),
            reflected=base.old.cpair(ct),
            relative_error=float(abs(ct-c.conjugate())/max(abs(ct),abs(c))),
            absolute_error=float(abs(ct-c.conjugate()))))
    return dict(mode=mode,lambda_mass=lam,smallest_N_singular=minimum,
        weight=base.old.cpair(q['weight']),reflected_weight=base.old.cpair(qt['weight']),
        scalar_reflection_relative=float(abs(qt['weight']/q['weight'].conjugate()-1)),
        source_rows=sources)


def run():
    links,e,phi=old.fixture()
    rng=np.random.default_rng(67911)
    z=rng.normal(size=(256,4))+1j*rng.normal(size=(256,4))
    z/=np.linalg.norm(z,axis=0)
    free=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    calibration=evaluate(free,e,phi,'exact',0.,z)
    assert max(r['relative_error'] for r in calibration['source_rows'])<2e-9
    rows=[evaluate(links,e,phi,mode,.37,z) for mode in ('exact','soft','near')]
    return dict(date='2026-10-02',entry_round=679,not_formal_round=True,
        mature_free_Weyl_reflection_calibration=calibration,
        same_original_dynamic_background=rows,
        local_Weyl_Theta_from661_not_fitted_to_results=True,
        antilinearity_and_reversed_Grassmann_order_retained=True,
        no_reflection_positivity_claim=True,no_full_Haar_or_S9_average_evaluated=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_661.md','research_note_674.md',
                      'research_note_677.md','research_note_678.md',
                      'joint_rational_physical_limit.py')})


if __name__=='__main__':
    result=run()
    if '--check' in sys.argv:assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result))

