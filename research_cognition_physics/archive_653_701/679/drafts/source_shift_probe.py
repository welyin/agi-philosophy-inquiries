"""679 exact-algebra candidate: equation-of-motion source shift."""
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_rational_physical_limit as old
base=old.base
def er(x):return float(np.max(abs(x),initial=0))
def run():
    rng=np.random.default_rng(67922)
    e=rng.normal(size=(1,10));e/=np.linalg.norm(e,axis=1)[:,None]
    mat=base.fixed_matrices(e,base.mass.car.PHI[:1]);jm,jp,m,mb,_=mat;n=64;r=32
    g0=np.kron(base.internal.spin.GAMMA[3],np.eye(16));g5=jp@jp.T-jm@jm.T
    t=g5@g0;c=jp.T@g0@jm
    rp=np.block([[np.zeros((r,r)),c.T],[c,np.zeros((r,r))]])
    rb=np.block([[np.zeros((n,n)),g0.T],[g0,np.zeros((n,n))]])
    s=np.block([[-.5*jm.T@m,jm.T],[jp.T,.5*jp.T@m.conj().T]])
    z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));v,_=np.linalg.qr(z)
    rows=[]
    for label,ev in [('soft',np.linspace(-.8,.9,n)),('projector',np.r_[-np.ones(30),np.ones(34)]),
                     ('singular',np.ones(n))]:
        eps=(v*ev)@v.conj().T;epst=-t@eps@t.conj().T
        q=old.soft(eps,mat,0);qt=old.soft(epst,mat,0)
        phi=q['Phi'].conj();nn=q['N'].conj();tilde=rp@qt['Phi']@rb
        diff=tilde-phi
        contact=s@phi.T-phi@s.T-s@nn@s.T
        errors=dict(source_shift=er(diff-s@nn),source_contact=er(contact),
            bare_reflection=er(rb.T@qt['N']@rb+q['N'].conj()),
            mass_reflection=er(rp.T@mat[-1].conj()@rp+mat[-1]))
        rows.append(dict(kind=label,errors=errors));assert max(errors.values())<1e-12
    return dict(rows=rows,only_gamma5_Hermiticity_not_idempotence_used=True)
if __name__=='__main__':
    result=run()
    with (HERE/'source_shift_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result))

