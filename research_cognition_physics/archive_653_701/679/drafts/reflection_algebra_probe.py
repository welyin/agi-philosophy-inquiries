"""679 algebra diagnostic: isolate which reflection assumptions are needed."""
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_rational_physical_limit as old
base=old.base
def er(a):return float(np.max(abs(a),initial=0))
def run():
    rng=np.random.default_rng(67922)
    e=rng.normal(size=(1,10));e/=np.linalg.norm(e,axis=1)[:,None]
    mat=base.fixed_matrices(e,base.mass.car.PHI[:1])
    jm,jp,m,mb,pair=mat;n=64;r=32
    g0=np.kron(base.internal.spin.GAMMA[3],np.eye(16))
    g5=np.kron(base.internal.spin.G5,np.eye(16));t=g5@g0
    c=jp.conj().T@g0@jm
    rp=np.block([[np.zeros((r,r)),c.T],[c,np.zeros((r,r))]])
    bare=np.block([[np.zeros((n,n)),g0.T],[g0,np.zeros((n,n))]])
    z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));v,_=np.linalg.qr(z)
    zz=rng.normal(size=(n,4))+1j*rng.normal(size=(n,4));zz/=np.linalg.norm(zz,axis=0)
    rows=[]
    for label,ev in [('soft',np.linspace(-.8,.9,n)),('projector',np.r_[-np.ones(30),np.ones(34)])]:
        eps=(v*ev)@v.conj().T;epst=-t@eps@t.conj().T
        q=old.soft(eps,mat);qt=old.soft(epst,mat)
        q0=old.soft(eps,mat,0);qt0=old.soft(epst,mat,0)
        wf=base.pf(q['N']);wt=base.pf(qt['N'])
        sr=[]
        for k in (2,4):
            a=old.source_coeff(q,zz[:,:k])
            b=old.source_coeff(qt,rp.T@zz[:,:k].conj()[:,::-1])
            sr.append(float(abs(b-a.conjugate())/max(abs(a),abs(b))))
        rows.append(dict(kind=label,scalar=float(abs(wt/wf.conjugate()-1)),
            sources=sr,bare_N0_congruence=er(bare.T@qt0['N']@bare+q0['N'].conj()),
            bare_N_congruence=er(bare.T@qt['N']@bare+q['N'].conj()),
            physical_bare_map=er(qt['Phi']@bare-rp@q['Phi'].conj())))
    return dict(B=base.internal.B.tolist() if False else None,
        identities=dict(Mbar_minus_Mstar=er(mb-m.conj()),MbarM_plus_I=er(mb@m+np.eye(n)),
        Mbar_minus_Mdagger=er(mb-m.conj().T),M_reflect=er(t.T@m@t-m)),
        rows=rows)
if __name__=='__main__':
    result=run()
    with (HERE/'reflection_algebra_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result))

