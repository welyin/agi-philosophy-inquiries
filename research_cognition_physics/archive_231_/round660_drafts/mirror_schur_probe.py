"""660 entry, not a completed round: actual mirror Schur matching.

Inspired by the readable expression on page19 of the original2026 slides.
Couplings y,z are defined by the Nambu matrix below; relative action factors
and reflection/Hamiltonian assumptions remain to be checked independently.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_spatial_auxiliary_geometry as old


def op(a):return float(np.linalg.norm(a,2))


def matrices(kind,value):
    if kind=='spatial':
        f=old.frame(value,2,2);count=4
        u=np.kron(f['V'],np.eye(16))
        g=np.kron(np.kron(np.eye(count),old.old.spin.G5),np.eye(16))
        j=np.kron(np.kron(np.eye(count),old.old.VM),np.eye(16))
        sign=(f['vec']*np.sign(f['ev']))@f['vec'].conj().T
        spin_d=(np.eye(4*count)+np.kron(np.eye(count),old.old.spin.G5)@sign)/2
        d=np.kron(spin_d,np.eye(16))
        e=np.random.default_rng(66001).normal(size=(count,10))*.08;e[:,0]=1
        e/=np.linalg.norm(e,axis=1)[:,None]
        m=np.zeros((64*count,64*count),complex);bar=m.copy()
        for i,ei in enumerate(e):
            t=sum(a*b for a,b in zip(ei,old.old.T));sl=slice(64*i,64*i+64)
            m[sl,sl]=np.kron(old.old.B,t);bar[sl,sl]=np.kron(old.old.B,t.conj().T)
    else:
        # Original615 charged flat holonomy, internal-first spin ordering.
        u,_,d,_=old.old.frames(value);g=np.kron(np.eye(16),old.old.spin.G5)
        j=np.kron(np.eye(16),old.old.VM)
        e=np.random.default_rng(66002).normal(size=10);e/=np.linalg.norm(e)
        t=sum(a*b for a,b in zip(e,old.old.T))
        m=np.kron(t,old.old.B);bar=np.kron(t.conj().T,old.old.B)
    qp=(np.eye(len(g))+g)/2;qm=np.eye(len(g))-qp
    k=j.conj().T@d@u;ap=u.T@m@qp@u;am=u.T@m@qm@u
    bb=j.conj().T@bar@j.conj()
    return u,d,qp,qm,k,ap,am,bb


def run():
    rows=[]
    for kind,value in (('spatial',1.),('holonomy',.17),('holonomy',.31)):
        u,d,qp,qm,k,ap,am,bb=matrices(kind,value)
        schur=k.T@np.linalg.solve(bb,k)
        error=max(op(d@u-qm@u),op(schur-am))
        assert error<3e-12
        pairs=[]
        for y,z in ((1.,1.),(.9,.6)):
            full=np.block([[y*ap,-z*k.T],[z*k,y*bb]])
            eff=y*ap+z*z/y*am
            actual=old.old.pfaffian(full)
            expected=old.old.pfaffian(y*bb)*old.old.pfaffian(eff)
            pf_error=float(abs(actual-expected));assert pf_error<2e-12
            # Preserving b-source insertions is stronger than a determinant identity.
            covariance=np.linalg.inv(full)[:len(ap),:len(ap)]
            covariance_error=op(covariance-np.linalg.inv(eff))
            assert covariance_error<2e-10
            pairs.append(dict(y=y,z=z,full_Pfaffian=[float(actual.real),float(actual.imag)],
                Schur_Pfaffian_error=pf_error,b_source_covariance_error=covariance_error,
                original_pairing_difference=op(eff-ap-am)))
        assert pairs[0]['original_pairing_difference']<1e-12
        assert pairs[1]['original_pairing_difference']>.05
        rows.append(dict(branch=kind,parameter=value,original_Weyl_and_Schur_error=error,
            bar_Pfaffian=[float(old.old.pfaffian(bb).real),float(old.old.pfaffian(bb).imag)],pairs=pairs))
    return dict(date='2026-10-02',status='660 entry; not a complete round',rows=rows,
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
            ('joint_spatial_auxiliary_geometry.py','joint_subgroup_measure_source.py','research_note_613.md')},
        source_receipt_sha256=hashlib.sha256((HERE/'source_receipt.json').read_bytes()).hexdigest(),
        scope='Exact original finite Nambu/Schur weight and retained b-source identity at matched couplings. No full physical-time reflection theorem, original32CAR map or gauge-interacting Hamiltonian established.')


if __name__=='__main__':
    r=run()
    with (HERE/'mirror_schur_probe_results.json').open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=r['status'],branches=len(r['rows']))))
