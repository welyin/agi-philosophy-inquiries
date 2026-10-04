"""671 entry: actual static non-Abelian spatial curvature and physical Gram.

Necessary finite Gaussian diagnostic only, not full gauge/S9 reflection proof.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_nonflat_mass_measure as prior
mass=prior.mass;old=prior.old;internal=prior.internal


def run():
    sites=[(t,x,y) for t in range(2) for x in range(2) for y in range(2)]
    count=len(sites);n=64*count;r=n//2
    ux=prior.rep(*prior.group(67101,.13));uy=prior.rep(*prior.group(67102,.12))
    plaquette=ux@uy@ux.conj().T@uy.conj().T
    defect=old.norm(plaquette-np.eye(16));assert defect>.01
    dw=np.zeros((n,n),complex)
    for mu in range(3):
        shift=np.zeros((n,n),complex)
        for i,(t,x,y) in enumerate(sites):
            dest=((t+1)%2,x,y) if mu==2 else (t,(x+1)%2,y) if mu==0 else (t,x,(y+1)%2)
            j=sites.index(dest);sign=-1 if mu==2 and t==1 else 1
            link=np.eye(16) if mu==2 else ux if mu==0 else uy
            shift[64*i:64*i+64,64*j:64*j+64]=sign*np.kron(np.eye(4),link)
        gamma=np.kron(np.kron(np.eye(count),internal.spin.GAMMA[3 if mu==2 else mu]),np.eye(16))
        dw+=np.eye(n)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    g5=np.kron(np.kron(np.eye(count),internal.spin.G5),np.eye(16))
    h=g5@(dw-np.eye(n));ev,vec=np.linalg.eigh(h)
    assert old.err(h-h.conj().T)<2e-13 and min(abs(ev))>.5
    v=vec[:,ev>0];u=vec[:,ev<0];assert v.shape==u.shape==(n,r)
    d=(np.eye(n)+g5@(vec*np.sign(ev))@vec.conj().T)/2
    jm=np.kron(np.kron(np.eye(count),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(count),internal.VP),np.eye(16))
    kl=jp.conj().T@d@v;w=jm.conj().T@v;l=old.diag(w,np.eye(r))
    z=np.zeros_like(kl);n0=np.block([[z,-kl.T],[kl,z]])
    pair,_=mass.mass_pairing(count,mass.car.PHI[0]);pm=l.T@pair@l
    reflection=np.eye(count)[[sites.index((1-t,x,y)) for t,x,y in sites]]
    tr=np.kron(reflection,np.eye(32));theta=np.block([[np.zeros_like(tr),tr],[tr,np.zeros_like(tr)]])
    positive=np.flatnonzero(np.tile(np.repeat([t==1 for t,x,y in sites],32),2))
    rows=[]
    for lam in (0.,.37,1.):
        nw=n0+lam*pm;g=-l@np.linalg.solve(nw,l.T)
        gram=(theta@g)[np.ix_(positive,positive)]
        herm=old.err(gram-gram.conj().T)
        minimum=float(min(np.linalg.eigvalsh((gram+gram.conj().T)/2)))
        weight=internal.pfaffian(nw)/internal.pfaffian(n0)
        rows.append(dict(mass_scale=lam,physical_reflection_hermiticity_error=herm,
            physical_one_field_Gram_minimum=minimum,conditional_mass_weight_ratio=old.cpair(weight),
            physical_Nambu_minimum=float(min(np.linalg.svd(nw,compute_uv=False)))))
        assert herm<4e-12 and minimum>-4e-12
        assert weight.real>0 and abs(weight.imag)<3e-10
    # The compatible external field is fixed, not a dynamic gauge path integral.
    return dict(date='2026-10-02',entry_for_round=671,formal_round_complete=False,
        lattice=dict(nx=2,ny=2,nt=2),internal_channels=16,
        static_nonabelian_spatial_plaquette_defect=defect,Wilson_gap=float(min(abs(ev))),
        minimum_massless_Kl=float(min(np.linalg.svd(kl,compute_uv=False))),
        time_reflection_compatible_by_static_links_and_identity_temporal_transport=True,
        rows=rows,all_finite_diagnostics_passed=True,
        auxiliary_S9_joint_or_dynamic_gauge_RP_not_proved=True,
        original_Hamiltonian_or_instrument_identity_not_proved=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in
            ('joint_nonflat_mass_measure.py','joint_mass_auxiliary_reflection.py','research_note_670.md')})


if __name__=='__main__':
    result=run();target=HERE/'magnetic_reflection_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_finite_diagnostics_passed')}))
