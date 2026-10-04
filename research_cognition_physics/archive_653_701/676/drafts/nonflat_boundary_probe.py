"""676 entry: test659's actual boundary dictionary on original nonflat fields.

Fifth-direction transfer is not physical time. No RP conclusion.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
old=base.old
TARGET=HERE/'nonflat_boundary_probe_results.json'

def run():
    mat=base.fixed_matrices(np.tile(np.eye(10)[0],(4,1)),base.mass.car.PHI)
    jm,jp=mat[:2];c=np.column_stack((jp,jm));r=128;n=256
    g5=np.diag(np.r_[np.ones(r),-np.ones(r)])
    seed=np.r_[np.zeros((r,r)),np.eye(r)]
    rows=[]
    for number,scale in ((67531,.28),(67351,.9),(67350,.9)):
        links=np.array([[base.prior.rep(*base.prior.group(number+8*mu+i,scale)) for i in range(4)] for mu in range(2)])
        u,v,d,h,gap=base.kernel(links)
        hc=c.conj().T@h@c;x=g5@hc
        target=(c.conj().T@u)@(c.conj().T@u).conj().T
        for a in (.15,.05):
            aa=np.eye(n)+a*x;wp=aa[:r,:r];wm=aa[r:,r:];z=aa[:r,r:]
            assert old.err(aa[r:,:r]+z.conj().T)<1e-13
            winv=np.linalg.inv(wp)
            m=np.block([[winv,-winv@z],[np.zeros((r,r)),np.eye(r)]])
            t=np.block([[winv,-winv@z],[-z.conj().T@winv,wm+z.conj().T@winv@z]])
            ha=g5@(a*x)@np.linalg.inv(2*np.eye(n)+a*x)
            he,hu=np.linalg.eigh((ha+ha.conj().T)/2)
            minus=hu[:,he<0]
            transverse=np.linalg.svd(minus.conj().T@seed,compute_uv=False)
            reference=np.linalg.qr(np.linalg.solve(2*np.eye(n)+a*x,minus))[0]
            reference=reference@reference.conj().T
            errors=dict(Shamir_hermiticity=old.err(ha-ha.conj().T),
                boundary_identity=old.err((2*np.eye(n)+a*x)@m-np.eye(n)-t),
                chiral_identity=old.err(g5@(a*x)@m-np.eye(n)+t),
                Cayley_identity=old.err(ha-(np.eye(n)-t)@np.linalg.inv(np.eye(n)+t)))
            assert max(errors.values())<3e-13
            minwp=float(np.linalg.eigvalsh(wp)[0]);minwm=float(np.linalg.eigvalsh(wm)[0])
            assert min(minwp,minwm)>=1-a-2e-13
            mineig=float(np.linalg.eigvalsh(t)[0]);assert mineig>0
            frame=seed.copy();levels=[]
            for level in range(1,97):
                if level>1:frame=np.linalg.qr(t@frame)[0]
                if level in (16,48,96):
                    boundary=np.linalg.qr(m@frame)[0]
                    p=boundary@boundary.conj().T
                    levels.append(dict(L=level,operator_error_to_mapped_limit=float(np.linalg.norm(p-reference,2)),
                        operator_error_to_original_Wilson=float(np.linalg.norm(p-target,2))))
            rows.append(dict(seed=number,gauge_scale=scale,fifth_step=a,Wilson_gap=gap,
                negative_Shamir_rank=minus.shape[1],minimum_Wplus=minwp,minimum_Wminus=minwm,
                minimum_fifth_transfer_eigenvalue=mineig,
                seed_transversality_minimum=float(transverse[-1]),errors=errors,
                mapped_limit_original_error=float(np.linalg.norm(reference-target,2)),levels=levels))
    deps=('research_note_659.md','research_note_673.md','research_note_675.md','joint_gauss_boundary_functional.py')
    return dict(date='2026-10-02',entry_round=676,rows=rows,
        full_original_16_channels_and_nonflat_links=True,fifth_direction_not_physical_time=True,
        fixed_E_reference_normalization_not_used=True,no_RP_or_original_HF_identity=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})
if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=676,rows=len(result['rows']),
        min_seed_transversality=min(q['seed_transversality_minimum'] for q in result['rows']),
        fifth_direction_not_physical_time=True)))

