"""676 follow-up: unused-spin symmetry and the fixed boundary seed rank."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
old=base.old
TARGET=HERE/'seed_sector_probe_results.json'

def run():
    sp=base.internal.spin
    ss=1j*sp.GAMMA[1]@sp.GAMMA[2]
    se,sv=np.linalg.eigh(ss)
    s=np.kron(np.kron(np.eye(4),ss),np.eye(16))
    frames=[np.kron(np.kron(np.eye(4),sv[:,se*sign>0]),np.eye(16)) for sign in (-1,1)]
    jm=np.kron(np.kron(np.eye(4),base.internal.VM),np.eye(16))
    g5=np.kron(np.kron(np.eye(4),sp.G5),np.eye(16))
    rows=[]
    for seed,scale in ((67531,.28),(67351,.9),(67350,.9)):
        links=np.array([[base.prior.rep(*base.prior.group(seed+8*mu+i,scale)) for i in range(4)] for mu in range(2)])
        u,v,d,h,gap=base.kernel(links)
        comm=old.err(s@h-h@s);assert comm<2e-13
        sectors=[]
        for sign,b in zip((-1,1),frames):
            hh=b.conj().T@h@b
            ev,evec=np.linalg.eigh(hh)
            residual=float(np.linalg.norm(hh@evec-evec*ev,2))
            projector=b@b.conj().T
            seed_dim=int(round(np.trace(jm.conj().T@projector@jm).real))
            sectors.append(dict(spin_sign=sign,sector_dimension=len(ev),
                negative_Wilson_count=int(sum(ev<0)),seed_dimension=seed_dim,
                sector_gap=float(min(abs(ev))),diagonalization_residual=residual))
        forced_loss=sum(max(q['negative_Wilson_count']-q['seed_dimension'],0) for q in sectors)
        transverse=np.linalg.svd(u.conj().T@jm,compute_uv=False)
        rows.append(dict(seed=seed,scale=scale,commutator_error=comm,sectors=sectors,
            algebraic_seed_rank_loss_if_verified_counts=forced_loss,
            measured_seed_nullity=int(sum(transverse<1e-10)),
            smallest_seed_singular_values=transverse[-4:].tolist(),
            counts_are_numeric_not_interval_certified=True))
    return dict(date='2026-10-02',entry_round=676,rows=rows,
        theorem='If S commutes with H and the fixed seed projector, then rank(Uminus^* Jminus) <= sum_s min(negative_rank_s,seed_rank_s). A sector rank excess is an exact seed obstruction.',
        existence_of_original_counterexample_still_requires_certified_counts=True,
        whole_domainwall_route_or_original_HF_not_refuted=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
          for p in ('research_note_659.md','research_note_673.md','joint_gauss_boundary_functional.py')})

if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps([dict(seed=q['seed'],counts=[s['negative_Wilson_count'] for s in q['sectors']],
        forced_seed_loss=q['algebraic_seed_rank_loss_if_verified_counts']) for q in result['rows']]))

