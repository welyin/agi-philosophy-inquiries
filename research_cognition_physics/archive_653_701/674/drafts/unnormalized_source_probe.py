"""674 entry: prepare Gauss-integrand sources without dividing near-zero weights.

Reuse673 identity and actual strong-group configurations. This is not a new
formal round, an exact-rank proof, or evaluation of the boundary integral.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_gauss_boundary_functional as m


def run():
    e=np.random.default_rng(67321).normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    phis=m.mass.car.PHI.copy();mat=m.fixed_matrices(e,phis)
    rng=np.random.default_rng(67401)
    source=rng.normal(size=(256,4))+1j*rng.normal(size=(256,4))
    source/=np.linalg.norm(source,axis=0)
    rows=[]
    for seed in (67350,67351,67352):
        links=np.array([[m.prior.rep(*m.prior.group(seed+8*mu+i,.9)) for i in range(4)] for mu in range(2)])
        u,v,d,h,gap=m.kernel(links);q=m.regular(u,v,d,mat)
        _,values,vh=np.linalg.svd(q['N'])
        near=values<1e-10;null=vh.conj().T[:,near]
        mapped=m.old.norm(q['Phi']@null)
        a=u.T@mat[2]@u
        asv=np.linalg.svd(a,compute_uv=False)
        coeff=[]
        for k in (2,4):
            c=q['Phi'].T@source[:,:k]
            aug=np.block([[q['N'],c],[-c.T,np.zeros((k,k))]])
            coeff.append(m.old.cpair(m.pf(aug)))
        rows.append(dict(seed=seed,Wilson_gap=gap,negative_rank=u.shape[1],positive_rank=v.shape[1],
            near_null_full_count=int(sum(near)),near_null_auxiliary_count=int(sum(asv<1e-10)),
            physical_image_of_near_null_subspace_norm=mapped,
            auxiliary_smallest_singular_values=asv[-6:].tolist(),
            original_weight=m.old.cpair(q['weight']),unnormalized_two_and_four_source_coefficients=coeff))
    return dict(date='2026-10-02',entry_for_round=674,formal_round_complete=False,
        rows=rows,normalized_inverse_or_weight_division_not_used=True,
        exact_rank_or_exact_zero_not_claimed_from_svd=True,
        complete_Haar_integral_not_computed=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in
            ('joint_gauss_boundary_functional.py','research_note_673.md')})


if __name__=='__main__':
    result=run();target=HERE/'unnormalized_source_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
