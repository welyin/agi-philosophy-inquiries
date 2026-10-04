"""680 executed entry: original half-mass algebra and nilpotent source insertions.
No reflection-positivity conclusion and no sampling of the full Haar integral.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_physical_source_reflection as previous
old=previous.old;base=previous.base
TARGET=HERE/'mass_congruence_probe_results.json'

def error(a):return float(np.max(abs(a),initial=0))

def run():
    links,e,phis=old.fixture();mat=base.fixed_matrices(e,phis)
    u,v,_,h,_=base.kernel(links)
    eps=v@v.conj().T-u@u.conj().T
    q=old.soft(eps,mat,.37)
    part=np.r_[np.repeat([0,0,1,1],32),np.repeat([0,0,1,1],32)]
    pos=np.diag(part);neg=np.eye(256)-pos
    pp=pos@mat[-1]@pos;pm=neg@mat[-1]@neg
    perm=np.eye(4)[[2,3,0,1]]
    _,ry,_=previous.reflection_matrices(mat,perm)
    matt=base.fixed_matrices(e[[2,3,0,1]],phis[[2,3,0,1]])
    ppt=pos@matt[-1]@pos
    groups=[base.prior.group(68010+i,.27) for i in range(4)]
    _,et,pt,rs=base.transform(links,e,phis,groups)
    mt=base.fixed_matrices(et,pt)
    blocks=[np.kron(np.eye(2),r) for r in rs]
    r2=base.old.diag(base.old.diag(*blocks[:2]),base.old.diag(*blocks[2:]))
    gy=base.old.diag(r2,r2.conj())
    errors=dict(full_mass_split=error(pp+pm-mat[-1]),
        physical_half_gauge_support=error(gy@pos-pos@gy),
        half_mass_gauge_covariance=error(gy.T@(pos@mt[-1]@pos)@gy-pp),
        half_mass_reflection=error(ppt+ry.T@pm.conj()@ry))
    assert max(errors.values())<3e-12
    rng=np.random.default_rng(68021)
    z=rng.normal(size=(256,2))+1j*rng.normal(size=(256,2))
    z/=np.linalg.norm(z,axis=0)
    rows=[]
    for side,pair in [('positive',pp),('negative',pm)]:
        i,j=np.unravel_index(np.argmax(abs(np.triu(pair,1))),pair.shape)
        alpha=.37*pair[i,j]
        select=np.eye(256,dtype=complex)[:,[i,j]]
        local=np.outer(select[:,0],select[:,1])-np.outer(select[:,1],select[:,0])
        delta=alpha*q['Phi'].T@local@q['Phi']
        before=dict(Phi=q['Phi'],N=q['N']-delta)
        after=dict(Phi=q['Phi'],N=before['N']+delta)
        assert abs(alpha)>1e-5
        smin=float(min(np.linalg.svd(before['N'],compute_uv=False)[-1],
                       np.linalg.svd(after['N'],compute_uv=False)[-1]))
        assert smin>1e-6
        checks=[]
        for k in (0,2):
            zz=z[:,:k]
            a=base.pf(before['N']) if k==0 else old.source_coeff(before,zz)
            b=base.pf(after['N']) if k==0 else old.source_coeff(after,zz)
            extra=np.column_stack((zz,select))
            c0=old.source_coeff(before,extra);c1=old.source_coeff(after,extra)
            predicted=a-alpha*c0
            inverse=b+alpha*c1
            scale=max(abs(a),abs(b),abs(alpha*c0),1e-280)
            errs=dict(forward=float(abs(b-predicted)/scale),
                inverse=float(abs(a-inverse)/scale),
                nilpotent_insertion_unchanged=float(abs(alpha*(c1-c0))/scale))
            assert max(errs.values())<1e-9,errs
            checks.append(dict(external_sources=k,errors=errs))
        rows.append(dict(side=side,physical_indices=[int(i),int(j)],
            original_mass_coefficient=base.old.cpair(alpha),
            minimum_kernel_singular=smin,checks=checks))
    return dict(date='2026-10-02',entry_round=680,not_formal_round=True,
        same_original_nonflat_full_group=True,original_16_channels=True,
        physical_half_mass_contract=errors,original_bilinear_rows=rows,
        augmented_Pfaffian_insertion_sign_minus=True,
        original_full_mass_not_replaced_by_selected_terms=True,
        selected_terms_only_test_exact_factorization_steps=True,
        complete_RP_and_nonzero_normalization_not_established=True,
        full_Haar_and_S9_not_sampled=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_668.md','research_note_669.md','research_note_673.md',
                      'research_note_679.md','joint_physical_source_reflection.py')})

if __name__=='__main__':
    result=run()
    if '--check' in sys.argv:assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result))
