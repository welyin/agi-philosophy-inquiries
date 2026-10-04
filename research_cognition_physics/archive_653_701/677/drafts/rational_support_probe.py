"""677 entry: whole-spectrum bounded Shamir sign on676's exact failed-seed case."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
old=base.old
TARGET=HERE/'rational_support_probe_results.json'

def finite(eps,mat,lam=.37):
    jm,jp,m,mb,pair=mat;n=len(eps);r=n//2
    g5=jp@jp.conj().T-jm@jm.conj().T
    pu=(np.eye(n)-eps)/2;pv=(np.eye(n)+eps)/2
    d=(np.eye(n)+g5@eps)/2
    phi=np.block([[jm.conj().T@pv,np.zeros((r,n))],
                  [-jp.T@m@(pu+.5*pv),jp.T]])
    nc=np.block([[m@(jp@jp.conj().T),-d.T],[d,mb@(jm@jm.conj().T)]])
    full=nc+lam*phi.T@pair@phi
    return dict(N=full,Phi=phi,weight=base.pf(full))

def run():
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
            links[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
    u,v,d,h,gap=base.kernel(links)
    rng=np.random.default_rng(67721);e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    mat=base.fixed_matrices(e,base.mass.car.PHI)
    exact_sign=v@v.conj().T-u@u.conj().T;exact=finite(exact_sign,mat)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    x=g5@h;n=len(h);rows=[]
    for a,L in ((.15,32),(.15,128),(.15,512),(.05,512),(.02,1024),(.005,4096)):
        ha=g5@(a*x)@np.linalg.inv(2*np.eye(n)+a*x)
        he,hv=np.linalg.eigh((ha+ha.conj().T)/2)
        assert max(abs(he))<1
        coeff=np.tanh(L*np.arctanh(he))
        eps=(hv*coeff)@hv.conj().T
        limit=(hv*np.sign(he))@hv.conj().T
        error_to_shamir=float(np.linalg.norm(eps-limit,2))
        bound=float(2*np.exp(-2*L*min(abs(np.arctanh(he)))))
        assert error_to_shamir<=bound+2e-13
        q=finite(eps,mat)
        rows.append(dict(a=a,L=L,epsilon_norm=float(max(abs(coeff))),
            error_to_Shamir_sign=error_to_shamir,finite_L_bound=bound,
            error_to_original_sign=float(np.linalg.norm(eps-exact_sign,2)),
            physical_map_error=float(np.linalg.norm(q['Phi']-exact['Phi'],2)),
            full_N_error=float(np.linalg.norm(q['N']-exact['N'],2)),
            raw_finite_weight=old.cpair(q['weight']),raw_exact_weight_numerical=old.cpair(exact['weight']),
            no_ratio_to_exact_zero_weight=True,
            finite_N_nullity_diagnostic=int(sum(np.linalg.svd(q['N'],compute_uv=False)<1e-10))))
    return dict(date='2026-10-02',entry_round=677,exact_failed_seed_flux=True,
        original_Wilson_gap=gap,rows=rows,
        finite_soft_projectors_not_claimed_exact_chiral_projectors=True,
        no_finite_RP_or_local_bulk_identity=True,no_physical_time_identification=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_673.md','research_note_676.md','joint_gauss_boundary_functional.py')})
if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps([dict(a=r['a'],L=r['L'],sign_error=r['error_to_original_sign'],
        Phi_error=r['physical_map_error']) for r in result['rows']]))

