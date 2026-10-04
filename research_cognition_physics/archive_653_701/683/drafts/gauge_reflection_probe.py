"""Actual683 gauge dictionary and a specified compensation reflection test."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_source_lift as body
import joint_physical_source_reflection as reflected
prior=body.prior;base=body.base;err=body.err
TARGET=HERE/'gauge_reflection_probe_results.json'

def blockdiag(*blocks):
    size=sum(len(b) for b in blocks);out=np.zeros((size,size),complex);i=0
    for b in blocks:out[i:i+len(b),i:i+len(b)]=b;i+=len(b)
    return out

def observation(mat):
    jm,jp,M,_,_=mat;n=len(M);r=n//2
    out=np.zeros((n,4*n),complex)
    out[:r,:n]=.5*jm.T;out[r:,:n]=-.75*jp.T@M;out[r:,n:2*n]=jp.T
    t=np.vstack((.5*jm.T,.25*jp.T@M))
    out[:,2*n:3*n]=2*t;out[:,3*n:]=-2*t
    return out

def effective(h,mat,a,layers,lam):
    jm,jp,M,_,pair=mat;n=len(M);g5=jp@jp.T-jm@jm.T
    eps,_,_=prior.regulate(h,g5,a,layers);b=2*np.eye(n)+a*g5@h
    pull=np.zeros((4*n,2*n),complex);pull[:2*n]=np.eye(2*n)
    pull[2*n:3*n,:n]=np.linalg.solve(b,(np.eye(n)+eps)/2)
    pull[3*n:,:n]=np.linalg.solve(b,(np.eye(n)-eps)/2)
    phi=observation(mat)@pull;q=prior.soft(eps,mat,lam)
    nm=prior.soft(eps,mat,0)['N']+lam*phi.T@pair@phi
    return nm,phi,q

def run():
    links,e,phis=prior.fixture();mat=base.fixed_matrices(e,phis);n=len(mat[2])
    _,_,_,h,_=base.kernel(links);g5=mat[1]@mat[1].T-mat[0]@mat[0].T
    groups=[base.prior.group(68340+i,.37) for i in range(4)]
    lt,et,pt,rs=base.transform(links,e,phis,groups);mt=base.fixed_matrices(et,pt)
    _,_,_,ht,_=base.kernel(lt)
    u=blockdiag(*[np.kron(np.eye(4),r) for r in rs])
    v=blockdiag(*[np.kron(np.eye(2),r) for r in rs])
    uy=blockdiag(v,v.conj());ux=blockdiag(u,u.conj())
    ub=blockdiag(u,u.conj(),u,u)
    nm,p,q=effective(h,mat,.23,3,.37);nmt,ptt,qt=effective(ht,mt,.23,3,.37)
    gauge=dict(local_observation=err(observation(mt)@ub-uy@observation(mat)),
        integrated_observation=err(ptt@ux-uy@p),
        full_mass_action=err(ux.T@nmt@ux-nm),
        original_kernel=err(ux.T@qt['N']@ux-q['N']))
    assert max(gauge.values())<3e-12,gauge
    # Original time link reversal, not simply exchanging unrelated backgrounds.
    order=[2,3,0,1];perm=np.eye(4)[order]
    linkst=links[:,order].copy();linkst[1]=np.array([z.conj().T for z in links[1]])
    matt=base.fixed_matrices(e[order],phis[order]);_,_,_,hr,_=base.kernel(linkst)
    ell,ry,rx=reflected.reflection_matrices(mat,perm)
    scalar_p=np.kron(perm,np.eye(64));pm=np.kron(np.eye(2),scalar_p)
    rows=[]
    for a in (.2,.1):
        layers=int(np.ceil(4/a**2));nn,pp,qq=effective(h,mat,a,layers,.37)
        nt,ptt,qtt=effective(hr,matt,a,layers,.37)
        # No claim that either finite modified functional is reflection invariant.
        z=np.eye(n)[:,[66,93]];zt=ry.T@z.conj()[:,::-1]
        def coeff(N,Phi,Z):return base.pf(body.augmented(N,Phi,Z))
        original=coeff(qq['N'],qq['Phi'],z)
        original_t=coeff(qtt['N'],qtt['Phi'],zt)
        modified=coeff(nn,pp,z);modified_t=coeff(nt,ptt,zt)
        defect=abs(modified_t-modified.conjugate())
        triangle=abs(modified-original)+abs(modified_t-original_t)+abs(original_t-original.conjugate())
        assert defect<=triangle+1e-30
        # Actual L=1 compensator, kept separate from the large-L endpoint run.
        k,_,_=body.blocks(g5@h,g5,a,1);kr,_,_=body.blocks(g5@hr,g5,a,1)
        gram=k.conj().T@k;gramr=kr.conj().T@kr
        defect_matrix=pm.T@gramr@pm-gram
        eigen,vec=np.linalg.eigh(defect_matrix)
        i=int(np.argmax(abs(eigen)));w=vec[:,i]
        witness=np.vdot(w,defect_matrix@w)
        assert abs(witness)>.01 and abs(witness.imag)<1e-12
        assert np.linalg.eigvalsh(gram)[0]>0
        # L1 exact algebra of the actual blocks, rather than assuming Wilson-PV identity.
        b=2*np.eye(n)+a*g5@h
        common=2*b.conj().T@b+a*a*h@h
        predicted=np.block([[common-2*a*h@b,a*a*h@h],
                            [a*a*h@h,common+2*a*h@b]])
        assert err(predicted-gram)<3e-13
        rows.append(dict(a=a,effective_L=layers,original_source=base.old.cpair(original),
            original_reflected=base.old.cpair(original_t),modified_source=base.old.cpair(modified),
            modified_reflected=base.old.cpair(modified_t),finite_modified_reflection_defect=float(defect),
            reflection_triangle_bound=float(triangle),
            compensation_L=1,compensation_Gram_formula_error=err(predicted-gram),
            compensation_Gram_minimum_eigenvalue=float(np.linalg.eigvalsh(gram)[0]),
            declared_scalar_reflection_Gram_defect=float(np.linalg.norm(defect_matrix,2)),
            normalized_vector_action_difference=base.old.cpair(witness),
            witness_vector=[[float(z.real),float(z.imag)] for z in w]))
    return dict(entry_round=683,gauge_errors=gauge,rows=rows,
        scalar_reflection_is_declared_choice_not_all_possible_reflections=True,
        failure_is_unintegrated_compensation_action_not_full_Gauss_counterexample=True,
        no_actual_original_RP_verdict=True)

if __name__=='__main__':
    data=run()
    if TARGET.exists():assert data==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(all_checks_passed=True,gauge=data['gauge_errors'],
        reflection=[(r['finite_modified_reflection_defect'],r['declared_scalar_reflection_Gram_defect']) for r in data['rows']])) )
