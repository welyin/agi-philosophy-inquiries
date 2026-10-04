"""682 probe on the ACTUAL678 finite auxiliary body.
Sigma is an antisymmetric Berezin inverse block, NOT a positive quantum noise.
The proposed eta-eta term is explicit and distinct from the original model.
"""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_source_lift as body
base=body.base
TARGET=HERE/'body_fluctuation_probe_results.json'
err=body.err

def small_aug(n,phi,contact,z):
    c=phi.T@z
    return np.block([[n,c],[-c.T,z.T@contact@z]])

def corrected_large_coefficient(matrix,k):
    phase,log=body.logpf(matrix);kd,kl=np.linalg.slogdet(k)
    return phase/kd*np.exp(log-kl)

def run():
    links,e,phis=body.prior.fixture();mat=base.fixed_matrices(e,phis)
    jm,jp,M,Mbar,pair=mat;n=len(M);r=n//2
    g5=jp@jp.T-jm@jm.T
    _,_,_,h,gap=base.kernel(links)
    T=np.vstack((.5*jm.T,.25*jp.T@M))
    sv=np.linalg.svd(T,compute_uv=False)
    assert err(sv[:r]-.5)<2e-14 and err(sv[r:]-.25)<2e-14
    rows=[]
    for label,lam in (('invisible',0.),('visible',0.),('invisible',.37)):
        lifted=body.lift(g5@h,mat,lam=lam,layers=1)
        k=lifted['K'];O=lifted['out'];m=len(k);delta=.09
        V=np.vstack((np.eye(n),np.eye(n)))
        if label=='invisible':sigma=delta*V@Mbar@V.T
        else:
            sigma=np.zeros((m,m),complex);sigma[:n,:n]=delta*Mbar
        assert err(sigma+sigma.T)<1e-13
        R=k@sigma@k.T
        full=lifted['N'].copy();full[2*n+m:,2*n+m:]+=R
        Y=lifted['Y'];reference=lifted['reference']
        assert err(Y[:,2*n:2*n+m]-T@O)<1e-13
        # General massive Schur identity, no assumed zero zz block.
        Q=full[2*n:,2*n:];Qi=np.linalg.inv(Q)
        mix=full[:2*n,2*n:];Yq=Y[:,2*n:]
        effective=full[:2*n,:2*n]+mix@Qi@mix.T
        phieff=Y[:,:2*n]+Yq@Qi@mix.T
        contact=Yq@Qi@Yq.T
        rphase,rlog=body.logpf(Q);kd,kl=np.linalg.slogdet(k)
        auxfactor=rphase/kd*np.exp(rlog-kl)
        assert abs(auxfactor-1)<2e-11
        omega=O@sigma@O.T
        control=dict(visible_inverse_block_norm=float(np.linalg.norm(omega,2)),
            physical_contact_norm=float(np.linalg.norm(contact,2)),
            original_action_difference=err(effective-reference['N']),
            original_source_difference=err(phieff-reference['Phi']))
        if label=='invisible':
            assert max(control.values())<3e-11
        else:
            # For lambda=0, exact action/source/contact changes in shifted variables.
            B=np.vstack((np.zeros((n,m)),.5*g5@O))
            U=T@O
            assert err(effective-reference['N']-B@sigma@B.T)<3e-12
            assert err(phieff-reference['Phi']-U@sigma@B.T)<3e-12
            assert err(contact-U@sigma@U.T)<3e-12
            assert min(control.values())>1e-4
        # Choose actual physical source coordinates; they are charged test
        # sources, not declared Gauss-invariant observables by themselves.
        if label=='visible':
            aa=np.triu(abs(contact),1);i,j=np.unravel_index(np.argmax(aa),aa.shape)
        else:i,j=4,129
        z=np.eye(n)[:,[i,j]]
        coeff=[]
        for count in (0,2):
            zz=z[:,:count]
            target=base.pf(small_aug(effective,phieff,contact,zz))*auxfactor
            correct=corrected_large_coefficient(body.augmented(full,Y,zz),k)
            old=base.pf(body.augmented(reference['N'],reference['Phi'],zz))
            omitted=base.pf(small_aug(effective,phieff,np.zeros_like(contact),zz))*auxfactor
            scale=max(abs(target),abs(correct),1e-250)
            residual=float(abs(correct-target)/scale)
            assert residual<3e-9
            if label=='invisible':assert abs(correct-old)/max(abs(correct),abs(old),1e-250)<3e-9
            if label=='visible' and count==2:
                assert abs(target-omitted)/max(abs(target),abs(omitted),1e-250)>1e-4
            coeff.append(dict(sources=count,source_rows=[int(i),int(j)][:count],
                full_with_old_bulk_compensation=base.old.cpair(correct),
                complete_effective=base.old.cpair(target),original=base.old.cpair(old),
                dropping_source_contact=base.old.cpair(omitted),
                full_effective_relative_error=residual))
        rows.append(dict(kind=label,physical_mass=lam,Sigma_operator_norm=float(np.linalg.norm(sigma,2)),
            auxiliary_pairing_norm=float(np.linalg.norm(R,2)),
            old_bulk_factor_ratio=base.old.cpair(auxfactor),**control,source_checks=coeff))
    return dict(entry_round=682,not_formal_round=True,original_16_internal_channels=True,
        original_full_nonflat_links=True,Wilson_gap=gap,a=.23,auxiliary_layers=1,
        physical_source_T_singular_values=[float(sv[0]),float(sv[-1])],rows=rows,
        conditional_all_source_matching_not_full_average_identity_necessity=True,
        no_auxiliary_RP_claim=True,no_original_Gauss_RP_verdict=True)

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=682,all_checks_passed=True,
        differences=[{k:x[k] for k in ('kind','physical_mass','physical_contact_norm','original_action_difference')}
                     for x in result['rows']])))

