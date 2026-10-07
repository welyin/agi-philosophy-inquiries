"""916: continuum Euler residuals of the SAME finite spacetime interpolants.
Independent of canonical evolution RHS; all bosonic sectors, reduced Einstein
operator and off-shell Noether identities retained. Diagnostics, not error bounds.
"""
from pathlib import Path
import sys,json,hashlib,argparse,time
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'915'))
import analytic_material_jets as aj
ex=aj.ex;ev=ex.ms.ev;REP=ex.REP;bracket=ex.ms.bracket
TARGET=HERE/'covariant_joint_residual_results.json'

def maximum(x):return float(np.max(abs(x)))
def parameters():
    with ev.ResearchRuntime(ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        m=ev.Model(17,old)
    return SimpleNamespace(K=m.K,L=m.L,u=m.u)

def geometry(p):
    g=p['g'];d=p['dg'];dd=p['ddg'];ig=np.linalg.inv(g)
    di=-np.einsum('bij,brjk,bkl->bril',ig,d,ig)
    Z=d+d.swapaxes(1,2)-np.moveaxis(d,1,-1)
    dZ=dd+dd.swapaxes(2,3)-np.moveaxis(dd,2,-1)
    Gamma=.5*np.einsum('brs,bmns->brmn',ig,Z)
    dGamma=.5*np.einsum('bars,bmns->barmn',di,Z)+.5*np.einsum('brs,bamns->barmn',ig,dZ)
    Ric=np.einsum('brrmn->bmn',dGamma)-np.einsum('bnrmr->bmn',dGamma)
    Ric+=np.einsum('brrl,blmn->bmn',Gamma,Gamma)-np.einsum('brnl,blmr->bmn',Gamma,Gamma)
    R=np.einsum('bmn,bmn->b',ig,Ric)
    C=np.einsum('bij,bijm->bm',ig,d)-.5*np.einsum('bij,bmij->bm',ig,d)
    dC=np.einsum('brij,bijm->brm',di,d)+np.einsum('bij,brijm->brm',ig,dd)
    dC-=.5*np.einsum('brij,bmij->brm',di,d)+.5*np.einsum('bij,brmij->brm',ig,dd)
    nablaC=dC-np.einsum('blmn,bl->bmn',Gamma,C)
    constraint_tensor=.5*(nablaC+nablaC.swapaxes(1,2))-.5*g*np.einsum('bmn,bmn->b',ig,nablaC)[:,None,None]
    gamma=g[:,1:,1:];invgamma=np.linalg.inv(gamma);vol=np.sqrt(np.linalg.det(gamma))
    alpha=1/np.sqrt(-ig[:,0,0]);beta=-ig[:,0,1:]/ig[:,0,0,None]
    normal=np.concatenate((np.ones((len(g),1)),-beta),axis=1)/alpha[:,None]
    return dict(ig=ig,di=di,Gamma=Gamma,dGamma=dGamma,Ric=Ric,R=R,C=C,nablaC=nablaC,
        constraint_tensor=constraint_tensor,gamma=gamma,invgamma=invgamma,vol=vol,
        alpha=alpha,beta=beta,normal=normal,mu=alpha*vol,d=d)

def target(phi,par):
    B=len(phi);dtype=phi.dtype;ph=phi[:,:5];F=2-np.sum(ph*ph,axis=-1)/6
    G=np.zeros((B,6,6),dtype=dtype);G[:,:5,:5]=np.eye(5)/F[:,None,None]+ph[:,:,None]*ph[:,None,:]/(6*F[:,None,None]**2);G[:,5,5]=1
    dG=np.zeros((B,6,6,6),dtype=dtype)
    for a in range(5):
        e=np.eye(5)[a]
        dG[:,a,:5,:5]=ph[:,a,None,None]*np.eye(5)/(3*F[:,None,None]**2)+(e[None,:,None]*ph[:,None,:]+ph[:,:,None]*e[None,None,:])/(6*F[:,None,None]**2)+ph[:,a,None,None]*ph[:,:,None]*ph[:,None,:]/(9*F[:,None,None]**3)
    delta=np.stack((np.sum(ph[:,:4]**2,axis=-1)-par.u[0],ph[:,4]**2-par.u[1]),axis=1)
    Ld=np.einsum('ab,nb->na',par.L,delta);V=.25*np.sum(delta*Ld,axis=1)
    U=V/F**2+.5*phi[:,5]**2
    Vp=np.concatenate((Ld[:,0,None]*ph[:,:4],Ld[:,1,None]*ph[:,4,None]),axis=1)
    Up=np.concatenate((Vp/F[:,None]**2+2*V[:,None]*ph/(3*F[:,None]**3),phi[:,5,None]),axis=1)
    return G,dG,U,Up

def euler(p,par):
    z=geometry(p);g=p['g'];ig=z['ig'];di=z['di'];Gamma=z['Gamma']
    phi=p['phi'];A=p['A'];rp=np.einsum('aij,bj->bai',REP,phi)
    D=p['dphi']+np.einsum('bma,bai->bmi',A,rp)
    dD=p['ddphi']+np.einsum('brma,bai->brmi',p['dA'],rp)+np.einsum('bma,aij,brj->brmi',A,REP,p['dphi'])
    F=p['dA']-p['dA'].swapaxes(1,2)+bracket(A[:,:,None,:],A[:,None,:,:])
    dF=p['ddA']-p['ddA'].swapaxes(2,3)+bracket(p['dA'][:,:,:,None,:],A[:,None,None,:,:])+bracket(A[:,None,:,None,:],p['dA'][:,:,None,:,:])
    Fup=np.einsum('bmr,bns,brsa->bmna',ig,ig,F)
    dFup=np.einsum('bkmr,bns,brsa->bkmna',di,ig,F)+np.einsum('bmr,bkns,brsa->bkmna',ig,di,F)+np.einsum('bmr,bns,bkrsa->bkmna',ig,ig,dF)
    G,Gphi,U,Up=target(phi,par)
    dG=np.einsum('bAIJ,brA->brIJ',Gphi,p['dphi'])
    flux=np.einsum('bmn,bAB,bnB->bmA',ig,G,D)
    dflux=np.einsum('brmn,bAB,bnB->brmA',di,G,D)+np.einsum('bmn,brAB,bnB->brmA',ig,dG,D)+np.einsum('bmn,bAB,brnB->brmA',ig,G,dD)
    density_derivative=np.einsum('bmmr->br',Gamma)
    scalar=np.einsum('bmmA->bA',dflux)+np.einsum('bm,bmA->bA',density_derivative,flux)+np.einsum('bma,aAB,bmB->bA',A,REP,flux)
    scalar-=.5*np.einsum('bABC,bmn,bmB,bnC->bA',Gphi,ig,D,D)+Up
    ym=(np.einsum('bmmna->bna',dFup)+np.einsum('bm,bmna->bna',density_derivative,Fup))*par.K
    ym+=sum(bracket(A[:,m,None,:],Fup[:,m,:,:]*par.K) for m in range(4))
    ym-=np.einsum('bmA,baA->bma',flux,rp)
    norm=np.einsum('bmn,bmA,bAB,bnB->b',ig,D,G,D)
    T=np.einsum('bmA,bAB,bnB->bmn',D,G,D)-g*(.5*norm+U)[:,None,None]
    T+=np.einsum('bmra,brs,bnsa,a->bmn',F,ig,F,par.K)-.25*g*np.einsum('bmna,bmna,a->b',F,Fup,par.K)[:,None,None]
    grav=z['Ric']-.5*g*z['R'][:,None,None]-T
    return dict(gravity=grav,scalar=scalar,YM=ym,reduced_gravity=grav-z['constraint_tensor'],stress=T,D=D,F=F,z=z)

def noether(p,dp,par):
    r=euler(p,par);h=1e-24
    derivatives={k:[] for k in ('gravity','scalar','YM','stress')}
    for mu in range(4):
        rr=euler({k:p[k].astype(complex)+1j*h*dp[mu][k] for k in p},par)
        for k in derivatives:derivatives[k].append(rr[k].imag/h)
    der={k:np.stack(v,axis=1) for k,v in derivatives.items()}
    z=r['z'];G=z['ig'];Gam=z['Gamma'];E=r['gravity']
    cov=der['gravity']-np.einsum('brma,brn->bman',Gam,E)-np.einsum('brmn,bar->bman',Gam,E)
    div=np.einsum('bma,bman->bn',G,cov)
    exchange=np.einsum('bA,bnA->bn',r['scalar'],r['D'])+np.einsum('bma,bnma->bn',r['YM'],r['F'])
    divY=np.einsum('bmma->ba',der['YM'])+np.einsum('bmmr,bra->ba',Gam,r['YM'])+sum(bracket(p['A'][:,m,:],r['YM'][:,m,:]) for m in range(4))
    charge=np.einsum('bA,aAB,bB->ba',r['scalar'],REP,p['phi'])
    return dict(diffeomorphism_noether_error=maximum(div+exchange),internal_noether_error=maximum(divY+charge),
        uncancelled_gravity_divergence=maximum(div),uncancelled_matter_exchange=maximum(exchange),
        uncancelled_YM_divergence=maximum(divY),uncancelled_scalar_charge=maximum(charge))

def run():
    start=time.time();par=parameters();bg,field,_=ex.previous.rebuilt.build_pair()
    src=ex.loop.LoopSource(bg,2,8);x=src.points
    bc=aj.Cache(bg,x);vc=aj.Cache(field,x);p=bc.jets();v=vc.jets();base=euler(p,par)
    h=1e-24;shifted=euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    tangent={k:shifted[k].imag/h for k in ('gravity','scalar','YM','reduced_gravity','stress')}
    weighted=aj.prior.weighted;u,_,_=weighted.receiver_modes()
    uu=u.evaluate(x);du=np.stack([u.evaluate(x,axis=i) for i in range(4)],axis=1)
    receiver=ex.previous.rebuilt.r
    frame=receiver.shared.frame_data(dict(g=p['g']),base['z'])
    data=receiver.bilinear_jet(dict(g=p['g']),base['z'],frame,uu,du)
    residual=dict(gravity=tangent['gravity']-data['stress'],scalar=tangent['scalar'],YM=tangent['YM'],reduced_gravity=tangent['reduced_gravity']-data['stress'])
    small=x[::64];cache=aj.Cache(bg,small);ps=cache.jets();dps=[cache.jets((i,)) for i in range(4)]
    identities=noether(ps,dps,par)
    assert identities['diffeomorphism_noether_error']<1e-9,identities
    assert identities['internal_noether_error']<1e-9,identities
    n=base['z']['normal'];dn=shifted['z']['normal'].imag/h;E=base['gravity'];R=residual['gravity']
    # Constraint derivatives include the moving normal on a nonzero residual background.
    normal_only=2*np.einsum('bm,bmn,bn->b',n,R,n)
    H=normal_only+4*np.einsum('bm,bmn,bn->b',dn,E,n)
    M=-np.einsum('bm,bmi->bi',n,R[:,:,1:])-np.einsum('bm,bmi->bi',dn,E[:,:,1:])
    gauss=base['z']['mu'][:,None]*residual['YM'][:,0,:]+shifted['z']['mu'].imag[:,None]/h*base['YM'][:,0,:]
    result=dict(round=916,status='working',date='2026-10-06',N=17,samples=len(x),
        original_background_response_receiver_retained=True,
        base_Euler_residual_sample_max={k:maximum(base[k]) for k in ('gravity','scalar','YM','reduced_gravity')},
        linear_joint_source_residual_sample_max={k:maximum(a) for k,a in residual.items()},
        base_harmonic_sample_max=maximum(base['z']['C']),
        receiver_Dirac_residual_sample_max=maximum(data['dirac_residual']),receiver_stress_sample_max=maximum(data['stress']),
        covariant_constraint_variation_sample_max=dict(H=maximum(H),M=maximum(M),Gauss=maximum(gauss)),
        moving_normal_H_correction_sample_max=maximum(H-normal_only),
        offshell_noether_checks=identities,
        source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),STAGE/'915/analytic_material_jets.py',STAGE/'913/independent_record_rebuild.py',STAGE/'906/compact_receiver_propagation.py')},
        sample_residual_is_not_uniform_bound=True,physical_solution_error_certified=False,
        full872_response_computed=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3))
    return result
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--write',action='store_true');a=pa.parse_args()
    if a.write:assert not TARGET.exists()
    r=run()
    if a.write:TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
