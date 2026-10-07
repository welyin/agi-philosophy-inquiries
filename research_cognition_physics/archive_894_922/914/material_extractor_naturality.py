"""914: naturality-reduced original773 extractor on actual912 source support.
Only the extractor and necessary scalar/metric projected values are computed.
No full weighted Dirac response, physical error certificate, or872 sum is claimed.
"""
from pathlib import Path
import sys,json,time,argparse,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
for n in ('911','913'):sys.path.insert(0,str(STAGE/n))
import discrete_record_tangent as previous
import future_test_projection as oldref
mb=previous.mb;ms=previous.ms;loop=previous.loop;REP=ms.ev.REP
TARGET=HERE/'material_extractor_naturality_results.json'
def maximum(v):return float(np.max(np.abs(v)))
def second_jets(field,points):
    p=field.jets(points)
    cache={(i,j):previous.derivative(field,points,(i,j)) for i in range(4) for j in range(i,4)}
    z=np.stack([np.stack([cache[tuple(sorted((i,j)))] for j in range(4)],axis=1) for i in range(4)],axis=1)
    p['ddg'],p['ddphi'],p['ddA']=mb.unpack(z)
    return p

def jacobian_from_jets(p):
    return np.stack([oldref.reference_delta(p,dict(g=p['dg'][:,i],phi=p['dphi'][:,i],dphi=p['ddphi'][:,i])) for i in range(4)],axis=-1)

def diff_extract(bg,field,points,p=None):
    p=second_jets(bg,points) if p is None else p
    v=field.jets(points);J=jacobian_from_jets(p);dX=oldref.reference_delta(p,v)
    xi=np.linalg.solve(J,dX[...,None])[...,0]
    return xi,J,p,v

def diff_jet(bg,field,points,p=None):
    xi,J,p,v=diff_extract(bg,field,points,p)
    h=1e-24;dxi=[]
    for mu in range(4):
        z=points.astype(complex);z[:,mu]+=1j*h
        dxi.append(diff_extract(bg,field,z)[0].imag/h)
    return xi,np.stack(dxi,axis=1),J,p,v

def features(p):
    D=p['dphi']+np.einsum('bma,aij,bj->bmi',p['A'],REP,p['phi'])
    F=p['dA']-p['dA'].swapaxes(1,2)+ms.bracket(p['A'][:,:,None,:],p['A'][:,None,:,:])
    return dict(H=p['phi'][:,:4],D=D[:,:,:4],F=F[...,:8])

def feature_delta(p,v):
    dD=v['dphi']+np.einsum('bma,aij,bj->bmi',v['A'],REP,p['phi'])+np.einsum('bma,aij,bj->bmi',p['A'],REP,v['phi'])
    dF=v['dA']-v['dA'].swapaxes(1,2)+ms.bracket(v['A'][:,:,None,:],p['A'][:,None,:,:])+ms.bracket(p['A'][:,:,None,:],v['A'][:,None,:,:])
    return dict(H=v['phi'][:,:4],D=dD[:,:,:4],F=dF[...,:8])

def lie_feature(p,xi,dxi):
    z=features(p)
    dD=p['ddphi']+np.einsum('brma,aij,bj->brmi',p['dA'],REP,p['phi'])+np.einsum('bma,aij,brj->brmi',p['A'],REP,p['dphi'])
    dF=p['ddA']-p['ddA'].swapaxes(2,3)+ms.bracket(p['dA'][:,:,:,None,:],p['A'][:,None,None,:,:])+ms.bracket(p['A'][:,None,:,None,:],p['dA'][:,:,None,:,:])
    h=np.einsum('br,bri->bi',xi,p['dphi'][:,:,:4])
    d=np.einsum('br,brmi->bmi',xi,dD[...,:4])+np.einsum('bmr,bri->bmi',dxi,z['D'])
    f=np.einsum('br,brmna->bmna',xi,dF[...,:8])+np.einsum('bmr,brna->bmna',dxi,z['F'])+np.einsum('bnr,bmra->bmna',dxi,z['F'])
    return dict(H=h,D=d,F=f)

def gauge_jets(p,xi,dxi,ddxi,alpha,dalpha,ddalpha):
    # Independent ordinary field-generator prolongation, includes all second
    # gauge-parameter jets; the natural-feature route does not use them.
    v={}
    v['g']=np.einsum('br,brmn->bmn',xi,p['dg'])+np.einsum('bmr,brn->bmn',dxi,p['g'])+np.einsum('bnr,bmr->bmn',dxi,p['g'])
    v['phi']=np.einsum('br,bri->bi',xi,p['dphi'])+np.einsum('ba,aij,bj->bi',alpha,REP,p['phi'])
    v['dphi']=np.einsum('bmr,bri->bmi',dxi,p['dphi'])+np.einsum('br,bmri->bmi',xi,p['ddphi'])+np.einsum('bma,aij,bj->bmi',dalpha,REP,p['phi'])+np.einsum('ba,aij,bmj->bmi',alpha,REP,p['dphi'])
    v['A']=np.einsum('br,brma->bma',xi,p['dA'])+np.einsum('bmr,bra->bma',dxi,p['A'])+ms.bracket(alpha[:,None,:],p['A'])-dalpha
    v['dA']=np.einsum('bnr,brma->bnma',dxi,p['dA'])+np.einsum('br,bnrma->bnma',xi,p['ddA'])+np.einsum('bnmr,bra->bnma',ddxi,p['A'])+np.einsum('bmr,bnra->bnma',dxi,p['dA'])+ms.bracket(dalpha[:,:,None,:],p['A'][:,None,:,:])+ms.bracket(alpha[:,None,None,:],p['dA'])-ddalpha
    return v

def internal_data(p):
    z=features(p);old=oldref.reference_data(p);q=np.sum(old['dh']*old['vh'],axis=-1)
    # Specific member of773's allowed positive auxiliary pairings: the clock
    # induces positive inverse metric G -2 vv/q; all component weights fixed1
    # in inherited model units. This choice must be retained in later vertices.
    Q=old['G']-2*old['vh'][:,:,None]*old['vh'][:,None,:]/q[:,None,None]
    rho=dict(H=np.einsum('aij,bj->bai',REP[:,:4,:4],z['H']),D=np.einsum('aij,bmj->bami',REP[:,:4,:4],z['D']))
    rho['F']=ms.bracket(np.eye(12)[None,:,None,None,:],np.pad(z['F'],((0,0),(0,0),(0,0),(0,4)))[:,None])[...,:8]
    # Real-jet bilinear pairing: analytic extension without conjugation.
    def pair(a,b):
        return (np.einsum('bai,bci->bac',a['H'],b['H'])+
                np.einsum('bmn,bami,bcni->bac',Q,a['D'],b['D'])+
                .5*np.einsum('bmr,bns,bamni,bcrsi->bac',Q,Q,a['F'],b['F']))
    gram=pair(rho,rho)
    def extract(dz):
        rhs=pair(rho,{k:v[:,None] for k,v in dz.items()})[...,0]
        return np.linalg.solve(gram,rhs[...,None])[...,0]
    return z,rho,gram,Q,extract

def actual(N=17):
    start=time.time();bg,field,init=previous.rebuilt.build_pair(N)
    src=loop.LoopSource(bg,2,8);points=src.points
    xi,dxi,J,p,v=diff_jet(bg,field,points)
    z,rho,gram,Q,extract=internal_data(p)
    dz=feature_delta(p,v);lie=lie_feature(p,xi,dxi)
    residual={k:dz[k]-lie[k] for k in dz};alpha=extract(residual)
    # For scalar/metric components Pi only needs xi, dxi, alpha (no dalpha).
    v0=dict(g=v['g']-np.einsum('br,brmn->bmn',xi,p['dg'])-np.einsum('bmr,brn->bmn',dxi,p['g'])-np.einsum('bnr,bmr->bmn',dxi,p['g']),
       phi=v['phi']-np.einsum('br,bri->bi',xi,p['dphi']),
       dphi=v['dphi']-np.einsum('bmr,bri->bmi',dxi,p['dphi'])-np.einsum('br,bmri->bmi',xi,p['ddphi']))
    pphi=v0['phi']-np.einsum('ba,aij,bj->bi',alpha,REP,p['phi'])
    # Natural-feature identity versus explicit second-jet generator, including
    # nonconstant internal transformations, on all original source samples.
    rng=np.random.default_rng(914)
    xig=rng.normal(size=xi.shape)*.02;dxig=rng.normal(size=dxi.shape)*.03
    ddxig=rng.normal(size=(len(points),4,4,4))*.04;ddxig=(ddxig+ddxig.swapaxes(1,2))/2
    ag=rng.normal(size=alpha.shape)*.05;dag=rng.normal(size=(len(points),4,12))*.04
    ddag=rng.normal(size=(len(points),4,4,12))*.07;ddag=(ddag+ddag.swapaxes(1,2))/2
    gauge=gauge_jets(p,xig,dxig,ddxig,ag,dag,ddag)
    dgz=feature_delta(p,gauge);lgz=lie_feature(p,xig,dxig)
    expected={k:np.einsum('ba,ba...->b...',ag,rho[k]) for k in rho}
    natur=max(maximum(dgz[k]-lgz[k]-expected[k]) for k in dz)
    extracted=extract({k:dgz[k]-lgz[k] for k in dz})
    diff_g=np.linalg.solve(J,oldref.reference_delta(p,gauge)[...,None])[...,0]
    # A deliberately wrong tensor rule drops derivative-of-xi terms.
    wronglie=lie_feature(p,xig,np.zeros_like(dxig))
    wrongalpha=extract({k:dgz[k]-wronglie[k] for k in dz})
    # Reevaluate old Xi derivatives by a separate real central difference.
    ids=np.arange(0,len(points),max(1,len(points)//16))[:16];central=[]
    for step in (2e-6,1e-6):
        errors=[]
        for mu in range(4):
            plus=points[ids].copy();minus=plus.copy();plus[:,mu]+=step;minus[:,mu]-=step
            fd=(diff_extract(bg,field,plus)[0]-diff_extract(bg,field,minus)[0])/(2*step)
            errors.append(maximum(fd-dxi[ids,mu]))
        central.append(dict(step=step,component_max_errors=errors))
    result=dict(N=N,points=len(points),original773_reference_used=True,original912_response_used=True,
       auxiliary_pairing='full H,D_muH,F_color with clock-induced positive inverse metric and fixed unit component weights',
       maximum_diff_parameter=maximum(xi),maximum_diff_parameter_derivative=maximum(dxi),
       maximum_internal_parameter=maximum(alpha),maximum_internal_color=maximum(alpha[:,:8]),maximum_internal_electroweak=maximum(alpha[:,8:]),
       max_projected_metric=maximum(v0['g']),max_projected_scalar=maximum(pphi),max_projected_neutral_p=maximum(pphi[:,5]),
       original_neutral_p_max=maximum(v['phi'][:,5]),
       minimum_sample_abs_old_Jdet=float(np.min(abs(np.linalg.det(J)))),
       maximum_sample_old_J_condition=float(np.max(np.linalg.cond(J))),
       minimum_sample_Gram_eigenvalue=float(np.min(np.linalg.eigvalsh(gram))),
       minimum_sample_positive_metric_eigenvalue=float(np.min(np.linalg.eigvalsh(Q))),
       reference_projected_residual=maximum(oldref.reference_delta(p,v0)),
       gauge_diff_extraction_error=maximum(diff_g-xig),
       gauge_internal_extraction_error=maximum(extracted-ag),
       feature_naturality_error=natur,
       omitted_tensor_index_terms_error=maximum(wrongalpha-ag),
       independent_real_derivative_diagnostics=central,
       full_projected_connection_computed=False,full_weighted_Dirac_response_computed=False,
       uniform_source_support_certified=False,physical_error_certified=False,
       elapsed_seconds=round(time.time()-start,3))
    assert natur<2e-12 and result['gauge_diff_extraction_error']<1e-8
    assert result['gauge_internal_extraction_error']<1e-10
    assert result['reference_projected_residual']<1e-13
    assert result['minimum_sample_Gram_eigenvalue']>0
    assert result['omitted_tensor_index_terms_error']>1e-5
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return result

def run():
    rows=[actual()]
    return dict(round=914,date='2026-10-06',status='working',rows=rows,
       derivative_order_bound=dict(Ldiff=1,Lint=2,Pi=3),
       naturality_reuses773_no_new_reference_field=True,
       full872_response_computed=False,actual_finite_observable_error_certified=False,full_goal_completed=False,
       source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),STAGE/'911/future_test_projection.py',STAGE/'913/discrete_record_tangent.py',STAGE/'913/independent_record_rebuild.py',STAGE/'908/magnetic_reference_source.py')})
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
