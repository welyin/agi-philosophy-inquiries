"""918: actual Euler-density Hessian residual, retaining background contacts.
Raw tensor components and variational density components are not interchangeable.
"""
from pathlib import Path
import sys,argparse,json,hashlib,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'917'))
import joint_constraint_driving as joint
stress=joint.stress;c=joint.c;validation=joint.validation
TARGET=HERE/'euler_density_contacts_results.json'

def density(E,z):
    mu=z['mu'];G=z['ig']
    return dict(gravity=-.5*mu[:,None,None]*np.einsum('bma,ban,bnr->bmr',G,E['gravity'],G),
        scalar=mu[:,None]*E['scalar'],YM=mu[:,None,None]*E['YM'])

def density_ward(Q,dQ,p,base):
    Gam=base['z']['Gamma'];g=p['g']
    div=np.einsum('bmml->bl',dQ['gravity'])+np.einsum('blmr,bmr->bl',Gam,Q['gravity'])
    diff=-2*np.einsum('bnl,bl->bn',g,div)+np.einsum('bA,bnA->bn',Q['scalar'],base['D'])+np.einsum('bma,bnma->bn',Q['YM'],base['F'])
    internal=np.einsum('bmma->ba',dQ['YM'])+sum(c.bracket(p['A'][:,mu,:],Q['YM'][:,mu,:]) for mu in range(4))
    internal+=np.einsum('bA,aAB,bB->ba',Q['scalar'],c.REP,p['phi'])
    return diff,internal

def actual_density(bg,field,u,x,par):
    p=c.aj.Cache(bg,x).jets();v=c.aj.Cache(field,x).jets();base=c.euler(p,par);h=1e-24
    varied=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    derivative={k:a.imag/h for k,a in density(varied,varied['z']).items()}
    T=stress.direct_stress(bg,u,x)
    source=density(dict(gravity=-T,scalar=np.zeros_like(base['scalar']),YM=np.zeros_like(base['YM'])),base['z'])
    return {k:derivative[k]+source[k] for k in derivative}

def run():
    start=time.time();par=c.parameters();bg,field,_=c.ex.previous.rebuilt.build_pair()
    src=c.ex.loop.LoopSource(bg,2,8);x=src.points[::64];bc=c.aj.Cache(bg,x);vc=c.aj.Cache(field,x)
    p=bc.jets();v=vc.jets();dp=[bc.jets((i,)) for i in range(4)]
    base,dE=validation.coord_derivatives(p,dp,par);z=base['z'];G=z['ig'];E=base['gravity']
    h=1e-24;varied=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    dGamma=varied['z']['Gamma'].imag/h;deltaG=varied['z']['ig'].imag/h
    delta={k:varied[k].imag/h for k in ('gravity','scalar','YM','D','F')}
    sigma=.5*np.einsum('bmn,bmn->b',G,v['g'])
    dsigma=.5*(np.einsum('brmn,bmn->br',z['di'],v['g'])+np.einsum('bmn,brmn->br',G,v['dg']))
    trace_dGamma=np.einsum('bmmr->br',dGamma)
    assert c.maximum(dsigma-trace_dGamma)<1e-14
    V=v['g'];dV=v['dg'];di=z['di']
    H=np.einsum('bma,ban->bmn',V,G)
    dH=np.einsum('brma,ban->brmn',dV,G)+np.einsum('bma,bran->brmn',V,di)
    corr={k:sigma.reshape((-1,)+(1,)*(base[k].ndim-1))*base[k] for k in ('gravity','scalar','YM')}
    corr['gravity']-=H@E+E@H.swapaxes(1,2)
    dcorr={}
    for k in corr:
        shape=(-1,4)+(1,)*(base[k].ndim-1)
        dcorr[k]=dsigma.reshape(shape)*base[k][:,None]+sigma.reshape((-1,)+(1,)*(dE[k].ndim-1))*dE[k]
    dcorr['gravity']-=dH@E[:,None]+H[:,None]@dE['gravity']+dE['gravity']@H[:,None].swapaxes(-1,-2)+E[:,None]@dH.swapaxes(-1,-2)
    u,_,_=stress.weighted.receiver_modes();T=stress.stress_jet(bg,u,x)
    rA={k:delta[k].copy() for k in corr};rA['gravity']-=T['stress']
    tilde={k:rA[k]+corr[k] for k in corr};Q=density(tilde,z)
    direct=actual_density(bg,field,u,x,par)
    conversion={k:c.maximum(Q[k]-direct[k]) for k in Q};assert max(conversion.values())<1e-12,conversion
    covE=joint.covariant_gravity_derivative(E,dE['gravity'],z)
    bgC=np.einsum('bma,bman->bn',deltaG,covE)
    bgC-=np.einsum('bma,brma,brn->bn',G,dGamma,E)+np.einsum('bma,brmn,bar->bn',G,dGamma,E)
    bgC+=np.einsum('bA,bnA->bn',base['scalar'],delta['D'])+np.einsum('bma,bnma->bn',base['YM'],delta['F'])
    bgIvol=np.einsum('br,bra->ba',trace_dGamma,base['YM'])
    bgIother=sum(c.bracket(v['A'][:,m,:],base['YM'][:,m,:]) for m in range(4))+np.einsum('bA,aAB,bB->ba',base['scalar'],c.REP,v['phi'])
    Wcorr,WcorrI=joint.ward(corr,dcorr,p,base)
    predicted=z['mu'][:,None]*(-T['divergence']-bgC+Wcorr)
    predictedI=z['mu'][:,None]*(-bgIvol-bgIother+WcorrI)
    pure_volume=z['mu'][:,None]*np.einsum('br,bra->ba',dsigma,base['YM'])
    cancellation=c.maximum(pure_volume-z['mu'][:,None]*bgIvol)
    reducedI=-z['mu'][:,None]*bgIother
    internal_reduction_error=c.maximum(predictedI-reducedI)
    assert internal_reduction_error<1e-14
    rows=[]
    for step in (2e-6,1e-6):
        dQ={k:[] for k in Q}
        for mu in range(4):
            xp=x.copy();xm=x.copy();xp[:,mu]+=step;xm[:,mu]-=step
            plus=actual_density(bg,field,u,xp,par);minus=actual_density(bg,field,u,xm,par)
            for k in dQ:dQ[k].append((plus[k]-minus[k])/(2*step))
        dQ={k:np.stack(a,axis=1) for k,a in dQ.items()};wd,wi=density_ward(Q,dQ,p,base)
        rows.append(dict(step=step,diffeomorphism_density_identity_error=c.maximum(wd-predicted),internal_density_identity_error=c.maximum(wi-predictedI),
                         directly_differenced_internal_Ward_max=c.maximum(wi)))
    assert max(r['diffeomorphism_density_identity_error'] for r in rows)<2e-9,rows
    assert max(r['internal_density_identity_error'] for r in rows)<2e-9,rows
    # This canonical Gauss quantity is unchanged by representing the SAME derivative.
    direct_gauss=z['mu'][:,None]*delta['YM'][:,0,:]+(varied['z']['mu'].imag/h)[:,None]*base['YM'][:,0,:]
    gauss_identity=c.maximum(Q['YM'][:,0,:]-direct_gauss)
    assert gauss_identity<1e-13
    return dict(round=918,status='working',date='2026-10-06',samples=len(x),N=17,
       same_action_background_response_receiver=True,direct_variation_density_conversion_error=conversion,
       volume_log_derivative_identity_error=c.maximum(dsigma-trace_dGamma),internal_volume_cancellation_error=cancellation,
       internal_density_reduction_error=internal_reduction_error,
       raw_tensor_internal_driver_times_volume_max=c.maximum(z['mu'][:,None]*(bgIvol+bgIother)),
       corrected_density_internal_driver_max=c.maximum(predictedI),corrected_density_diffeomorphism_driver_max=c.maximum(predicted),
       density_background_contact_correction_max={k:c.maximum(a) for k,a in corr.items()},
       same_Gauss_variation_identity_error=gauss_identity,same_Gauss_variation_sample_max=c.maximum(direct_gauss),
       independent_coordinate_checks=rows,small_internal_driver_resolved_by_coordinate_difference=False,
       physical_error_not_removed_by_representation=True,actual_finite_observable_error_certified=False,
       full872_response_computed=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3),
       source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),STAGE/'917/joint_constraint_driving.py',STAGE/'917/receiver_stress_derivative.py',STAGE/'916/covariant_joint_residual.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))
