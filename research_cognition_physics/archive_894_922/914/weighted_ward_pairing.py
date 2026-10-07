"""914: enforce the proven closed-source Ward identity before quadrature.
No post-hoc subtraction of an inferred physical signal: the identity and the
independently evaluated pure-gauge quadrature defect are both retained.
"""
from pathlib import Path
import json,hashlib,argparse,time
import numpy as np
import weighted_receiver_pairing as weighted
ex=weighted.ex;previous=ex.previous;loop=ex.loop;ms=ex.ms
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;TARGET=HERE/'weighted_ward_pairing_results.json'

def ward_pair(src,base,dalpha,w,dw):
    xi,dxi,p,v,alpha,rho,residual=base
    eta=np.einsum('br,bra->ba',xi,p['A'])-alpha
    deta=np.einsum('bmr,bra->bma',dxi,p['A'])+np.einsum('br,bmra->bma',xi,p['dA'])-dalpha
    Deta=deta+ms.bracket(p['A'],eta[:,None,:])
    lowxi=np.einsum('bmn,bn->bm',p['g'],xi)
    qg=dw[:,:,None]*lowxi[:,None,:]+dw[:,None,:]*lowxi[:,:,None]
    qa=dw[:,:,None]*eta[:,None,:]
    qf=dw[:,None,:,None]*Deta[:,:,None,:]-dw[:,:,None,None]*Deta[:,None,:,:]
    # q has zero scalar component: delta h,delta dh,delta s,delta p are zero.
    qm=np.sum(src.data['metric']*qg,axis=(1,2))+np.sum(src.data['B']*qf[...,:8],axis=(1,2,3))
    correction=np.array([np.sum(C*qa)+np.sum(k[:,3]*qm) for C,k in zip(src.currents,src.ks)])
    raw=ex.oldref.multiply(v,w,dw)
    rawpair=np.array([np.sum(ms.pair(s,raw)) for s in src.sources])
    # Independent evaluation of j(R(w xi,w alpha)), using scalar covariance
    # of X_new. Only ordinary gauge generator values enter the loop current.
    xiw=w[:,None]*xi;ax=w[:,None]*alpha
    dxi_w=w[:,None,None]*dxi+dw[:,:,None]*xi[:,None,:]
    da_w=w[:,None,None]*dalpha+dw[:,:,None]*alpha[:,None,:]
    Rconn=np.einsum('br,brma->bma',xiw,p['dA'])+np.einsum('bmr,bra->bma',dxi_w,p['A'])+ms.bracket(ax[:,None,:],p['A'])-da_w
    Jnew=np.linalg.inv(src.Jinv);RX=np.einsum('bam,bm->ba',Jnew,xiw)
    gauge=np.array([np.sum(C*Rconn)+np.sum(k*RX) for C,k in zip(src.currents,src.ks)])
    return rawpair+correction,gauge,correction

def run():
    start=time.time();bg,field,_=previous.rebuilt.build_pair();u,z,mode=weighted.receiver_modes()
    src=loop.LoopSource(bg,2,8);src.kernels();base=weighted.extract_values(bg,field,src.points)
    w,dw,stats=weighted.weight_jets(u,z,src.points);step=1e-6;da=[]
    for mu in range(4):
        xp=src.points.copy();xm=xp.copy();xp[:,mu]+=step;xm[:,mu]-=step
        da.append((weighted.extract_values(bg,field,xp)[4]-weighted.extract_values(bg,field,xm)[4])/(2*step))
    da=np.stack(da,axis=1);pv=weighted.projected_inputs(base,da)
    direct=weighted.weighted_pair(src,pv,w,dw);ward,gauge,corr=ward_pair(src,base,da,w,dw)
    ones=np.ones_like(w);zeros=np.zeros_like(dw)
    constant,gauge_constant,_=ward_pair(src,base,da,ones,zeros)
    direct_constant=weighted.weighted_pair(src,pv,ones,zeros)
    original=np.array(src.response(lambda x,p:field.jets(x))['total'])
    error=ex.maximum(ward-direct-gauge);constant_error=ex.maximum(constant-original)
    constant_decomp=ex.maximum(constant-direct_constant-gauge_constant)
    assert error<2e-12 and constant_error<1e-16 and constant_decomp<2e-12
    old=json.loads((HERE/'weighted_receiver_pairing_results.json').read_text('utf-8'))['rows'][-1]
    oldvalue=np.array(old['weighted_source']['real'])+1j*np.array(old['weighted_source']['imag'])
    assert ex.maximum(direct-oldvalue)<1e-14
    result=dict(round=914,date='2026-10-06',N=17,segments=8,samples=len(src.points),
       auxiliary_pairing='773 full feature tuple with clock-positive metric, component weights1 in inherited units',
       auxiliary_mode_initial_data='z=(p/.02)u at t=0, evolved by same906 free Dirac equation',
       weight=stats,mode=mode,alpha_derivative_step=step,
       direct_projected_weighted_source=weighted.complex_list(direct),
       Ward_reduced_weighted_source=weighted.complex_list(ward),
       independent_pure_gauge_quadrature_defect=weighted.complex_list(gauge),
       Ward_principal_correction=weighted.complex_list(corr),
       conditional_non_diagonal_vA_matrix_element=weighted.complex_list(1j*np.array([6.,6.,.75])*ward),
       weighted_Ward_decomposition_error=error,
       constant_weight_original_source_error=constant_error,
       constant_weight_Ward_decomposition_error=constant_decomp,
       constant_weight_quadrature_defect=weighted.complex_list(gauge_constant),
       previous_direct_calculation_reproduced_error=ex.maximum(direct-oldvalue),
       direct_gauge_defect_not_misidentified_as_physical_signal=True,
       uniform_reference_support_certified=False,finite_source_error_certified=False,
       full872_response_computed=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3),
       source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),HERE/'weighted_receiver_pairing.py',HERE/'material_extractor_naturality.py',STAGE/'908/relational_loop_source.py')})
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
