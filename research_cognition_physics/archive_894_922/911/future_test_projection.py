"""911: actual old-reference cutoff commutator in the neutral p component.
This component of full Pi is independent of the internal gauge extractor because
p is an internal singlet. It does not compute the complete Pi array or872 source.
"""
from pathlib import Path
import sys,itertools,json,time,argparse,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
for n in ('908','910'):sys.path.insert(0,str(STAGE/n))
import physical_direction_dual_bridge as bridge
mb=bridge.mb
TARGET=HERE/'future_test_projection_results.json'

def reference_data(p):
    g=p['g'];phi=p['phi'];dphi=p['dphi'];G=np.linalg.inv(g)
    h=np.sqrt(np.sum(phi[...,:4]**2,axis=-1));r=phi[...,:4]/h[...,None]
    dh=np.einsum('...A,...mA->...m',r,dphi[...,:,:4]);ds=dphi[..., :,4]
    vh=np.einsum('...mn,...n->...m',G,dh);vs=np.einsum('...mn,...n->...m',G,ds)
    X=np.stack((h,phi[...,4],np.sum(dh*vh,axis=-1),np.sum(ds*vs,axis=-1)),axis=-1)
    return dict(G=G,h=h,r=r,dh=dh,ds=ds,vh=vh,vs=vs,X=X)

def reference_delta(p,v):
    d=reference_data(p);dh=np.sum(d['r']*v['phi'][...,:4],axis=-1)
    dr=(v['phi'][...,:4]-d['r']*dh[...,None])/d['h'][...,None]
    dn=np.einsum('...A,...mA->...m',dr,p['dphi'][...,:,:4])+np.einsum('...A,...mA->...m',d['r'],v['dphi'][...,:,:4])
    da=2*np.sum(d['vh']*dn,axis=-1)-np.einsum('...m,...mn,...n->...',d['vh'],v['g'],d['vh'])
    db=2*np.sum(d['vs']*v['dphi'][..., :,4],axis=-1)-np.einsum('...m,...mn,...n->...',d['vs'],v['g'],d['vs'])
    return np.stack((dh,v['phi'][...,4],da,db),axis=-1)

def old_J(bg,points):
    columns=[];h=1e-24
    for mu in range(4):
        z=points.astype(complex);z[:,mu]+=1j*h
        columns.append(reference_data(bg.jets(z))['X'].imag/h)
    return np.stack(columns,axis=-1)

def smooth_cutoff(t,start=5e-5,end=1.5e-4):
    z=(t-start)/(end-start);chi=np.zeros_like(z);prime=np.zeros_like(z)
    chi[z>=1]=1;inside=(z>0)&(z<1);a=z[inside]
    chi[inside]=1/(1+np.exp(1/a-1/(1-a)))
    prime[inside]=chi[inside]*(1-chi[inside])*(1/a**2+1/(1-a)**2)/(end-start)
    return chi,prime

def multiply(v,chi,dchi):
    out={k:a*chi.reshape((-1,)+(1,)*(a.ndim-1)) for k,a in v.items()}
    out['dphi']+=dchi[:, :,None]*v['phi'][:,None,:]
    out['dA']+=dchi[:,:,None,None]*v['A'][:,None,:,:]
    if 'dg' in out:out['dg']+=dchi[:,:,None,None]*v['g'][:,None,:,:]
    return out

def actual_check(N):
    bg=mb.Background(N);u=bridge.PhysicalDirection(N)
    times=np.array([2.5e-5,7.5e-5,1e-4,1.25e-4,1.75e-4])
    shifts=np.array(list(itertools.product((-.005,0.,.005),repeat=3)))
    points=np.concatenate([np.column_stack((np.full(len(shifts),t),mb.CENTER[None,1:]+shifts)) for t in times])
    p=bg.jets(points);v=u.jets(points);J=old_J(bg,points);Xv=reference_delta(p,v)
    chi,prime=smooth_cutoff(points[:,0]);dc=np.zeros((len(points),4));dc[:,0]=prime
    scaled=multiply(v,chi,dc);Xscaled=reference_delta(p,scaled)
    principal=np.zeros_like(Xv);d=reference_data(p)
    principal[:,2]=2*np.sum(d['vh']*dc,axis=-1)*Xv[:,0]
    principal[:,3]=2*np.sum(d['vs']*dc,axis=-1)*Xv[:,1]
    comm_actual=np.linalg.solve(J,(Xscaled-chi[:,None]*Xv)[...,None])[...,0]
    comm_formula=np.linalg.solve(J,principal[...,None])[...,0]
    xi=np.linalg.solve(J,Xv[...,None])[...,0];xic=np.linalg.solve(J,Xscaled[...,None])[...,0]
    projected_p=v['phi'][:,5]-np.einsum('bm,bm->b',p['dphi'][:,:,5],xi)
    projected_cut_p=scaled['phi'][:,5]-np.einsum('bm,bm->b',p['dphi'][:,:,5],xic)
    comm_p=projected_cut_p-chi*projected_p
    formula_p=-np.einsum('bm,bm->b',p['dphi'][:,:,5],comm_formula)
    rows=[]
    for t in times:
        sel=points[:,0]==t
        rows.append(dict(time=float(t),chi=float(chi[sel][0]),chi_prime=float(prime[sel][0]),
            max_abs_diff_extractor_commutator=float(np.max(abs(comm_actual[sel]))),
            max_abs_full_projector_p_commutator=float(np.max(abs(comm_p[sel]))),
            center_full_projector_p_commutator=float(comm_p[sel][13]),
            max_abs_projected_p=float(np.max(abs(projected_p[sel]))),
            minimum_sample_abs_old_Jdet=float(np.min(abs(np.linalg.det(J[sel]))))))
    # Independent nonlinear field-jet directional difference for X'[u].
    eps=1e-24;z={k:np.asarray(p[k],complex)+1j*eps*v[k] for k in p}
    dx=reference_data(z)['X'].imag/eps
    return dict(N=N,rows=rows,independent_reference_variation_error=float(np.max(abs(dx-Xv))),
        full_principal_reference_identity_error=float(np.max(abs(Xscaled-chi[:,None]*Xv-principal))),
        extraction_commutator_error=float(np.max(abs(comm_actual-comm_formula))),
        neutral_p_projector_commutator_error=float(np.max(abs(comm_p-formula_p))),
        outside_transition_commutator=float(np.max(abs(comm_p[prime==0]))),
        all_other_projector_components_computed=False,uniform_reference_support_certified=False)

def run():
    rows=[]
    for N in (16,24):
        start=time.time();row=actual_check(N);rows.append(row)
        print(json.dumps(row,ensure_ascii=False),flush=True);print('seconds',round(time.time()-start,2),flush=True)
    return dict(round=911,status='working',date='2026-10-06',rows=rows,
        original910_constrained_direction_used=True,original773_reference_used=True,
        full872_response_computed=False,actual_net_tail_error_evaluated=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      (Path(__file__),STAGE/'910/physical_direction_dual_bridge.py',STAGE/'908/material_background.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    v=run()
    if a.write:TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
