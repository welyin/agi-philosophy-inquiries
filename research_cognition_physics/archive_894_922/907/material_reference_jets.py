"""907 working: original old/new material reference jets on the actual902 trajectory.
Not the full Pi/Pi_adjoint, smeared Wilson source, or a certified reference patch.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'904'))
import coupled_boson_tangent as tangent
ev=tangent.base;TARGET=HERE/'material_reference_jets_results.json'

def references(model,y):
    F,z,m=model.rhs(y,aux=True);phi=y['phi'];h=np.sqrt(np.sum(phi[...,:4]**2,axis=-1));s=phi[...,4]
    dphi=np.concatenate((F['phi'][...,None,:],model.grad(phi)),axis=3)
    dh=np.einsum('...A,...mA->...m',phi[...,:4],dphi[...,:,:4])/h[...,None];ds=dphi[..., :,4]
    aa=np.einsum('...m,...mn,...n->...',dh,z['ig'],dh);bb=np.einsum('...m,...mn,...n->...',ds,z['ig'],ds)
    vh=np.einsum('...mn,...n->...m',z['ig'],dh)
    P=z['ig']-vh[..., :,None]*vh[...,None,:]/aa[...,None,None]
    field=np.zeros(h.shape+(4,4,8),dtype=np.result_type(phi,complex));field[...,1:,1:,:]=m['curv'][...,:,:,:8]
    F0=z['alpha'][...,None,None]*m['e']+np.einsum('...j,...jia->...ia',z['beta'],m['curv'])
    field[...,0,1:,:]=F0[...,:8];field[...,1:,0,:]=-F0[...,:8]
    Mh=.5*np.einsum('...mr,...ns,...mna,...rsa->...',P,P,field,field)
    return dict(old=np.stack((h,s,aa,bb),axis=-1),new=np.stack((h,s,phi[...,5],Mh),axis=-1))

def derivative(model,y,v):
    eps=1e-24;out=references(model,{k:np.asarray(y[k],complex)+1j*eps*v[k] for k in y})
    return {k:a.imag/eps for k,a in out.items()}

def diffeo(model,y):
    N=model.N;x=2*np.pi*np.arange(N)/N;X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    xi=np.stack((.03*np.sin(X)+.02*np.cos(Y),.04*np.sin(Z)+.015*np.cos(Y),.02*np.cos(X)+.013*np.sin(Z)),axis=-1)
    dxi=model.grad(xi);div=np.trace(dxi,axis1=-2,axis2=-1)
    D=np.zeros(xi.shape[:-1]+(4,4));D[...,1:,1:]=dxi
    v={}
    for k,a in y.items():
        grads=model.grad(a);v[k]=sum(xi[...,i].reshape(xi.shape[:3]+(1,)*(a.ndim-3))*np.take(grads,i,axis=3) for i in range(3))
    for k in ('g','gd'):v[k]+=np.einsum('...ml,...ln->...mn',D,y[k])+np.einsum('...ml,...nl->...mn',y[k],D)
    v['pi']+=div[...,None]*y['pi']
    v['A']+=np.einsum('...ja,...ij->...ia',y['A'],dxi)
    v['E']+=div[...,None,None]*y['E']-np.einsum('...ja,...ji->...ia',y['E'],dxi)
    return xi,v

def run():
    with ev.ResearchRuntime(ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        N=24;model=ev.Model(N,old);analytic=tangent.AnalyticModel(N,old);y,_=model.initial();point=(0,6,3);rows=[]
        for t in (0.,.005):
            if t:y=ev.rk4(model,y,.005)
            F=model.rhs(y);ref=references(analytic,y);directions=[F,*[{k:np.take(model.grad(a),i,axis=3) for k,a in y.items()} for i in range(3)]]
            dref=[derivative(analytic,y,v) for v in directions]
            J={k:np.stack([d[k][point] for d in dref],axis=-1) for k in ref}
            xi,v=diffeo(model,y);gauge=derivative(analytic,y,v)
            fd=[]
            for eps in (1e-3,5e-4,2.5e-4):
                rp=references(analytic,ev.add(y,F,eps));rm=references(analytic,ev.add(y,F,-eps))
                fd.append({k:float(np.max(abs((rp[k][point]-rm[k][point]).real/(2*eps)-J[k][:,0]))) for k in ref})
            row=dict(time=t,point=[0.,float(np.pi/2),float(np.pi/4)],clock_margin=float(-ref['old'][point][2].real),data={})
            for k in ref:
                residual=gauge[k][point]-J[k][:,1:]@xi[point]
                row['data'][k]=dict(value=ref[k][point].real.tolist(),Jacobian=J[k].tolist(),determinant=float(np.linalg.det(J[k])),
                    spatial_diffeomorphism_scalar_residual=residual.tolist(),recovered_spatial_generator_error=(np.linalg.solve(J[k],gauge[k][point])-np.r_[0.,xi[point]]).tolist(),
                    time_direction_central_difference_errors=[q[k] for q in fd])
            rows.append(row)
    return dict(round=907,status='working',formal_previous=906,cumulative_previous=3691,N=N,rows=rows,
        no_frozen_reference_replaced=True,actual_old_and_new_references_both_used=True,
        auxiliary_Jacobian_not_a_new_physical_field=True,full_physical_projection_implemented=False,
        full_smeared_loop_source_implemented=False,patch_uniform_bound_certified=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();v=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert v==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(v,ensure_ascii=False,indent=2))
