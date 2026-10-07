"""919: separate continuum reconstruction Gauss from canonical finite Gauss.
No change of state, source, evolution grid or physics. Identify the nonlinear
interpolation/derivative commutator using the actual 912 initial data and flow.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse,json,hashlib,time
import numpy as np
import weak_gauss_transport as weak
c=weak.c;mb=weak.mb;density=weak.density
response=c.ex.previous.rebuilt.response;r=response.r
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;TARGET=HERE/'canonical_gauss_commutator_results.json'

def momenta(p,par):
    z=c.geometry(p);ig=z['ig'];di=z['di'];mu=z['mu'];A=p['A']
    F=p['dA']-p['dA'].swapaxes(1,2)+c.bracket(A[:,:,None,:],A[:,None,:,:])
    dF=p['ddA']-p['ddA'].swapaxes(2,3)+c.bracket(p['dA'][:,:,:,None,:],A[:,None,None,:,:])+c.bracket(A[:,None,:,None,:],p['dA'][:,:,None,:,:])
    Fu=np.einsum('bmr,bns,brsa->bmna',ig,ig,F)
    dFu=np.einsum('bkmr,bns,brsa->bkmna',di,ig,F)+np.einsum('bmr,bkns,brsa->bkmna',ig,di,F)+np.einsum('bmr,bns,bkrsa->bkmna',ig,ig,dF)
    lm=.5*np.einsum('bmn,brmn->br',ig,p['dg'])
    E=mu[:,None,None]*par.K*Fu[:,1:,0,:]
    dE=mu[:,None,None,None]*par.K*(dFu[:,:,1:,0,:]+lm[:,:,None,None]*Fu[:,None,1:,0,:])
    G,*_=c.target(p['phi'],par)
    rp=np.einsum('aij,bj->bai',c.REP,p['phi'])
    D=p['dphi']+np.einsum('bma,bai->bmi',A,rp)
    pi=-mu[:,None]*np.einsum('bm,bAB,bmB->bA',ig[:,0,:],G,D)
    gauss=np.einsum('biia->ba',dE[:,1:,:,:])+sum(c.bracket(A[:,i+1,:],E[:,i,:]) for i in range(3))+np.einsum('bA,baA->ba',pi,rp)
    return dict(E=E,dE=dE,pi=pi,gauss=gauss)

def snapshot_field(y,dy,k,T):
    val=mb.pack(y);vel=mb.pack(dy);N=val.shape[0]
    coeff=np.stack((val,T*vel,np.zeros_like(val),np.zeros_like(val)),axis=-2)
    return SimpleNamespace(k=k.reshape((-1,3)),T=T,c=np.fft.fftn(coeff,axes=(0,1,2)).reshape((-1,4,58))/N**3)

def inspect(model,analytic,y,u,A,par,time_value,bg,field):
    dy,du,dA=response.rhs(model,analytic,y,u,A)
    snap=snapshot_field(y,dy,model.k,bg.T);vsnap=snapshot_field(A,dA,model.k,bg.T)
    sg=weak.GridField(snap,0.,model.N);vg=weak.GridField(vsnap,0.,model.N)
    hg=weak.GridField(bg,time_value,model.N);hv=weak.GridField(field,time_value,model.N)
    flat=lambda a:a.reshape((-1,)+a.shape[3:])
    yy={k:flat(v) for k,v in y.items()};aa={k:flat(v) for k,v in A.items()}
    gdE=flat(model.grad(A['E']));gE=flat(model.grad(y['E']))
    discrete=np.einsum('biia->ba',gdE)
    discrete+=sum(c.bracket(aa['A'][:,i,:],yy['E'][:,i,:])+c.bracket(yy['A'][:,i,:],aa['E'][:,i,:]) for i in range(3))
    discrete+=np.einsum('bA,aAB,bB->ba',aa['pi'],c.REP,yy['phi'])+np.einsum('bA,aAB,bB->ba',yy['pi'],c.REP,aa['phi'])
    check=response.tangent.constraints(analytic,{k:v.astype(complex)+1j*1e-24*A[k] for k,v in y.items()})['Gauss'].imag/1e-24
    direct_discrete_error=c.maximum(discrete-flat(check))
    arrays={k:[] for k in ('canonical_Gauss','reconstructed_Gauss','derivative_commutator','algebraic_momentum_mismatch','decomposition_error','Hermite_Gauss','Hermite_minus_snapshot')}
    maxerr={k:0. for k in ('base_momentum_at_nodes','tangent_momentum_at_nodes','Euler_density_dictionary','nodal_fields','nodal_tangent_fields')}
    for lo in range(0,model.N**3,256):
        hi=min(lo+256,model.N**3);p=sg.jets(lo,hi);v=vg.jets(lo,hi);h=1e-24
        b=momenta(p,par);var=momenta({k:p[k].astype(complex)+1j*h*v[k] for k in p},par);d={k:z.imag/h for k,z in var.items()}
        Eb=b['E']-yy['E'][lo:hi];dEb=d['E']-aa['E'][lo:hi]
        pib=b['pi']-yy['pi'][lo:hi];dpib=d['pi']-aa['pi'][lo:hi]
        comm=np.einsum('biia->ba',d['dE'][:,1:,:,:]-gdE[lo:hi])
        algebra=sum(c.bracket(aa['A'][lo:hi,i,:],Eb[:,i,:])+c.bracket(yy['A'][lo:hi,i,:],dEb[:,i,:]) for i in range(3))
        algebra+=np.einsum('bA,aAB,bB->ba',dpib,c.REP,yy['phi'][lo:hi])+np.einsum('bA,aAB,bB->ba',pib,c.REP,aa['phi'][lo:hi])
        closure=d['gauss']-discrete[lo:hi]-comm-algebra
        hp=hg.jets(lo,hi);vv=hv.jets(lo,hi)
        shifted=c.euler({k:hp[k].astype(complex)+1j*h*vv[k] for k in hp},par)
        H=density.density(shifted,shifted['z'])['YM'][:,0,:].imag/h
        shifted_snapshot=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
        canonical_density=density.density(shifted_snapshot,shifted_snapshot['z'])['YM'][:,0,:].imag/h
        for key,value in dict(canonical_Gauss=discrete[lo:hi],reconstructed_Gauss=d['gauss'],derivative_commutator=comm,
          algebraic_momentum_mismatch=algebra,decomposition_error=closure,Hermite_Gauss=H,Hermite_minus_snapshot=H-d['gauss']).items(): arrays[key].append(value)
        maxerr['base_momentum_at_nodes']=max(maxerr['base_momentum_at_nodes'],c.maximum(Eb),c.maximum(pib))
        maxerr['tangent_momentum_at_nodes']=max(maxerr['tangent_momentum_at_nodes'],c.maximum(dEb),c.maximum(dpib))
        maxerr['Euler_density_dictionary']=max(maxerr['Euler_density_dictionary'],c.maximum(canonical_density-d['gauss']))
        maxerr['nodal_fields']=max(maxerr['nodal_fields'],c.maximum(p['g']-yy['g'][lo:hi]),c.maximum(p['phi']-yy['phi'][lo:hi]),c.maximum(p['A'][:,1:,:]-yy['A'][lo:hi]))
        maxerr['nodal_tangent_fields']=max(maxerr['nodal_tangent_fields'],c.maximum(v['g']-aa['g'][lo:hi]),c.maximum(v['phi']-aa['phi'][lo:hi]),c.maximum(v['A'][:,1:,:]-aa['A'][lo:hi]))
    arrays={k:np.concatenate(v) for k,v in arrays.items()}
    norms={k:dict(max=c.maximum(v),RMS=float(np.sqrt(np.mean(v*v)))) for k,v in arrays.items()}
    assert direct_discrete_error<1e-12 and maxerr['Euler_density_dictionary']<1e-12
    assert norms['decomposition_error']['max']<1e-12
    return dict(time=time_value,N=model.N,derivative_commutator_decomposition=norms,nodal_dictionary_errors=maxerr,
      independent_canonical_constraint_error=direct_discrete_error,physical_error_certified=False)

def run():
    start=time.time();bg,field,_=c.ex.previous.rebuilt.build_pair();rows=[]
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(17,old);analytic=response.tangent.AnalyticModel(17,old)
        par=SimpleNamespace(K=model.K,L=model.L,u=model.u)
        y0,u0,A0,_=response.initial(model,analytic,old)
        rows.append(inspect(model,analytic,y0,u0,A0,par,0.,bg,field));print(json.dumps(rows[-1]),flush=True)
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();A={k:v.copy() for k,v in A0.items()}
            for _ in range(2):y,u,A=response.step(model,analytic,y,u,A,sign*bg.T/2)
            rows.append(inspect(model,analytic,y,u,A,par,sign*bg.T,bg,field));print(json.dumps(rows[-1]),flush=True)
    return dict(round=919,date='2026-10-06',same912_preparation_and_flow=True,rows=rows,
      physical_initial_surface_time=0.,Hermite_left_endpoint_is_not_preparation_time=True,
      diagnostic_identity_not_physical_error_certificate=True,full872_response_computed=False,full_goal_completed=False,
      elapsed_seconds=round(time.time()-start,3),source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
       (Path(__file__),HERE/'weak_gauss_transport.py',STAGE/'912/receiver_forced_response.py',STAGE/'904/coupled_boson_tangent.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'saved':a.write,'elapsed_seconds':out['elapsed_seconds']}))
