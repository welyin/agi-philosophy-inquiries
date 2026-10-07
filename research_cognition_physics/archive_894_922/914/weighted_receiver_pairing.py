"""914: homogeneous Dirac test-mode pairing with original weighted source.
Uses natural curvature/clock variations to avoid a full projected connection jet.
The stored values are numerical interface diagnostics, not continuum bounds.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import material_extractor_naturality as ex
previous=ex.previous;work=previous.work;response=work.response;r=work.r;mb=ex.mb;loop=ex.loop;ms=ex.ms
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;TARGET=HERE/'weighted_receiver_pairing_results.json'

class SpinorInterpolant:
    def __init__(self,N,T,minus,plus,waves):
        ym,dm=minus;yp,dp=plus;self.T=T;self.k=waves.reshape((-1,3))
        co=np.stack(((yp+ym)/2-T*(dp-dm)/4,3*(yp-ym)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(yp-ym))/4),axis=-2)
        self.co=np.fft.fftn(co,axes=(0,1,2)).reshape((-1,4,4))/N**3
    def evaluate(self,x,axis=None):
        t=x[:,0]/self.T
        powers=np.stack((np.ones_like(t),t,t*t,t*t*t),axis=1)
        co=self.co
        if axis==0:powers=np.stack((t*0,t*0+1,2*t,3*t*t),axis=1)/self.T
        elif axis is not None:co=co*(1j*self.k[:,axis-1,None,None])
        values=[]
        for a in range(0,len(x),128):
            q=np.exp(1j*x[a:a+128,1:]@self.k.T)@co.reshape((-1,16))
            values.append(np.einsum('bp,bpf->bf',powers[a:a+128],q.reshape((-1,4,4))))
        return np.concatenate(values)

def receiver_modes(N=17,T=.00025,steps=2):
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(N,old);y0,_=model.initial();u0=r.initial(N)
        # Auxiliary homogeneous test solution, not a new physical preparation:
        # initial datum z=(original neutral p/.02)u=sin(x)u.
        z0=(y0['phi'][...,5]/.02)[...,None]*u0
        def rhs(y,u,z):
            dy,g,_=model.rhs(y,aux=True);fr=r.shared.frame_data(y,g)
            return dy,-1j*r.apply(model,y,g,fr,u),-1j*r.apply(model,y,g,fr,z)
        ends=[];norms=[];overlap0=(2*np.pi)**3*np.mean(np.sum(u0.conj()*z0,axis=-1))
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();z=z0.copy()
            dt=sign*T/steps
            for _ in range(steps):
                stages=[]
                for c,i in ((0,None),(.5,0),(.5,1),(1,2)):
                    yy=y if i is None else r.ev.add(y,stages[i][0],c*dt)
                    uu=u if i is None else u+c*dt*stages[i][1]
                    zz=z if i is None else z+c*dt*stages[i][2]
                    stages.append(rhs(yy,uu,zz))
                w=(1,2,2,1)
                y={k:y[k]+dt/6*sum(a*s[0][k] for a,s in zip(w,stages)) for k in y}
                u=u+dt/6*sum(a*s[1] for a,s in zip(w,stages));z=z+dt/6*sum(a*s[2] for a,s in zip(w,stages))
            dy,du,dz=rhs(y,u,z);ends.append(((u,du),(z,dz)))
            overlap=(2*np.pi)**3*np.mean(np.sum(u.conj()*z,axis=-1))
            norms.append(dict(sign=sign,u_norm_drift=abs(r.norm(u)-r.norm(u0)),z_norm_drift=abs(r.norm(z)-r.norm(z0)),overlap_drift=float(abs(overlap-overlap0))))
        return SpinorInterpolant(N,T,ends[0][0],ends[1][0],model.k),SpinorInterpolant(N,T,ends[0][1],ends[1][1],model.k),dict(test_initial_norm=r.norm(z0),propagation=norms)

def weight_jets(u,z,points):
    uu=u.evaluate(points);zz=z.evaluate(points)
    bil=lambda a,b:np.einsum('bi,ij,bj->b',a.conj(),r.BETA,b)
    den=bil(uu,uu);num=bil(zz,uu);w=num/den;dw=[]
    for mu in range(4):
        du=u.evaluate(points,mu);dz=z.evaluate(points,mu)
        dn=bil(dz,uu)+bil(zz,du);dd=bil(du,uu)+bil(uu,du)
        dw.append((dn-w*dd)/den)
    assert den.real.min()>0 and ex.maximum(den.imag)<1e-14
    return w,np.stack(dw,axis=1),dict(minimum_sample_calibration_density=float(den.real.min()),maximum_weight_abs=ex.maximum(w),maximum_weight_derivative=ex.maximum(np.stack(dw,axis=1)))

def extract_values(bg,field,x):
    xi,dxi,J,p,v=ex.diff_jet(bg,field,x)
    _,rho,gram,Q,extract=ex.internal_data(p)
    dz=ex.feature_delta(p,v);lie=ex.lie_feature(p,xi,dxi)
    residual={k:dz[k]-lie[k] for k in dz};alpha=extract(residual)
    return xi,dxi,p,v,alpha,rho,residual

def projected_inputs(base,dalpha):
    xi,dxi,p,v,alpha,rho,residual=base;ref=ex.oldref.reference_data(p)
    pg=v['g']-np.einsum('br,brmn->bmn',xi,p['dg'])-np.einsum('bmr,brn->bmn',dxi,p['g'])-np.einsum('bnr,bmr->bmn',dxi,p['g'])
    pa=v['A']-np.einsum('br,brma->bma',xi,p['dA'])-np.einsum('bmr,bra->bma',dxi,p['A'])-ms.bracket(alpha[:,None,:],p['A'])+dalpha
    pphi=v['phi']-np.einsum('br,bri->bi',xi,p['dphi'])-np.einsum('ba,aij,bj->bi',alpha,ex.REP,p['phi'])
    dh=np.sum(ref['r']*v['phi'][:,:4],axis=-1)
    dr=(v['phi'][:,:4]-ref['r']*dh[:,None])/ref['h'][:,None]
    ddh=np.einsum('bi,bmi->bm',dr,p['dphi'][:,:,:4])+np.einsum('bi,bmi->bm',ref['r'],v['dphi'][:,:,:4])
    hh=(np.einsum('bmi,bni->bmn',p['dphi'][:,:,:4],p['dphi'][:,:,:4])+np.einsum('bi,bmni->bmn',p['phi'][:,:4],p['ddphi'][...,:4])-ref['dh'][:,:,None]*ref['dh'][:,None,:])/ref['h'][:,None,None]
    ph=dh-np.einsum('br,br->b',xi,ref['dh'])
    pdh=ddh-np.einsum('bmr,br->bm',dxi,ref['dh'])-np.einsum('br,bmr->bm',xi,hh)
    pf=residual['F']-np.einsum('ba,bamni->bmni',alpha,rho['F'])
    return dict(g=pg,A=pa,phi=pphi,h=ph,dh=pdh,F=pf)

def weighted_pair(src,pv,w,dw):
    # j is a first-jet functional. Multiplication by w occurs before its
    # derivatives; the dw terms are essential, even though w is a scalar.
    data=src.data
    wF=w[:,None,None,None]*pv['F']+dw[:,:,None,None]*pv['A'][:,None,:,:8]-dw[:,None,:,None]*pv['A'][:,:,None,:8]
    dclock=w[:,None]*pv['dh']+dw*pv['h'][:,None]
    dM=(np.sum(data['metric']*(w[:,None,None]*pv['g']),axis=(1,2))+
        np.sum(data['clock']*dclock,axis=1)+np.sum(data['B']*wF,axis=(1,2,3)))
    dX=np.stack((w*pv['h'],w*pv['phi'][:,4],w*pv['phi'][:,5],dM),axis=-1)
    return np.array([np.sum(C*(w[:,None,None]*pv['A']))+np.sum(k*dX) for C,k in zip(src.currents,src.ks)])

def complex_list(z):return {'real':np.asarray(z).real.tolist(),'imag':np.asarray(z).imag.tolist()}

def run():
    start=time.time();bg,field,_=previous.rebuilt.build_pair();u,z,mode=receiver_modes()
    src=loop.LoopSource(bg,2,8);src.kernels();base=extract_values(bg,field,src.points)
    w,dw,stats=weight_jets(u,z,src.points);rows=[]
    # Algebraic cross-check of the reduced source formula without projection.
    p=base[2];v=base[3];ref=ex.oldref.reference_data(p)
    dh=np.sum(ref['r']*v['phi'][:,:4],axis=-1);dr=(v['phi'][:,:4]-ref['r']*dh[:,None])/ref['h'][:,None]
    ddh=np.einsum('bi,bmi->bm',dr,p['dphi'][:,:,:4])+np.einsum('bi,bmi->bm',ref['r'],v['dphi'][:,:,:4])
    raw=dict(g=v['g'],phi=v['phi'],A=v['A'],h=dh,dh=ddh,F=ex.feature_delta(p,v)['F'])
    rawreduced=weighted_pair(src,raw,w,dw)
    weighted=ex.oldref.multiply(v,w,dw)
    rawjet=np.array([np.sum(ms.pair(s,weighted)) for s in src.sources])
    assert ex.maximum(rawreduced-rawjet)<1e-15
    # Constant test recovers the previous prepared-mode identity algebraically.
    wc,dwc,constantstats=weight_jets(u,u,src.points)
    assert ex.maximum(wc-1)<1e-13 and ex.maximum(dwc)<1e-10
    for step in (2e-6,1e-6):
        da=[]
        for mu in range(4):
            xp=src.points.copy();xm=xp.copy();xp[:,mu]+=step;xm[:,mu]-=step
            ap=extract_values(bg,field,xp)[4];am=extract_values(bg,field,xm)[4]
            da.append((ap-am)/(2*step))
        da=np.stack(da,axis=1);pv=projected_inputs(base,da)
        value=weighted_pair(src,pv,w,dw)
        omitdw=weighted_pair(src,pv,w,np.zeros_like(dw))
        constant=weighted_pair(src,pv,wc,dwc)
        original=np.array(src.response(lambda x,p:field.jets(x))['total'])
        row=dict(alpha_derivative_step=step,weighted_source=complex_list(value),
          conditional_non_diagonal_vA_matrix_element=complex_list(1j*np.array([6.,6.,.75])*value),
          omit_projection=complex_list(rawjet),omit_weight_derivatives=complex_list(omitdw),
          projection_effect_diagnostic=complex_list(value-rawjet),weight_derivative_effect_diagnostic=complex_list(value-omitdw),
          constant_weight_pairing=complex_list(constant),constant_weight_minus_original_source=complex_list(constant-original),
          maximum_internal_parameter_derivative=ex.maximum(da),maximum_projected_connection=ex.maximum(pv['A']),elapsed_seconds=round(time.time()-start,3))
        rows.append(row);print(json.dumps(row,ensure_ascii=False),flush=True)
    return dict(round=914,date='2026-10-06',status='working',N=17,source_samples=len(src.points),mode=mode,weight=stats,
       rows=rows,source_identity_error=ex.maximum(rawreduced-rawjet),constant_test_weight_error=ex.maximum(wc-1),constant_test_weight_derivative_error=ex.maximum(dwc),
       auxiliary_pairing='773 full feature tuple; clock-positive inverse metric; fixed unit component weights',
       original_preparation_and_source_unchanged=True,test_mode_is_auxiliary_not_new_physical_preparation=True,
       full_projected_connection_jet_not_required=True,full_vA_field_not_required_for_this_matrix_element=True,
       full872_response_computed=False,physical_error_certified=False,full_goal_completed=False,
       source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'material_extractor_naturality.py',STAGE/'913/independent_record_rebuild.py',STAGE/'911/future_test_projection.py')})
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');a=parser.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
