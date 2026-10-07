"""826 working: constrained causal reconstruction and cutoff collars.

An exactly local discrete wave recurrence calibrates the operator identities;
it is not a discretization of the original Einstein-matter background.
"""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'causal_region_extension_probe_results.json'

def run():
    nt,nx=52,101;x=np.arange(nx)-nx//2
    def plus(a):return np.concatenate((a[1:],np.zeros_like(a[:1])),axis=0)
    def minus(a):return np.concatenate((np.zeros_like(a[:1]),a[:-1]),axis=0)
    dt=lambda a:(plus(a)-minus(a))/2
    dx=lambda a:(np.roll(a,-1,axis=1)-np.roll(a,1,axis=1))/2
    lap=lambda a:np.roll(a,-1,axis=1)+np.roll(a,1,axis=1)-2*a
    D=lambda a:dt(a)+.3*dx(a)
    T=lambda a:.4*dt(a)+.2*dx(a)
    F=lambda v:v[0]-D(v[1])
    S=lambda h:np.array([h-D(T(h)),-T(h)])
    Fadj=lambda h:np.array([h,D(h)])
    Sadj=lambda v:v[0]-T(D(v[0]))+T(v[1])
    K=lambda a:np.array([D(a),a])
    Kadj=lambda v:-D(v[0])+v[1]
    Pi=lambda v:S(F(v))
    Piadj=lambda v:Fadj(Sadj(v))
    L=lambda v:v[1]+T(F(v))
    W=lambda h:plus(h)+minus(h)-2*h-.4**2*lap(h)+.15**2*h
    P=lambda v:Fadj(W(F(v)))
    def ret(source):
        u=np.zeros_like(source)
        for t in range(1,nt-1):
            u[t+1]=source[t]+(2-.15**2)*u[t]-u[t-1]+.4**2*lap(u)[t]
        return u
    h=np.zeros((nt,nx))
    h[0]=.3*np.sin(2*np.pi*x/nx)+.2*np.cos(6*np.pi*x/nx)
    h[1]=h[0]+.07*np.cos(4*np.pi*x/nx)
    for t in range(1,nt-1):
        h[t+1]=(2-.15**2)*h[t]-h[t-1]+.4**2*lap(h)[t]
    physical=S(h)
    time=np.arange(nt);chi=np.zeros(nt);chi[time>=14]=1
    transition=(time>8)&(time<14)
    chi[transition]=.5*(1-np.cos(np.pi*(time[transition]-8)/6))
    chi=chi[:,None]
    def collar(inner,outer):
        a=abs(x);rho=np.zeros(nx);rho[a<=inner]=1
        m=(a>inner)&(a<outer)
        rho[m]=.5*(1+np.cos(np.pi*(a[m]-inner)/(outer-inner)))
        return rho
    target=(slice(22,27),abs(x)<=3)
    rows=[]
    for name,rho in [('safe',collar(30,37)),('too_narrow',collar(2,5))]:
        u=Pi(physical*rho[None,None,:])
        raw=P(u*chi[None])-P(u)*chi[None]
        source=Piadj(raw)
        out=S(ret(Sadj(source)))
        defect=out-physical
        local_error=float(np.max(abs(defect[:,target[0],:][:,:,target[1]])))
        ward=float(np.max(abs(Kadj(source))))
        raw_ward=float(np.max(abs(Kadj(raw))))
        slice_error=float(np.max(abs(L(u))))
        unprojected_slice=float(np.max(abs(L(physical*rho[None,None,:]))))
        assert ward<1e-12 and slice_error<1e-12
        assert raw_ward>1e-8 and unprojected_slice>1e-8
        if name=='safe':assert local_error<1e-11
        else:assert local_error>1e-3
        rows.append(dict(collar=name,all_target_region_field_error=local_error,
            physical_source_Ward_error=ward,unprojected_source_Ward_defect=raw_ward,
            projected_cutoff_slice_error=slice_error,
            unprojected_cutoff_slice_defect=unprojected_slice,
            source_time_support=np.where(np.max(abs(source),axis=(0,2))>1e-13)[0].tolist()))
    rng=np.random.default_rng(826);v=rng.normal(size=(2,nt,nx));a=rng.normal(size=(nt,nx))
    identities=dict(Pi_K=float(np.max(abs(Pi(K(a))))),Kadj_P=float(np.max(abs(Kadj(P(v))))),
        Pi_idempotence=float(np.max(abs(Pi(Pi(v))-Pi(v)))),
        P_Pi=float(np.max(abs(P(Pi(v))-P(v)))))
    assert max(identities.values())<1e-11
    return dict(round=826,status='working_not_formal',all_checks_passed=True,
        local_discrete_wave_with_differential_oblique_gauge=True,
        operator_identity_errors=identities,rows=rows,
        original_continuum_region_extension_proven=False,
        original_interacting_state_scope_checked=False,
        formal_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
