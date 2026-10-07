"""Independent finite-amplitude ODE check of921; same source/field identities."""
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
import kinematic_defect_readout as test
c=test.c;r=test.r;mb=test.mb;response=test.response
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;TARGET=HERE/'kinematic_response_validation_results.json'

def full_tangent(fam,forces):
    u={k:np.zeros_like(v) for k,v in fam.y(0.).items()};t=0.;dt=fam.T/2
    def rhs(t,u):
        L=fam.analytic.jvp(fam.y(t),u);f=test.interpolate_force(forces,t,fam.T)
        return {k:L[k]+f[k] for k in u}
    for _ in range(2):
        s=[]
        for a,j in ((0.,None),(.5,0),(.5,1),(1.,2)):
            v=u if j is None else {k:u[k]+a*dt*s[j][k] for k in u}
            s.append(rhs(t+a*dt,v))
        u={k:u[k]+dt/6*sum(w*q[k] for w,q in zip((1,2,2,1),s)) for k in u};t+=dt
    return u

def finite_ode(fam,forces,amplitude):
    # W'=F(y_hat(t)+W)-F(y_hat(t))+amplitude*f(t) has the921 variational
    # equation exactly at amplitude0, without pretending y_hat solves F.
    u={k:np.zeros_like(v) for k,v in fam.y(0.).items()};t=0.;dt=fam.T/2
    def rhs(t,u):
        y=fam.y(t);plus=fam.model.rhs({k:y[k]+u[k] for k in y});base=fam.model.rhs(y)
        f=test.interpolate_force(forces,t,fam.T)
        return {k:plus[k]-base[k]+amplitude*f[k] for k in u}
    for _ in range(2):
        s=[]
        for a,j in ((0.,None),(.5,0),(.5,1),(1.,2)):
            v=u if j is None else {k:u[k]+a*dt*s[j][k] for k in u}
            s.append(rhs(t+a*dt,v))
        u={k:u[k]+dt/6*sum(w*q[k] for w,q in zip((1,2,2,1),s)) for k in u};t+=dt
    return u

def run():
    start=time.time();bg,original,_=c.ex.previous.rebuilt.build_pair()
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        fam=test.Family(bg.T,old);forces=[test.projected_residual(fam,bg,original,t,25)[0] for t in (-bg.T,0.,bg.T)]
        u=full_tangent(fam,forces);rows=[]
        for amplitude in (1024.,512.):
            plus=finite_ode(fam,forces,amplitude);minus=finite_ode(fam,forces,-amplitude)
            err={k:c.maximum((plus[k]-minus[k])/(2*amplitude)-u[k]) for k in u}
            rows.append(dict(amplitude=amplitude,finite_parameter_max_deformation=max(c.maximum(a) for a in plus.values()),
              maximum_error_by_field=err,maximum_error=max(err.values())))
        correction,_=test.integrate(fam,forces)
        grid=test.prior.weak.GridField(correction,bg.T,17).get()
        field_error=c.maximum(grid-mb.pack(u).reshape((-1,58)))
    assert max(row['maximum_error'] for row in rows)<1e-15,rows
    assert field_error<1e-18
    src=c.ex.loop.LoopSource(bg,2,8);src.kernels();weighted=c.aj.prior.weighted
    uu,zz,_=weighted.receiver_modes();w,dw,_=weighted.weight_jets(uu,zz,src.points)
    b,da,_=c.aj.extractor(bg,correction,src.points)
    xi,dxi,J,p,v,alpha,rho,residual=b;base=(xi,dxi,p,v,alpha,rho,residual)
    weighted_value,gauge,principal=c.aj.prior.ward_pair(src,base,da,w,dw)
    constant,_,_=c.aj.prior.ward_pair(src,base,da,np.ones_like(w),np.zeros_like(dw))
    bare=np.array(src.response(lambda x,p:correction.jets(x))['total'])
    assert c.maximum(constant-bare)<1e-20
    old=json.loads((HERE/'kinematic_defect_readout_results.json').read_text('utf-8'))
    oldv=np.array(old['original_weighted_Ward_source']['real'])+1j*np.array(old['original_weighted_Ward_source']['imag'])
    same=c.maximum(weighted_value-oldv);assert same<1e-16
    raw=np.array([np.sum(c.ex.ms.pair(s,c.ex.oldref.multiply(v,w,dw))) for s in src.sources])
    assert c.maximum(raw+principal-weighted_value)<1e-18
    show=weighted.complex_list
    return dict(round=921,date='2026-10-06',all_checks_passed=True,
      independent_finite_amplitude_rows=rows,integrated_field_final_value_error=field_error,
      constant_weight_original_source_identity_error=c.maximum(constant-bare),weighted_original_source_reproduction_error=same,
      weighted_raw_source=show(raw),weighted_projection_principal_correction=show(principal),
      actual_extractor_sample_max={k:c.maximum(a) for k,a in dict(xi=xi,dxi=dxi,alpha=alpha,dalpha=da).items()},
      full_observable_error_certified=False,physical_GP_inverse_certified=False,
      large_validation_amplitude_is_numerical_check_not_physical_preparation=True,
      source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'kinematic_defect_readout.py')},
      elapsed_seconds=round(time.time()-start,3))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))
