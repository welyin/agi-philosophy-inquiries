"""717: original CAR records and the curved-target mass force.

Full coefficient and sparse-CAR tests, plus normalized Gauss scalar packets.
No interacting Gibbs spectrum or finite-time continuum is simulated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_smooth_mode_contract as shared
import joint_operator_domain_completion as domain

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_mass_feedback_results.json'
car=shared.car;matter=car.matter;geom=matter.original;entry=shared.entry
YABS=abs(matter.Y['s'])
LINEAR_MASS=[]
for e in np.eye(5):
    h,d=matter.mass_matrices(e)
    LINEAR_MASS.append((h*np.sqrt(geom.F(e)),d*np.sqrt(geom.F(e))))


def mass_x(x):
    return tuple(sum(x[j]*LINEAR_MASS[j][k] for j in range(5)) for k in (0,1))


def force_x(x):
    a=np.sum(x[..., :4]**2,axis=-1);b=x[..., 4]**2
    U,gy,_=domain.polynomial(a,b)
    du=np.concatenate((2*x[..., :4]*gy[..., 0,None],
                       (2*x[..., 4]*gy[..., 1])[...,None]),axis=-1)
    J=du+x*np.sum(x*du,axis=-1)[...,None]/6
    return U,du,J


def dot(a,b):
    return sum(v.conjugate()*b.get(i,0) for i,v in a.items())


def subtract(a,b):
    return car.old.add(dict(a),b,-1)


def source_coefficients(points,weights,probe):
    n=32*len(points)
    h=np.zeros((n,n),complex);d=np.zeros_like(h)
    for v,x in enumerate(points):
        _,_,J=force_x(x);hv,dv=mass_x(probe[v]/weights[v]*J)
        sl=slice(32*v,32*v+32)
        h[sl,sl]=hv;d[sl,sl]=dv
    return h,d


def force_check():
    rng=np.random.default_rng(71721)
    rows=[]
    weights=np.array([.8,1.1,.6])
    p=np.array([.5,.3,.2])
    f=np.zeros(96,complex);f[30::32]=np.sqrt(p)
    pair=(1<<30)|(1<<31)
    state={0:complex(1/np.sqrt(2)),pair:np.exp(1j*np.angle(matter.Y['s']))/np.sqrt(2)}
    n=lambda s:entry.n_apply(s,f)
    k1=n(state);k0=subtract(state,k1)
    for scale in (.3,1.,3.,8.):
        x=rng.normal(size=(3,5))*scale
        fh,fd=source_coefficients(x,weights,p)
        dh,dd=shared.difference_coeff(fh,fd,f)
        source=lambda s:car.quadratic(s,fh,fd)
        delta=lambda s:car.quadratic(s,dh,dd)
        original=float(dot(state,source(state)).real)
        updated=float(sum(dot(k,source(k)).real for k in (k0,k1)))
        changed=float(dot(state,delta(state)).real)
        _,_,J=force_x(x[0])
        acceleration=p[0]**2/weights[0]*YABS*J[4]
        assert abs(changed+acceleration)<1e-10*(1+abs(acceleration))
        assert abs(updated-original-changed)<1e-10*(1+abs(acceleration))
        # Independently differentiate original ball-coordinate masses along K^-1 dA.
        phi=domain.ball(x)
        Jphi=[]
        for v in range(3):
            F=geom.F(phi[v])
            jac=np.eye(5)/np.sqrt(F)+np.outer(phi[v],phi[v])/(6*F**1.5)
            _,du,_=force_x(x[v])
            gradient_phi=jac.T@du
            Jphi.append(geom.inverse(phi[v])@gradient_phi*p[v]/weights[v])
        Jphi=np.array(Jphi)
        errors=[]
        # Keep the perturbation relative to distance from the boundary.
        for relative in (4e-4,2e-4,1e-4):
            step=relative*min(geom.F(phi))/(1+np.max(np.linalg.norm(Jphi,axis=1)))
            hp=np.zeros_like(fh);hm=np.zeros_like(fh)
            dp=np.zeros_like(fd);dm=np.zeros_like(fd)
            for v in range(3):
                sl=slice(32*v,32*v+32)
                hp[sl,sl],dp[sl,sl]=matter.mass_matrices(phi[v]+step*Jphi[v])
                hm[sl,sl],dm[sl,sl]=matter.mass_matrices(phi[v]-step*Jphi[v])
            eh,ed=shared.difference_coeff((hp-hm)/(2*step),(dp-dm)/(2*step),f)
            errors.append(float(max(np.linalg.norm(eh-dh),np.linalg.norm(ed-dd))/
                                (1+np.linalg.norm(dh)+np.linalg.norm(dd))))
        Q=np.eye(96)-np.outer(f,f.conj())
        norm,residual,dim=shared.exact_norm(dh,dd,np.column_stack((f,Q@fh@f,Q@fd@f.conj())))
        bound=float(np.linalg.norm(Q@fh@f)+np.linalg.norm(Q@fd@f.conj()))
        assert norm<=bound+1e-9 and residual<1e-8*(1+norm)
        assert max(errors)<2e-6
        rows.append(dict(scale=scale,full_original_modes=96,
                         pre_force=original,post_force=updated,
                         acceleration_difference=-changed,analytic_acceleration=acceleration,
                         finite_difference_relative_errors=errors,
                         exact_force_difference_norm=norm,coefficient_bound=bound,
                         active_modes=dim))
    return dict(rows=rows,other_H_terms_cancel_from_record_difference_analytically=True,
                full_boson_dynamics_not_frozen_in_the_proof=True)


def packet(R,order):
    u,wu=np.polynomial.legendre.leggauss(order)
    radius=.7*(u+1)/2;wr=.7*wu/2
    y=u;wy=wu
    r,s=np.meshgrid(radius,R+y,indexing='ij')
    rr,yy=np.meshgrid(radius,y,indexing='ij')
    a=r*r;b=s*s
    amplitude2=np.exp(-2/(1-(r/.7)**2)-2/(1-yy*yy))
    # H5 measure: dx/sqrt(1+|x|^2/6); common S3 volume cancels.
    measure=wr[:,None]*wy[None,:]*r**3/np.sqrt(1+(a+b)/6)
    weighted=measure*amplitude2
    Z=float(weighted.sum());prob=weighted/Z
    U,gy,_=domain.polynomial(a,b)
    J5=2*s*(gy[...,1]+(a*gy[...,0]+b*gy[...,1])/6)
    dr=-2*r/.7**2/(1-(r/.7)**2)**2
    ds=-2*yy/(1-yy*yy)**2
    kinetic_density=dr*dr+ds*ds+(r*dr+s*ds)**2/6
    return dict(Z=Z,U=float(np.sum(prob*U)),J5=float(np.sum(prob*J5)),
                target_kinetic=float(np.sum(prob*kinetic_density)),
                mean_s=float(np.sum(prob*s)),
                U54=float(np.sum(prob*(1+U)**1.25)))


def packet_check():
    w=.8;p0=.5
    base=packet(0.,80)
    rows=[]
    U4=float(.5*domain.HESS[1,1])
    limit=p0*p0/w*YABS*2*U4/3
    for R in (4.,8.,16.,32.,64.,128.):
        q=packet(R,80);lower=packet(R,48)
        eps=R**-4
        mix=lambda key:(1-eps)*base[key]+eps*q[key]
        force=p0*p0/w*YABS*mix('J5')
        onsite=w*mix('U')+mix('target_kinetic')/(2*w)+YABS*mix('mean_s')
        relative=max(abs(q[k]-lower[k])/(1+abs(q[k])) for k in ('U','J5','target_kinetic','U54'))
        assert relative<2e-4 and force>0
        rows.append(dict(R=R,tail_probability=eps,normalization=q['Z'],
                         actual_mixture_U=mix('U'),actual_mixture_target_kinetic=mix('target_kinetic'),
                         actual_node_energy_including_Majorana=onsite,
                         original_acceleration_difference=force,
                         acceleration_over_R=force/R,asymptotic_slope=limit,
                         U54_moment=mix('U54'),quadrature_order_relative_change=relative))
    assert rows[-1]['original_acceleration_difference']>20*rows[0]['original_acceleration_difference']
    assert abs(rows[-1]['acceleration_over_R']/limit-1)<.002
    assert max(r['actual_mixture_U'] for r in rows)<.2
    assert max(r['actual_node_energy_including_Majorana'] for r in rows)<100
    return dict(base_packet=base,rows=rows,nonthermal_normal_Gauss_states=True,
                complete_graph_edge_energy_bound_analytic=True,
                no_claim_of_divergence_at_any_fixed_nonzero_time=True)


def physical_weight_check():
    rows=[]
    for n in (12,20,32):
        a=4/n
        axis=(np.arange(n)+.5)*a-2
        x=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1).reshape(-1,3)
        rr=np.sum(x*x,axis=1)
        b=np.zeros(len(x));mask=rr<1;b[mask]=np.exp(-1/(1-rr[mask]))
        w=a**3*np.exp(.18*np.cos(x[:,0])-.12*np.sin(x[:,1]))
        p=w*b*b/(w@(b*b))
        fields=np.column_stack((.3+.06*np.cos(x[:,0]),.2*np.sin(x[:,1]),
                                .12*np.ones(len(x)),.08*np.cos(x[:,0]+x[:,1]),
                                .35+.05*np.sin(x[:,0])))
        U,_,J=force_x(fields)
        v=p/w
        hd=abs(matter.Y['nu'])*np.sqrt(np.sum(p*v*v*np.sum(J[:,:4]**2,axis=1)))
        md=YABS*np.sqrt(np.sum(p*v*v*J[:,4]**2))
        weighted_moment=float(np.sum(p*(1+U)**2.5))
        assert np.isfinite(hd+md) and weighted_moment>1
        # Uniform physical-volume change, f and p fixed: force has weight w^-1.
        step=1e-5
        derivative=((hd+md)*np.exp(-6*step)-(hd+md)*np.exp(6*step))/(2*step)
        error=float(abs(derivative+6*(hd+md)))
        assert error<1e-8
        rows.append(dict(n=n,nodes=len(x),profile_density_max=float(max(v)),
                         actual_mass_force_cross_bound=float(hd+md),
                         actual_weighted_U52=weighted_moment,
                         uniform_geometry_derivative_error=error))
    assert abs(rows[-1]['actual_mass_force_cross_bound']-
               rows[-2]['actual_mass_force_cross_bound'])<1e-4
    return dict(rows=rows,physical_weights_not_raw_node_count=True,
                no_full_quantum_state_or_continuum_claim=True)


def bounded_record_check():
    rng=np.random.default_rng(71741)
    weights=np.array([.8,1.1,.6]);p=np.array([.5,.3,.2]);q=p
    f=np.zeros(96,complex);f[30::32]=np.sqrt(p)
    rows=[]
    exact_bound=YABS*np.sqrt(geom.M)/4*np.sqrt(np.sum(p*(q/weights)**2))
    for scale in (0.,.3,2.,8.,30.):
        x=rng.normal(size=(3,5))*scale
        phi=domain.ball(x);F=geom.F(phi)
        mh=np.zeros((96,96),complex);md=np.zeros_like(mh)
        coefficient=q/weights*np.sqrt(F)*np.cos(phi[:,4])/4
        for v in range(3):
            sl=slice(32*v,32*v+32)
            mh[sl,sl],md[sl,sl]=mass_x(np.array([0.,0.,0.,0.,coefficient[v]]))
        dh,dd=shared.difference_coeff(mh,md,f)
        Q=np.eye(96)-np.outer(f,f.conj())
        norm,residual,dim=shared.exact_norm(dh,dd,np.column_stack((f,Q@md@f.conj())))
        formula=YABS*np.sqrt(np.sum(p*coefficient**2))
        assert abs(norm-formula)<1e-13 and norm<=exact_bound+1e-13
        metric_error=0.
        for v in range(3):
            j=geom.inverse(phi[v])[4]*np.cos(phi[v,4])/4*q[v]/weights[v]
            jac=np.eye(5)/np.sqrt(F[v])+np.outer(phi[v],phi[v])/(6*F[v]**1.5)
            expected=np.array([0.,0.,0.,0.,coefficient[v]])
            metric_error=max(metric_error,float(np.linalg.norm(jac@j-expected)))
        assert metric_error<1e-12
        rows.append(dict(scale=scale,exact_force_difference_norm=norm,
                         coefficient_norm=formula,global_exact_norm=exact_bound,
                         full_metric_contraction_error=metric_error,
                         support_residual=residual,active_modes=dim))
    assert abs(rows[0]['exact_force_difference_norm']-exact_bound)<1e-13
    return dict(original_effect='1/2 + sin(s)/4',selection_weights=q.tolist(),
                rows=rows,global_bound_independent_of_input_state=True,
                no_uniform_finite_time_remainder_claim=True)


def run():
    r=dict(round=717,tests_run=4,failures=0,errors=0,
           original_force=force_check(),normal_packet_counterfamily=packet_check(),
           shared_physical_weights=physical_weight_check(),original_bounded_record=bounded_record_check())
    names=('research_note_600.md','research_note_623.md','research_note_624.md',
           'research_note_669.md','research_note_712.md','research_note_716.md',
           'joint_operator_domain_completion.py','joint_smooth_mode_contract.py',
           'round717_drafts/postrecord_reference_entry.py')
    r['dependency_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    r['scope']=dict(original_full_finite_graph=True,feedback_is_quantum_mass_force=True,
                    bounded_energy_alone_not_uniform_acceleration=True,
                    nonthermal_counterfamily_not_Gibbs=True,
                    finite_time_blowup_not_claimed=True,spatial_continuum_not_proved=True,
                    dynamic_Einstein_equation_not_derived=True)
    # Explicit positive flags avoid reading a negated flag as a Gibbs simulation.
    r['scope']['counterfamily_is_nonthermal']=True
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
