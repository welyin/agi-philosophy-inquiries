"""707: one nonlinear scalar spatial block of the full original Gauss model.
Analytic theorem keeps all links/CAR and the relative S4 carrier, tracing only
a gauge-trivial radial factor. Numerics check actual H5, original full mass
matrices and the new midpoint instrument; no full Gibbs spectrum is computed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_tree_matter_source as tree
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_geodesic_spatial_block_results.json'
original=tree.original;matter=tree.matter
R=np.sqrt(6.);M=original.M
ETA=np.diag([-1.,1.,1.,1.,1.,1.])

def hyper(phi):
    F=original.F(phi)
    return np.r_[np.sqrt(6*M/F),phi/np.sqrt(F)]

def klein(X):
    return np.sqrt(6*M)*X[1:]/X[0]

def frame(m):
    X=hyper(m)
    return np.vstack((X[1:]/R,np.eye(5)+np.outer(X[1:],X[1:])/(R*(X[0]+R))))

def endpoints(m,z):
    X=hyper(m);a=np.linalg.norm(z)/R
    sinc=1. if a<1e-12 else np.sinh(a)/a
    tangent=frame(m)@z
    return klein(np.cosh(a)*X-sinc*tangent),klein(np.cosh(a)*X+sinc*tangent)

def midpoint(x,y):
    ax=1/np.sqrt(original.F(x));ay=1/np.sqrt(original.F(y))
    return (ax*x+ay*y)/(ax+ay)

def inverse_pair(x,y):
    m=midpoint(x,y);X=hyper(m);Y=hyper(y)
    a=np.arccosh(max(1.,float(-X@ETA@Y/6)))
    factor=1. if a<1e-9 else a/np.sinh(a)
    z=frame(m).T@ETA@(factor*(Y-np.cosh(a)*X))
    return m,z

def density(x):
    return np.sqrt(M)/original.F(x)**3

def jacobian(r):
    a=2*r/R
    return 32*(1. if abs(a)<1e-12 else np.sinh(a)/a)**4

def grad(f,x,step=2e-6):
    eye=np.eye(len(x))*step
    return np.array([(f(x+e)-f(x-e))/(2*step) for e in eye])

def geometry_check():
    rng=np.random.default_rng(707)
    rows=[];inverse=equiv=edge=frame_err=0.
    for radius in (0.,.35,1.2,2.4):
        m=.18*rng.normal(size=5);z=rng.normal(size=5);z*=radius/np.linalg.norm(z)
        x,y=endpoints(m,z);mm,zz=inverse_pair(x,y)
        inverse=max(inverse,float(np.linalg.norm(mm-m)),float(np.linalg.norm(zz-z)))
        edge=max(edge,abs(float(original.distance_squared(x,y))-4*radius*radius))
        frame_err=max(frame_err,float(np.linalg.norm(frame(m).T@ETA@frame(m)-np.eye(5))))
        q=np.r_[m,z]
        def fun(q):return np.concatenate(endpoints(q[:5],q[5:]))
        errors=[];metric_errors=[]
        for step in (2e-5,1e-5):
            E=np.eye(10)*step
            D=np.stack([(fun(q+e)-fun(q-e))/(2*step) for e in E],axis=1)
            actual=abs(np.linalg.det(D))*density(x)*density(y)/density(m)
            errors.append(abs(actual/jacobian(radius)-1))
            # Independent horizontal/vertical Jacobi-field metric, including frame connection.
            dm=np.eye(5)*step;Ef=frame(m)
            DX=np.stack([(hyper(m+e)-hyper(m-e))/(2*step) for e in dm],axis=1)
            DE=[(frame(m+e)-frame(m-e))/(2*step) for e in dm]
            H=np.c_[Ef.T@ETA@DX,np.zeros((5,5))]
            V=np.c_[np.stack([Ef.T@ETA@d@z for d in DE],axis=1),np.eye(5)]
            unit=z/radius if radius else np.eye(5)[0]
            parallel=np.outer(unit,unit);perp=np.eye(5)-parallel;a=radius/R
            av=1. if radius==0 else np.sinh(a)/a
            expected=2*H.T@(parallel+np.cosh(a)**2*perp)@H+2*V.T@(parallel+av**2*perp)@V
            G=np.zeros((10,10));G[:5,:5]=original.metric(x);G[5:,5:]=original.metric(y)
            pulled=D.T@G@D
            metric_errors.append(float(np.linalg.norm(pulled-expected)/np.linalg.norm(expected)))
        g=tree.old.sample(rng)
        gm,gz=tree.scalar_action(g,m),tree.scalar_action(g,z)
        xx,yy=endpoints(gm,gz)
        equiv=max(equiv,float(np.linalg.norm(xx-tree.scalar_action(g,x))),
                  float(np.linalg.norm(yy-tree.scalar_action(g,y))))
        rows.append(dict(relative_radius=radius,jacobian=jacobian(radius),
                         independent_volume_relative_errors=errors,
                         full_pullback_metric_relative_errors=metric_errors,
                         minimum_F=float(min(original.F(x),original.F(y)))))
    assert max(inverse,equiv,edge,frame_err)<2e-12
    assert max(v for r in rows for v in r['independent_volume_relative_errors'])<2e-7
    assert max(v for r in rows for v in r['full_pullback_metric_relative_errors'])<2e-7
    radial=[]
    for r in (.2,.7,1.4):
        step=1e-4
        logroot=lambda t:np.log(jacobian(t)*t**4)/2
        first=(logroot(r+step)-logroot(r-step))/(2*step)
        second=(logroot(r+step)+logroot(r-step)-2*logroot(r))/step**2
        prediction=16/R**2+8/(R*np.sinh(2*r/R))**2
        error=abs(first*first+second-prediction)
        assert error<2e-5
        radial.append(dict(radius=r,flat_radial_measure_potential=prediction,independent_difference_error=error))
    return dict(rows=rows,inverse_error=inverse,full_quotient_equivariance_error=equiv,
                original_edge_distance_error=edge,parallel_frame_error=frame_err,
                radial_measure_checks=radial)

def instrument_check():
    rng=np.random.default_rng(7071);w1=.8;w2=1.1;hbar=.7
    bound=hbar*hbar*M/128*(1/w1+1/w2)
    rows=[];electric=identity=gradient_bound=0.
    samples=[(np.zeros(5),np.zeros(5))]
    for radius in (.25,.8,1.5,2.5):
        m=.2*rng.normal(size=5);u=rng.normal(size=5);u/=np.linalg.norm(u)
        samples.append(endpoints(m,radius*u))
    for x,y in samples:
        S=float(midpoint(x,y)[4])
        gx=grad(lambda xx:midpoint(xx,y)[4],x)
        gy=grad(lambda yy:midpoint(x,yy)[4],y)
        nx=float(gx@original.inverse(x)@gx);ny=float(gy@original.inverse(y)@gy)
        mm,zz=inverse_pair(x,y);rr=np.linalg.norm(zz);uu=zz/rr if rr>1e-10 else np.eye(5)[0]
        XX=hyper(mm);direction=frame(mm)@uu
        along=np.sqrt(6*M)*(direction[5]*XX[0]-XX[5]*direction[0])/XX[0]**2
        whole=original.inverse(mm)[4,4]
        exact_grad=(along**2+(whole-along**2)/np.cosh(rr/R)**2)/4
        assert max(abs(nx-exact_grad),abs(ny-exact_grad))<2e-9
        # Exact scalar action leaves endpoint F and singlet separately invariant.
        g=tree.old.sample(rng)
        electric=max(electric,abs(float(midpoint(x,tree.scalar_action(g,y))[4])-S))
        L=np.sqrt(.5+np.array([1,-1])*.25*np.sin(S))
        dL=np.array([1,-1])*np.cos(S)/(8*L)
        ell=float(dL@dL)
        identity=max(identity,abs(ell-np.cos(S)**2/(4*(4-np.sin(S)**2))),abs(float(L@L)-1))
        inj=hbar*hbar*ell/2*(nx/w1+ny/w2)
        gradient_bound=max(gradient_bound,nx-M/4,ny-M/4)
        assert inj<=bound+1e-10
        rows.append(dict(midpoint_s=S,arithmetic_s=float((x[4]+y[4])/2),
            energy_injection=inj,exact_operator_norm=bound,ratio=inj/bound,
            endpoint_gradient_norms=[nx,ny],independent_geodesic_gradient_squared=float(exact_grad)))
    assert max(electric,identity,gradient_bound)<1e-9
    assert abs(rows[0]['ratio']-1)<1e-9
    assert max(abs(r['midpoint_s']-r['arithmetic_s']) for r in rows)>1e-4
    return dict(rows=rows,original_link_derivative_zero_error=electric,
                completeness_and_derivative_identity_error=identity,
                endpoint_CAT0_bound_excess=gradient_bound,
                exact_norm_attained_at_coincident_origin_analytically=True,
                no_electric_injection_even_though_full_midpoint_is_link_covariant=True)

def mass_check():
    unit=[]
    for e in np.eye(5)*.2:
        h,d=matter.mass_matrices(e)
        unit.append((h*np.sqrt(original.F(e))/.2,d*np.sqrt(original.F(e))/.2))
    def linear(y):
        return tuple(sum(y[j]*unit[j][k] for j in range(5)) for k in range(2))
    m=np.array([.25,-.18,.09,.23,.41])
    u=np.array([.4,.1,-.6,.2,.65]);u/=np.linalg.norm(u)
    E=frame(m)@u;hm,dm=matter.mass_matrices(m);hd,dd=linear(E[1:])
    O=np.block([[np.eye(32),np.eye(32)],[np.eye(32),-np.eye(32)]])/np.sqrt(2)
    rows=[];error=0.;herm=0.;pairing=0.
    for radius in (.2,.7,1.4):
        x,y=endpoints(m,radius*u);hx,dx=matter.mass_matrices(x);hy,dy=matter.mass_matrices(y)
        a=radius/R
        h_avg=np.cosh(a)*hm;d_avg=np.cosh(a)*dm
        h_diff=-R*np.sinh(a)*hd;d_diff=-R*np.sinh(a)*dd
        hpair=np.block([[hx,np.zeros((32,32))],[np.zeros((32,32)),hy]])
        dpair=np.block([[dx,np.zeros((32,32))],[np.zeros((32,32)),dy]])
        expected_h=np.block([[h_avg,h_diff],[h_diff,h_avg]])
        expected_d=np.block([[d_avg,d_diff],[d_diff,d_avg]])
        residual=max(np.linalg.norm(O@hpair@O.T-expected_h),
                     np.linalg.norm(O@dpair@O.T-expected_d))
        error=max(error,float(residual))
        herm=max(herm,float(np.linalg.norm(expected_h-expected_h.conj().T)))
        pairing=max(pairing,float(np.linalg.norm(expected_d+expected_d.T)))
        naive_h=np.kron(np.eye(2),hm);naive_d=np.kron(np.eye(2),dm)
        rows.append(dict(relative_radius=radius,shared_mass_factor=float(np.cosh(a)),
            complete64mode_identity_error=float(residual),
            contrast_Dirac_norm=float(np.linalg.norm(h_diff)),
            contrast_Majorana_norm=float(np.linalg.norm(d_diff)),
            naive_midpoint_only_mass_error=float(np.sqrt(np.linalg.norm(hpair-naive_h)**2+np.linalg.norm(dpair-naive_d)**2))))
    assert max(error,herm,pairing)<5e-13
    assert all(r['contrast_Dirac_norm']>0 and r['contrast_Majorana_norm']>0 for r in rows)
    return dict(rows=rows,complete_mass_identity_error=error,hermiticity_error=herm,
                pairing_antisymmetry_error=pairing,
                original_complex_Y_and_both32mode_sites_preserved=True,
                no_full_H_thermal_spectrum_numerically_claimed=True)

def run():
    deps=('research_note_574.md','research_note_598.md','research_note_639.md',
          'research_note_640.md','research_note_642.md','research_note_704.md',
          'research_note_706.md','round707_drafts/spatial_map_entry.md',
          'joint_curved_quantum_source.py','joint_fermion_gauss_completion.py',
          'joint_tree_matter_source.py')
    return dict(date='2026-10-02',round=707,tests_run=3,failures=0,errors=0,
        geometry=geometry_check(),instrument=instrument_check(),mass=mass_check(),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(full_original_Gauss_scalar_block_analytic=True,
            actual_nonlinear_H5_not_Gaussian_replacement=True,
            retains_relative_S4_and_all_CAR_and_loop_variables=True,
            traces_gauge_trivial_radial_factor_only=True,
            original_Gibbs_and_geometric_preparation_derivatives_transport=True,
            new_declared_midpoint_instrument_not_coarse_grained_old_records=True,
            memory_in_real_histories_retained=True,
            no_bare_one_node_H_or_continuum_or_quantum_GR_claim=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=707,tests=3,geometry=result['geometry'],instrument=result['instrument'],
                         mass=result['mass'],all_checks_passed=True)))
