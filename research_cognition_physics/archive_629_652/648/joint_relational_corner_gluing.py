"""648: original nonminimal corner charge and the boost gluing interface.

Classical original573 reference/geometry and an explicitly declared regular
boost representation test. No gravitational Hilbert space or entropy derived.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gravity_material_coordinates as original
import joint_geometric_boundary_matching as boundary

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_relational_corner_gluing_results.json'


def curvature_charge(phi,g,n,r,tangents):
    F,_=boundary.data(phi);gi=np.linalg.inv(g)
    nc=g@n;rc=g@r;eps=np.outer(nc,rc)-np.outer(rc,nc)
    E=F/4*(np.einsum('ac,bd->abcd',gi,gi)-np.einsum('ad,bc->abcd',gi,gi))
    contraction=float(np.einsum('abcd,ab,cd',E,eps,eps))
    area=float(np.sqrt(np.linalg.det(tangents.T@g@tangents)))
    return -area*contraction,area,F,contraction


def corner_frame_check():
    rng=np.random.default_rng(648);eta4=np.diag([-1.,1.,1.,1.]);rows=[]
    for radius in (.1,.5,.85):
        phi=rng.normal(size=5);phi*=np.sqrt(12)*radius/np.linalg.norm(phi)
        T=np.eye(4)+.1*rng.normal(size=(4,4));Ti=np.linalg.inv(T)
        g=T.T@eta4@T;n=Ti[:,0];r=Ti[:,1]
        tangents=Ti[:,2:]@np.array([[1.1,.13],[0.,.9]])
        c,area,F,contraction=curvature_charge(phi,g,n,r,tangents)
        areaE=float(np.sqrt(np.linalg.det(tangents.T@(F*g)@tangents)))
        assert abs(contraction+F)<3e-14 and abs(c-areaE)<3e-14
        # An off-shell joint variation: independent boundary metric, scalar, angle.
        q=tangents.T@g@tangents
        dq=np.array([[.07,.02],[.02,-.03]]);dphi=rng.normal(size=5)*.025
        eta=.17;deta=-.03;dF=-float(phi@dphi)/3
        dc=area*dF+.5*c*np.trace(np.linalg.solve(q,dq))
        analytic=eta*dc+c*deta
        steps=[]
        for step in (1e-3,5e-4,2.5e-4):
            def action(t):
                Ft,_=boundary.data(phi+t*dphi)
                return float(Ft*np.sqrt(np.linalg.det(q+t*dq))*(eta+t*deta))
            fd=(action(step)-action(-step))/(2*step)
            steps.append(dict(step=step,derivative=fd,error=abs(fd-analytic)))
        assert steps[-1]['error']<2e-10
        rows.append(dict(radius_fraction=radius,F=F,curvature_contraction=contraction,
                         boost_charge_density=c,einstein_area_density=areaE,
                         scalar_corner_source_norm=float(eta*area*np.linalg.norm(phi)/3),
                         independent_action_variation=steps))
    return dict(rows=rows,all_five_original_scalars_retained=True,
                geometric_Wald_density_not_statistical_entropy=True)


def original_material_corner_check():
    rows=[];lam=.005
    for N in (24,32,48):
        f=original.fields(N);i=(0,N//4,N//8)
        psi=float(f['psi'][i]);phi=f['q']['phi'][i];F=float(f['q']['F'][i])
        vh=float(f['v'][i][1]);vs=float(f['v'][i][4])
        sh=float(.06*f['sstar']*psi**-2)
        A=float(f['rh'][i]);B=float(f['rs'][i]);C=-vh*vs
        K=A+2*lam*C+lam*lam*B
        assert A<0 and K<0 and C*C-A*B>0 and vh+lam*vs>0
        # Two independent normal-angle evaluations of the same original fields.
        sinh=lam*np.sqrt(C*C-A*B)/np.sqrt((-A)*(-K))
        eta=float(np.arcsinh(sinh))
        normal=np.array([vh+lam*vs,lam*sh]);normal/=np.sqrt(-K)
        eta2=float(np.arccosh(normal[0]))
        eta3=float(np.arctanh(lam*sh/(vh+lam*vs)))
        assert abs(eta-eta2)<2e-12 and abs(eta-eta3)<2e-13
        # At this original line dh has no spatial components and ds only y.
        # Hence the h,s intersection has tangent basis (partial_x,partial_z).
        areaE=psi**4;areaJ=areaE/F
        gJ=np.diag([-1.,psi**4,psi**4,psi**4])/F
        n=np.array([np.sqrt(F),0.,0.,0.]);r=np.array([0.,0.,np.sqrt(F)/psi**2,0.])
        tangents=np.eye(4)[:,[1,3]]
        charge,_,_,_=curvature_charge(phi,gJ,n,r,tangents)
        assert abs(charge-areaE)<2e-14
        dphi=np.zeros(5);dphi[1]=.07;dphi[4]=-.02
        dF=-float(phi@dphi)/3;dc=areaJ*dF
        step=1e-4
        plus=curvature_charge(phi+step*dphi,gJ,n,r,tangents)[0]
        minus=curvature_charge(phi-step*dphi,gJ,n,r,tangents)[0]
        fd=(plus-minus)/(2*step)
        assert abs(dc-fd)<1e-11 and abs(dc)>1e-3
        rows.append(dict(N=N,psi=psi,F=F,normal_h=vh,normal_s=vs,s_spatial_norm=sh,
            clock_gram=[A,C,B],second_clock_norm=K,boost_angle=eta,
            normal_angle_discrepancy=abs(eta-eta2),einstein_area_density=areaE,
            jordan_area_density=areaJ,corner_charge_from_curvature=charge,
            bare_M_area_relative_error=boundary.M/F-1,
            same_corner_action_density=eta*charge,
            fixed_Jordan_metric_scalar_charge_derivative=dc,
            independent_charge_difference_error=abs(dc-fd)))
    return dict(reference_clock_combination_lambda=lam,rows=rows,
                lambda_is_boundary_choice_not_new_dynamics=True,
                original_on_shell_background_but_variations_may_be_off_shell=True,
                local_corner_not_a_horizon_or_global_chart=True)


def integrate_interval(fn,lo,hi,order):
    if hi<=lo:return 0.
    x,w=np.polynomial.legendre.leggauss(order)
    return float((hi-lo)/2*np.sum(w*fn((lo+hi)/2+(hi-lo)/2*x)))


def regular_boost_check():
    # u=a+b is unchanged by the matching action (a,b)->(a+t,b-t).
    # Gaussian f(u) is only a normalized diagnostic, not a gravity state.
    rows=[]
    for order in (48,64):
        norm_u=integrate_interval(lambda u:np.exp(-u*u)/np.sqrt(np.pi),-9.,9.,order)
        assert abs(norm_u-1)<1e-10
        for L in (2.,8.,32.):
            t=.7
            cuts=sorted(set([-L,L,-L-t,L-t]))
            def square(a):
                left=(abs(a)<=L).astype(float)/np.sqrt(2*L)
                right=(abs(a+t)<=L).astype(float)/np.sqrt(2*L)
                return (left-right)**2
            normdiff=sum(integrate_interval(square,a,b,order) for a,b in zip(cuts[:-1],cuts[1:]))*norm_u
            inner=integrate_interval(lambda a:np.ones_like(a)/np.sqrt(2*L*8*L),-L,L,order)*norm_u
            assert abs(normdiff-t/L)<2e-10 and abs(inner-.5)<2e-10
            rows.append(dict(quadrature_order=order,L=L,unit_norm=norm_u,
                fixed_boost_defect_squared=normdiff,overlap_JL_J4L=inner,
                distance_JL_J4L_squared=2*norm_u-2*inner))
    return dict(rows=rows,representation='L2(R da) regular boost frame; explicit candidate only',
                invariant_L2_R2_subspace_zero_proved_analytically=True,
                cutoff_isometry_no_strong_limit_proved_analytically=True,
                refined_algebraic_quantization_and_other_representations_not_excluded=True)


def run():
    deps=('research_note_573.md','research_note_617.md','research_note_618.md',
          'research_note_647.md','joint_gravity_material_coordinates.py',
          'joint_gravity_material_coordinates_results.json','joint_geometric_boundary_matching.py')
    return dict(round=648,tests_run=3,failures=0,errors=0,
        original_corner_and_frame=corner_frame_check(),
        original_reference_and_corner=original_material_corner_check(),
        specified_quantum_gluing_obstruction=regular_boost_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_two_derivative_classical_corner_interface=True,
                   quantum_result_limited_to_declared_regular_boost_gluing=True,
                   no_cognition_design_or_new_physical_cutoff=True,
                   no_full_quantum_gravity_Hilbert_space_or_entropy_derivation=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
