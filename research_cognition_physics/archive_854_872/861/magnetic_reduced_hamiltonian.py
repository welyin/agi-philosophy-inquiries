"""861: original classical clock reduction and full source pullback.
Exact magnetic metric-momentum calibration plus the original 859 constrained
background. These checks do not quantize the reduced Hamiltonian or claim a
fixed-physical-hbar graph continuum limit.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
import numpy as np
import magnetic_coordinate_canonical_probe as canonical
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'magnetic_reduced_hamiltonian_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def metric_transport():
    gen=canonical.generators();A=[F(i+2,5)*gen[i]+F(1,7)*gen[(i+1)%3] for i in range(3)]
    curvature=[[canonical.comm(A[i],A[j]) for j in range(3)] for i in range(3)]
    gram=np.array([[sum(canonical.inner(curvature[i][k],curvature[j][k]) for k in range(3)) for j in range(3)] for i in range(3)],dtype=object)
    magnetic=sum(gram[i,i] for i in range(3))/2
    pi_r=F(7,11);oldP=np.diag([F(1,3),F(1,5),-F(8,15)])
    grad=-gram;grad_tf=grad-np.eye(3,dtype=object)*sum(grad[i,i] for i in range(3))/3
    newP=oldP+F(3,4)*pi_r/magnetic*grad_tf
    directions=[]
    for diagonal in ((1,-1,0),(1,0,-1)):
        directions.append(np.diag([F(v) for v in diagonal]))
    for i,j in ((0,1),(0,2),(1,2)):
        d=canonical.zero((3,3));d[i,j]=d[j,i]=F(1);directions.append(d)
    rows=[]
    for d in directions:
        db=sum((grad*d).flat);dr=F(3,4)*db/magnetic
        theta_old=pi_r*dr+sum((oldP*d).flat);theta_new=sum((newP*d).flat)
        assert theta_old==theta_new
        # At fixed M and material labels, V=exp(r)/eta, U=eta^2 exp(-r).
        volume_log_variation=dr;counterterm_log_variation=-dr
        rows.append(dict(tracefree_direction=[[str(x) for x in row] for row in d],
            magnetic_variation=str(db),symplectic_potential_residual='0',
            volume_log_variation=str(volume_log_variation),potential_log_variation=str(counterterm_log_variation)))
    assert any(F(r['volume_log_variation']) for r in rows)
    return dict(all_five_unimodular_directions_tested=True,rows=rows)

def exact_source():
    h,s=F(1),F(1,2);vol=F(3,2);targetF=2-(h*h+s*s)/6
    invhh=targetF*(1-h*h/12);invhs=-targetF*h*s/12;invss=targetF*(1-s*s/12)
    a=invhh/(2*vol);P_s=F(-2,5);pi_h=F(3,7);b=invhs*P_s/vol
    speed=2*a*pi_h+b;assert speed>0
    c=-a*pi_h*pi_h-b*pi_h
    # Declared compatible-family differential after solving the spatial
    # constraints: P_s(j)=P_s+j*P_s_prime. This is a chain-rule check, not
    # a new covariant source model or a solution of the Einstein PDE.
    P_s_prime=F(2,9);U=F(5,13)
    b_prime=invhs*P_s_prime/vol
    c_prime=U+invss*P_s*P_s_prime/vol
    correct=(b_prime*pi_h+c_prime)/speed
    naive=U/speed
    assert correct!=naive
    pi_h_prime=-correct
    differentiated_constraint=speed*pi_h_prime+b_prime*pi_h+c_prime
    assert differentiated_constraint==0
    # When the spatial constraint is unaffected, H''(0)=2a U^2/v_h^3.
    second=2*a*U*U/speed**3;assert second>0
    return dict(a=str(a),b=str(b),clock_speed=str(speed),
        compatible_source_chain_derivative=str(correct),frozen_spatial_momentum_derivative=str(naive),
        omitted_source_reaction=str(correct-naive),differentiated_constraint_residual='0',
        potential_second_derivative=str(second),calibration_not_a_covariant_source_construction=True)

def original_background():
    rows=[]
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        geo=old.geo
        for N in (16,24):
            q,psi0,tensor0,_,color=old.completed(N,1.)
            x=q['grid'][...,0];eps=.02;probe=eps*np.sin(x)
            new=dict(q);new['B']=q['B']+eps**2*np.cos(x)**2
            new['U']=q['U']+.5*probe**2;new['C']=q['C']-probe**2
            psi,tensor,stats=geo.solve_hamiltonian(new,initial=psi0)
            idx=(0,N//4,N//8);vol=float(psi[idx]**6)
            phi=q['phi'][idx];momentum=q['p'][idx]
            h=float(phi[1]);s=float(phi[4]);ph=float(momentum[1]);ps=float(momentum[4])
            F0=float(q['F'][idx]);kinv=F0*(np.eye(5)-np.outer(phi,phi)/12)
            assert np.max(abs(phi[[0,2,3]]))<1e-14
            a=kinv[1,1]/(2*vol);b=kinv[1,4]*ps/vol
            vh=float((kinv@momentum)[1]/vol);assert vh>0
            assert abs(vh-(2*a*ph+b))<1e-15
            scalar_curvature=-8*psi**-5*geo.laplace(psi)
            extrinsic=psi**-12*np.sum(tensor*tensor,axis=(-1,-2))-2*q['tau2']/3
            Cgrav=.5*vol*float((-scalar_curvature+extrinsic)[idx])
            Ckin=.5*float(momentum@kinv@momentum)/vol
            Cother=vol*float((.5*psi**-4*new['B']+new['U']+psi**-8*q['Y'])[idx])
            Ctotal=Cgrav+Ckin+Cother
            c_full=Ctotal-a*ph**2-b*ph
            discriminant=b*b-4*a*c_full
            assert discriminant>0 and abs(Ctotal)<2e-10
            reconstructed=(-b+np.sqrt(discriminant))/(2*a)
            root_error=abs(reconstructed-ph);assert root_error<2e-7
            # Exact on-shell identity is used for the source diagnostic;
            # actual independently evaluated constraint residual is kept above.
            c=-a*ph**2-b*ph;U=1/vol
            source_radius=vh**2/(4*a*U)
            def shift(t):
                # Stable difference H(lambda)-H(0), lambda=t*source_radius.
                return 2*t*source_radius*U/(vh+np.sqrt(vh**2-4*a*t*source_radius*U))
            e=1e-4
            central=(shift(e)-shift(-e))/(2*e*source_radius)
            predicted=U/vh
            err=abs(central/predicted-1);assert err<2e-8
            branch=[]
            for t in (-.5,.1,.5):
                lam=t*source_radius;Hshift=shift(t);newph=ph-Hshift
                resid=a*newph**2+b*newph+c+lam*U
                assert abs(resid)<1e-16
                new_speed=2*a*newph+b
                assert abs(new_speed/vh-np.sqrt(1-t))<1e-12
                branch.append(dict(source_fraction_of_positive_branch_endpoint=t,
                    constraint_residual=float(resid),clock_speed_ratio=float(new_speed/vh),physical_H_shift=float(Hshift)))
            rows.append(dict(grid_N=N,probe_amplitude=eps,volume_density=vol,
                original_constraint_density_residual=float(Ctotal),a=float(a),b=float(b),
                reconstructed_clock_momentum_error=float(root_error),actual_clock_speed=vh,
                positive_discriminant=float(discriminant),physical_lapse=float(1/vh),
                unit_eta_potential=U,source_derivative=predicted,finite_difference_relative_error=float(err),
                unit_eta_positive_branch_endpoint=float(source_radius),source_branch=branch,
                numerical_point_is_h_normal_at_point_only=True,full_material_chart_grid_not_solved=True))
    return rows

def run():
    assert canonical.run()==json.loads(canonical.TARGET.read_text('utf-8'))
    scales=[dict(label_cell_size=f'1/{n}',inverse_volume_potential_multiplier=n**6,
                 analytic_coupling_radius_multiplier=f'1/{n**6}') for n in (1,2,4,8)]
    return dict(round=861,fresh_test_groups=1,all_checks_passed=True,
        canonical_color_transport_reproduced=True,metric_transport=metric_transport(),
        full_source_chain=exact_source(),original_859_background_clock_reduction=original_background(),
        fixed_geometry_material_refinement_scaling=scales,
        analytic_result='A local classical reduction using h,s,p and the original color magnetic scalar yields an explicit physical Hamiltonian and source/volume pullbacks. It uses the original Einstein constraints and the full canonical momentum transformation.',
        exact_high_derivative_EFT_Hamiltonian_claimed=False,
        arbitrary_added_potential_preserves_full_off_gauge_constraint_algebra_claimed=False,
        original_graph_quantization_or_continuum_equivalence_proved=False,
        residual_579_divergence_cancelled=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
