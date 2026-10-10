"""1097: finite checks for the explicitly conditional task-cone interface.
The code does not infer all-resource closure or affine motion from samples.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
ETA=np.diag([1.,-1.,-1.,-1.])


def boost(beta,direction):
    n=np.asarray(direction,float);n=n/np.linalg.norm(n)
    g=1/math.sqrt(1-beta*beta)
    a=np.eye(4);a[0,0]=g;a[0,1:]=-g*beta*n;a[1:,0]=-g*beta*n
    a[1:,1:]+=(g-1)*np.outer(n,n)
    return a


def q(v):return float(v@ETA@v)


def calculate():
    groups=[]
    def add(name,ok,**data):
        assert bool(ok),name
        groups.append(dict(name=name,status='PASS',**data))

    tasks=[dict(n=n,time=str(F(1,n)),cost_a=1,cost_b=n) for n in range(1,5)]
    a_budget=min(F(x['time']) for x in tasks if x['cost_a']<=1)
    b_budget=min(F(x['time']) for x in tasks if x['cost_b']<=1)
    add('resource_union_versus_fixed_budget',a_budget==F(1,4) and b_budget==1,
        finite_menu=tasks,all_resource_event_menus_equal=True,
        budget_one_optimum_a=str(a_budget),budget_one_optimum_b=str(b_budget),
        analytic_extension_to_all_finite_n_requires_no_uniform_cost_bound=True)

    menu=[dict(n=n,duration=1,speed=str(1-F(1,n))) for n in [2,3,5,10,100]]
    two_leg_distance=2*(1-F(1,2))
    # Every primitive has length < 1 and duration 1; any serial duration is integer.
    add('cone_envelope_is_not_the_event_set',two_leg_distance==1,
        finite_samples=menu,baseline_one_minimum_duration=2,
        actual_half_duration_task=False,finite_primitive_attains_speed_one=False,
        analytical_speed_supremum=1,finite_scale_handoff_A=1,
        inherited_1093_upper_bound=2*1*1/2,
        discrete_task_menu_is_not_a_microtime_axiom=True)

    one_way=np.eye(4);one_way[0,0]=2;one_way[1,0]=.5
    target=np.array([1.,-1.,0.,0.]);pre=np.linalg.solve(one_way,target)
    qs=q(pre)
    residual=float(np.linalg.norm(one_way.T@ETA@one_way-(one_way.T@ETA@one_way)[0,0]*ETA))
    add('forward_simulation_not_bidirectional_symmetry',qs==float(-F(21,16)) and residual>0,
        matrix=one_way.tolist(),forward_spatial_bound=1.5,forward_time=2.,
        inverse_witness=pre.tolist(),inverse_quadratic_form=qs,
        conformal_residual=residual)

    source=F(4,5);u=F(1,4);transformed=source+u
    add('finite_galilean_witness_and_background',transformed>1,
        source_event=[1,str(-source)],target_event=[1,str(-transformed)],
        source_resource_index=5,boost=str(u),declared_same_menu_bound=1,
        same_menu_symmetry_fails=True,time_order_preserved=True,
        transformed_background_velocity_ball_center=str(-u),
        transformed_background_radius=1,background_covariance_not_excluded=True)

    metric_error=0.;inverse_error=0.;minimum_future=math.inf;rows=[]
    for b in [.1,.3,.7]:
        for n in [[1,0,0],[1,2,3]]:
            for scale in [.5,1.,2.]:
                a=scale*boost(b,n)
                metric_error=max(metric_error,float(np.linalg.norm(a.T@ETA@a-scale**2*ETA)))
                inverse_error=max(inverse_error,float(np.linalg.norm(np.linalg.inv(a)@a-np.eye(4))))
                # Analytical positivity is lambda*gamma*(1-|beta|)>0 for every unit n.
                minimum_future=min(minimum_future,scale*(1-b)/math.sqrt(1-b*b))
                rows.append(dict(beta=b,scale=scale))
    dilation=2*boost(.3,[1,0,0])
    add('conditional_conformal_class_and_free_units',max(metric_error,inverse_error)<1e-11 and minimum_future>0,
        cases=len(rows),maximum_metric_error=metric_error,maximum_inverse_error=inverse_error,
        minimum_null_future_lower_bound=minimum_future,
        scaled_lorentz_metric_factor=float((dilation.T@ETA@dilation)[0,0]),
        physical_scale_eliminated=False)

    delta=.07
    bell=np.array([1,0,0,1],complex)/math.sqrt(2)
    phase=np.array([1,0,0,-1],complex)/math.sqrt(2)
    ideal=np.outer(bell,bell.conj())
    changed=(1-delta)*ideal+delta*np.outer(phase,phase.conj())
    half_trace=float(np.sum(np.abs(np.linalg.eigvalsh(changed-ideal)))/2)
    plus=np.array([1,1],complex)/math.sqrt(2);minus=np.array([1,-1],complex)/math.sqrt(2)
    out=(1-delta)*np.outer(plus,plus.conj())+delta*np.outer(minus,minus.conj())
    local_error=float(1-np.real(plus.conj()@out@plus))
    add('transport_error_does_not_preserve_exact_error_class',abs(half_trace-delta)<1e-12 and local_error>0,
        phase_flip_probability=delta,normalized_choi_half_trace_distance=half_trace,
        plus_input_error=local_error,exact_identity_menu_qualification_lost=True,
        no_claim_about_unrestricted_error_corrected_capacity=True)

    rng=np.random.default_rng(1097)
    a=1.2*boost(.3,[1,2,3])+.001*rng.normal(size=(4,4))
    b=a.T@ETA@a;aa=float(b[0,0]);m=float(np.linalg.norm(a,ord=2));r=.998
    directions=[]
    for e in np.eye(3):directions.extend([e,-e])
    for i,j in [(0,1),(0,2),(1,2)]:directions.append((np.eye(3)[i]+np.eye(3)[j])/math.sqrt(2))
    near=[];null=[]
    for n in directions:
        z=np.r_[1.,r*n];ell=np.r_[1.,n]
        near.append(abs(float(z@b@z)));null.append(abs(float(ell@b@ell)))
    eps=max(near);zeta=eps+3*m*m*(1-r)
    actual=float(np.linalg.norm(b-aa*ETA,ord='fro'))
    add('finite_nine_direction_algebra_certificate',max(null)<=zeta and actual<=9*zeta and aa>0,
        probes=len(directions),near_boundary_radius=r,matrix_norm_bound=m,
        measured_quadratic_residual_upper=eps,inferred_null_residual_upper=zeta,
        actual_null_residual=max(null),actual_conformal_residual=actual,
        certified_residual_bound=9*zeta,positive_a=aa,
        all_resource_task_coverage_certified=False,empirical_measurements=False)
    return dict(schema='round1097_complete_task_cone_bridge_v1',round=1097,status='PASS',
        numpy_version=np.__version__,groups=groups,
        scope=dict(conditional_conformal_bridge=True,adopted_R=False,
            cognitive_source_of_finite_envelope_proved=False,
            full_influence_cone_certified=False,physical_clock_scale_fixed=False,
            all_matter_motion_transport_certified=False,entire_conjecture_decided=False,
            new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(path.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1097,status='PASS',groups=len(out['groups']),saved_result_matches=True),ensure_ascii=False))


if __name__=='__main__':main()
