"""866: exact finite joint instrument on the original quotient graph.
Sparse rational Fourier proof certifies the finite cubature, separately from
floating-point matrix/Dirichlet calibration on the original constrained family.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json
import numpy as np
import invariant_joint_loop_probe as probe
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_loop_instrument_results.json'
ZERO=(0,0,0)
def add(a,b,scale=Q(1)):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,Q(0))+scale*v
    return {k:v for k,v in out.items() if v}
def mul(a,b):
    out={}
    for k,v in a.items():
        for l,w in b.items():
            m=tuple(x+y for x,y in zip(k,l));out[m]=out.get(m,Q(0))+v*w
    return {k:v for k,v in out.items() if v}
def conj(a):return {tuple(-x for x in k):v for k,v in a.items()}
def mono(k):return {k:Q(1)}

def exact_certificate():
    one=mono(ZERO);u=mono((1,0,0));v=mono((0,1,0));w=mono((-1,-1,0))
    color=add(add(u,v),w);R=mul(color,mono((0,0,-2)));Rbar=conj(R)
    A=add(mul(R,Rbar),one,-1)
    delta=mul(mul(add(u,v,-1),add(u,w,-1)),add(v,w,-1))
    haar={k:value/6 for k,value in mul(delta,conj(delta)).items()}
    bases=[one,R,Rbar,A];observables=[one,R,Rbar,A,mul(R,Rbar)]
    expected=[[1,0,0,0,1],[0,0,1,0,0],[0,1,0,0,0],[0,0,0,1,1]]
    maximum=0;checked=0;table=[]
    for i,b in enumerate(bases):
        row=[]
        for j,f in enumerate(observables):
            poly=mul(haar,mul(b,f))
            maximum=max(maximum,max(abs(x) for k in poly for x in k))
            actual=poly.get(ZERO,Q(0))
            grid=sum((value for key,value in poly.items() if all(x%12==0 for x in key)),Q(0))
            assert grid==actual==expected[i][j]
            row.append(str(actual));checked+=1
        table.append(row)
    assert maximum<12 and maximum==8
    etaR=Q(1,30);etaA=Q(1,320)
    lower=1-18*etaR-64*etaA
    variance_intercept=etaR**-2-1
    variance_slope=etaA/etaR**2-1
    assert lower==Q(1,5) and variance_intercept==899 and variance_slope==Q(29,16)
    # Independent upper estimates using CF*d_R^2=12 and CA*d_A^2=192.
    color_bound=2*((2*etaR*3)**2*12+(etaA*8)**2*192)/(4*lower)
    u1_bound=(2*etaR*3)**2*36/(4*lower)
    assert color_bound==Q(3,2) and u1_bound==Q(9,5)
    return dict(exact_cubature_products=checked,max_absolute_Fourier_index=maximum,
        chosen_N=12,orthogonality_table=table,uniform_kernel_lower_bound=str(lower),
        variance_intercept=str(variance_intercept),variance_slope=str(variance_slope),
        color_Dirichlet_bound_per_simple_link=str(color_bound),
        U1_Dirichlet_bound_per_simple_link=str(u1_bound))

def run():
    prior=probe.run()
    assert prior==json.loads(probe.TARGET.read_text('utf-8'))
    exact=exact_certificate()
    indices=np.arange(12)
    a,b,c=np.meshgrid(indices,indices,indices,indexing='ij')
    w=(-a-b)%12
    mask=(a!=b)&(a!=w)&(b!=w)
    rows=[]
    with probe.ResearchRuntime(probe.Layout()).installed():
        import joint_reference_constraint_strata as old
        for L in (.4,.2):
            U,dU,_,_=probe.original.matrices(old,probe.original.inherited.constants(),L)
            rows.append(dict(coordinate_side=L,finite_instrument=probe.integrate(U,dU,old.T,12)))
    return dict(round=866,date='2026-10-06',formal_reports=866,
        cumulative_numbered_groups=3651,fresh_numbered_groups=1,all_checks_passed=True,
        exact_certificate=exact,finite_outcome_grid_size=12**3,
        nonzero_weight_outcomes=int(np.count_nonzero(mask)),rows=rows,
        input_parameters='eta_R=1/30, eta_A=1/320, finite class-function cubature',
        branchwise_Gauss_preservation=True,normal_CP_instrument_on_original_finite_graph=True,
        exact_first_and_selected_joint_moments_for_all_graph_states=True,
        bare_fusion_invariant_means_joint_estimator_not_product_of_unbiased_single_estimates=True,
        energy_scope='original fixed finite graph electric quadratic form; off-diagonal coefficients retained analytically',
        continuum_correspondence='classical material symbols and their full first physical source; not a continuum quantum instrument',
        continuum_CP_lift_proved=False,autonomous_original_action_implementation_proved=False,
        complete_graph_continuum_quantum_dynamics_match=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
