"""864 working: common-anchor cone margin and nonlinear color Ward source.
Cone box is a mathematical sufficient-condition calibration, not the original
spacetime reference Jacobian. Matrix Ward probe uses the original color field;
it is not a computation of the full relational gravitational current.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'higher_loop_source_probe_results.json'
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(HERE.parent/'863'))
from research_layout import Layout,ResearchRuntime
from physical_relational_loop_bridge import exp_with_tangent


def anchor_bound():
    delta=Q(1,16);theta=Q(1,8);pairings=[]
    # Future cube |k_i| <= k_0 contains the local Minkowski future cone.
    # V lies in a componentwise delta-box around (1,0,0,0).
    for signs_k in product((-1,1),repeat=3):
        k=[Q(1)]+list(map(Q,signs_k))
        for signs_V in product((-1,1),repeat=4):
            V=[Q(1)+delta*signs_V[0]]+[delta*s for s in signs_V[1:]]
            pairings.append(sum(a*b for a,b in zip(k,V)))
    assert min(pairings)==1-4*delta==Q(3,4)
    assert Q(1,4)>theta>delta
    return dict(componentwise_velocity_deviation=str(delta),closed_cone_threshold=str(theta),
        future_covector_vertex_minimum_pairing=str(min(pairings)),
        future_cone_l1_margin=str(Q(1,4)-theta),outside_closed_cone_phase_margin=str(theta-delta),
        corner_cases=len(pairings),
        extension_to_any_number_of_insertions='sum of the same inequality, including zero legs',
        original_reference_chart_numerically_reconstructed=False)


def derivative(A,B,L=.8):
    V=A[0]-A[1];Z=A[2];dV=B[0]-B[1];dZ=B[2]
    U=np.eye(3,dtype=complex);dU=np.zeros_like(U)
    for a,b in ((-L*V,-L*dV),(-L*Z,-L*dZ),(L*V,L*dV),(L*Z,L*dZ)):
        C,dC=exp_with_tangent(a,b);dU=dC@U+C@dU;U=C@U
    return float(np.trace(dU).real)


def color_higher_ward():
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        A=-1j*old.color(1.)['A'];T=-1j*(.3*old.T[1]+.7*old.T[3]+.2*old.T[6])
        B=-1j*np.array([.23*old.T[5]+.07*old.T[2],.13*old.T[4],.17*old.T[0]+.11*old.T[6]])
    K=np.array([T@a-a@T for a in A]);KB=np.array([T@b-b@T for b in B])
    first=derivative(A,K);change=derivative(A,KB)
    assert abs(first)<1e-14 and abs(change)>1e-6
    rows=[]
    for step in (1e-3,5e-4,1e-4):
        mixed=(derivative(A+step*B,K)-derivative(A-step*B,K))/(2*step)
        rows.append(dict(step=step,hessian_with_frozen_generator=mixed,
            generator_change_source=change,full_second_Ward_residual=mixed+change,
            omitted_source_defect=abs(mixed)))
    errors=[abs(r['full_second_Ward_residual']) for r in rows]
    assert all(a>b for a,b in zip(errors,errors[1:])) and errors[-1]<2e-11
    return dict(first_Ward_residual=first,rows=rows,
        scope='Original color factor with a fixed path; B is a diagnostic variation, not a proved full constrained relational solution. Full second Ward identity is derived analytically in the draft.')


def run():
    return dict(kind='round_864_working_higher_loop_sources',formal_reports=863,new_numbered_scientific_groups=0,
        all_checks_passed=True,common_anchor_sufficient_cone_bound=anchor_bound(),
        original_color_nonlinear_Ward=color_higher_ward(),
        full_relational_equicausal_membership_published=False,
        interacting_quantum_Ward_or_BV_completion_proved=False,
        original_graph_quantum_matching_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
