"""864: finite-order source calibration for a local equicausal loop germ.
Uniform infinite-dimensional estimates are analytic in research_note_864.md.
This code reproduces the old cone test and removes finite differencing from
the original-color second Ward calibration by using two-variable matrix jets.
"""
from pathlib import Path
import argparse,json,math,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'equicausal_relational_loop_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
import higher_loop_source_probe as working
INDICES=((0,0),(1,0),(0,1),(1,1))

def zero():return {i:np.zeros((3,3),complex) for i in INDICES}
def identity():
    r=zero();r[(0,0)]=np.eye(3,dtype=complex);return r
def mul(A,B):
    out=zero()
    for i,a in A.items():
        for j,b in B.items():
            k=(i[0]+j[0],i[1]+j[1])
            if k in out:out[k]+=a@b
    return out
def exponential(A,order=36):
    out=identity();term=identity()
    for n in range(1,order+1):
        term={i:a/n for i,a in mul(term,A).items()}
        out={i:out[i]+term[i] for i in INDICES}
    rho=sum(float(np.linalg.norm(a,2)) for a in A.values())
    tail=math.exp(rho)*rho**(order+1)/math.factorial(order+1)
    return out,tail

def color_two_jet():
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        A=-1j*old.color(1.)['A']
        T=-1j*(.3*old.T[1]+.7*old.T[3]+.2*old.T[6])
        B=-1j*np.array([.23*old.T[5]+.07*old.T[2],.13*old.T[4],.17*old.T[0]+.11*old.T[6]])
    K=np.array([T@a-a@T for a in A]);KB=np.array([T@b-b@T for b in B])
    rows=[]
    for include_change in (False,True):
        fields={ (0,0):A,(1,0):B,(0,1):K,(1,1):KB if include_change else np.zeros_like(KB)}
        U=identity();bounds=[]
        for sign,axis in ((-1,'v'),(-1,'z'),(1,'v'),(1,'z')):
            block={i:sign*.8*(a[0]-a[1] if axis=='v' else a[2]) for i,a in fields.items()}
            factor,bound=exponential(block);bounds.append(bound);U=mul(factor,U)
        rows.append(dict(generator_change_included=include_change,
            background_trace=float(np.trace(U[(0,0)]).real),
            first_gauge_coefficient=float(np.trace(U[(0,1)]).real),
            mixed_coefficient=float(np.trace(U[(1,1)]).real),
            largest_single_exponential_series_tail_bound=max(bounds)))
    source=working.derivative(A,KB)
    assert abs(rows[0]['mixed_coefficient'])>1e-6
    assert abs(rows[1]['mixed_coefficient'])<2e-17
    assert abs(rows[0]['mixed_coefficient']+source)<2e-17
    assert max(r['largest_single_exponential_series_tail_bound'] for r in rows)<1e-40
    # Ordinary finite-difference check is retained as an independent method.
    old=working.color_higher_ward()
    assert abs(old['rows'][-1]['hessian_with_frozen_generator']-rows[0]['mixed_coefficient'])<2e-11
    return dict(rows=rows,generator_change_source=source,
        frozen_generator_hessian_plus_source_residual=abs(rows[0]['mixed_coefficient']+source),
        finite_difference_crosscheck_error=abs(old['rows'][-1]['hessian_with_frozen_generator']-rows[0]['mixed_coefficient']),
        calibration_scope='Original 753 color factor and fixed path, diagnostic B; not a full gravitational relational source calculation.',
        analytic_series_tail_bound_excludes_floating_point_roundoff=True)

def run():
    prior=working.run()
    assert prior==json.loads(working.TARGET.read_text('utf-8'))
    return dict(round=864,date='2026-10-06',formal_reports=864,fresh_numbered_groups=1,
        cumulative_numbered_groups=3649,all_checks_passed=True,
        common_anchor_cone_margin=prior['common_anchor_sufficient_cone_bound'],
        original_color_second_source_jet=color_two_jet(),
        same_relational_loop_family_as_863=True,
        uniform_equicausal_germ='analytic finite pushforward and fixed-cone seminorm proof',
        same_fixed_mixed_Hadamard_star_coefficients='well-defined and smooth at each finite order',
        unchanged_original_physical_linear_variance_from_863=True,
        full_nonlinear_quantum_BV_Ward_completed=False,
        time_ordered_extension_with_loop_insertions_completed=False,
        finite_hbar_convergence_or_bounded_operator_proved=False,
        original_graph_quantum_matching_proved=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
