"""782: full-domain direct normalization diagnostics.

These finite sign, source, and linear-kernel checks are not computations of
the original continuum anomaly. No nonlinear quantum equivalence is assumed.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('bv777',HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv=importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
add,mul,scale,diff=bv.add,bv.mul,bv.scale,bv.diff


def full_nonlinear_input():
    v={name:bv.var(name) for name in bv.NAMES}
    q,r,b,c,h,p,t,a=[v[x] for x in ('q','r','b','c','h','p','t','a')]
    alpha,coupling=Q(2,5),Q(3,5)
    invariant=add(q,scale(mul(r,r),alpha))
    frame_anti=add(t,scale(mul(r,p),-2*alpha))
    dictionary={'q':invariant,'t':frame_anti}
    base=add(scale(mul(q,q),Q(1,2)),scale(bv.prod(q,q,q),coupling/6),
             mul(t,c),mul(b,r),scale(mul(a,b),-1),scale(mul(h,c),-1))
    full=bv.substitute(base,dictionary)
    free={key:value for key,value in full.items() if sum(key)==2}
    interaction=add(full,scale(free,-1))
    assert not bv.bracket(full,full)
    assert not bv.bracket(free,free)
    assert not bv.bracket(full,invariant)
    canonical=0
    for x in bv.NAMES:
        for y in bv.NAMES:
            assert bv.bracket(bv.substitute(v[x],dictionary),bv.substitute(v[y],dictionary))==bv.bracket(v[x],v[y])
            canonical+=1
    s0=lambda f:bv.bracket(free,f)
    contraction=lambda f:add(scale(diff(diff(f,'q'),'q'),Q(1,2)),
                              diff(diff(f,'r'),'b'),scale(diff(diff(f,'c'),'h'),-1))
    laplacian=lambda f:add(*(scale(diff(diff(f,anti),field),-1 if field in bv.ODD else 1)
                            for field,anti in bv.PAIRS))
    checked=0
    for degree in range(4):
        for names in combinations_with_replacement(bv.NAMES,degree):
            f=bv.prod(*(v[name] for name in names))
            if not f: continue
            assert not s0(s0(f))
            assert not add(s0(contraction(f)),scale(contraction(s0(f)),-1),laplacian(f))
            checked+=1
    cubic={key:value for key,value in interaction.items() if sum(key)==3}
    assert not s0(cubic)
    # Same smooth physical Wick contraction, no auxiliary state covariance.
    smooth=Q(2,3)
    source=scale(diff(diff(cubic,'q'),'q'),smooth/2)
    assert source==scale(q,coupling*smooth/2)
    base_cubic={key:value for key,value in base.items() if sum(key)==3}
    assert source==scale(diff(diff(base_cubic,'q'),'q'),smooth/2)
    # Nonlinear vertices stay as inputs; replacing the whole action by Q fails.
    assert interaction and len(interaction)>1
    # Dropping the transformed antifield vertex destroys the CME.
    wrong=add(full,scale(bv.prod(r,p,c),2*alpha))
    defect=bv.bracket(wrong,wrong)
    assert defect
    # A second local chart has a nonzero pure-frame tadpole. Work to degree 6
    # in an analytic chart; do not assert an exact CME for a truncated series.
    order=6
    inverse_factor={}
    for power in range(order):
        inverse_factor=add(inverse_factor,scale(bv.prod(*([r]*power)),(-alpha)**power))
    scaled=mul(q,add(bv.ONE,scale(r,alpha)))
    scaled_t=add(t,scale(bv.prod(q,p,inverse_factor),-alpha))
    scaled_action=add(scale(mul(scaled,scaled),Q(1,2)),
                      scale(bv.prod(scaled,scaled,scaled),coupling/6),
                      mul(scaled_t,c),mul(b,r),scale(mul(a,b),-1),scale(mul(h,c),-1))
    scaled_action={key:value for key,value in scaled_action.items() if sum(key)<=order}
    master=bv.bracket(scaled_action,scaled_action)
    assert not {key:value for key,value in master.items() if sum(key)<=order}
    scaled_cubic={key:value for key,value in scaled_action.items() if sum(key)==3}
    scaled_source=scale(diff(diff(scaled_cubic,'q'),'q'),smooth/2)
    frame_source=scale(r,alpha*smooth)
    assert scaled_source==add(source,frame_source)
    first_frame_counterterm=scale(frame_source,-1)
    assert not add(s0(scaled_source),s0(first_frame_counterterm))
    return dict(canonical_generator_pairs=canonical,free_bv_monomials=checked,
                full_nonlinear_master_equation_exact=True,
                nonlinear_vertices_retained=bv.display(interaction),
                same_cubic_smooth_source=bv.display(source),
                analytic_frame_chart_cme_checked_through_degree=order,
                second_chart_frame_source=bv.display(frame_source),
                first_frame_counterterm=bv.display(first_frame_counterterm),
                frame_source_repaired_without_changing_physical_source=True,
                omitted_antifield_interaction_master_defect=bv.display(defect),
                same_free_contraction_suffices_for_complete_input_domain=True)


def linear_ghost_kernel_transport():
    n=3; eye=np.eye(n); zero=np.zeros((n,n))
    k=np.array([[3.,-1,0],[-1,4,-1],[0,-1,3]])
    derivative=np.array([[-2.,2.,0.],[0.,-3.,3.],[1.,0.,-1.]])
    p=np.block([[zero,k,zero,zero],[-k.T,zero,zero,zero],
                [zero,zero,zero,eye],[zero,zero,-eye,zero]])
    # Old ordering (c,bar c,lambda,bar lambda), adapted (c,bar c,theta,bar theta).
    m=np.eye(4*n); m[2*n:3*n,:n]=-derivative
    inverse=np.linalg.inv(m)
    old_p=inverse.T@p@inverse
    g=np.linalg.inv(p)
    old_g=m@g@m.T
    u=np.array([[.8,.1,0],[.2,.7,.1],[0,.3,.9]])
    w=np.block([[zero,u,zero,zero],[-u.T,zero,zero,zero],
                [zero,zero,zero,zero],[zero,zero,zero,zero]])
    old_w=m@w@m.T
    assert np.array_equal(old_w[:2*n,:2*n],w[:2*n,:2*n])
    errors=dict(green_identity=float(np.linalg.norm(old_p@old_g-np.eye(4*n))),
                covariance_round_trip=float(np.linalg.norm(inverse@old_w@inverse.T-w)),
                cotangent_pairing=float(np.linalg.norm(m.T@inverse.T-np.eye(4*n))))
    wrong=old_g.copy()
    wrong[2*n:3*n,n:2*n]=0
    wrong[n:2*n,2*n:3*n]=0
    defect=float(np.linalg.norm(old_p@wrong-np.eye(4*n)))
    assert max(errors.values())<1e-12 and defect>1
    assert np.linalg.norm(old_w[2*n:3*n,n:2*n])>1
    return dict(sites=n,maximum_residuals=errors,
                original_propagating_covariance_unchanged=True,
                original_ghost_mixed_block_required=True,
                omit_mixed_green_block_defect=defect,
                fixed_background_difference_operator=derivative.tolist(),
                finite_w_matrix_claimed_as_continuum_bisolution=False)


def nonlinear_star_pullback_boundary():
    # Two free linear variables with [q,r]=i*hbar. The ordinary symbol
    # Poisson coefficient exactly fixes the commutator in this quadratic-linear case.
    q,r=bv.var('q'),bv.var('r')
    kappa=Q(2,7)
    f=add(q,scale(mul(q,q),kappa))
    poisson=add(mul(diff(f,'q'),diff(r,'r')),scale(mul(diff(f,'r'),diff(r,'q')),-1))
    defect=add(poisson,scale(bv.ONE,-1))
    assert defect==scale(q,2*kappa)
    # A different, source-map error occurs even for an ordinary Gaussian.
    # <f(q)> - f(<q>) = kappa*hbar*C, with C=3/5 and centered q.
    source_shift=kappa*Q(3,5)
    assert source_shift==Q(6,35)
    return dict(exact_commutator_over_i_hbar=bv.display(poisson),
                naive_same_free_star_pullback_defect=bv.display(defect),
                quantum_correction_starting_at_hbar_cannot_change_leading_poisson=True,
                centered_gaussian_composite_source_shift_over_hbar=str(source_shift),
                counterexample_scope='Naive unchanged free-star pullback only; not a no-go for interacting equivalence or full-domain direct normalization.')


def run():
    return dict(round=782,groups={
        'full_nonlinear_bv_input_in_fixed_free_family':full_nonlinear_input(),
        'background_linear_full_kernel_dictionary':linear_ghost_kernel_transport(),
        'naive_nonlinear_star_transport_boundary':nonlinear_star_pullback_boundary()},
        all_checks_passed=True,
        analytical_claim='Original N1/N2 existence in the declared full local field chart; stronger independent-chart equivalence is not inferred.',
        finite_checks_prove_continuum_theorem=False,
        original_anomaly_coefficients_computed=False,
        all_loop_quantum_master_equation_proven=False,
        interacting_positive_state_proven=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true'); parser.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=run(); destination=HERE/'direct_full_normalization_results.json'
    if args.write: destination.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if args.check: assert result==json.loads(destination.read_text(encoding='utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
