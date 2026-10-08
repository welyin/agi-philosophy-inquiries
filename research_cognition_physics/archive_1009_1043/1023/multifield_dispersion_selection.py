"""1023: a two-species crossing cone and its controlled-error enlargement.

This calibrates an adopted improved forward-amplitude contract, not positivity
of bare Wilson coefficients or existence of a UV theory/physical instrument.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "multifield_dispersion_selection_results.json"
ORDER = [(0,0),(1,1),(0,1),(1,0)]
HISTORICAL = (
    "archive_1009_/1022/research_round_1022_checks.json",
    "archive_1009_/1022/NEXT.md",
    "archive_935_955/research_note_952.md",
    "archive_990_1008/research_note_993.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
)


def norm(a):
    return float(np.linalg.norm(a))


def fractions(values):
    return tuple(F(x) for x in values)


def cone_exact(values, epsilon=F(0)):
    """No square-root tolerance at the cone's equality or zero-coefficient faces."""
    a,b,c,d=fractions(values); a+=epsilon; b+=epsilon; c+=epsilon
    if min(a,b,c)<0:
        return False
    excess=abs(d)-c
    return excess<=0 or excess*excess<=a*b


def matrix(values):
    a,b,c,d=map(float,values)
    return np.array([[a,d/2,0,0],[d/2,b,0,0],[0,0,c,d/2],[0,0,d/2,c]])


def tensor_from_coefficients(values):
    """Direct four-index assignment, independently of the block-matrix helper."""
    a,b,c,d=map(float,values); T=np.zeros((2,2,2,2))
    T[0,0,0,0]=a; T[1,1,1,1]=b
    T[0,1,0,1]=T[1,0,1,0]=c
    for index in ((0,0,1,1),(0,1,1,0),(1,0,0,1),(1,1,0,0)):
        T[index]=d/2
    return T


def flatten(T):
    return np.array([[T[i,j,k,l] for k,l in ORDER] for i,j in ORDER])


def pair_vector(u,v):
    return np.array([u[i]*v[j] for i,j in ORDER])


def q_tensor(T,u,v):
    result=np.einsum("i,j,ijkl,k,l->",u,v,T,np.conj(u),np.conj(v))
    assert abs(result.imag)<2e-10
    return float(result.real)


def q_polynomial(values,u,v):
    a,b,c,d=map(float,values)
    return float(a*u[0]**2*v[0]**2+b*u[1]**2*v[1]**2
        +c*(u[0]**2*v[1]**2+u[1]**2*v[0]**2)+2*d*u[0]*u[1]*v[0]*v[1])


def generator(m):
    return np.einsum("ij,kl->ijkl",m,m)+np.einsum("il,kj->ijkl",m,m)


def cone_generators(values):
    assert cone_exact(values)
    a,b,c,d=map(float,values); s=math.sqrt(a*b); scale=s+c
    Dplus=np.diag([math.sqrt(a/2),math.sqrt(b/2)])
    Dminus=np.diag([math.sqrt(a/2),-math.sqrt(b/2)])
    Oplus=math.sqrt(c/2)*np.array([[0.,1.],[1.,0.]])
    Ominus=math.sqrt(c/2)*np.array([[0.,1.],[-1.,0.]])
    if scale==0:
        assert d==0
        return [Dplus]
    alpha=(1+d/scale)/2
    assert -1e-14<=alpha<=1+1e-14
    alpha=max(0.,min(1.,alpha))
    return [math.sqrt(alpha)*Dplus,math.sqrt(1-alpha)*Dminus,
            math.sqrt(alpha)*Oplus,math.sqrt(1-alpha)*Ominus]


def violation_witness(values):
    """A finite product witness for every failing case used below."""
    a,b,c,d=fractions(values)
    if a<0:
        return np.array([1.,0.]),np.array([1.,0.]),a,None
    if b<0:
        return np.array([0.,1.]),np.array([0.,1.]),b,None
    if c<0:
        return np.array([1.,0.]),np.array([0.,1.]),c,None
    k=abs(d)-c
    assert k>0 and k*k>a*b
    x=(a+k)/(b+k); gap=(k*k-a*b)/(a+b+2*k)
    sign=1 if d>0 else -1
    u=np.array([1.,math.sqrt(float(x))])/math.sqrt(float(1+x))
    v=np.array([1.,-sign*math.sqrt(float(x))])/math.sqrt(float(1+x))
    exact=(a+b*x*x-2*k*x)/(1+x)**2
    assert exact==-gap<0
    return u,v,exact,x


def direction_calibrations(values,seed):
    T=tensor_from_coefficients(values); M=matrix(values)
    rng=np.random.default_rng(seed)
    real_pairs=[(np.array([1.,0.]),np.array([1.,0.])),
                (np.array([0.,1.]),np.array([0.,1.])),
                (np.array([1.,0.]),np.array([0.,1.])),
                (np.ones(2)/math.sqrt(2),np.array([1.,-1.])/math.sqrt(2))]
    for _ in range(6):
        u=rng.normal(size=2); v=rng.normal(size=2)
        real_pairs.append((u/np.linalg.norm(u),v/np.linalg.norm(v)))
    complex_pairs=[(np.array([1.,1j])/math.sqrt(2),np.array([1.,1.])/math.sqrt(2)),
                   (np.array([1.,1j])/math.sqrt(2),np.array([1.,-1j])/math.sqrt(2)),
                   (np.array([1.,2j])/math.sqrt(5),np.array([2j,1.])/math.sqrt(5)),
                   (np.array([1j,0.]),np.array([0.,1.]))]
    for _ in range(4):
        u=rng.normal(size=2)+1j*rng.normal(size=2)
        v=rng.normal(size=2)+1j*rng.normal(size=2)
        complex_pairs.append((u/np.linalg.norm(u),v/np.linalg.norm(v)))
    errors=dict(tensor_matrix=norm(flatten(T)-M),
                crossing_i_k=norm(T-T.transpose(2,1,0,3)),
                crossing_j_l=norm(T-T.transpose(0,3,2,1)),
                simultaneous_pair_reversal=norm(T-T.transpose(1,0,3,2)),
                pair_matrix_hermiticity=norm(M-M.T),
                direct_polynomial=0.,product_matrix=0.,complex_real_decomposition=0.)
    real_values=[]; complex_values=[]
    for u,v in real_pairs:
        q=q_tensor(T,u,v); w=pair_vector(u,v)
        errors["direct_polynomial"]=max(errors["direct_polynomial"],abs(q-q_polynomial(values,u,v)))
        errors["product_matrix"]=max(errors["product_matrix"],abs(q-float(w@M@w)))
        real_values.append(q)
    for u,v in complex_pairs:
        q=q_tensor(T,u,v); w=pair_vector(u,v)
        decomposed=sum(q_tensor(T,x,y) for x in (u.real,u.imag) for y in (v.real,v.imag))
        errors["complex_real_decomposition"]=max(errors["complex_real_decomposition"],abs(q-decomposed))
        errors["product_matrix"]=max(errors["product_matrix"],abs(q-np.vdot(w,M@w)))
        complex_values.append(q)
    if cone_exact(values):
        assert min(real_values+complex_values)>=-2e-12
    assert max(errors.values())<3e-11
    return dict(real_direction_count=len(real_pairs),complex_direction_count=len(complex_pairs),
                smallest_sampled_real_product=float(min(real_values)),
                smallest_sampled_complex_product=float(min(complex_values)),errors=errors,
                sampled_directions_do_not_prove_cone=True)


def cone_cases():
    cases=[("interior_PSD",(1,2,1,1)),("non_PSD_legal",(1,1,4,3)),
           ("positive_surface",(1,1,1,2)),("negative_surface",(4,9,2,-8)),
           ("zero_c_surface",(1,4,0,2)),("zero_a_surface",(0,1,2,2)),
           ("zero_b_surface",(1,0,2,-2)),("zero_a_b_surface",(0,0,1,1)),
           ("zero_origin",(0,0,0,0)),("a_axis",(1,0,0,0)),
           ("b_axis",(0,3,0,0)),("c_axis",(0,0,2,0)),
           ("nonsquare_ab_interior",(2,3,1,3)),
           ("fixed_flavors_miss_violation",(1,1,1,3)),
           ("negative_a",(-1,2,1,0)),("negative_b",(2,-1,1,0)),
           ("negative_c",(1,2,-1,0)),("zero_a_violation",(0,1,1,2)),
           ("zero_b_violation",(1,0,1,-2)),("zero_a_b_violation",(0,0,1,2))]
    records=[]
    for index,(name,raw) in enumerate(cases):
        values=fractions(raw); allowed=cone_exact(values); T=tensor_from_coefficients(values)
        record=dict(name=name,parameters=[str(x) for x in values],in_exact_cone=allowed,
                    pair_matrix_eigenvalues=np.linalg.eigvalsh(matrix(values)).tolist(),
                    direction_calibration=direction_calibrations(values,20261023+index))
        if allowed:
            generators=cone_generators(values)
            generated=sum((generator(m) for m in generators),np.zeros((2,2,2,2)))
            error=norm(generated-T); assert error<3e-12
            record.update(generator_matrices=[m.tolist() for m in generators],
                          generator_tensor_reconstruction_error=error)
        else:
            u,v,q,x=violation_witness(values)
            actual=q_tensor(T,u,v)
            assert math.isclose(actual,float(q),rel_tol=2e-12,abs_tol=2e-12)
            assert q<0
            record.update(witness_u=u.tolist(),witness_v=v.tolist(),
                          witness_value_exact=str(q),witness_value=actual,
                          witness_squared_component_ratio=None if x is None else str(x))
        records.append(record)
    return records


def non_PSD_generator_witness():
    generators=[np.diag([1.,1.])/math.sqrt(2),
                math.sqrt(1.5)*np.array([[0.,1.],[1.,0.]]),
                np.array([[0.,1.],[-1.,0.]])/math.sqrt(2)]
    T=sum((generator(m) for m in generators),np.zeros((2,2,2,2)))
    target=tensor_from_coefficients((1,1,4,3))
    error=norm(T-target); assert error<3e-12
    z=np.array([1.,-1.,0.,0.])/math.sqrt(2)
    assert np.linalg.matrix_rank(np.array([[z[0],z[2]],[z[3],z[1]]]))==2
    entangled_value=float(z@flatten(T)@z)
    assert math.isclose(entangled_value,-.5,abs_tol=3e-12)
    u=np.array([2.,-1.])/math.sqrt(5); v=np.array([1.,3.])/math.sqrt(10)
    squares=[2*float(u@m@v)**2 for m in generators]
    assert abs(sum(squares)-q_tensor(T,u,v))<3e-12
    return dict(parameters=[1,1,4,3],generator_matrices=[m.tolist() for m in generators],
                reconstruction_error=error,minimum_matrix_eigenvalue=float(np.linalg.eigvalsh(flatten(T))[0]),
                entangled_vector=z.tolist(),entangled_coefficient_value=entangled_value,
                entangled_vector_Schmidt_rank=2,positive_real_product_squares=squares,
                generator_decomposition_proves_product_positivity=True,
                coefficient_tensor_is_not_a_Choi_operator=True,
                positive_generators_do_not_prove_UV_completion=True)


def minimal_repair_cases():
    inputs=[(1,1,1,3),(4,9,2,10),(0,1,1,2),(1,0,1,-2),
            (0,0,1,2),(F(1,3),F(2,5),F(1,7),2)]
    records=[]
    for raw in inputs:
        a,b,c,d=values=fractions(raw); k=abs(d)-c
        assert min(a,b,c)>=0 and k>0 and k*k>a*b
        gap=(k*k-a*b)/(a+b+2*k)
        u,v,q,x=violation_witness(values)
        assert q==-gap
        assert not cone_exact(values,gap/2)
        assert cone_exact(values,gap) and cone_exact(values,3*gap/2)
        shifted=(a+gap,b+gap,c+gap,d)
        assert (abs(d)-(c+gap))**2==(a+gap)*(b+gap)
        M=matrix(values); C=matrix(shifted)
        shift_error=norm(C-M-float(gap)*np.eye(4))
        repair_norm=float(np.linalg.norm(C-M,ord=2))
        boundary_value=q_tensor(tensor_from_coefficients(shifted),u,v)
        assert shift_error<3e-12 and abs(boundary_value)<3e-12
        assert math.isclose(repair_norm,float(gap),abs_tol=3e-12)
        generators=cone_generators(shifted)
        generator_error=norm(sum((generator(m) for m in generators),np.zeros((2,2,2,2)))-tensor_from_coefficients(shifted))
        assert generator_error<3e-12
        records.append(dict(parameters=[str(z) for z in values],k_exact=str(k),
            squared_component_ratio_exact=str(x),minimal_repair_exact=str(gap),
            witness_value_exact=str(q),witness_value=q_tensor(tensor_from_coefficients(values),u,v),
            witness_u=u.tolist(),witness_v=v.tolist(),shifted_parameters=[str(z) for z in shifted],
            half_gap_repair_allowed=False,exact_gap_repair_allowed=True,one_and_half_gap_repair_allowed=True,
            shift_identity_error=shift_error,operator_norm_repair=repair_norm,
            shifted_witness_value=boundary_value,shifted_generator_reconstruction_error=generator_error))
    general=[]
    for raw,eps,expected in [((-1,2,1,0),F(1),True),((-1,2,1,0),F(1,2),False),
                              ((0,1,-1,0),F(1),True),((0,0,0,1),F(1,2),True),
                              ((0,0,0,1),F(1,3),False)]:
        allowed=cone_exact(raw,eps); assert allowed==expected
        general.append(dict(parameters=[str(z) for z in raw],epsilon_exact=str(eps),in_enlarged_cone=allowed))
    return dict(nonnegative_diagonal_exact_distance_cases=records,
                general_shifted_cone_zero_and_negative_cases=general,
                product_uniform_distance_equals_operator_norm_distance_in_this_slice=True,
                distance_to_spectral_cone_not_distance_to_PSD_cone=True,
                analytic_characterization_not_numeric_optimization=True)


def finite_uncertainty():
    eta=F(1,10); coefficient_bound=3*eta/2
    u=v=np.ones(2)/math.sqrt(2)
    perturbation=matrix((eta,eta,eta,eta))
    actual=float(pair_vector(u,v)@perturbation@pair_vector(u,v))
    assert math.isclose(actual,float(coefficient_bound),abs_tol=2e-12)
    assert math.isclose(float(np.linalg.norm(perturbation,ord=2)),float(coefficient_bound),abs_tol=2e-12)
    budgets=[]; gap=F(1,2)
    for eta,outer in ((F(1,10),F(1,10)),(F(1,5),F(1,5)),(F(1,5),F(1,4))):
        total=3*eta/2+outer
        budgets.append(dict(each_coefficient_error_bound_exact=str(eta),
            outer_contour_remainder_bound_exact=str(outer),total_product_error_exact=str(total),
            bad_example_gap_exact=str(gap),strict_exclusion=gap>total,
            residual_exclusion_margin_exact=str(gap-total),
            enlarged_cone_membership=cone_exact((1,1,1,3),total)))
        assert (gap>total)==(not cone_exact((1,1,1,3),total))
    return dict(sharp_coefficient_to_product_error_factor="3/2",
                sharp_example_eta_exact="1/10",sharp_uniform_bound_exact=str(coefficient_bound),
                sharp_example_product_error=actual,budgets=budgets,
                outer_contour_bound_is_independent_physical_input=True,
                outer_bound_inferred_from_low_energy_coefficients=False,
                fixed_finite_energy_contour_is_not_removed_without_bound=True,
                full_same_order_improved_amplitude_required=True,
                bare_Wilson_coefficient_signs_not_the_claim=True,
                error_budgets_are_not_claimed_experimental_certification=True)


def Choi_type_control():
    records=[]; duration=.23; dim=4
    omega=np.eye(dim).reshape(-1)/math.sqrt(dim)
    for name,values in (("violating_coefficient",(1,1,1,3)),("allowed_non_PSD_coefficient",(1,1,4,3))):
        M=matrix(values); eig,vec=np.linalg.eigh(M)
        U=(vec*np.exp(-1j*duration*eig))@vec.conj().T
        psi=np.kron(np.eye(dim),U)@omega
        J=np.outer(psi,psi.conj())
        partial=np.einsum("iaja->ij",J.reshape(dim,dim,dim,dim))
        eigen=np.linalg.eigvalsh(J)
        errors=dict(unitarity=norm(U.conj().T@U-np.eye(dim)),
                    trace_preservation=norm(partial-np.eye(dim)/dim),
                    normalization=abs(float(np.trace(J).real)-1),
                    pure_normalized_Choi=norm(J@J-J))
        assert max(errors.values())<3e-12 and min(eigen)>-3e-12
        records.append(dict(name=name,parameters=list(values),duration=duration,
            coefficient_in_dispersion_cone=cone_exact(values),
            coefficient_matrix_minimum_eigenvalue=float(min(eig)),
            normalized_Choi_minimum_eigenvalue=float(min(eigen)),
            normalized_Choi_eigenvalues=eigen.tolist(),errors=errors,
            auxiliary_unitary_channel_CP_TP=True))
    return dict(cases=records,auxiliary_channel_dimension=4,
                actual_channel_Choi_not_amplitude_tensor=True,
                passing_CP_does_not_imply_dispersion_admissibility=True,
                these_channels_are_not_claimed_to_be_scattering_UV_completions=True)


def run():
    cases=cone_cases(); nonpsd=non_PSD_generator_witness(); repair=minimal_repair_cases()
    uncertainty=finite_uncertainty(); choi=Choi_type_control()
    return dict(round=1023,date="2026-10-08",all_scientific_calibrations_passed=True,
                new_calibration_groups=1,cumulative_test_groups=3801,new_cognitive_axioms=0,
                scope="Two equal positive-mass real scalar species with a global Z2 x Z2 selection rule, not an imposed flavor superselection rule, in the four-parameter crossing slice of the complete same-order light-subtracted forward amplitude; no massless/long-range exchange, analyticity, crossing, unitarity and controlled outer contour are adopted. Actual coherent preparations remain an implementation obligation.",
                pair_index_order=["11","22","12","21"],
                exact_cone="a,b,c >= 0 and |d| <= c + sqrt(a*b)",
                enlarged_cone="a+epsilon,b+epsilon,c+epsilon >= 0 and |d| <= c+epsilon+sqrt((a+epsilon)*(b+epsilon))",
                arbitrary_real_and_complex_product_directions_proved_analytically=True,
                entangled_coefficient_PSD_required=False,
                all_multispecies_EFT_cones_classified=False,
                finite_sampling_proves_global_cone=False,
                quantum_CP_alone_derives_dispersion=False,
                full_UV_completion_constructed=False,
                gravity_forward_pole_problem_solved=False,
                cognition_principles_derive_analyticity_or_mass_gap=False,
                physical_instrument_or_experimental_budget_certified=False,
                four_generator_sufficiency_is_algebraic_not_UV_completion=True,
                cone_calibrations=cases,legal_non_PSD_example=nonpsd,
                tolerance_cone=repair,finite_error_contract=uncertainty,Choi_object_type_control=choi,
                historical_source_sha256={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                                          for name in HISTORICAL})


def compare(fresh,saved,path="result"):
    if isinstance(fresh,dict):
        assert isinstance(saved,dict) and fresh.keys()==saved.keys(),path
        for key in fresh:
            compare(fresh[key],saved[key],path+"."+key)
    elif isinstance(fresh,list):
        assert isinstance(saved,list) and len(fresh)==len(saved),path
        for i,(a,b) in enumerate(zip(fresh,saved)):
            compare(a,b,path+f"[{i}]")
    elif isinstance(fresh,float):
        assert math.isclose(fresh,saved,rel_tol=5e-10,abs_tol=3e-11),(path,fresh,saved)
    else:
        assert fresh==saved,(path,fresh,saved)


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); start=time.perf_counter(); result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        compare(result,json.loads(OUT.read_text("utf8")))
    print(json.dumps(dict(round=1023,passed=True,elapsed_seconds=time.perf_counter()-start,
        cone_cases=len(result["cone_calibrations"]),
        exact_repair_cases=len(result["tolerance_cone"]["nonnegative_diagonal_exact_distance_cases"]),
        bad_example_exact_gap="1/2",cumulative_test_groups=3801),ensure_ascii=False,indent=2))
