"""Limited 1023 delivery verification; never freeze live navigation."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import multifield_dispersion_selection as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
NOTE=HERE.parent/"research_note_1023.md"
OUT=HERE/"research_round_1023_checks.json"
OWN=["multifield_dispersion_selection.py","multifield_dispersion_selection_results.json",
     "review.md","selection_audit.md","input_dependency_update_v0_12.md","NEXT.md","verify_round1023.py"]
HISTORICAL={"archive_1009_/1022/research_round_1022_checks.json","archive_1009_/1022/NEXT.md",
            "archive_935_955/research_note_952.md","archive_990_1008/research_note_993.md",
            "archive_1009_/1009/input_dependency_ledger_v0_1.md",
            "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md"}
LINK=re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$",re.M)
EXCLUDED=re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$",re.S|re.M)
ORDER=[(0,0),(1,1),(0,1),(1,0)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a,b):
    assert math.isclose(a,b,rel_tol=5e-10,abs_tol=3e-11),(a,b)


def cone(params,epsilon=F(0)):
    a,b,c,d=map(F,params); a+=epsilon; b+=epsilon; c+=epsilon
    return min(a,b,c)>=0 and (abs(d)<=c or (abs(d)-c)**2<=a*b)


def expected_matrix(params):
    a,b,c,d=map(float,params)
    return np.array([[a,d/2,0,0],[d/2,b,0,0],[0,0,c,d/2],[0,0,d/2,c]])


def generator_matrix(generators):
    # Deliberately use scalar index products rather than the science einsum path.
    result=np.zeros((4,4))
    for m in generators:
        for row,(i,j) in enumerate(ORDER):
            for col,(k,l) in enumerate(ORDER):
                result[row,col]+=m[i][j]*m[k][l]+m[i][l]*m[k][j]
    return result


def local_links(text):
    masked=EXCLUDED.sub(lambda match:" "*len(match[0]),text)
    for match in LINK.finditer(masked):
        value=match[1] if match[1] is not None else match[2]; stripped=value.strip()
        if stripped.startswith("<") and ">" in stripped:
            target=stripped[1:stripped.index(">")]
        else:
            target=re.split(r"\s+[\"']",stripped,1)[0]
        local=unquote(target.split("#",1)[0])
        if not local or re.match(r"^[a-zA-Z]+:",local) or local.startswith("//"):
            continue
        yield target,local


def verify(prospective=False):
    fresh=science.run()
    science.compare(fresh,json.loads((HERE/"multifield_dispersion_selection_results.json").read_text("utf8")))
    assert fresh["round"]==1023 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"]==1 and fresh["cumulative_test_groups"]==3801
    assert fresh["new_cognitive_axioms"]==0
    assert fresh["pair_index_order"]==["11","22","12","21"]
    assert "equal positive-mass" in fresh["scope"] and "not an imposed flavor superselection" in fresh["scope"]
    for key in ("entangled_coefficient_PSD_required","all_multispecies_EFT_cones_classified",
                "finite_sampling_proves_global_cone","quantum_CP_alone_derives_dispersion",
                "full_UV_completion_constructed","gravity_forward_pole_problem_solved",
                "cognition_principles_derive_analyticity_or_mass_gap","physical_instrument_or_experimental_budget_certified"):
        assert fresh[key] is False,key
    assert fresh["arbitrary_real_and_complex_product_directions_proved_analytically"]
    assert fresh["four_generator_sufficiency_is_algebraic_not_UV_completion"]

    cases=fresh["cone_calibrations"]; assert len(cases)==20
    rows={row["name"]:row for row in cases}
    allowed_count=0; max_direction=0.; max_generator=0.
    for row in cases:
        p=list(map(F,row["parameters"])); a,b,c,d=map(float,p)
        allowed=cone(p); assert allowed==row["in_exact_cone"]
        expected=sorted([(a+b-math.sqrt((a-b)**2+d*d))/2,
                         (a+b+math.sqrt((a-b)**2+d*d))/2,c-d/2,c+d/2])
        for actual,target in zip(row["pair_matrix_eigenvalues"],expected):
            close(actual,target)
        directions=row["direction_calibration"]
        assert directions["real_direction_count"]==10 and directions["complex_direction_count"]==8
        assert directions["sampled_directions_do_not_prove_cone"]
        max_direction=max(max_direction,max(directions["errors"].values()))
        assert max(directions["errors"].values())<3e-11
        if allowed:
            allowed_count+=1
            error=float(np.linalg.norm(generator_matrix(row["generator_matrices"])-expected_matrix(p)))
            assert error<3e-11
            close(error,row["generator_tensor_reconstruction_error"])
            max_generator=max(max_generator,error)
            assert directions["smallest_sampled_real_product"]>-3e-11
            assert directions["smallest_sampled_complex_product"]>-3e-11
        else:
            u=np.array(row["witness_u"]); v=np.array(row["witness_v"])
            close(float(u@u),1); close(float(v@v),1)
            w=np.array([u[i]*v[j] for i,j in ORDER])
            q=float(w@expected_matrix(p)@w)
            close(q,float(F(row["witness_value_exact"]))); close(q,row["witness_value"])
            assert q<0
    assert allowed_count==13
    bad=rows["fixed_flavors_miss_violation"]
    assert bad["parameters"]==["1","1","1","3"] and F(bad["witness_value_exact"])==-F(1,2)
    assert np.array_equal(np.diag(expected_matrix((1,1,1,3))),np.ones(4))
    legal=rows["non_PSD_legal"]
    assert legal["in_exact_cone"] and legal["parameters"]==["1","1","4","3"]
    close(legal["pair_matrix_eigenvalues"][0],-.5)

    example=fresh["legal_non_PSD_example"]
    M=generator_matrix(example["generator_matrices"])
    assert np.linalg.norm(M-expected_matrix((1,1,4,3)))<3e-11
    alternative=[(np.diag([1.,-1.])/math.sqrt(2)).tolist(),
                 (math.sqrt(2)*np.array([[0.,1.],[1.,0.]])).tolist()]
    assert np.linalg.norm(generator_matrix(alternative)-M)<3e-11
    close(example["minimum_matrix_eigenvalue"],-.5)
    z=np.array(example["entangled_vector"]); close(float(z@M@z),-.5)
    assert np.linalg.matrix_rank(np.array([[z[0],z[2]],[z[3],z[1]]]))==example["entangled_vector_Schmidt_rank"]==2
    assert example["generator_decomposition_proves_product_positivity"]
    assert example["coefficient_tensor_is_not_a_Choi_operator"] and example["positive_generators_do_not_prove_UV_completion"]

    tolerance=fresh["tolerance_cone"]
    repair=tolerance["nonnegative_diagonal_exact_distance_cases"]; assert len(repair)==6
    max_repair=0.
    for row in repair:
        a,b,c,d=map(F,row["parameters"]); k=abs(d)-c
        assert min(a,b,c)>=0 and k>0 and k*k>a*b
        gap=(k*k-a*b)/(a+b+2*k); x=(a+k)/(b+k)
        assert gap==F(row["minimal_repair_exact"])>0
        assert x==F(row["squared_component_ratio_exact"])
        assert (a+b*x*x-2*k*x)/(1+x)**2==-gap==F(row["witness_value_exact"])
        assert not cone((a,b,c,d),gap/2) and cone((a,b,c,d),gap) and cone((a,b,c,d),3*gap/2)
        assert [a+gap,b+gap,c+gap,d]==list(map(F,row["shifted_parameters"]))
        assert (abs(d)-c-gap)**2==(a+gap)*(b+gap)
        close(row["operator_norm_repair"],float(gap))
        close(row["witness_value"],-float(gap)); close(row["shifted_witness_value"],0)
        assert not row["half_gap_repair_allowed"] and row["exact_gap_repair_allowed"] and row["one_and_half_gap_repair_allowed"]
        max_repair=max(max_repair,row["shift_identity_error"],row["shifted_generator_reconstruction_error"],abs(row["shifted_witness_value"]))
    assert F(repair[0]["minimal_repair_exact"])==F(1,2)
    general=tolerance["general_shifted_cone_zero_and_negative_cases"]; assert len(general)==5
    assert [row["in_enlarged_cone"] for row in general]==[True,False,True,True,False]
    for row in general:
        assert cone(row["parameters"],F(row["epsilon_exact"]))==row["in_enlarged_cone"]
    assert tolerance["product_uniform_distance_equals_operator_norm_distance_in_this_slice"]
    assert tolerance["distance_to_spectral_cone_not_distance_to_PSD_cone"]
    assert tolerance["analytic_characterization_not_numeric_optimization"]

    errors=fresh["finite_error_contract"]
    assert F(errors["sharp_coefficient_to_product_error_factor"])==F(3,2)
    eta=F(errors["sharp_example_eta_exact"]); sharp=F(errors["sharp_uniform_bound_exact"])
    assert sharp==3*eta/2==F(3,20); close(errors["sharp_example_product_error"],float(sharp))
    # Independent two-block norm bound and its saturating normalized vector.
    E=expected_matrix((eta,eta,eta,eta)); w=np.ones(4)/2
    close(float(w@E@w),float(sharp)); close(float(np.linalg.norm(E,ord=2)),float(sharp))
    budgets=errors["budgets"]; assert len(budgets)==3
    for row in budgets:
        eta=F(row["each_coefficient_error_bound_exact"]); outer=F(row["outer_contour_remainder_bound_exact"])
        total=3*eta/2+outer; gap=F(row["bad_example_gap_exact"])
        assert total==F(row["total_product_error_exact"]) and gap==F(1,2)
        assert gap-total==F(row["residual_exclusion_margin_exact"])
        assert row["strict_exclusion"]==(gap>total)
        assert row["enlarged_cone_membership"]==cone((1,1,1,3),total)
    assert [row["strict_exclusion"] for row in budgets]==[True,False,False]
    assert errors["outer_contour_bound_is_independent_physical_input"]
    assert not errors["outer_bound_inferred_from_low_energy_coefficients"]
    assert errors["fixed_finite_energy_contour_is_not_removed_without_bound"]
    assert errors["full_same_order_improved_amplitude_required"] and errors["bare_Wilson_coefficient_signs_not_the_claim"]
    assert errors["error_budgets_are_not_claimed_experimental_certification"]

    controls=fresh["Choi_object_type_control"]; assert controls["auxiliary_channel_dimension"]==4
    assert len(controls["cases"])==2
    for row in controls["cases"]:
        assert row["coefficient_in_dispersion_cone"]==cone(row["parameters"])
        close(row["coefficient_matrix_minimum_eigenvalue"],-.5)
        eigen=row["normalized_Choi_eigenvalues"]
        assert len(eigen)==16 and min(eigen)>-3e-11
        close(sum(eigen),1); close(eigen[-1],1)
        assert max(row["errors"].values())<3e-11 and row["auxiliary_unitary_channel_CP_TP"]
    assert not controls["cases"][0]["coefficient_in_dispersion_cone"]
    assert controls["cases"][1]["coefficient_in_dispersion_cone"]
    assert controls["actual_channel_Choi_not_amplitude_tensor"]
    assert controls["passing_CP_does_not_imply_dispersion_admissibility"]
    assert controls["these_channels_are_not_claimed_to_be_scattering_UV_completions"]

    assert set(fresh["historical_source_sha256"])==HISTORICAL
    for name,digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name)==digest,name
    assets=[NOTE]+[HERE/name for name in OWN]; assert len(assets)==8
    assert all(path.name not in {"README.md","RESEARCH_STATE.md","research_direction.md"} for path in assets)
    count=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix==".md":
            text=path.read_text("utf-8-sig"); assert text.count("$$")%2==0,path
            for target,local in local_links(text):
                resolved=(path.parent/local.replace("\\","/")).resolve()
                assert resolved.exists() or (prospective and resolved==OUT),(path,target)
                count+=1
    assert all(f"## {i}." in NOTE.read_text("utf8") for i in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1023,date="2026-10-08",all_delivery_checks_passed=True,
        scientific_result_reproduced=True,new_calibration_groups=1,cumulative_research_groups=3801,
        local_links_checked=count,frozen_current_files=len(assets),historical_input_files=len(HISTORICAL),
        live_navigation_frozen=False,neighboring_round_frozen=False,cone_cases=20,
        allowed_cone_calibrations=allowed_count,real_direction_checks=200,complex_direction_checks=160,
        exact_minimal_repair_cases=6,general_tolerance_boundary_cases=5,
        fixed_flavor_false_positive_gap="1/2",legal_non_PSD_minimum_eigenvalue=-.5,
        coefficient_error_factor="3/2",outer_contour_budgets=3,
        auxiliary_Choi_CP_TP_controls=2,maximum_direction_identity_error=max_direction,
        maximum_generator_reconstruction_error=max_generator,maximum_repair_identity_error=max_repair,
        actual_instrument_or_UV_completion_constructed=False,gravity_forward_pole_solved=False,
        new_cognitive_axiom=False,goal_completed=False,visual_checks_performed=False,
        analytic_scope=fresh["scope"],historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(path.relative_to(ROOT)).replace("\\","/"):sha(path) for path in assets})


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); result=verify(prospective=args.write)
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        assert result==json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},ensure_ascii=False,indent=2))
