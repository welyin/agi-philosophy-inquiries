"""Recompute round 1022 and freeze only its scientific delivery and old inputs."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import chiral_higgs_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1022.md"
OUT = HERE / "research_round_1022_checks.json"
OWN = ["chiral_higgs_selection.py", "chiral_higgs_selection_results.json", "review.md",
       "selection_audit.md", "input_dependency_update_v0_11.md", "NEXT.md",
       "verify_round1022.py"]
HISTORICAL = {
    "archive_531_553/research_note_531.md",
    "archive_585_628/research_note_627.md",
    "archive_1009_/research_note_1012.md",
    "archive_1009_/1021/research_round_1021_checks.json",
    "archive_1009_/1021/input_dependency_update_v0_10.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=4e-10, abs_tol=3e-11), (a, b)


def local_links(text):
    masked = EXCLUDED.sub(lambda match: " "*len(match[0]), text)
    for match in LINK.finditer(masked):
        value = match[1] if match[1] is not None else match[2]
        stripped = value.strip()
        if stripped.startswith("<") and ">" in stripped:
            target = stripped[1:stripped.index(">")]
        else:
            target = re.split(r"\s+[\"']", stripped, 1)[0]
        local = unquote(target.split("#", 1)[0])
        if not local or re.match(r"^[a-zA-Z]+:", local) or local.startswith("//"):
            continue
        yield target, local


def check_charges(record, linear=0, cubic=0, chiral=True):
    q = [F(x) for x in record["charges"]]
    assert sum(q) == F(record["linear_anomaly"]) == linear
    assert sum(x**3 for x in q) == F(record["cubic_anomaly"]) == cubic
    assert sum(x*x for x in q) == F(record["quadratic_sum"])
    assert record["no_zero_charge"] == all(qi != 0 for qi in q)
    assert record["no_opposite_charges"] == all(-qi not in q for qi in q)
    if chiral:
        assert record["no_zero_charge"] and record["no_opposite_charges"]
    return q


def check_mass(q, mass):
    n = len(q); matrix = np.zeros((n,n))
    for i,j,value in mass["mass_terms"]:
        assert q[i]+q[j] in (mass["higgs_charge"],-mass["higgs_charge"])
        matrix[i,j] += value
        if i != j:
            matrix[j,i] += value
    assert np.array_equal(matrix,matrix.T)
    assert np.linalg.matrix_rank(matrix) == mass["numerical_rank"] == n
    actual = np.linalg.eigvalsh(matrix)
    expected = sorted(mass["Takagi_singular_values"])
    for a,b in zip(sorted(abs(actual)),expected):
        close(a,b)
    close(float(np.linalg.det(matrix)),float(F(mass["exact_determinant"])))
    assert max(mass["errors"].values()) < 3e-10
    assert mass["arbitrary_same_charge_unitary_covariance_is_analytic"]
    assert mass["arbitrary_coherent_mass_matrix_selection_is_analytic"]


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh,json.loads((HERE/"chiral_higgs_selection_results.json").read_text("utf8")))
    assert fresh["round"] == 1022 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"] == 1 and fresh["cumulative_test_groups"] == 3800
    assert fresh["new_cognitive_axioms"] == 0
    for key in ("spatial_dimension_or_SM_group_selected", "unique_chiral_charge_spectrum_selected",
                "quantum_gravity_or_nonperturbative_UV_completion_proved",
                "actual_Wilson_instrument_or_preparation_built"):
        assert fresh[key] is False, key
    for key in ("all_mass_mixings_covered_by_analytic_cycle_cover_argument",
                "arbitrary_N_existence_from_8_and_12_template_union",
                "templates_have_no_cross_opposite_charges",
                "charged_Higgs_VEV_assumed_not_derived",
                "single_Higgs_renormalizable_mass_rule_is_extra_input",
                "mass_generation_does_not_claim_discrete_gauge_remnant_removed"):
        assert fresh[key],key

    examples = fresh["existence_calibrations"]
    assert [row["number"] for row in examples] == [8,12,16,20,24,28,32]
    maximum_covariance_error = 0.
    for row in examples:
        n = row["number"]; q = check_charges(row)
        assert len(q) == n == 8*row["copies_of_8"]+12*row["copies_of_12"]
        assert all(s in (-2,2) for s in row["pair_sums"])
        assert row["pair_sums"].count(2) == row["pair_sums"].count(-2) == n//4
        assert row["primitive_integer_charge_gcd"] == 1
        check_mass(q,row["mass"])
        assert F(row["mass"]["exact_determinant"]) == math.factorial(n//2)**2
        for a,b in zip(sorted(row["mass"]["Takagi_singular_values"]),
                       [float(j) for j in range(1,n//2+1) for _ in range(2)]):
            close(a,b)
        maximum_covariance_error=max(maximum_covariance_error,max(row["mass"]["errors"].values()))
    # The two full templates have compatible charge signs under direct sums.
    joint=examples[0]["charges"]+examples[1]["charges"]
    assert all(q and -q not in joint for q in joint)

    walks=fresh["odd_cycle_lemma_calibration"]
    assert walks["total_odd_sign_words"] == 682
    assert walks["cases"] == [dict(length=n,sign_words_checked=2**n) for n in (1,3,5,7,9)]
    assert walks["every_walk_has_a_self_pairable_charge"]
    assert not walks["both_loop_signs_in_one_orbit"]
    assert walks["all_lengths_proved_by_charge_orbit_graph_not_enumeration"]
    # Independent exact identity behind the N=4 exclusion; finite values only calibrate it.
    for a,b,c in ((1,2,3),(-2,4,5),(F(2,3),F(7,5),F(-9,4)),(2,-2,7)):
        d=-a-b-c
        assert a**3+b**3+c**3+d**3 == -3*(a+b)*(a+c)*(b+c)

    control=fresh["removed_assumption_controls"]
    higher=control["higher_dimension"]; q=check_charges(higher)
    assert higher["number"] == higher["numerical_rank"] == 6
    assert higher["signed_Higgs_powers"] == [-2,3,-1]
    assert higher["operator_dimensions"] == [5,6,4]
    assert F(higher["exact_determinant"]) == -36
    for (i,j,_),k in zip(higher["mass_terms"],higher["signed_Higgs_powers"]):
        assert q[i]+q[j]+k*higher["higgs_charge"] == 0
    two=control["second_Higgs"]; q=check_charges(two)
    assert two["number"] == two["numerical_rank"] == 5
    assert two["higgs_charges"] == [10,15] and two["operator_dimensions"] == [4,4,4]
    assert F(two["exact_determinant"]) == 18
    for (i,j,_),h in zip(two["mass_terms"],two["signed_Higgs_charges"]):
        assert q[i]+q[j]+h == 0
    other=control["other_controls"]
    assert [row["number"] for row in other] == [6,4,2]
    for row,lin,cub,chiral in zip(other,(6,0,0),(0,-162,0),(True,True,False)):
        q=check_charges(row,lin,cub,chiral); check_mass(q,row["mass"])
        maximum_covariance_error=max(maximum_covariance_error,max(row["mass"]["errors"].values()))
    assert not other[-1]["no_opposite_charges"]

    freedom=fresh["residual_charge_freedom"]
    candidates=freedom["candidates"]; assert len(candidates)==2
    raw=[]
    for row in candidates:
        q=check_charges(row); check_mass(q,row["mass"])
        maximum_covariance_error=max(maximum_covariance_error,max(row["mass"]["errors"].values()))
        assert row["mass"]["higgs_charge"]==26
        assert F(row["quadratic_sum"])==18252 and row["primitive_integer_charge_gcd"]==1
        raw.append(q)
    assert sorted(raw[0]) != sorted(raw[1]) and sorted(raw[0]) != sorted(-q for q in raw[1])
    scaled=[[q/13 for q in qs] for qs in raw]; x=F(1,10)
    fourth=sum(q**4 for q in scaled[0])-sum(q**4 for q in scaled[1])
    sixth=sum(abs(q)**6 for qs in scaled for q in qs)
    assert fourth==F(freedom["fourth_moment_difference"])==F(9596160,28561)
    assert sixth==F(freedom["sixth_absolute_moment_sum"])==F(3140425416,28561)
    signal=fourth*x**4/384-sixth*x**6/11520
    assert signal==F(freedom["strict_effect_gap_lower_exact"])==F(1068668941,13709280000000)>0
    close(float(signal),freedom["strict_effect_gap_lower"])
    effects=[.5+sum(math.cos(float(q*x)) for q in qs)/16 for qs in scaled]
    for a,b in zip(effects,freedom["effect_probabilities"]):
        close(a,b); assert 0<=a<=1
    close(effects[0]-effects[1],freedom["numerical_effect_gap"])
    assert effects[0]-effects[1]>float(signal)
    lipschitz=sum(abs(q) for qs in scaled for q in qs)/16
    assert lipschitz==F(freedom["difference_Lipschitz_bound_exact"])==F(175,52)
    width=F(freedom["smooth_packet_support_halfwidth_x"]); assert width==F(1,10**6)
    packet_lower=signal-lipschitz*width
    assert packet_lower==F(freedom["packet_averaged_gap_lower_exact"])==F(1022531941,13709280000000)>0
    close(float(packet_lower),freedom["packet_averaged_gap_lower"])
    assert freedom["any_probability_distribution_in_support_satisfies_bound"]
    assert freedom["compact_smooth_wavepacket_has_finite_electric_quadratic_form_by_analysis"]
    assert not freedom["packet_or_physical_instrument_constructed"]
    assert set(fresh["historical_source_sha256"])==HISTORICAL
    for name,digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name)==digest,name

    assets=[NOTE]+[HERE/name for name in OWN]
    assert len(assets)==8
    assert all(path.name not in {"README.md","RESEARCH_STATE.md","research_direction.md"} for path in assets)
    count=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix==".md":
            text=path.read_text("utf-8-sig")
            assert text.count("$$")%2==0,path
            for target,local in local_links(text):
                resolved=(path.parent/local.replace("\\","/")).resolve()
                assert resolved.exists() or (prospective and resolved==OUT),(path,target)
                count+=1
    assert all(f"## {i}." in NOTE.read_text("utf8") for i in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1022,date="2026-10-08",all_delivery_checks_passed=True,
                scientific_result_reproduced=True,new_calibration_groups=1,cumulative_research_groups=3800,
                local_links_checked=count,frozen_current_files=len(assets),historical_input_files=len(HISTORICAL),
                live_navigation_frozen=False,neighboring_round_frozen=False,
                arbitrary_coherent_mass_selection_is_analytic=True,
                species_number_classification="N in 4*positive integers and N >= 8, within the stated candidate class",
                existence_calibrations=7,odd_cycle_sign_words=682,
                removed_assumption_controls=5,maximum_matrix_covariance_error=maximum_covariance_error,
                same_quadratic_charge_sum_does_not_identify_spectrum=True,
                exact_Wilson_gap_lower=str(signal),exact_packet_gap_lower=str(packet_lower),
                finite_holonomy_width_certified=True,actual_Wilson_instrument_or_preparation_built=False,
                SM_selected=False,spacetime_dimension_selected=False,
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
