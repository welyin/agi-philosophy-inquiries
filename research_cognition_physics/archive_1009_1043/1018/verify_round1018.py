"""Reproduce round 1018 and freeze only its explicit delivery/input files.

The local-link checker is inlined from the historical verification convention,
so verification does not import a changing migration/navigation inventory.
No README, RESEARCH_STATE, research_direction, neighboring delivery or goal
state is included in this round's current-file hash manifest.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import re
from urllib.parse import unquote

import gravity_ideal_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1018.md"
OUT = HERE / "research_round_1018_checks.json"
OWN = ["gravity_ideal_selection.py", "gravity_ideal_selection_results.json",
       "review.md", "NEXT.md", "input_dependency_update_v0_7.md",
       "verify_round1018.py"]
HISTORICAL = {
    "archive_301_341/research_note_326.md",
    "archive_342_369/research_note_358.md",
    "archive_1009_/research_note_1017.md",
    "archive_1009_/1017/stress_source_commutant_results.json",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def all_zero(value):
    if isinstance(value, list):
        return all(all_zero(v) for v in value)
    return F(value) == 0


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE / OWN[1]).read_text("utf8")))
    assert fresh["round"] == 1018
    assert fresh["new_calibration_groups"] == 1
    assert fresh["cumulative_test_groups"] == 3796
    assert fresh["new_cognitive_axioms"] == 0
    assert fresh["all_scientific_calibrations_passed"]
    for name in ("first_order_minimal_data_adopted", "regular_open_algebra_allowed",
                 "first_order_recoil_exists_by_analytic_lemma",
                 "all_parameter_classification_is_analytic_not_sampling"):
        assert fresh[name] is True, name
    for name in ("higher_local_corrections_arbitrarily_assumed_zero",
                 "M1_assumed_zero", "nonlinear_completion_unique",
                 "arbitrary_R2_or_open_terms_numerically_verified",
                 "full_backreaction_solution_numerically_produced",
                 "Einstein_action_derived_by_this_code", "gravity_generated",
                 "number_of_spin2_fields_selected"):
        assert fresh[name] is False, name

    killing = fresh["killing_jet_calibration"]
    assert killing["bracket"] == ["0", "0", "1", "0"]
    assert killing["leading_commutator"] == ["2/7", "1/6"]
    assert killing["delta_eta_after_delta_xi"] == ["2/7", "1/6"]
    assert all_zero(killing["delta_xi_after_delta_eta"])
    assert all_zero(killing["D_xi"]) and all_zero(killing["D_eta"])
    assert killing["parameter_jets_retained"]
    assert killing["all_higher_parameter_derivatives_zero_by_affinity"]
    assert killing["all_derivatives_of_D_xi_D_eta_zero"]

    cases = fresh["exact_coupling_cases"]
    assert len(cases) == 7
    assert [case["obstruction_zero"] for case in cases] == [True, False, False, False, True, True, False]
    assert cases[1]["D"][0][1] == "6"
    assert cases[2]["D"][0][0] == "-1"
    assert cases[3]["D"][0][0] == "1"
    pair_count = 0
    for case in cases:
        n = len(case["q"])
        assert len(case["sequential_checks"]) == n*n
        for item in case["sequential_checks"]:
            a, b = item["a"], item["b"]
            coefficient = F(case["D"][a][b])
            assert list(map(F, item["residual"])) == [coefficient*F(2, 7), coefficient*F(1, 6)]
            assert list(map(F, item["residual"])) == [F(x)-F(y) for x, y in
                zip(item["commutator"], item["complete_B1_action"])]
            pair_count += 1
    assert pair_count == 45
    assert fresh["finite_support_classification_samples"] == 188

    ode = fresh["scalar_on_shell_local_jet_calibration"]
    assert ode["exact_commutant_dimension"] == 1
    assert ode["degree"] == 10 and ode["verified_equation_orders"] == list(range(9))
    assert len(ode["series_coefficients"]) == 2
    assert all(len(row) == 11 for row in ode["series_coefficients"])
    assert all_zero(ode["equation_residual_coefficients"])
    assert ode["finite_Taylor_polynomial_claimed_exact_solution"] is False
    assert ode["actual_local_analytic_solution_uses_ODE_existence"] is True

    rotation = fresh["source_alignment_control"]
    assert rotation["correct_defect_prime"] == [["0", "0"], ["0", "-2"]]
    assert rotation["q_prime_div_sqrt2"] == ["1", "0"]
    assert rotation["deleted_A_1_22_div_sqrt2"] == "1"
    assert all_zero(rotation["all_mixed_product_entries_wrongly_deleted_defect"])
    assert rotation["deleting_mixed_product_changes_the_theory"] is True
    assert rotation["source_alignment_alone_sufficient"] is False
    assert rotation["maximum_numpy_cross_check_residual"] < 3e-14
    rational = fresh["rational_orthogonal_rotation_controls"]
    assert [len(item["O"]) for item in rational] == [3, 4]
    assert all(item["covariance_exact"] for item in rational)

    recoil = fresh["first_order_recoil_initial_constraints"]
    assert len(recoil["samples"]) == 4 and recoil["exact_initial_constraints_pass"]
    assert recoil["numerical_PDE_solution_claimed"] is False
    assert recoil["source_is_not_claimed_to_be_a_complete_scalar_backreaction"] is True
    for sample in recoil["samples"]:
        assert sample["spatial_divergence_p"] == sample["laplacian_W"]
        assert all_zero(sample["H_nu"]) and all_zero(sample["dt_H_nu"])
        lapu = sum((F(sample["u_hessian"][i][i]) for i in range(3)), F(0))
        assert F(sample["J_0nu"][0]) == lapu
        assert list(map(F, sample["J_0nu"][1:])) == [-F(v) for v in sample["laplacian_W"]]
    separated = fresh["separated_matter_boundary"]
    assert separated["both_sources_nonzero"] and separated["distinct_decoupled_matter_sectors"]
    assert all_zero(separated["product_residuals"])

    assert set(fresh["historical_source_sha256"]) == HISTORICAL
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    assert len(assets) == 7
    assert all(path.name not in {"README.md", "RESEARCH_STATE.md", "research_direction.md"} for path in assets)
    count = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == ".md":
            text = path.read_text("utf-8-sig")
            assert text.count("$$") % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent / local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                count += 1
    note = NOTE.read_text("utf8")
    assert all(f"## {n}." in note for n in range(1, 11))
    review = (HERE / "review.md").read_text("utf8")
    assert "独立科学签审" in review

    return dict(round=1018, date="2026-10-08", all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3796, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(HISTORICAL),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        exact_coupling_cases=7, explicit_killing_pair_checks=pair_count,
        finite_support_classification_samples=188,
        rational_orthogonal_rotation_controls=2, recoil_initial_data_samples=4,
        exact_scalar_ODE_residual_orders=list(range(9)),
        finite_samples_not_general_proof=True,
        first_order_minimal_data_adopted=True,
        arbitrary_higher_local_corrections_not_set_zero=True,
        regular_open_algebra_allowed=True, M1_assumed_zero=False,
        first_order_recoil_exists_by_analytic_lemma=True,
        rotated_source_alone_not_sufficient=True,
        separated_matter_counterexample_boundary_preserved=True,
        nonlinear_completion_unique=False, gravity_generated=False,
        new_cognitive_axiom=False, goal_completed=False,
        visual_checks_performed=False, analytic_scope=fresh["scope"],
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(path.relative_to(ROOT)).replace("\\", "/"):sha(path)
                       for path in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({key:value for key, value in result.items()
                      if not key.endswith("sha256")}, ensure_ascii=False, indent=2))
