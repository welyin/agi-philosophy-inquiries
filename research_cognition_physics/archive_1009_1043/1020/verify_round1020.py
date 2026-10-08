"""Reproduce and verify the limited round 1020 delivery, excluding live navigation."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import passivity_reference_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent/"research_note_1020.md"
OUT = HERE/"research_round_1020_checks.json"
OWN = ["passivity_reference_selection.py", "passivity_reference_selection_results.json",
       "review.md", "input_dependency_update_v0_9.md", "NEXT.md",
       "next_selection_audit.md", "verify_round1020.py"]
HISTORICAL = {
    "archive_1009_/research_note_1019.md",
    "archive_1009_/1019/spin_statistics_locality_results.json",
    "archive_301_341/research_note_305.md",
    "archive_370_428/research_note_422.md",
    "archive_923_934/research_note_929.md",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-14), (a, b)


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


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE/"passivity_reference_selection_results.json").read_text("utf8")))
    assert fresh["round"] == 1020 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"] == 1 and fresh["cumulative_test_groups"] == 3798
    assert fresh["new_cognitive_axioms"] == 0
    for key in ("complete_passivity_generated_from_internal_resource_accounting",
                "mature_complete_passivity_theorem_claimed_as_new",
                "finite_copy_samples_claimed_to_prove_complete_passivity",
                "full_FUCP_nonimplication_proved"):
        assert fresh[key] is False, key
    assert fresh["physical_H_positive_in_all_matrix_examples"]
    assert fresh["unrestricted_operation_menu_is_an_extra_input"]
    assert fresh["all_independent_copies_is_stronger_than_fixed_finite_task_budget"]
    assert fresh["positive_reference_generator_requires_ground_support_in_stated_matrix_scope"]
    assert fresh["no_unique_gap_ratio_selected"]

    families = fresh["exact_finite_copy_families"]
    assert [f["guaranteed_finite_window"] for f in families] == [1,2,3,6]
    assert [f["eta"] for f in families] == ["2/3", "4/5", "6/7", "12/13"]
    assert [f["first_active_copy_count"] for f in families] == [5,9,11,19]
    all_families = families+[fresh["equality_boundary"]]
    scans = 0; classes = 0; comparisons = 0
    for family in all_families:
        eta = F(family["eta"])
        k = 1
        while eta**k >= F(1,2):
            k += 1
        first = 2*k+1
        assert family["k_first_strict_crossing"] == k
        assert family["first_active_copy_count"] == first
        assert F(family["exact_last_nonstrict_power"]) == eta**(k-1) >= F(1,2)
        assert F(family["exact_first_strict_power"]) == eta**k < F(1,2)
        assert len(family["scans"]) == first+1
        for row in family["scans"]:
            n = row["copies"]
            criterion = eta**((n-1)//2) >= F(1,2)
            assert row["passive"] == row["adjacent_energy_envelopes_passive"] == row["exact_criterion"] == criterion
            assert row["represented_basis_states"] == 3**n
            assert row["occupation_classes"] == (n+1)*(n+2)//2
            assert row["equal_energy_pairs_excluded"]
            assert (row["inversion_pairs"] == 0) == row["passive"]
            scans += 1; classes += row["occupation_classes"]; comparisons += row["distinct_energy_comparisons"]
        window = family["guaranteed_finite_window"]
        if window is not None:
            assert all(row["passive"] for row in family["scans"][:window])
            assert eta**window >= F(window+1, 2*window+1) > F(1,2)
        z = F(3,2)+eta/4
        witness = family["exact_swap_witness"]
        high, low = F(1,2)**first/z**first, (eta/4)**k/z**first
        assert F(witness["high_individual_population"]) == high
        assert F(witness["low_individual_population"]) == low
        assert F(witness["extracted_work"]) == high-low > 0
        close(witness["extracted_work_decimal"], float(high-low))
        assert witness["high_energy"]-witness["low_energy"] == 1
        assert witness["high_occupation"] == [0, first, 0]
        assert witness["low_occupation"] == [k+1, 0, k]
        assert witness["degeneracy_factor_not_multiplied_into_individual_population"]
        assert witness["physical_battery_or_cyclic_controller_implemented"] is False
    assert (scans, classes, comparisons) == (54,2675,143151)
    boundary = fresh["equality_boundary"]
    assert boundary["eta"] == "1/2" and boundary["first_active_copy_count"] == 5
    assert boundary["scans"][2]["passive"] and boundary["scans"][3]["passive"]
    assert not boundary["scans"][4]["passive"]

    gns = fresh["GNS_reference_spectrum"]
    assert gns["support_cases"] == len(gns["support_enumeration"]) == 32
    assert gns["mixed_ground_positive_examples"] == 5
    for row in gns["support_enumeration"]:
        e, support = row["energies"], row["support"]
        differences = sorted(e[i]-e[j] for j in support for i in range(len(e)))
        assert row["matrix_eigenvalues"] == differences == row["expected_energy_differences"]
        assert row["GNS_dimension"] == len(e)*len(support)
        ground = all(e[j] == min(e) for j in support)
        assert row["GNS_generator_positive"] == row["support_in_ground_eigenspace"] == ground
        assert row["stationary"] and row["physical_H_positive"]
    thermal = gns["thermal_qubit"]
    assert thermal["matrix_eigenvalues"] == [-1.,0.,0.,1.]
    assert thermal["physical_H_positive"] and not thermal["GNS_generator_positive"]
    assert gns["thermal_all_copy_passivity_is_analytic_Gibbs_identity"]
    assert gns["canonical_GNS_generator_not_same_as_physical_H"]
    assert gns["full_QFT_forward_cone_derived"] is False
    assert len(gns["finite_thermal_copy_calibrations"]) == 6

    tfd = fresh["purification_boundary"]
    close(tfd["purification_purity"], 1)
    close(tfd["physical_state_purity"], 5/9)
    close(tfd["plus_energy_variance"], 8/9)
    assert tfd["minus_generator_eigenvalues"] == [-1.,0.,0.,1.]
    assert tfd["plus_generator_eigenvalues"] == [0.,1.,1.,2.]
    assert tfd["minus_generator_fixes_omega"] and not tfd["plus_generator_fixes_omega"]
    assert tfd["vector_purity_does_not_imply_algebraic_state_purity"]
    assert tfd["auxiliary_copy_not_treated_as_free_physical_resource"]
    freedom = fresh["remaining_gap_freedom"]
    assert [row["lambda_gap_ratio"] for row in freedom] == ["2", "3", "5/2"]
    for row, expected in zip(freedom, [1.,0.,.5]):
        close(row["common_plus_read_probability"], expected)
        assert row["Gibbs_all_copy_passivity_by_analytic_identity"]
        assert row["all_FUCP_and_all_physical_recovery_claimed"] is False

    assert set(fresh["historical_source_sha256"]) == HISTORICAL
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name) == digest, name
    assets = [NOTE]+[HERE/name for name in OWN]
    assert len(assets) == 8
    assert all(path.name not in {"README.md", "RESEARCH_STATE.md", "research_direction.md"} for path in assets)
    count = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == ".md":
            text = path.read_text("utf-8-sig")
            assert text.count("$$") % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent/local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                count += 1
    assert all(f"## {i}." in NOTE.read_text("utf8") for i in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1020, date="2026-10-08", all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3798, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(HISTORICAL),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        exact_copy_scenarios=scans, occupation_classes=classes,
        exact_distinct_energy_comparisons=comparisons,
        first_active_copies=[5,9,11,19], equality_first_active=5,
        stationary_GNS_support_cases=32, mixed_ground_positive_examples=5,
        purification_counterexample_preserved=True,
        canonical_GNS_generator_distinguished_from_physical_H=True,
        finite_test_window_not_complete_passivity=True,
        unrestricted_cyclic_operations_explicit_input=True,
        dimensionless_gap_freedom_preserved=True,
        complete_passivity_generated_from_internal_resources=False,
        all_FUCP_nonimplication_proved=False, new_cognitive_axiom=False,
        goal_completed=False, visual_checks_performed=False, analytic_scope=fresh["scope"],
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(path.relative_to(ROOT)).replace("\\", "/"):sha(path) for path in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = verify(prospective=args.write)
    if args.write:
        serialized = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
        with OUT.open("x", encoding="utf8") as stream:
            stream.write(serialized)
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, ensure_ascii=False, indent=2))
