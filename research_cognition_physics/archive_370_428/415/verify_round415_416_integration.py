"""Integrate rounds 415 and 416 without changing or rerunning historical science."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round413_414_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / "round415_416_integration_checks.json"
CONFIG = {415: "bounded_control_budget_audit", 416: "finite_steering_budget_audit"}


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    saved_by_round = {}
    links = formulas = 0
    for number, stem in CONFIG.items():
        checked = core.read(HERE / f"research_round_{number}_checks.json")
        assert checked["all_reported_checks_passed"]
        assert checked["parallel_batch"] == [415, 416]
        assert checked["scientific_base_through_round"] == 414
        assert checked["batch_scientific_dependencies"] == []
        assert checked["fresh_tests"] == dict(run=8, failures=0, errors=0)
        assert len(checked["new_file_hashes"]) == 3
        for name, sha in checked["new_file_hashes"].items():
            assert core.digest(HERE / name) == sha, name
        saved = core.read(HERE / (stem + "_results.json"))
        assert saved["checks"] == checked["fresh_tests"]
        saved_by_round[number] = saved
        note = HERE / f"research_note_{number}.md"
        text_checks = core.text_checks(note)
        assert text_checks == checked["text_checks"]
        formulas += text_checks["display_formulas"]
        for link in core.link_parser()(note.read_text(encoding="utf-8")):
            assert (HERE / link).resolve().exists(), (note.name, link)
            links += 1
        ast.parse((HERE / (stem + ".py")).read_text(encoding="utf-8"))
    a, b = saved_by_round[415], saved_by_round[416]
    assert a["finite_bounded_generators_are_extra_input"]
    assert a["measurable_convex_control_access_is_extra_input"]
    assert a["optimal_cost_attained_in_declared_measurable_control_class"]
    for flag in ("compactness_requires_finite_hilbert_dimension", "finite_program_exact_optimality_claimed",
                 "endpoint_resource_cost_inferred_continuous", "scaling_continuity_inferred_from_amplitude_adjustment",
                 "single_device_autonomous_implementation_proved", "budget_identified_with_total_internal_resource_cost",
                 "actual_position_or_dimension_generated", "new_cognitive_axiom_adopted", "phase_closure_triggered"):
        assert not a[flag], flag
    assert b["euclidean_rotation_action_is_input"] and b["euclidean_sublevel_compactness_proved"]
    assert b["uniform_limit_is_on_declared_euclidean_windows"]
    for flag in ("input_position_dimension_is_derived", "arbitrary_rotation_is_a_free_primitive",
                 "steering_cost_omitted", "universal_rate_from_density_alone", "all_round385_conditions_satisfied",
                 "proper_pointed_gh_convergence_claimed", "cognition_to_space_completed"):
        assert not b[flag], flag
    for link in core.link_parser()((HERE / "spatial_premise_closure_audit.md").read_text(encoding="utf-8")):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    for name in ("verify_budget_origin_rounds.py", Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_414") == 553
    assert result.pop("unchanged_prior_science_hashes_231_412") == 547
    assert result["stage_saved_tests"] == 1883
    assert result["total_protected_evidence_hashes"] == 578
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)", (HERE.parent / "research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 416 and int(latest[2]) >= 1899
    result.update(date="2026-09-24", rounds=[415, 416],
        execution_mode="two independent budget-origin rounds from frozen 414",
        scientific_base_through_round_by_round={415: 414, 416: 414},
        fresh_tests_by_round={415: 8, 416: 8}, fresh_tests=16, stage_saved_tests=1899,
        science_hashes_verified_231_416=559, unchanged_prior_science_hashes_231_414=553,
        total_protected_evidence_hashes=584, unchanged_prior_evidence_hashes=578,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"] + links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"] + 3,
        bounded_control_budget_compactness_proved=True,endpoint_cost_continuity_proved=False,
        finite_paid_steering_large_scale_isotropy_proved=True,
        microscopic_control_cost_identified_with_spatial_topology=False,
        actual_endpoint_contract_still_input=True, actual_subject_position_generated=False,
        spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, imaging_branch_deferred=True,
        old_scientific_experiments_rerun=False, independent_final_agent_review_completed=True,
        all_reported_checks_passed=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-checks", action="store_true")
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
