"""Integrate round 406 without rerunning or rewriting prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round405_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/"round406_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_406_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 405
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"variable_completion_time_audit_results.json")
    assert saved["checks"] == {"run": 7, "failures": 0, "errors": 0}
    for flag in ("variable_completion_time_included",
                 "fixed_finite_processor_single_run_exact_universality_excluded",
                 "countably_many_finite_hardware_types_extension_proved"):
        assert saved[flag], flag
    for flag in ("continuous_program_state_coordinates_counted_as_external_settings",
                 "arbitrary_accuracy_at_unbounded_times_excluded",
                 "dense_reorientation_contract_excluded",
                 "finite_group_position_certificate_excluded",
                 "all_internal_resource_models_excluded",
                 "control_parameter_dimension_is_spatial_dimension",
                 "spatial_dimension_generated", "full_cognition_to_gr_refuted"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_406.md")["display_formulas"]
    assert formulas == 13
    links = 0
    for path in (HERE/"research_note_406.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path, link)
            links += 1
    for name in ("variable_completion_time_audit.py", "verify_variable_completion_round.py", Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_405") == 526
    assert result.pop("unchanged_prior_science_hashes_231_404") == 523
    assert result["stage_saved_tests"] == 1820
    assert result["total_protected_evidence_hashes"] == 549
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                       (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 406 and int(latest[2]) >= 1827
    result.update(
        rounds=[406], execution_mode="analytic dimension of exact variable-time implementation",
        scientific_base_through_round_by_round={406: 405},
        fresh_tests_by_round={406: 7}, fresh_tests=7, stage_saved_tests=1827,
        science_hashes_verified_231_406=529, unchanged_prior_science_hashes_231_405=526,
        total_protected_evidence_hashes=552, unchanged_prior_evidence_hashes=549,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        round345_variable_completion_time_gap_closed=True,
        finite_or_dense_position_protocols_invalidated=False,
        control_parameter_dimension_identified_as_space=False,
        old_scientific_experiments_rerun=False,
        actual_subject_position_generated=False, spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False,
        independent_final_agent_review_completed=False, all_reported_checks_passed=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-checks", action="store_true")
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

