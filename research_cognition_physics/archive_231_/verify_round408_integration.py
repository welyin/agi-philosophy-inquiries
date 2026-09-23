"""Integrate round 408, preserving previous scientific evidence and snapshots."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round407_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/"round408_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_408_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 407
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"calibrated_effect_anchor_audit_results.json")
    assert saved["checks"] == {"run": 8, "failures": 0, "errors": 0}
    for name, sha in saved["frozen_source_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    cert = saved["exact_commutant_certificate"]
    assert cert["constraint_shape"] == [4096, 63] and cert["rank"] == 63
    assert cert["prime"] == 1009 and cert["minor_determinant_mod_prime"] == 26
    assert len(cert["rows"]) == len(cert["columns"]) == 63
    for flag in ("scalar_joint_commutant_proved",
                 "original_h_commuting_unitary_freedom_removed_with_fixed_effect",
                 "robust_unknown_reference_bound_proved"):
        assert saved[flag], flag
    for flag in ("one_binary_effect_is_one_measurement_shot",
                 "anchor_preparation_derived", "actual_position_generated",
                 "global_unique_tps_proved", "full_unitary_controllability_proved",
                 "spatial_dimension_generated", "full_cognition_to_gr_refuted",
                 "infinite_resource_choice_required", "estimated_gap_is_certified_lower_bound"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_408.md")["display_formulas"]
    assert formulas == 12
    links = 0
    for path in (HERE/"research_note_408.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path, link)
            links += 1
    for name in ("calibrated_effect_anchor_audit.py", "verify_effect_anchor_round.py", Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_407") == 532
    assert result.pop("unchanged_prior_science_hashes_231_406") == 529
    assert result["stage_saved_tests"] == 1834
    assert result["total_protected_evidence_hashes"] == 555
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                      (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 408 and int(latest[2]) >= 1842
    result.update(
        rounds=[408], execution_mode="calibrated effect removes residual commuting unitary freedom",
        scientific_base_through_round_by_round={408: 407},
        fresh_tests_by_round={408: 8}, fresh_tests=8, stage_saved_tests=1842,
        science_hashes_verified_231_408=535, unchanged_prior_science_hashes_231_407=532,
        total_protected_evidence_hashes=558, unchanged_prior_evidence_hashes=555,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        exact_fixed_effect_commutant_certificate_checked=True,
        approximate_anchor_bound_preserves_arbitrary_reference=True,
        estimated_gap_promoted_to_certified_lower_bound=False,
        actual_anchor_preparation_derived=False, actual_subject_position_generated=False,
        spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
        infinite_resource_choice_required=False,
        old_scientific_experiments_rerun=False,
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

