"""Integrate round 410; historical science is hashed, not rerun."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round409_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / "round410_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_410_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 409
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"shared_record_topology_audit_results.json")
    assert saved["checks"] == dict(run=6,failures=0,errors=0)
    assert saved["exact_dyadic_fibers"]["exact_pairs"] == 63
    assert saved["exact_dyadic_fibers"]["sharp_projection_distance_from_continuous_coordinate_algebra"] == .5
    for flag in ("finite_experiments_prove_infinite_topology","all_quantum_states_assumed_classical",
                 "complete_record_topology_identified_as_physical_space",
                 "mathematical_coordinate_choices_are_distinct_physical_worlds",
                 "position_quotient_erases_actual_history","actual_subject_movement_implemented",
                 "spatial_dimension_generated","full_cognition_to_gr_refuted",
                 "initial_infinite_resource_adopted","new_cognitive_axiom_adopted"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_410.md")["display_formulas"]
    assert formulas == 12
    links = 0
    for path in (HERE/"research_note_410.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path,link)
            links += 1
    for name in ("shared_record_topology_audit.py","verify_record_topology_round.py",Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_409") == 538
    assert result.pop("unchanged_prior_science_hashes_231_408") == 535
    assert result["stage_saved_tests"] == 1848
    assert result["total_protected_evidence_hashes"] == 561
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                       (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 410 and int(latest[2]) >= 1854
    result.update(rounds=[410],execution_mode="current shared records and continuous position topology",
                  scientific_base_through_round_by_round={410:409},
                  fresh_tests_by_round={410:6},fresh_tests=6,stage_saved_tests=1854,
                  science_hashes_verified_231_410=541,unchanged_prior_science_hashes_231_409=538,
                  total_protected_evidence_hashes=564,unchanged_prior_evidence_hashes=561,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result["local_links_checked"]+links,
                  navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
                  full_sharp_record_position_candidate_excluded_under_explicit_topology=True,
                  continuous_position_quotient_selection_still_input=True,
                  finite_invertibility_not_extrapolated_to_limit=True,
                  actual_subject_position_generated=False,spatial_dimension_generated=False,
                  full_cognition_to_gr_refuted=False,infinite_resource_choice_required=False,
                  old_scientific_experiments_rerun=False,
                  independent_final_agent_review_completed=False,all_reported_checks_passed=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-checks",action="store_true")
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
