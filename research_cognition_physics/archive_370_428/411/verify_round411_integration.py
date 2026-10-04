"""Integrate round 411; historical science is hashed, not rerun."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round410_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / "round411_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_411_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 410
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"source_interaction_algebra_audit_results.json")
    assert saved["checks"] == dict(run=6,failures=0,errors=0)
    assert [r["analytic_algebra_dimension"] for r in saved["algebra_cases"]] == [4,16,10,16]
    assert [r["traceless_minor_determinant_mod_prime"] for r in saved["controllability"]] == [927,6]
    assert saved["source_factor_and_generator_still_inputs"]
    for flag in ("generated_algebra_equals_finite_cost_control_set",
                 "source_only_statistics_identify_global_hamiltonian",
                 "environmental_tensor_sites_given_in_general_theorem",
                 "unitary_covariance_counted_as_physical_nonuniqueness",
                 "all_time_silent_region_inferred_from_first_order",
                 "hilbert_dimension_identified_as_spatial_dimension",
                 "spatial_dimension_generated","full_cognition_to_gr_refuted",
                 "new_cognitive_axiom_adopted"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_411.md")["display_formulas"]
    assert formulas == 10
    links = 0
    for path in (HERE/"research_note_411.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path,link)
            links += 1
    for name in ("source_interaction_algebra_audit.py","verify_source_interaction_round.py",Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_410") == 541
    assert result.pop("unchanged_prior_science_hashes_231_409") == 538
    assert result["stage_saved_tests"] == 1854
    assert result["total_protected_evidence_hashes"] == 564
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                       (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 411 and int(latest[2]) >= 1860
    result.update(rounds=[411],execution_mode="minimal direct source interaction algebra",
                  scientific_base_through_round_by_round={411:410},
                  fresh_tests_by_round={411:6},fresh_tests=6,stage_saved_tests=1860,
                  science_hashes_verified_231_411=544,unchanged_prior_science_hashes_231_410=541,
                  total_protected_evidence_hashes=567,unchanged_prior_evidence_hashes=564,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result["local_links_checked"]+links,
                  navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
                  canonical_minimal_interaction_algebra_from_one_source=True,
                  exact_interaction_algebra_not_stable_spatial_neighbor=True,
                  source_and_generator_still_input=True,
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
