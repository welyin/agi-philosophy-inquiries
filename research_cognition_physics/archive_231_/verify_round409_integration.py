"""Integrate round 409 and preserve all historical scientific snapshots."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round408_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/"round409_integration_checks.json"


def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/"research_round_409_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"]==408
    assert checked["batch_scientific_dependencies"]==[]
    assert len(checked["new_file_hashes"])==3
    for name,sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/"separable_program_regularity_audit_results.json")
    assert saved["checks"]==dict(run=6,failures=0,errors=0)
    assert saved["bounded_generator_exact_channel_hausdorff_dimension_upper"]==1
    assert saved["finite_mean_energy_exact_channel_hausdorff_dimension_upper"]==2
    for flag in ("uniform_program_energy_cap_required","finite_program_dimension_required",
                 "program_selection_continuity_required","approximation_excluded",
                 "spatial_dimension_generated","full_cognition_to_gr_refuted",
                 "strengthened_exact_autonomous_control_adopted_as_axiom",
                 "initial_resource_amount_choice_required",
                 "infinite_dimension_theorem_established_by_finite_numerics"):
        assert not saved[flag],flag
    formulas=core.text_checks(HERE/"research_note_409.md")["display_formulas"]
    assert formulas==14
    links=0
    for path in (HERE/"research_note_409.md",HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest=(path.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),(path,link)
            links+=1
    for name in ("separable_program_regularity_audit.py","verify_separable_program_round.py",Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_408")==535
    assert result.pop("unchanged_prior_science_hashes_231_407")==532
    assert result["stage_saved_tests"]==1842
    assert result["total_protected_evidence_hashes"]==558
    latest=re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                     (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1])>=409 and int(latest[2])>=1848
    result.update(
        rounds=[409],execution_mode="exact autonomous scope for separable normal programs",
        scientific_base_through_round_by_round={409:408},
        fresh_tests_by_round={409:6},fresh_tests=6,stage_saved_tests=1848,
        science_hashes_verified_231_409=538,unchanged_prior_science_hashes_231_408=535,
        total_protected_evidence_hashes=561,unchanged_prior_evidence_hashes=558,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        finite_program_assumption_removed_under_explicit_regularities=True,
        infinite_program_approximation_boundary_preserved=True,
        all_normal_unbounded_energy_programs_excluded=False,
        exact_universal_autonomous_control_adopted_as_cognitive_axiom=False,
        dimension_of_control_targets_identified_as_space=False,
        actual_subject_position_generated=False,spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False,infinite_resource_choice_required=False,
        old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=False,all_reported_checks_passed=True)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-checks",action="store_true")
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
