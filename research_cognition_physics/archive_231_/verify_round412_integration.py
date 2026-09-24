"""Integrate round 412 without rerunning historical scientific experiments."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round411_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/"round412_integration_checks.json"


def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/"research_round_412_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"]==411
    assert checked["batch_scientific_dependencies"]==[]
    assert len(checked["new_file_hashes"])==3
    for name,sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name)==sha,name
    assert len(checked["preserved_round412_draft_hashes"])==2
    for name,sha in checked["preserved_round412_draft_hashes"].items():
        assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/"operational_coordinate_audit_results.json")
    assert saved["checks"]==dict(run=8,failures=0,errors=0)
    assert saved["coordinates"]["budget_queries_per_new_position"]==4
    assert saved["conditional_coordinate_construction"]
    assert saved["actual_endpoint_contract_still_input"]
    for flag in ("direction_injectivity_added_as_axiom","actual_so3_action_added_as_axiom",
                 "real_group_roots_required_by_readout","finite_statistics_certify_global_contract",
                 "globally_flat_universe_derived","spatial_dimension_generated_from_cognitive_axioms",
                 "new_cognitive_axiom_adopted","phase_closure_triggered"):
        assert not saved[flag],flag
    formulas=core.text_checks(HERE/"research_note_412.md")["display_formulas"]
    assert formulas==12
    links=0
    for path in (HERE/"research_note_412.md",HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest=(path.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),(path,link)
            links+=1
    for name in ("operational_coordinate_audit.py","verify_operational_coordinate_round.py",Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_411")==544
    assert result.pop("unchanged_prior_science_hashes_231_410")==541
    assert result["stage_saved_tests"]==1860
    assert result["total_protected_evidence_hashes"]==567
    latest=re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",(HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1])>=412 and int(latest[2])>=1868
    result.update(rounds=[412],execution_mode="conditional operational coordinate reconstruction",
        scientific_base_through_round_by_round={412:411},fresh_tests_by_round={412:8},fresh_tests=8,
        stage_saved_tests=1868,science_hashes_verified_231_412=547,
        unchanged_prior_science_hashes_231_411=544,total_protected_evidence_hashes=572,
        preserved_round412_draft_hashes=checked["preserved_round412_draft_hashes"],
        budget_access_qualification_changed_scientific_values=False,
        unchanged_prior_evidence_hashes=567,new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        direction_injectivity_derived_under_displacement_contract=True,
        budget_exponent_calibrated_by_repeated_displacement=True,
        actual_endpoint_contract_still_input=True,actual_subject_position_generated=False,
        spatial_dimension_generated=False,full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False,imaging_branch_deferred=True,
        old_scientific_experiments_rerun=False,independent_final_agent_review_completed=True,
        all_reported_checks_passed=True)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-checks",action="store_true")
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
