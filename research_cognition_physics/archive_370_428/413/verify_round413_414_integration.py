"""Integrate independent rounds 413 and 414; preserve all historical evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round412_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / "round413_414_integration_checks.json"
CONFIG = {413: ("direction_only_coordinate_audit", 8, 9),
          414: ("word_budget_limit_audit", 7, 12)}


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    saved_by_round = {}
    links = 0
    for number, (stem, count, formulas) in CONFIG.items():
        checked = core.read(HERE / f"research_round_{number}_checks.json")
        assert checked["all_reported_checks_passed"]
        assert checked["parallel_batch"] == [413, 414]
        assert checked["scientific_base_through_round"] == 412
        assert checked["batch_scientific_dependencies"] == []
        assert checked["fresh_tests"] == dict(run=count, failures=0, errors=0)
        assert len(checked["new_file_hashes"]) == 3
        for name, sha in checked["new_file_hashes"].items():
            assert core.digest(HERE / name) == sha, name
        saved = core.read(HERE / (stem + "_results.json"))
        assert saved["checks"] == checked["fresh_tests"]
        saved_by_round[number] = saved
        note = HERE / f"research_note_{number}.md"
        assert core.text_checks(note)["display_formulas"] == formulas
        for link in core.link_parser()(note.read_text(encoding="utf-8")):
            assert (HERE / link).resolve().exists(), (note.name, link)
            links += 1
        ast.parse((HERE / (stem + ".py")).read_text(encoding="utf-8"))
    a, b = saved_by_round[413], saved_by_round[414]
    assert a["coordinates"]["budget_queries"] == 0
    assert a["full_displacement_effect_access_is_additional_input"]
    for flag in ("budget_structure_removed", "global_nonzero_contrast_inferred_from_unit_shell",
                 "unknown_single_state_tomography_or_cloning_claimed", "position_topology_generated_from_cognitive_axioms",
                 "exact_group_contract_certified_by_finite_tests", "new_cognitive_axiom_adopted", "phase_closure_triggered"):
        assert not a[flag], flag
    assert b["fixed_finite_translation_alphabet_is_model_input"] and b["endpoint_group_identity_is_model_input"]
    assert b["limit_symmetry"]["linear_isometry_count"] == 48
    for flag in ("polynomial_growth_is_cognitive_axiom", "schreier_quotient_automatically_a_group",
                 "arbitrary_finite_control_alphabet_ruled_out", "random_walk_limit_identified_with_word_budget",
                 "finite_symmetry_count_selects_spatial_dimension", "three_dimensional_topology_ruled_out",
                 "cognition_to_gr_refuted", "new_cognitive_axiom_adopted"):
        assert not b[flag], flag
    for link in core.link_parser()((HERE / "spatial_premise_closure_audit.md").read_text(encoding="utf-8")):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    for name in ("verify_coordinate_source_rounds.py", Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_412") == 547
    assert result.pop("unchanged_prior_science_hashes_231_411") == 544
    assert result["stage_saved_tests"] == 1868
    assert result["total_protected_evidence_hashes"] == 572
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)", (HERE.parent / "research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 414 and int(latest[2]) >= 1883
    result.update(
        rounds=[413, 414], execution_mode="two independent coordinate-source rounds from frozen 412",
        scientific_base_through_round_by_round={413: 412, 414: 412},
        fresh_tests_by_round={413: 8, 414: 7}, fresh_tests=15, stage_saved_tests=1883,
        science_hashes_verified_231_414=553, unchanged_prior_science_hashes_231_412=547,
        total_protected_evidence_hashes=578, unchanged_prior_evidence_hashes=572,
        new_scientific_display_formulas_checked=21,
        local_links_checked=result["local_links_checked"] + links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"] + 3,
        numerical_budget_oracle_replaced_by_explicit_effect_access=True,
        fixed_translation_word_budget_has_polyhedral_limit=True,
        arbitrary_finite_control_alphabets_ruled_out=False,
        actual_endpoint_contract_still_input=True, actual_subject_position_generated=False,
        spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, imaging_branch_deferred=True,
        old_scientific_experiments_rerun=False, independent_final_agent_review_completed=True,
        all_reported_checks_passed=True,
    )
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
