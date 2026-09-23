"""Integrate round 405 without rerunning or rewriting prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round404_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/"round405_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_405_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 404
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"scattering_dictionary_audit_results.json")
    assert saved["checks"] == {"run": 7, "failures": 0, "errors": 0}
    for flag in ("exact_fixed_global_unitary_dictionary_proved",
                 "normal_program_future_bound_uniform_over_future_waits",
                 "passive_moving_readouts_covered"):
        assert saved[flag], flag
    for flag in ("dictionary_operator_norm_quasilocal",
                 "newly_inserted_active_interventions_covered",
                 "global_information_erased", "convergence_uniform_over_all_programs",
                 "actual_subject_position_generated", "spatial_dimension_generated",
                 "full_cognition_to_gr_refuted"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_405.md")["display_formulas"]
    assert formulas == 12
    links = 0
    for path in (HERE/"research_note_405.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path, link)
            links += 1
    for name in ("scattering_dictionary_audit.py", "verify_scattering_dictionary_round.py", Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_404") == 523
    assert result.pop("unchanged_prior_science_hashes_231_403") == 520
    assert result["stage_saved_tests"] == 1813
    assert result["total_protected_evidence_hashes"] == 546
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                       (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 405 and int(latest[2]) >= 1820
    result.update(
        rounds=[405], execution_mode="fixed global scattering dictionary versus physical local access",
        scientific_base_through_round_by_round={405: 404},
        fresh_tests_by_round={405: 7}, fresh_tests=7, stage_saved_tests=1820,
        science_hashes_verified_231_405=526, unchanged_prior_science_hashes_231_404=523,
        total_protected_evidence_hashes=549, unchanged_prior_evidence_hashes=546,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        exact_global_dictionary_distinguished_from_physical_access=True,
        moving_passive_readout_future_bound_proved=True,
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

