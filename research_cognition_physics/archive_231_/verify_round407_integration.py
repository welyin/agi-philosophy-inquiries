"""Integrate round 407; preserve all earlier scientific evidence and snapshots."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round406_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/"round407_integration_checks.json"


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/"research_round_407_checks.json")
    assert checked["all_reported_checks_passed"]
    assert checked["scientific_base_through_round"] == 406
    assert checked["batch_scientific_dependencies"] == []
    assert len(checked["new_file_hashes"]) == 3
    for name, sha in checked["new_file_hashes"].items():
        assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/"spectral_tps_audit_results.json")
    assert saved["checks"] == {"run": 7, "failures": 0, "errors": 0}
    rows = saved["exact_certificates"]
    assert [(r["qubits"], r["moment_jacobian_certificate"]["rank"],
             r["gauge_certificate"]["rank"], r["hankel_certificate"]["rank"],
             r["infinitesimal_quotient_dimension"]) for r in rows] == [
                 (5, 32, 15, 32, 5), (6, 46, 18, 64, 0)]
    for row in rows:
        for key in ("moment_jacobian_certificate", "gauge_certificate", "hankel_certificate"):
            cert = row[key]
            assert len(cert["rows"]) == len(cert["columns"]) == cert["rank"]
            assert 0 < cert["minor_determinant_mod_prime"] < 1009
    for flag in ("exact_integer_certificates", "six_qubit_local_rigidity_mod_symmetry",
                 "five_qubit_nontrivial_local_isospectral_family",
                 "generic_claim_only_within_specified_locality_class"):
        assert saved[flag], flag
    for flag in ("global_unique_tps_proved", "all_graph_classes_excluded",
                 "actual_access_algebras_identified", "spatial_dimension_generated",
                 "full_cognition_to_gr_refuted", "infinite_resource_choice_required"):
        assert not saved[flag], flag
    formulas = core.text_checks(HERE/"research_note_407.md")["display_formulas"]
    assert formulas == 12
    links = 0
    for path in (HERE/"research_note_407.md", HERE/"spatial_premise_closure_audit.md"):
        for link in core.link_parser()(path.read_text(encoding="utf-8")):
            dest = (path.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (path, link)
            links += 1
    for name in ("spectral_tps_audit.py", "verify_spectral_tps_round.py", Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding="utf-8"))
    assert result.pop("science_hashes_verified_231_406") == 529
    assert result.pop("unchanged_prior_science_hashes_231_405") == 526
    assert result["stage_saved_tests"] == 1827
    assert result["total_protected_evidence_hashes"] == 552
    latest = re.search(r"最新科学轮次与检查数为(\d+)／(\d+)",
                      (HERE.parent/"research_direction.md").read_text(encoding="utf-8"))
    assert latest and int(latest[1]) >= 407 and int(latest[2]) >= 1834
    result.update(
        rounds=[407], execution_mode="exact certificates for conditional local spectral identification",
        scientific_base_through_round_by_round={407: 406},
        fresh_tests_by_round={407: 7}, fresh_tests=7, stage_saved_tests=1834,
        science_hashes_verified_231_407=532, unchanged_prior_science_hashes_231_406=529,
        total_protected_evidence_hashes=555, unchanged_prior_evidence_hashes=552,
        new_scientific_display_formulas_checked=formulas,
        local_links_checked=result["local_links_checked"]+links,
        navigation_and_notes_checked=result["navigation_and_notes_checked"]+2,
        exact_finite_field_certificates_checked=True,
        spectral_identification_distinguished_from_actual_access=True,
        global_unique_tps_proved=False, actual_access_algebras_identified=False,
        actual_subject_position_generated=False, spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False, infinite_resource_choice_required=False,
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

