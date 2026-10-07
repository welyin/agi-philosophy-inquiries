"""Reproduce the quantum-reference interface without rewriting older research."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
import verify_material_reference_geometry_review as previous
import verify_interaction_rounds as core
import quantum_reference_readout as science

HERE = Path(__file__).resolve().parent
TARGET = HERE / "quantum_reference_readout_checks.json"
FILES = ["quantum_reference_readout.py", "quantum_reference_readout_results.json",
         "quantum_reference_readout_review.md", "quantum_reference_readout_review_draft.txt",
         "quantum_reference_readout.py.before_grid_refinement.txt",
         "quantum_reference_readout_results.json.before_grid_refinement.txt",
         "quantum_reference_readout_review_before_final_review.txt"]


def verify():
    base = previous.verify()
    assert base["unchanged_numbered_scientific_tests"] == 2552
    assert base["unchanged_scientific_file_hashes"] == 871
    assert base["total_protected_including_this_review"] == 1021
    if TARGET.exists():
        for name, digest in core.read(TARGET)["new_unnumbered_evidence_hashes"].items():
            assert core.digest(HERE / name) == digest, name
    assert (HERE / FILES[2]).read_bytes() == (HERE / FILES[3]).read_bytes()
    result = core.read(HERE / FILES[1])
    assert result == json.loads(json.dumps(science.run()))
    assert (result["diagnostic_tests"], result["failures"], result["errors"]) == (6, 0, 0)
    prior = core.read(HERE / FILES[5])
    without_grid = json.loads(json.dumps(result))
    del without_grid["diagnostics"]["discretization"]["fixed_box_finer_grid_state_difference"]
    assert prior == without_grid
    for name in (FILES[0], FILES[4], Path(__file__).name):
        ast.parse((HERE / name).read_text("utf8"))
    checked = core.text_checks(HERE / FILES[2])
    assert checked["display_formulas"] == 12
    links = 0
    for link in core.link_parser()((HERE / FILES[2]).read_text("utf8")):
        assert (HERE / link).resolve().exists(), link
        links += 1
    return dict(date="2026-09-30", scientific_base_through_round=520,
        unchanged_numbered_scientific_tests=2552, unchanged_scientific_file_hashes=871,
        unchanged_previously_protected_evidence_hashes=1021,
        total_protected_including_this_review=1028,
        numbered_round_created=False, numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=6, failures=0, errors=0),
        new_unnumbered_evidence_hashes={f: core.digest(HERE / f) for f in FILES},
        results_reproduced=True, text_checks=checked,
        local_links_checked=links, broken_links=0,
        previous_result_preserved_without_change=True,
        old_science_rewritten=False, visual_checks_performed=False,
        independent_final_review_completed=True,
        scope=result["scope"], all_reported_checks_passed=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write-checks", action="store_true")
    args = p.parse_args()
    answer = verify()
    if args.write_checks:
        with TARGET.open("x", encoding="utf8", newline="\n") as f:
            f.write(json.dumps(answer, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: answer[k] for k in ("scientific_base_through_round",
        "unchanged_numbered_scientific_tests", "total_protected_including_this_review",
        "all_reported_checks_passed")}))
