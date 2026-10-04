"""Verify independent rounds 413 and 414 against frozen round 412 evidence."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest

import verify_round412_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
CONFIG = {
    413: ("direction_only_coordinate_audit", 8, 9, [287, 288, 370, 371, 380, 385, 412]),
    414: ("word_budget_limit_audit", 7, 12, [285, 348, 378, 379, 384, 385, 412]),
}


def verify(number, pending=False):
    stem, count, formulas, dependencies = CONFIG[number]
    target = HERE / f"research_round_{number}_checks.json"
    old_target = previous.TARGET
    previous.TARGET = target
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen["science_hashes_verified_231_412"] == 547
    assert frozen["total_protected_evidence_hashes"] == 572
    names = [stem + ".py", stem + "_results.json", f"research_note_{number}.md"]
    if target.exists():
        for name, sha in core.read(target)["new_file_hashes"].items():
            assert core.digest(HERE / name) == sha, name
    source = HERE / names[0]
    ast.parse(source.read_text(encoding="utf-8"))
    spec = importlib.util.spec_from_file_location(stem + "_checked", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    saved = core.read(HERE / names[1])
    assert {k: v for k, v in saved.items() if k not in ("checks", "runtime")} == json.loads(json.dumps(module.report()))
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(module))
    assert tests.wasSuccessful() and tests.testsRun == count, output.getvalue()
    assert saved["checks"] == dict(run=count, failures=0, errors=0)
    checked = core.text_checks(HERE / names[2])
    assert checked["display_formulas"] == formulas
    links = 0
    for link in core.link_parser()((HERE / names[2]).read_text(encoding="utf-8")):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == target), link
        links += 1
    review = {
        413: "bearing_coordinate_review independently checked proofs, raw probability and position error bounds, blind-shell counterexample and eight checks; source version and numerical rejection thresholds clarified before freezing",
        414: "parent independently checked the full note and code: LP rounding proof, finite polytope symmetry, exact word length, all-point 4/m bound, and restriction to a fixed translation alphabet",
    }
    return dict(
        date="2026-09-24", round=number, parallel_batch=[413, 414],
        execution_mode="independent complete rounds from one frozen scientific baseline",
        scientific_base_through_round=412, additional_frozen_dependency_rounds=dependencies,
        batch_scientific_dependencies=[], fresh_tests=dict(run=count, failures=0, errors=0),
        saved_results_reproduced=True, scientific_results_rewritten=False,
        previous_scientific_file_hashes_verified=547, previous_protected_evidence_hashes_verified=572,
        text_checks=checked, local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in names},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review=review[number], actual_endpoint_contract_still_input=True,
        phase_closure_triggered=False, scope=saved["scope"], all_reported_checks_passed=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round", type=int, choices=CONFIG)
    parser.add_argument("--write-checks", action="store_true")
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        with (HERE / f"research_round_{args.round}_checks.json").open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
