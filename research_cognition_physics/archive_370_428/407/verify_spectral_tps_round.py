"""Read-only science verification for round 407, with the 406 evidence audit."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest
import verify_round406_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/"research_round_407_checks.json"


def verify(number=407, pending=False):
    assert number == 407
    # Verifies protected evidence and old links; does not execute old experiments.
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen["science_hashes_verified_231_406"] == 529
    assert frozen["total_protected_evidence_hashes"] == 552
    names = ["spectral_tps_audit.py", "spectral_tps_audit_results.json",
             "research_note_407.md"]
    if TARGET.exists():
        for name, sha in core.read(TARGET)["new_file_hashes"].items():
            assert core.digest(HERE/name) == sha, name
    source = HERE/names[0]
    ast.parse(source.read_text(encoding="utf-8"))
    spec = importlib.util.spec_from_file_location("spectral_tps_checked", source)
    science = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(science)
    saved = core.read(HERE/names[1])
    assert {k: v for k, v in saved.items() if k not in ("checks", "runtime")} == json.loads(json.dumps(science.report()))
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(science))
    assert tests.wasSuccessful() and tests.testsRun == 7, output.getvalue()
    assert saved["checks"] == {"run": 7, "failures": 0, "errors": 0}
    checked = core.text_checks(HERE/names[2])
    assert checked["display_formulas"] == 12
    links = 0
    for link in core.link_parser()((HERE/names[2]).read_text(encoding="utf-8")):
        dest = (HERE/link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(
        date="2026-09-23", round=407, parallel_batch=[407],
        execution_mode="exact modular certificates for local spectral rigidity",
        scientific_base_through_round=406, additional_frozen_dependency_rounds=[],
        batch_scientific_dependencies=[], fresh_tests=dict(run=7, failures=0, errors=0),
        saved_results_reproduced=True, scientific_results_rewritten=False,
        previous_scientific_file_hashes_verified=529,
        previous_protected_evidence_hashes_verified=552,
        text_checks=checked, local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in names},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        final_science_review="primary-agent source, exact-minor proof and executable review; no independent agent final review",
        scope="Specified five/six-qubit open-chain class; parameter-local spectral rigidity modulo local basis and Hamiltonian symmetries. No global TPS uniqueness or physical space generation.",
        all_reported_checks_passed=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round", type=int, choices=[407])
    parser.add_argument("--write-checks", action="store_true")
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
