"""Verify round 412; historical scientific evidence is hashed, not rerun."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest
import verify_round411_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/"research_round_412_checks.json"


def verify(number=412,pending=False):
    assert number == 412
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        frozen=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    assert frozen["science_hashes_verified_231_411"] == 544
    assert frozen["total_protected_evidence_hashes"] == 567
    draft_names=["round412_drafts/research_note_412_pre_budget_access.md",
                 "round412_drafts/research_round_412_checks_pre_budget_access.json"]
    draft_check=core.read(HERE/draft_names[1])
    assert core.digest(HERE/draft_names[0])==draft_check["new_file_hashes"]["research_note_412.md"]
    for name in ("operational_coordinate_audit.py","operational_coordinate_audit_results.json"):
        assert core.digest(HERE/name)==draft_check["new_file_hashes"][name]
    names=["operational_coordinate_audit.py","operational_coordinate_audit_results.json","research_note_412.md"]
    if TARGET.exists():
        for name,sha in core.read(TARGET)["new_file_hashes"].items():
            assert core.digest(HERE/name)==sha,name
    source=HERE/names[0]
    ast.parse(source.read_text(encoding="utf-8"))
    spec=importlib.util.spec_from_file_location("operational_coordinates_checked",source)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    saved=core.read(HERE/names[1])
    assert {k:v for k,v in saved.items() if k not in ("checks","runtime")} == module.report()
    output=io.StringIO()
    tests=unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(module))
    assert tests.wasSuccessful() and tests.testsRun==8,output.getvalue()
    assert saved["checks"] == dict(run=8,failures=0,errors=0)
    checked=core.text_checks(HERE/names[2])
    assert checked["display_formulas"]==12
    links=0
    for link in core.link_parser()((HERE/names[2]).read_text(encoding="utf-8")):
        dest=(HERE/link).resolve()
        assert dest.exists() or (pending and dest==TARGET),link
        links+=1
    return dict(date="2026-09-23",round=412,parallel_batch=[412],
        execution_mode="conditional operational coordinates from budget measurements",
        scientific_base_through_round=411,additional_frozen_dependency_rounds=[287,288,380,381,383,384,385],
        batch_scientific_dependencies=[],fresh_tests=dict(run=8,failures=0,errors=0),
        saved_results_reproduced=True,scientific_results_rewritten=False,
        previous_scientific_file_hashes_verified=544,previous_protected_evidence_hashes_verified=567,
        text_checks=checked,local_links_checked=links,broken_links=0,
        new_file_hashes={name:core.digest(HERE/name) for name in names},
        preserved_round412_draft_hashes={name:core.digest(HERE/name) for name in draft_names},
        budget_access_qualification_changed_scientific_values=False,
        visual_rendering_performed=False,legacy_science_tests_rerun=False,
        independent_review="transport_bridge_audit: compact-group extension, radial budget, four-query coordinates, chart transitions and code reviewed; rtol made explicit; raw measurement error scope retained",
        actual_endpoint_contract_still_input=True,phase_closure_triggered=False,
        scope=saved["scope"],all_reported_checks_passed=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round",type=int,choices=[412])
    parser.add_argument("--write-checks",action="store_true")
    args=parser.parse_args()
    result=verify(args.round,args.write_checks)
    if args.write_checks:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
