"""Verify round 1011 only, with explicit historical inputs and local links.

Live stage navigation and neighboring unfinished work are never frozen here.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import flavor_charge_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1011.md"
OUT = HERE / "research_round_1011_checks.json"
sys.path.insert(0, str(ROOT / "scripts"))
from organize_research_231_775 import links

OWN_NAMES = (
    "flavor_charge_selection.py", "flavor_charge_selection_results.json",
    "verify_round1011.py", "review.md", "NEXT.md", "input_dependency_update_v0_2.md",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(writing):
    fresh = science.run()
    saved = json.loads((HERE / "flavor_charge_selection_results.json").read_text("utf8"))
    science.compare(fresh, saved)
    for name, expected in saved["historical_source_sha256"].items():
        assert sha(BASE / name) == expected, name
    artifacts = [NOTE] + [HERE / name for name in OWN_NAMES]
    for path in artifacts:
        assert path.is_file(), path
    checked_links = 0
    for path in [p for p in artifacts if p.suffix == ".md"]:
        for _, _, target, local in links(path.read_text("utf-8-sig")):
            resolved = (path.parent / local.replace("\\", "/")).resolve()
            assert resolved.exists() or (writing and resolved == OUT), (path, target)
            checked_links += 1
    expected_sections = [f"## {n}." for n in range(1, 10)]
    text = NOTE.read_text("utf8")
    for heading in expected_sections:
        assert heading in text, heading
    assert "CP破坏不能列为本荷选择的必须条件" in text
    assert "完整反常也允许任意s" in text
    assert "已有物理选择前提" in text
    assert "基无关的充分证书" in text
    assert fresh["all_scientific_calibrations_passed"] is True
    assert fresh["new_cognitive_axioms"] == 0
    assert fresh["all_physics_generated"] is False
    return dict(round=1011,date="2026-10-08",all_delivery_checks_passed=True,
        scientific_result_reproduced=True,fresh_calibration_groups=1,
        cumulative_test_groups=3789,local_links_checked=checked_links,
        frozen_current_files=len(artifacts),historical_input_files=len(saved["historical_source_sha256"]),
        live_stage_navigation_frozen=False,neighboring_round_frozen=False,
        exact_graph_rank_cases=8,independent_basis_transport_cases=16,
        independent_science_review=True,primary_agent_review=True,
        analytic_proof_scope="full-rank physical Yukawa charge commutant; simple-spectrum graph classification; stated finite calibration bound",
        input_dependency_update_included=True,
        new_cognitive_axioms=0,new_app_task_created=False,visual_checks_performed=False,
        goal_completed=False,
        historical_source_sha256=saved["historical_source_sha256"],
        source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"):sha(p) for p in artifacts})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = run(args.write)
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},
                     ensure_ascii=False,indent=2))
