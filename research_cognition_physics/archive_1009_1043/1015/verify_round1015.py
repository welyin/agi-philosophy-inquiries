"""Round 1015 delivery and numerical reproduction; no navigation is frozen."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import vector_consistency_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1015.md"
OUT = HERE / "research_round_1015_checks.json"
sys.path.insert(0, str(ROOT / "scripts"))
from organize_research_231_775 import links

OWN = ["vector_consistency_selection.py", "vector_consistency_selection_results.json",
       "verify_round1015.py", "review.md", "NEXT.md", "input_dependency_update_v0_4.md"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(prospective=False):
    fresh = science.run()
    assert fresh == json.loads((HERE / OWN[1]).read_text("utf8"))
    assert fresh["all_scientific_calibrations_passed"]
    assert fresh["cumulative_research_groups"] == 3793
    assert fresh["cumulative_count_pending_primary_confirmation"] is False
    assert all(v is False for v in fresh["claim_boundaries"].values())
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE / name) == digest, name
    assert fresh["code_sha256"] == sha(HERE / OWN[0])
    assets = [NOTE] + [HERE / name for name in OWN]
    count = 0
    for p in assets:
        assert p.is_file(), p
        if p.suffix == ".md":
            text = p.read_text("utf-8-sig")
            assert text.count("$$") % 2 == 0, p
            for _, _, target, local in links(text):
                resolved = (p.parent / local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (p, target)
                count += 1
    note = NOTE.read_text("utf8")
    for n in range(1, 11):
        assert f"## {n}." in note
    assert "主代理审阅通过" in (HERE / "review.md").read_text("utf8")
    assert "独立终审通过" in (HERE / "review.md").read_text("utf8")
    return dict(round=1015, date="2026-10-08", all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3793, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(fresh["historical_source_sha256"]),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        first_jet_calibrations=fresh["first_jet_count"],
        classification_completeness_adopted_from_primary_theorem=True,
        finite_samples_not_general_classification_proof=True,
        local_total_derivative_not_global_physical_equivalence=True,
        independent_primary_math_and_code_review=True,
        independent_scope_review=True, new_cognitive_axiom=False,
        goal_completed=False, visual_checks_performed=False,
        analytic_scope="Canonical positive free Maxwell sector; four-dimensional local Poincare-invariant leading dimension-four regular deformation preserving gauge freedom; Yang-Mills class modulo declared local equivalence.",
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"):sha(p) for p in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, ensure_ascii=False, indent=2))
