"""Round-local delivery check; no navigation or neighboring round is frozen."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import multi_higgs_charge_transport as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1013.md"
OUT = HERE / "research_round_1013_checks.json"
sys.path.insert(0, str(ROOT / "scripts"))
from organize_research_231_775 import links

OWN = ["multi_higgs_charge_transport.py", "multi_higgs_charge_transport_results.json",
       "verify_round1013.py", "review.md", "NEXT.md"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(prospective=False, draft=False):
    fresh = science.run()
    saved = json.loads((HERE / "multi_higgs_charge_transport_results.json").read_text("utf8"))
    assert fresh == saved
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name) == digest, name
    assets = [NOTE] + [HERE/name for name in OWN]
    link_count = 0
    for p in assets:
        assert p.is_file(), p
        if p.suffix == ".md":
            for _, _, target, local in links(p.read_text("utf-8-sig")):
                resolved = (p.parent/local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (p, target)
                link_count += 1
    note = NOTE.read_text("utf8")
    for n in range(1,11):
        assert f"## {n}." in note
    for marker in ["相反超荷、相反方向可以兼容", "不要求每个Yukawa",
                   "封闭三维轻上分量质量块", "纯Y", "严格正定", "零VEV旁观者",
                   "没有生成质量机制或真空", "Higgs基只"]:
        # Last marker belongs to the scope review rather than the note.
        assert marker in note + (HERE/"review.md").read_text("utf8"), marker
    reviewed = "主代理审阅通过" in (HERE/"review.md").read_text("utf8")
    assert draft or reviewed, "Await primary-agent review before freezing the receipt."
    assert all(value is False for value in fresh["claim_boundaries"].values())
    assert fresh["all_scientific_calibrations_passed"] is True
    return dict(round=1013,date="2026-10-08",all_delivery_checks_passed=True,
        scientific_result_reproduced=True,new_calibration_groups=1,cumulative_research_groups=3791,
        local_links_checked=link_count,frozen_current_files=len(assets),
        historical_input_files=len(fresh["historical_source_sha256"]),
        live_navigation_frozen=False,neighboring_round_frozen=False,
        gram_cases=7,rotated_vacua=16,exact_aligned_determinants=4,
        independent_primary_agent_review=reviewed,
        new_cognitive_axiom=False,goal_completed=False,visual_checks_performed=False,
        analytic_scope="Positive kinetic Gram kernel; arbitrary doublet VEV orientation with mixed unbroken generator; conjugate-doublet dictionary; total Yukawa and stipulated Majorana mass transport.",
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"):sha(p) for p in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--draft", action="store_true")
    args = parser.parse_args()
    result = verify(prospective=args.write or args.draft, draft=args.draft)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    elif not args.draft:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},ensure_ascii=False,indent=2))
