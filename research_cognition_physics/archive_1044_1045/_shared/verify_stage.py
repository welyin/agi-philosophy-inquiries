"""Read-only evidence and link checks for the 1044--1045 phase.

This does not prove mathematics or replace the per-round scientific verifiers.
--record exclusively creates a receipt after actual goal completion.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / "stage_evidence_checks.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check():
    first = json.loads((PHASE / "1044/research_round_1044_checks.json").read_text(encoding="utf-8"))
    second = json.loads((PHASE / "1045/research_round_1045_checks.json").read_text(encoding="utf-8"))
    accepted = json.loads((PHASE / "1045/mainline_acceptance.json").read_text(encoding="utf-8"))
    verified = []
    for base, table in (
        (ROOT, first["scientific_assets"]),
        (ROOT, first["history_inputs"]),
        (PHASE / "1045", second["owned_sha256"]),
        (ROOT, second["historical_sha256"]),
    ):
        for path, expected in table.items():
            full = (base / path).resolve()
            assert sha(full) == expected, str(full)
            verified.append(str(full.relative_to(ROOT.resolve())))
    for key in ("author_receipt", "independent_review"):
        item = accepted[key]
        assert sha(PHASE / "1045" / item["path"]) == item["sha256"]
    assert accepted["status"] == "accepted"
    assert accepted["cumulative_scientific_groups"] == 3821
    assert first["cumulative_after"] == accepted["previous_cumulative_scientific_groups"]
    assert accepted["new_cognitive_axioms"] == first["new_cognitive_axioms"] == 0
    paths = [
        HERE / "阶段综合与输入账.md", HERE / "goal_completion_audit.md",
        HERE / "independent_completion_review.md",
        PHASE / "README.md", PHASE / "文件索引.md", ROOT / "ROADMAP.md",
        ROOT / "README.md", ROOT / "RESEARCH_STATE.md", ROOT / "research_direction.md",
    ]
    links = 0
    for path in paths:
        content = path.read_text(encoding="utf-8")
        # Root historical snapshots remain as written, outside current navigation.
        if path == ROOT / "README.md":
            content = content.split("## 1042：", 1)[0]
        elif path == ROOT / "RESEARCH_STATE.md":
            content = content.split("## 1042时点状态：", 1)[0]
        elif path == ROOT / "research_direction.md":
            content = content.split("## 本阶段原目标正文：", 1)[0]
        content = re.sub(r"\x60{3}[\s\S]*?\x60{3}|\x60[^\x60\n]*\x60", "", content)
        for target in re.findall(r"(?<![\w\\])\[[^\]\n]*\]\(([^)\n]+)\)", content):
            target = target.strip().strip("<>")
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith("#"):
                continue
            full = (path.parent / unquote(target.split("#", 1)[0])).resolve()
            assert full.exists(), (str(path), target)
            links += 1
    completion = HERE / "goal_completion.json"
    completed = completion.exists()
    if completed:
        saved = json.loads(completion.read_text(encoding="utf-8"))
        assert saved["goal"]["status"] == "complete"
        start = json.loads((PHASE / "1044/goal_start.json").read_text(encoding="utf-8"))
        assert saved["goal"]["objective"] == start["application_result"]["goal"]["objective"]
    evidence = [
        PHASE / "1044/goal_start.json", PHASE / "1044/research_round_1044_checks.json",
        PHASE / "1044/verification_entry_correction.md",
        PHASE / "1045/research_round_1045_checks.json", PHASE / "1045/mainline_acceptance.json",
        PHASE / "1045/independent_review.md",
        HERE / "阶段综合与输入账.md", HERE / "goal_completion_audit.md",
        HERE / "independent_completion_review.md", Path(__file__),
    ] + ([completion] if completed else [])
    return {
        "schema": "phase1044_1045_evidence_v1", "all_checks_passed": True,
        "frozen_asset_references_checked": len(verified),
        "unique_frozen_assets_checked": len(set(verified)),
        "additional_acceptance_hashes_checked": 2, "current_local_links_checked": links,
        "latest_completed_round": 1045, "cumulative_scientific_groups": 3821,
        "new_cognitive_axioms": 0, "goal_completion_receipt_verified": completed,
        "evidence_sha256": {str(p.relative_to(ROOT)): sha(p) for p in evidence},
        "scope": "Evidence integrity and current links only; per-round mathematics and science are separately reviewed.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    result = check()
    if args.record:
        assert result["goal_completion_receipt_verified"]
        with RECEIPT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif RECEIPT.exists():
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        assert result["evidence_sha256"] == saved["evidence_sha256"]
    print(json.dumps({k: v for k, v in result.items() if k != "evidence_sha256"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
