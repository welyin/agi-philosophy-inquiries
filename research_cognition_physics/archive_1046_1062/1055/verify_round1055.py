"""Author science/history/link replay. Only --freeze creates the first receipt."""
from pathlib import Path
import argparse
import json
import re
import moving_optical_interface as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "research_round_1055_checks.json"
OWN = ["../research_note_1055.md", "proof.md", "moving_optical_interface.py",
       "results.json", "dependency_update.md", "NEXT.md", "review.md", "sources.md",
       "verify_round1055.py"]


def verify(freeze=False):
    current = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(current, saved)
    for rel, sha in saved["historical_sha256"].items():
        assert science.digest(ROOT/rel) == sha, rel
    count = 0
    for rel in OWN:
        path = HERE/rel
        assert path.is_file(), rel
        if path.suffix != ".md":
            continue
        # TeX such as [channel](rho) is not a Markdown link.
        body = path.read_text(encoding="utf-8")
        body = re.sub(r"\$\$[\s\S]*?\$\$", "", body)
        body = re.sub(r"(?<!\\)\$(?:\\.|[^$\n])*(?<!\\)\$", "", body)
        for target in re.findall(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)", body):
            target = target.strip().strip("<>")
            if target.startswith(("https://", "http://", "#")):
                continue
            destination = (path.parent/target.split("#", 1)[0]).resolve()
            if not (freeze and destination == RECEIPT):
                assert destination.is_file(), (rel, target)
            count += 1
    return {
        "round": 1055, "date": "2026-10-09", "author_checks_passed": True,
        "independent_final_review_completed_at_author_freeze": False,
        "new_project_conditional_physical_dictionary_groups": 1,
        "new_cognitive_axioms": 0, "whole_roadmap_completed_here": False,
        "owned_sha256": {rel: science.digest(HERE/rel) for rel in OWN},
        "historical_sha256": saved["historical_sha256"],
        "local_links_checked": count,
        "co_lipschitz_lower": saved["rational_certificate"]["exact"]["new_co_lipschitz_lower"],
        "direction_gap_lower": saved["rational_certificate"]["exact"]["direction_antipodal_gap_lower"],
        "scope": saved["scope"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    result = verify(args.freeze)
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        assert result == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({"round": 1055, "author_checks_passed": True,
                      "mode": "exclusive_freeze" if args.freeze else "read_only",
                      "owned": len(OWN), "history": len(result["historical_sha256"]),
                      "local_links": result["local_links_checked"]}))


if __name__ == "__main__":
    main()
