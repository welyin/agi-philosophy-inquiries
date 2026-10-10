"""Read-only replay and fixed-scope evidence verification for round 1050."""
from pathlib import Path
import argparse
import json
import re
import coulomb_natural_instrument as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RECEIPT = HERE / "research_round_1050_checks.json"
OWN = ["../research_note_1050.md", "proof.md", "coulomb_natural_instrument.py",
       "coulomb_natural_instrument_results.json", "dependency_update.md",
       "review.md", "NEXT.md", "sources.md", "verify_round1050.py"]


def verify(record=False):
    current = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(current, saved)
    for rel, sha in saved["historical_sha256"].items():
        assert science.digest(ROOT / rel) == sha, rel
    links = 0
    for rel in OWN:
        path = HERE / rel
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            target = target.strip().strip("<>")
            if target.startswith(("https://", "http://", "#")):
                continue
            dest = (path.parent / target.split("#", 1)[0]).resolve()
            if record and dest == RECEIPT:
                pass
            else:
                assert dest.exists(), (rel, target)
            links += 1
    return {
        "round": 1050, "date": "2026-10-08", "all_passed": True,
        "new_scientific_groups": 1, "new_cognitive_axioms": 0,
        "full_roadmap_completed": False,
        "owned_sha256": {rel: science.digest(HERE / rel) for rel in OWN},
        "historical_sha256": saved["historical_sha256"],
        "local_links_checked": links,
        "scientific_scope": saved["scope"],
        "finite_window_error_exact": saved["finite_window"]["instrument_error_strict_upper_bound_exact"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    result = verify(args.freeze)
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert result == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({"round": 1050, "all_passed": True,
                      "mode": "exclusive_freeze" if args.freeze else "read_only",
                      "owned": len(OWN), "history": len(result["historical_sha256"]),
                      "local_links": result["local_links_checked"]}))


if __name__ == "__main__":
    main()
