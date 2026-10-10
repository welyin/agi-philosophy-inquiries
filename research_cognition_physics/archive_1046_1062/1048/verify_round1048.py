"""Read-only scientific replay, historical hashes and local evidence links."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import material_matching as science

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / "research_round_1048_checks.json"
OWN = [
    "../research_note_1048.md", "proof.md", "material_matching.py",
    "material_matching_results.json", "NEXT.md", "verify_round1048.py",
]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify():
    result = science.run()
    saved = json.loads(science.RESULT.read_text(encoding="utf-8"))
    science.compare(result, saved)
    history = saved["historical_sha256"]
    for rel, value in history.items():
        assert digest(ROOT / rel) == value, rel
    link_count = 0
    for rel in OWN:
        path = HERE / rel
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if target.startswith(("https://", "http://", "#")):
                continue
            dest = target.split("#", 1)[0].strip("<>")
            assert (path.parent / dest).resolve().exists(), (rel, dest)
            link_count += 1
    return {
        "round": 1048,
        "all_scientific_checks_passed": True,
        "new_scientific_groups": 1,
        "new_cognitive_axioms": 0,
        "owned_sha256": {rel: digest(HERE / rel) for rel in OWN},
        "historical_sha256": history,
        "local_links_checked": link_count,
        "full_roadmap_completed": False,
        "scope": saved["scope"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    receipt = verify()
    if args.freeze:
        with RECEIPT.open("x", encoding="utf-8") as f:
            json.dump(receipt, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        assert receipt == json.loads(RECEIPT.read_text(encoding="utf-8"))
    print(json.dumps({
        "round": 1048, "all_passed": True,
        "mode": "exclusive_freeze" if args.freeze else "read_only",
        "owned": len(receipt["owned_sha256"]),
        "history": len(receipt["historical_sha256"]),
        "local_links": receipt["local_links_checked"],
    }))


if __name__ == "__main__":
    main()
