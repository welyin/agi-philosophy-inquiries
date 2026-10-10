"""Read-only verification of count-zero mature adoptions after round 1055.

--freeze creates the first receipt exclusively.  Navigation hashes describe the
snapshot at adoption, not a prohibition on future authorized navigation edits.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
RECEIPT = HERE / "adoptions_after1055_checks.json"
OWN = [
    "pion_anomaly_adoption.md", "pion_anomaly_adoption_review.md",
    "acoustic_transport_adoption.md", "acoustic_transport_adoption_review.md",
    "acoustic_transport_adoption/check.py", "acoustic_transport_adoption/results.json",
    "integration_adoptions_after1055.md", "verify_adoptions_after1055.py",
    "../_admission/physical_subsystems_after1055/selection.md",
    "../_admission/anomaly_matching_after1055/selection.md",
    "../_admission/acoustic_transport_after1055/selection.md",
]
HISTORY = [
    "可组合认知结构与复量子状态空间_阶段论文.md",
    "可组合认知结构与有限维量子理论_阶段论文.md",
    "archive_342_369/research_note_360.md", "archive_342_369/research_note_365.md",
    "archive_585_628/research_note_617.md", "archive_702_741/research_note_705.md",
    "archive_702_741/research_note_706.md", "archive_819_853/research_note_829.md",
    "archive_990_1008/research_note_1000.md", "archive_990_1008/research_note_1005.md",
    "archive_1009_1043/research_note_1035.md", "archive_1046_/research_note_1049.md",
    "archive_1046_/1055/mainline_acceptance.json",
    "archive_1046_/_shared/navigation_1055_checks.json",
    "archive_1046_/_shared/goal_start.json",
]
NAV = [ROOT / n for n in ("README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md")]
NAV += [PHASE / "README.md", PHASE / "文件索引.md"]
TAILS = {
    "README.md": ("## 当前入口：1044—1045阶段已结项", "69ac60ee0679e7587adec3e75f884de91b79700d8a24b13ecdef4ab44363ce7f"),
    "research_direction.md": ("以下保存上一阶段结项及历史方向，不覆盖当前目标。", "b06d7acbb9e2308e491ac9cbdea1f48400083abb454ff359031b2737d17231c2"),
    "RESEARCH_STATE.md": ("以下为1044—1045结项快照和更早历史，不覆盖本栏。", "15b1d1d926f6c4de40e891a7712308c904b4a4a4fd4a7c99578d0d7c3db68838"),
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def local_links(path, allow_pending=False, text=None):
    text = path.read_text(encoding="utf-8") if text is None else text
    count = 0
    for target in re.findall(r"\]\(([^)]+)\)", text):
        target = target.strip("<>").split("#", 1)[0]
        if not target or re.match(r"[a-zA-Z][\w+.-]*://", target):
            continue
        dest = (path.parent / unquote(target)).resolve()
        assert dest.exists() or (allow_pending and dest == RECEIPT.resolve()), (path, target)
        count += 1
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    if args.freeze:
        assert not RECEIPT.exists(), "Preserve the existing frozen receipt."

    # Mature-relation checks only: no cosmological solver or data fit.
    run = subprocess.run(
        [sys.executable, "-B", "-X", "utf8", str(HERE / "acoustic_transport_adoption/check.py")],
        capture_output=True, text=True, encoding="utf-8", check=True,
    )
    algebra = json.loads(run.stdout)
    assert algebra["all_checks_passed"]
    assert algebra["counting"] == {"new_science_groups": 0, "new_cognitive_axioms": 0}
    assert algebra["polarization"]["shear_coefficient"]["exact"] == "16/15"
    assert algebra["momentum"]["weighted_collision_sum"]["exact"] == "0"

    # Reviews must explicitly accept the exact author snapshots being frozen.
    for topic in ("pion_anomaly", "acoustic_transport"):
        report = HERE / f"{topic}_adoption.md"
        review = (HERE / f"{topic}_adoption_review.md").read_text(encoding="utf-8")
        assert "通过" in review and digest(report) in review.lower(), topic

    assets = []
    for rel in OWN:
        path = (HERE / rel).resolve()
        links = local_links(path, args.freeze) if path.suffix == ".md" else 0
        assets.append({"file": rel, "sha256": digest(path), "local_links": links})
    history = {rel: digest(ROOT / rel) for rel in HISTORY}
    tails = {}
    for name, (marker, expected) in TAILS.items():
        data = (ROOT / name).read_bytes()
        tail = data[data.index(marker.encode("utf-8")):]
        actual = sha256(tail).hexdigest()
        assert actual == expected, (name, "historical tail changed")
        tails[name] = actual

    nav_links = 0
    for path in NAV:
        text = path.read_text(encoding="utf-8")
        if path.parent == ROOT and path.name in TAILS:
            text = text.split(TAILS[path.name][0], 1)[0]
        nav_links += local_links(path, args.freeze, text)
    stable = {
        "date": "2026-10-09",
        "kind": "mature_physics_adoption_and_historical_nonduplication",
        "new_scientific_groups": 0, "new_cognitive_axioms": 0,
        "latest_accepted_round": 1055, "cumulative_scientific_groups": 3831,
        "stage_scientific_groups": 10,
        "assets": assets, "history_sha256": history,
        "historical_tails_sha256": tails,
        "independent_reviews_passed": True,
        "mature_algebra_readonly_reproduction_passed": True,
        "raw_data_fit_or_large_physics_solver_run": False,
        "scope": "Pion anomaly and photon-baryon transport are mature adoptions. Gauge/SSR composition reuses existing contracts. No original round or new empirical significance.",
        "all_passed": True, "roadmap_complete": False, "goal_status": "active",
    }
    if args.freeze:
        receipt = dict(stable)
        receipt["navigation_snapshot"] = {
            "sha256": {str(p.relative_to(ROOT)): digest(p) for p in NAV},
            "local_link_occurrences": nav_links,
        }
        with RECEIPT.open("x", encoding="utf-8") as stream:
            json.dump(receipt, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        saved.pop("navigation_snapshot")
        assert stable == saved, "Frozen adoption asset or history mismatch"
    print(json.dumps({"all_passed": True, "assets": len(assets), "history": len(history),
                      "adoption_local_links": sum(x["local_links"] for x in assets),
                      "current_navigation_local_links": nav_links,
                      "science_count_added": 0, "latest_round": 1055}, ensure_ascii=False))


if __name__ == "__main__":
    main()
