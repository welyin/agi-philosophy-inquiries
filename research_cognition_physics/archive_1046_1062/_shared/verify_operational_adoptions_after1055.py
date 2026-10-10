"""Freeze/read-only checks for the operational mature-adoption batch.

This verifies assets, links and algebra, not the truth of a unified theory.
Prior frozen assets and the original ROADMAP requirements remain unchanged.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BASE = runpy.run_path(str(HERE / "verify_adoptions_after1055.py"))
ROOT, PHASE = BASE["ROOT"], BASE["PHASE"]
RECEIPT = HERE / "operational_adoptions_after1055_checks.json"
OWN = [
    "neutrino_propagation_adoption.md", "neutrino_propagation_adoption_review.md",
    "muon_decay_readout_adoption.md", "muon_decay_readout_adoption_review.md",
    "muon_decay_readout_adoption/check.py", "muon_decay_readout_adoption/results.json",
    "optomechanical_readout_adoption.md", "optomechanical_readout_adoption_review.md",
    "early_late_expansion_audit.md", "early_late_expansion_audit_review.md",
    "../_admission/neutrino_propagation_after1055/selection.md",
    "../_admission/record_feedback_after1055/selection.md",
    "../_admission/common_history_after1055/selection.md",
    "integration_operational_adoptions_after1055.md", "verify_operational_adoptions_after1055.py",
]
HISTORY = [
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_956_989/research_note_985.md", "archive_956_989/research_note_986.md",
    "archive_990_1008/research_note_1002.md", "archive_1044_/research_note_1045.md",
    "archive_990_1008/1007/neutrino_observation_adoption_v1.md",
    "archive_990_1008/research_note_1008.md",
    "archive_1046_/research_note_1052.md", "archive_1046_/1055/mainline_acceptance.json",
    "archive_1046_/_shared/adoptions_after1055_checks.json",
]
ROADMAP_TAIL_SHA = "0572fd4469c75ffa86d4e08162ca11784dfdbaa983762c520076a5e152f10fa4"
ROADMAP_NAVIGATION_EDITS = [
    ("必需组织任务量词尚缺，统一一维自主旧家族未有全认知合同映射；自然选择仍开放",
     "必需组织任务与空间维数谓词尚缺；独立A合同不包含E_nat现实恢复规格，一维旧家族未有完整同义映射；自然选择仍开放"),
    ("候选全控制器问题不作全路线门槛",
     "弱衰变与LIGO实际读口另按各自任务采用；候选全控制器问题不作全路线门槛"),
    ("核子／核弱、中微子退耦及π⁰异常按原有效阶采用",
     "核子／核弱、中微子退耦／太阳传播、缪子形状及π⁰异常按原有效阶采用"),
    ("1045真实热后态；1050有限Coulomb仪器与材料能源界，实际准备／控制／记录总资源未闭合",
     "1045真实热后态；1050有限Coulomb仪器；旧有限共同历史直接复用，破坏性弱读数与光力反作用按域接入，跨接口资源接续未闭合"),
    ("宇宙碰撞史和钟的成熟桥按范围复用，整体未完",
     "宇宙碰撞史、味传播和钟的成熟桥按范围复用；早晚绝对标定保留已发表张力，整体未完"),
    ("认知特有预测与整体覆盖仍开放",
     "新增真实弱／光力读口采用及早晚宇宙分支审计；认知特有预测与整体覆盖仍开放"),
]


def prior_roadmap_hash():
    prior = json.loads((HERE / "adoptions_after1055_checks.json").read_text(encoding="utf-8"))
    return prior["navigation_snapshot"]["sha256"]["ROADMAP.md"]


def verify_scope_transition(raw):
    """Reverse only this batch's prefix edits in memory; never alter the file.

    The initial tail digest transcribed from a conversation summary was wrong.
    This full-file comparison against the prior frozen receipt establishes the
    actual unchanged tail without treating that summary as authoritative.
    """
    start = raw.index("### 实际操作与早晚宇宙重叠验收".encode("utf-8"))
    end = raw.index("## 研究方法补充".encode("utf-8"), start)
    restored = raw[:start] + raw[end:]
    for old, new in ROADMAP_NAVIGATION_EDITS:
        assert restored.count(new.encode("utf-8")) == 1
        restored = restored.replace(new.encode("utf-8"), old.encode("utf-8"), 1)
    assert sha256(restored).hexdigest() == prior_roadmap_hash()
    marker = "## M1 检验认知要求的独立选择力".encode("utf-8")
    assert restored[restored.index(marker):] == raw[raw.index(marker):]


def checked_links(path, pending):
    # The base checker permits only its own prior receipt to be pending.
    content = path.read_text(encoding="utf-8")
    if pending:
        content = content.replace("(operational_adoptions_after1055_checks.json)", "(#pending-receipt)")
        content = content.replace("(archive_1046_/_shared/operational_adoptions_after1055_checks.json)", "(#pending-receipt)")
        content = content.replace("(_shared/operational_adoptions_after1055_checks.json)", "(#pending-receipt)")
    if path.parent == ROOT and path.name in BASE["TAILS"]:
        content = content.split(BASE["TAILS"][path.name][0], 1)[0]
    return BASE["local_links"](path, False, content)


def run_check(script):
    run = subprocess.run([sys.executable, "-B", "-X", "utf8", str(script)],
                         capture_output=True, text=True, encoding="utf-8", check=True)
    return json.loads(run.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    if args.freeze:
        assert not RECEIPT.exists(), "Keep the existing receipt frozen."
    # The previous checker scans live navigation; create our receipt only after
    # checking its assets, then run that previous entry separately once recorded.
    algebra = run_check(HERE / "muon_decay_readout_adoption/check.py")
    assert algebra["all_checks_passed"] and algebra["new_science_groups"] == 0
    assert algebra["exact"]["hemisphere_spin_coefficient"] == "1/12"
    assets = []
    for rel in OWN:
        path = (HERE / rel).resolve()
        if path.suffix == ".md":
            checked_links(path, args.freeze)
        assets.append({"file": rel, "sha256": BASE["digest"](path)})
    for topic in ("neutrino_propagation", "muon_decay_readout", "optomechanical_readout"):
        review = (HERE / (topic + "_adoption_review.md")).read_text(encoding="utf-8")
        assert "通过" in review
        assert BASE["digest"](HERE / (topic + "_adoption.md")) in review.lower()
    early_late_review = (HERE / "early_late_expansion_audit_review.md").read_text(encoding="utf-8")
    assert "通过" in early_late_review
    assert BASE["digest"](HERE / "early_late_expansion_audit.md") in early_late_review.lower()
    history = {p: BASE["digest"](ROOT / p) for p in HISTORY}
    for name, (marker, expected) in BASE["TAILS"].items():
        raw = (ROOT / name).read_bytes()
        assert sha256(raw[raw.index(marker.encode("utf-8")):]).hexdigest() == expected
    roadmap = (ROOT / "ROADMAP.md").read_bytes()
    marker = "## M1 检验认知要求的独立选择力".encode("utf-8")
    assert sha256(roadmap[roadmap.index(marker):]).hexdigest() == ROADMAP_TAIL_SHA
    if args.freeze:
        verify_scope_transition(roadmap)
    nav_links = sum(checked_links(path, args.freeze) for path in BASE["NAV"])
    stable = {
        "date": "2026-10-09", "kind": "operational_mature_adoptions_and_scope_audit",
        "assets": assets, "history_sha256": history,
        "roadmap_original_requirements_sha256": ROADMAP_TAIL_SHA,
        "roadmap_previous_full_snapshot_sha256": prior_roadmap_hash(),
        "roadmap_prefix_only_transition_verified_at_freeze": True,
        "new_scientific_groups": 0, "new_cognitive_axioms": 0,
        "latest_accepted_round": 1055, "cumulative_scientific_groups": 3831,
        "stage_scientific_groups": 10, "independent_reviews_passed": True,
        "muon_algebra_readonly_reproduction_passed": True,
        "new_experimental_fit": False, "full_quantum_instruments_certified": False,
        "roadmap_complete": False, "goal_status": "active", "all_passed": True,
    }
    if args.freeze:
        saved = dict(stable)
        saved["navigation_snapshot"] = {
            "sha256": {str(p.relative_to(ROOT)): BASE["digest"](p) for p in BASE["NAV"]},
            "local_link_occurrences": nav_links,
        }
        with RECEIPT.open("x", encoding="utf-8") as out:
            json.dump(saved, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        saved.pop("navigation_snapshot")
        assert stable == saved, "Frozen asset/history mismatch"
    print(json.dumps({"all_passed": True, "assets": len(assets), "history": len(history),
                      "current_navigation_local_links": nav_links,
                      "latest_accepted_round": 1055, "new_science_groups": 0}))


if __name__ == "__main__":
    main()
