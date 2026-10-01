"""Reproduce the historical index and publish the forward-bridge review.

Documentation audit only: no historical experiment is rerun and no scientific
round is added. With no arguments, verify the frozen index and publication.
"""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import os
import re

import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
STEM = "cognition_forward_bridge_review_585"
REVIEW = HERE / (STEM + ".md")
INVENTORY = HERE / (STEM + "_inventory.json")
CHECKS = HERE / (STEM + "_checks.json")
SNAPSHOT = HERE / "navigation_before_forward_review_585_20261001"
NAVIGATION = [ROOT / "README.md", RESEARCH / "README.md",
              RESEARCH / "research_direction.md", RESEARCH / "RESEARCH_STATE.md",
              HERE / "README.md"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_new(path, value):
    with path.open("x", encoding="utf8", newline="\n") as f:
        f.write(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def index():
    links = core.link_parser()
    notes = []
    for number in range(231, 586):
        p = HERE / f"research_note_{number}.md"
        body = p.read_text(encoding="utf-8-sig")
        heading = next(line for line in body.splitlines() if line.startswith("# "))
        assert str(number) in heading, p
        notes.append(dict(round=number, path=p.name, title=heading[2:],
                          sha256=sha(p.read_bytes())))
    cited = {}
    dependencies = {}
    for link in links(REVIEW.read_text(encoding="utf8")):
        p = (HERE / link).resolve()
        if p.name.startswith(STEM):
            continue  # Audit products are hashed separately, without a cycle.
        assert p.is_relative_to(ROOT) and p.exists(), p
        cited[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
        if re.fullmatch(r"research_note_\d+\.md", p.name):
            for child in links(p.read_text(encoding="utf-8-sig")):
                q = (p.parent / child).resolve()
                if q.suffix in (".py", ".json"):
                    assert q.is_relative_to(ROOT) and q.exists(), q
                    dependencies[q.relative_to(ROOT).as_posix()] = sha(q.read_bytes())
    assert len(notes) == 355
    return dict(date="2026-10-01",scope="historical_documentation_audit_only",
                note_range=[231, 585], title_inventory_count=len(notes), notes=notes,
                review_direct_local_sources=dict(sorted(cited.items())),
                cited_notes_code_and_result_links=dict(sorted(dependencies.items())),
                coverage="All titles indexed; selected text, formulas and premises reviewed, not all full texts.",
                historical_experiments_rerun=False, new_scientific_round=False)


def check_links(path, raw, allow_pending=False):
    count = 0
    for link in core.link_parser()(raw.decode("utf-8-sig")):
        q = (path.parent / link).resolve()
        assert q.exists() or (allow_pending and q in (INVENTORY, CHECKS)), (path, link)
        count += 1
    return count


def validate_review(allow_pending=False):
    body = REVIEW.read_text(encoding="utf8")
    assert len(re.findall(r"^\$\$\s*$", body, re.M)) == 4
    assert not any(ord(c) < 32 and c not in "\r\n\t" for c in body)
    ast.parse(Path(__file__).read_text(encoding="utf8"))
    return check_links(REVIEW, REVIEW.read_bytes(), allow_pending)


def publication():
    assert not CHECKS.exists() and not INVENTORY.exists() and not SNAPSHOT.exists(), "already published"
    latest = read(HERE / "research_round_585_checks.json")
    assert latest["round"] == 585 and latest["all_reported_checks_passed"]
    original_support = read(HERE / "round585_navigation_checks.json")["supporting_document_hashes"]
    for name, digest in original_support.items():
        assert sha((HERE / name).read_bytes()) == digest, name
    inventory = index()
    review_links = validate_review(True)
    originals = {p: p.read_bytes() for p in NAVIGATION}
    planned = {}
    for p, raw in originals.items():
        encoding = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
        newline = "\r\n" if b"\r\n" in raw else "\n"
        body = raw.decode(encoding).replace("\r\n", "\n")
        assert STEM not in body, p
        prefix = "research_cognition_physics/archive_231_/" if p.parent == ROOT else "archive_231_/" if p.parent == RESEARCH else ""
        block = ("**231—585轮正向成果回顾（2026-10-01）：** [认知结果如何约束统一模型](" + prefix + STEM + ".md)"
                 "已整理355份报告目录，并重点复核共同记录、跨层过程、有限反作用、规范拼接、动态局域性及物理条件接口。"
                 "量子结构已进入当前模型；未完成的是同一过程与动态几何、记录和尺度的连接。本次为回顾，不增加研究轮次，科学基线仍为585／2984。\n\n"
                 "**当前执行顺序（优先于下方历史安排）：** 先将旧正向成果落实为当前物质模型的对象、仪器、组合、参考与尺度映射；"
                 "再检验局部时间、能源流和几何响应。357／361的有限反作用极限作为待核候选，原586能源流入口保留；"
                 "不把旧模型的Hamiltonian或额外前提直接移入新模型。目标不改，不重做已有实验。")
        first, rest = body.split("\n\n", 1)
        planned[p] = (first + "\n\n" + block + "\n\n" + rest).replace("\n", newline).encode(encoding)
    nav_links = sum(check_links(p, data, True) for p, data in planned.items())
    # Preserve every byte before touching navigation; reject concurrent changes.
    assert all(p.read_bytes() == raw for p, raw in originals.items()), "concurrent change"
    SNAPSHOT.mkdir(exist_ok=False)
    manifest = {}
    for i, (p, raw) in enumerate(originals.items()):
        before = f"{i}_before_{p.name}"
        after = f"{i}_published_{p.name}"
        (SNAPSHOT / before).write_bytes(raw)
        (SNAPSHOT / after).write_bytes(planned[p])
        manifest[p.relative_to(ROOT).as_posix()] = dict(before=before, before_sha256=sha(raw),
                                                       published=after, published_sha256=sha(planned[p]))
    write_new(SNAPSHOT / "manifest.json", manifest)
    write_new(INVENTORY, inventory)
    for p, data in planned.items():
        assert p.read_bytes() == originals[p], ("concurrent change", p)
        temp = p.with_name(p.name + ".forward_review_585.tmp")
        with temp.open("xb") as f:
            f.write(data)
        os.replace(temp, p)
    report = dict(date="2026-10-01", latest_scientific_round=585, cumulative_numbered_tests=2984,
                  new_scientific_round=False, title_inventory_count=355,
                  directly_cited_documents=len(inventory["review_direct_local_sources"]),
                  cited_code_and_result_files=len(inventory["cited_notes_code_and_result_links"]),
                  review_links=review_links, navigation_files=len(planned), navigation_links=nav_links,
                  broken_links=0, snapshot_manifest_sha256=sha((SNAPSHOT / "manifest.json").read_bytes()),
                  review_sha256=sha(REVIEW.read_bytes()), inventory_sha256=sha(INVENTORY.read_bytes()),
                  audit_script_sha256=sha(Path(__file__).read_bytes()),
                  preserved_round585_support_hashes=original_support,
                  full_historical_text_review=False, historical_experiments_rerun=False,
                  navigation_history_preserved=True, active_goal_unchanged=True,
                  no_new_automation_or_task=True, all_documentation_checks_passed=True)
    write_new(CHECKS, report)
    return verify()


def verify():
    report = read(CHECKS)
    assert read(INVENTORY) == index(), "source/index changed"
    assert report["review_sha256"] == sha(REVIEW.read_bytes())
    assert report["inventory_sha256"] == sha(INVENTORY.read_bytes())
    assert report["audit_script_sha256"] == sha(Path(__file__).read_bytes())
    assert report["review_links"] == validate_review()
    manifest_path = SNAPSHOT / "manifest.json"
    assert sha(manifest_path.read_bytes()) == report["snapshot_manifest_sha256"]
    count = 0
    for name, entry in read(manifest_path).items():
        for label in ("before", "published"):
            assert sha((SNAPSHOT / entry[label]).read_bytes()) == entry[label + "_sha256"]
        count += check_links(ROOT / name, (SNAPSHOT / entry["published"]).read_bytes())
    assert count == report["navigation_links"]
    for name, digest in report["preserved_round585_support_hashes"].items():
        assert sha((HERE / name).read_bytes()) == digest, name
    return {k: report[k] for k in ("latest_scientific_round", "title_inventory_count",
            "directly_cited_documents", "cited_code_and_result_files", "navigation_links",
            "broken_links", "new_scientific_round", "all_documentation_checks_passed")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    print(json.dumps(publication() if args.publish else verify(), ensure_ascii=False))
