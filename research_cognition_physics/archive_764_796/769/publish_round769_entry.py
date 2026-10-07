"""Publish a verified working entry; do not increment scientific round counters."""
import ast
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
TARGET = HERE/"round769_drafts/absolute_source_entry_checks.json"
assert not TARGET.exists()
history = dict(core.read(HERE/"round584_drafts/historical_evidence_manifest.json")["evidence_hashes"])
for n in range(584, 769):
    receipt = core.read(HERE/f"research_round_{n}_checks.json")
    for key in ("new_file_hashes", "preserved_draft_hashes"):
        for name, digest in receipt[key].items():
            assert name not in history or history[name] == digest
            history[name] = digest
history.update(core.read(HERE/"cognitive_foundation_bridge_605_navigation.json")["supplementary_artifact_hashes"])
assert len(history) == 3726
for name, digest in history.items():
    assert core.digest(HERE/name) == digest, name
artifacts = (
    "round769_drafts/research_note_769_working.md",
    "round769_drafts/bv_source_entry_calibration.py",
    "round769_drafts/bv_source_entry_calibration_results.json",
    "round769_drafts/entry_condition_ledger.md",
    "round769_drafts/entry_scope_review.json",
)
spec = importlib.util.spec_from_file_location("entry769calibration", HERE/artifacts[1])
calibration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calibration)
assert calibration.run() == core.read(HERE/artifacts[2])
ast.parse((HERE/artifacts[1]).read_text("utf8"))
report_text = (HERE/artifacts[0]).read_text("utf8")
assert report_text.count("$$") == 14
for link in core.link_parser()(report_text):
    assert ((HERE/artifacts[0]).parent/link).resolve().exists(), link
text_checks = core.text_checks(HERE/artifacts[0])
assert text_checks["display_formulas"] == 7
paths = [ROOT/"README.md", RESEARCH/"README.md", RESEARCH/"research_direction.md",
         RESEARCH/"RESEARCH_STATE.md", HERE/"README.md",
         HERE/"spatial_premise_closure_audit.md", HERE/"three_dimensional_four_conditions_review.md"]
before = {p: p.read_bytes() for p in paths}
snapshot = HERE/"navigation_before_round769_entry_20261004"
assert not snapshot.exists()
summary = "**769绝对来源工作报告已保存，正式仍768／3491：** [研究报告]({p}round769_drafts/research_note_769_working.md)核成熟EFT/BV定理的原对象映射，区分规范反常、来源接触及背景比较；两组原势/规范固定局部校准已复算，未宣称实际量子反常或全部门Ward完成。下一项先接固定背景首阶绝对来源，不先要求全阶背景独立。[条件增量]({p}round769_drafts/entry_condition_ledger.md)、[入口核验]({p}round769_drafts/absolute_source_entry_checks.json)。不新增完成轮次，目标保持。"
planned = {}
for p, raw in before.items():
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
    nl = "\r\n" if b"\r\n" in raw else "\n"
    txt = raw.decode(enc).replace("\r\n", "\n")
    assert "**769绝对来源工作报告已保存" not in txt
    prefix = "research_cognition_physics/archive_231_/" if p.parent == ROOT else "archive_231_/" if p.parent == RESEARCH else ""
    addition = summary.format(p=prefix)
    if p in paths[:5]:
        head, tail = txt.split("\n\n", 1)
        txt = head+"\n\n"+addition+"\n\n"+tail
    else:
        txt += "\n\n"+addition+"\n"
    planned[p] = txt.replace("\n", nl).encode(enc)
links = 0
for p, raw in planned.items():
    for link in core.link_parser()(raw.decode("utf-8-sig")):
        target = (p.parent/link).resolve()
        assert target.exists() or target == TARGET.resolve(), (p, link)
        links += 1
assert all(p.read_bytes() == raw for p, raw in before.items())
snapshot.mkdir()
manifest = {}
for i, (p, raw) in enumerate(before.items()):
    name = f"{i}_{p.name}"
    with (snapshot/name).open("xb") as f:
        f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()] = dict(snapshot=name, sha256=hashlib.sha256(raw).hexdigest())
with (snapshot/"manifest.json").open("x", encoding="utf8") as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)
for p, raw in planned.items():
    assert p.read_bytes() == before[p]
    temp = p.with_name(p.name+".round769entry.tmp")
    with temp.open("xb") as f:
        f.write(raw)
    os.replace(temp, p)
report = dict(status="working entry only", latest_completed_round=768,
              cumulative_scientific_tests=3491, protected_scientific_evidence_files=3726,
              previous_evidence_unchanged=True, new_scientific_tests=0, entry_calibrations_reproduced=2, entry_text_checks=text_checks,
              working_artifact_hashes={n: core.digest(HERE/n) for n in artifacts},
              navigation_links=links, broken_links=0,
              published_navigation_hashes={p.relative_to(ROOT).as_posix(): core.digest(p) for p in paths},
              inherited_767_768_free_state_preserved=True, background_split_anomaly_computed=False, all_sector_Ward_claimed=False, active_goal_unchanged=True,
              all_checks_passed=True)
with TARGET.open("x", encoding="utf8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)
print(json.dumps({k: report[k] for k in ("latest_completed_round", "navigation_links", "all_checks_passed")}))
