"""Verify the active 1353 stage and supply pre-stage bytes to historical replay."""
import hashlib
import json
from pathlib import Path, PurePosixPath
from research_navigation_1364 import before_current_stage
from research_sm_1353_closure import CLOSED, original_1353_bytes, physical_1353

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research_cognition_physics"
STAGE = CLOSED
MANIFEST = STAGE / "_shared/navigation_manifest.json"
ADAPTER = ROOT / "scripts/research_closure_1352.py"
ROOT_NAV = {"README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md"}
BEGIN = b"<!-- RESEARCH_STAGE_1353_BEGIN -->"
END = b"<!-- RESEARCH_STAGE_1353_END -->"
# Prefix-free byte edits. A token, not a line break, ends each replacement.
ADAPTER_EDITS = (
    (b"from research_layout import ROOT, no_archive_or_research_writes",
     b"from research_layout import ROOT, no_archive_or_research_writes\nfrom research_navigation_1353 import before_stage"),
    (b"restore(path.read_bytes(), entry)", b"restore(before_stage(path), entry)"),
    (b"digest(bounded(ROOT, name).read_bytes()) == expected, 'Closure product",
     b"digest(before_stage(bounded(ROOT, name))) == expected, 'Closure product"),
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def check(test, message):
    if not test:
        raise ValueError(message)


def before_stage(path):
    raw = before_current_stage(path)
    if path == ADAPTER:
        for before, after in ADAPTER_EDITS:
            check(raw.count(after) == 1, "Unrecognized closure adapter change")
            raw = raw.replace(after, before)
        return raw
    if path.parent != BASE or path.name not in ROOT_NAV:
        return raw
    # The frozen closure independently verifies the recovered document digest.
    info = json.loads(MANIFEST.read_bytes())["root_navigation"][path.name]
    check(sha(raw) == info["current_sha256"], "Current navigation changed: " + path.name)
    cut = raw.index(b"\n") + 1
    newline = b"\r\n" if raw[cut-2:cut] == b"\r\n" else b"\n"
    start, end = newline + BEGIN + newline, END + newline
    check(raw[cut:].startswith(start) and raw.count(BEGIN) == raw.count(END) == 1,
          "Invalid stage navigation block: " + path.name)
    stop = raw.index(end, cut) + len(end)
    old = raw[:cut] + raw[stop:]
    check(sha(old) == info["closure_sha256"], "Navigation preimage changed: " + path.name)
    return old


def verify():
    data = json.loads(MANIFEST.read_bytes())
    check(data["schema"] == "research_navigation_1353_v1" and set(data["root_navigation"]) == ROOT_NAV,
          "Unexpected active stage navigation")
    old = json.loads((BASE / "_migration/closure_1323_1352/manifest.json").read_bytes())
    old_digests = {p["logical"]:p["current_sha256"] for p in old["changed_markdown"]}
    for name in ROOT_NAV:
        path = BASE / name
        check(sha(before_stage(path)) == old_digests[path.relative_to(ROOT).as_posix()], "Prior closure mismatch")
    for name, expected in data["support_sha256"].items():
        check(name in ("scripts/research_closure_1352.py", "scripts/research_navigation_1353.py")
              and sha(before_current_stage(ROOT / name)) == expected, "Active stage replay support changed")
    check(len(data["support_sha256"]) == 2, "Incomplete replay support")
    checked = 0
    def located(name):
        rel = PurePosixPath(name)
        check(not rel.is_absolute() and ".." not in rel.parts and ":" not in name and "\\" not in name,
              "Unsafe evidence name: " + name)
        path = physical_1353(ROOT / rel)
        check(path.resolve().is_relative_to(STAGE), "Evidence outside active stage")
        return path
    for name, expected in data["stage_files_sha256"].items():
        check(sha(original_1353_bytes(located(name))) == expected, "Changed stage evidence: " + name)
    for n in range(1353, data["accepted_through"] + 1):
        path = STAGE / str(n) / "acceptance.json"
        historical = f"research_cognition_physics/archive_1353_/{n}/acceptance.json"
        check(historical in data["stage_files_sha256"], "Unbound acceptance")
        item = json.loads(path.read_bytes())
        check(item["round"] == n and item["status"].startswith("ACCEPTED"), "Unexpected acceptance status")
        for name, expected in item["sha256"].items():
            check(sha(original_1353_bytes(located(name))) == expected, "Acceptance evidence changed: " + name)
            checked += 1
    print(json.dumps(dict(status="PASS",accepted_through=data["accepted_through"],root_navigation=4,
                          bound_evidence=checked,stage_files=len(data["stage_files_sha256"]),
                          stage_manifest_sha256=sha(MANIFEST.read_bytes()))))


if __name__ == "__main__":
    verify()

