"""Verify the twelve closed SM archives and reproduce the preceding navigation."""
import hashlib
import json
from pathlib import Path
from research_sm_1353_closure import original_1353_bytes
from research_navigation_1818 import before_stage_1818
from run_research_sm_archive import history, verify as verify_registrations

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research_cognition_physics"
CLOSURE = BASE / "_migration/closure_1364_1815"
MANIFEST = BASE / "archive_1798_1815/_shared/navigation_manifest.json"
PREVIOUS = BASE / "archive_1353_1363/_shared/navigation_manifest.json"
PREVIOUS_RUNNER = ROOT / "scripts/research_navigation_1353.py"
ROOT_NAV = {"README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md"}
BEGIN = b"<!-- RESEARCH_STAGE_1364_BEGIN -->"
END = b"<!-- RESEARCH_STAGE_1364_END -->"
RUNNER_EDITS = (
    (b"from pathlib import Path, PurePosixPath",
     b"from pathlib import Path, PurePosixPath\nfrom research_navigation_1364 import before_current_stage"),
    (b"    raw = path.read_bytes()", b"    raw = before_current_stage(path)"),
    (b"and sha((ROOT / name).read_bytes()) == expected",
     b"and sha(before_current_stage(ROOT / name)) == expected"),
)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, why):
    if not condition:
        raise ValueError(why)


def before_current_stage(path):
    raw = original_1353_bytes(path, before_stage_1818(path))
    if path == PREVIOUS_RUNNER:
        for old, new in RUNNER_EDITS:
            require(raw.count(new) == 1, "Unexpected replay runner change")
            raw = raw.replace(new,old)
        return raw
    if path.parent != BASE or path.name not in ROOT_NAV:
        return raw
    item = json.loads(MANIFEST.read_bytes())["root_navigation"][path.name]
    closure = json.loads((CLOSURE / "closure.json").read_bytes())
    products = closure["products"]
    require(digest(before_stage_1818(path)) == products[path.relative_to(ROOT).as_posix()],
            "Closed navigation changed: " + path.name)
    for old, new in closure.get("root_reference_edits", {}).get(path.name, []):
        require(raw.count(new.encode("utf8")) == 1, "Ambiguous navigation reference")
        raw = raw.replace(new.encode("utf8"), old.encode("utf8"))
    saved = json.loads((CLOSURE / "root_blocks_at_acceptance.json").read_bytes())
    require(raw.count(BEGIN) == raw.count(END) == 1, "Ambiguous navigation block")
    a = raw.index(BEGIN)
    b = raw.index(END, a) + len(END)
    raw = raw[:a] + saved[path.name].encode("utf8") + raw[b:]
    require(digest(raw) == item["current_sha256"], "Could not reproduce final acceptance: " + path.name)
    title_end = raw.index(b"\n") + 1
    newline = b"\r\n" if raw[:title_end].endswith(b"\r\n") else b"\n"
    begin, end = newline + BEGIN + newline, END + newline
    require(raw[title_end:].startswith(begin) and raw.count(BEGIN) == raw.count(END) == 1,
            "Invalid current navigation block: " + path.name)
    block_end = raw.index(end, title_end) + len(end)
    before = raw[:title_end] + raw[block_end:]
    require(digest(before) == item["previous_sha256"], "Could not restore prior navigation: " + path.name)
    return before


def verify():
    raw = MANIFEST.read_bytes()
    data = json.loads(raw)
    require(data["schema"] == "research_navigation_1364_v1" and set(data["root_navigation"]) == ROOT_NAV,
            "Wrong current manifest schema or missing navigation")
    prior_raw = PREVIOUS.read_bytes()
    require(digest(prior_raw) == data["previous_manifest_sha256"], "Prior closed manifest changed")
    prior = json.loads(prior_raw)
    for name in ROOT_NAV:
        require(digest(before_current_stage(BASE / name)) == prior["root_navigation"][name]["current_sha256"],
                "Failed prior stage boundary: " + name)
    support = {"scripts/research_navigation_1364.py": digest((CLOSURE / "navigation_1364_v1.py").read_bytes()),
               "scripts/research_navigation_1353.py": digest(original_1353_bytes(PREVIOUS_RUNNER))}
    require(support == data["support_sha256"], "Frozen verification support changed")
    closed = json.loads((CLOSURE / "closure.json").read_bytes())
    for name, expected in closed["products"].items():
        file = (ROOT / name).resolve()
        require(file.is_relative_to(ROOT) and digest(before_stage_1818(file)) == expected, "Changed closure product: " + name)
    files, manifest_hash = history(None)
    rounds = verify_registrations(files)
    checks = sum(len(v) for v in rounds.values())
    require(list(rounds) == list(range(1364,1816)) and checks == 930, "Incomplete round coverage")
    move = json.loads((CLOSURE / "manifest.json").read_bytes())
    require([n for p in move["phases"] for n in range(p["lo"],p["hi"]+1)] == list(rounds),
            "Physical phase coverage differs from the accepted rounds")
    source = json.loads(files["research_cognition_physics/archive_1364_/_shared/goal_start_receipt.json"])
    require(source["first_round"] == 1364 and source["status_at_start"] == data["goal_at_start"] == "active",
            "Original goal boundary changed")
    return dict(status="PASS", previous_stage=1363, accepted_through=1815, root_navigation=len(ROOT_NAV),
                phases=len(move["phases"]), round_acceptances=len(rounds), registered_checks=checks,
                migration_manifest_sha256=manifest_hash, navigation_manifest_sha256=digest(raw))


if __name__ == "__main__":
    print(json.dumps(verify(), ensure_ascii=False))

