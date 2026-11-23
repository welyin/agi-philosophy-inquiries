"""Current research navigation and reversible access to the preceding closed stage."""
import hashlib
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research_cognition_physics"
STAGE = BASE / "archive_1364_"
MANIFEST = STAGE / "_shared/navigation_manifest.json"
PREVIOUS = BASE / "archive_1353_/_shared/navigation_manifest.json"
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
    raw = path.read_bytes()
    if path == PREVIOUS_RUNNER:
        for old, new in RUNNER_EDITS:
            require(raw.count(new) == 1, "Unexpected replay runner change")
            raw = raw.replace(new,old)
        return raw
    if path.parent != BASE or path.name not in ROOT_NAV:
        return raw
    item = json.loads(MANIFEST.read_bytes())["root_navigation"][path.name]
    require(digest(raw) == item["current_sha256"], "Current navigation changed: " + path.name)
    title_end = raw.index(b"\n") + 1
    newline = b"\r\n" if raw[:title_end].endswith(b"\r\n") else b"\n"
    begin, end = newline + BEGIN + newline, END + newline
    require(raw[title_end:].startswith(begin) and raw.count(BEGIN) == raw.count(END) == 1,
            "Invalid current navigation block: " + path.name)
    block_end = raw.index(end, title_end) + len(end)
    before = raw[:title_end] + raw[block_end:]
    require(digest(before) == item["previous_sha256"], "Could not restore prior navigation: " + path.name)
    return before


def stage_file(name):
    rel = PurePosixPath(name)
    require(not rel.is_absolute() and '..' not in rel.parts and ':' not in name and '\\' not in name,
            "Unsafe evidence path")
    target = (ROOT / name).resolve()
    require(target.is_relative_to(STAGE), "Evidence outside the active archive")
    return target


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
    support = [Path(__file__).resolve(), PREVIOUS_RUNNER]
    require({p.relative_to(ROOT).as_posix(): digest(p.read_bytes()) for p in support} == data["support_sha256"],
            "Replay support changed")
    for name, expected in data["stage_files_sha256"].items():
        require(digest(stage_file(name).read_bytes()) == expected, "Stage file changed: " + name)
    round_count = 0
    evidence_count = 0
    for n in range(1364, data["accepted_through"] + 1):
        acc = STAGE / str(n) / "acceptance.json"
        require(acc.relative_to(ROOT).as_posix() in data["stage_files_sha256"], "Missing acceptance")
        item = json.loads(acc.read_bytes())
        require(item["round"] == n and item["status"].startswith("ACCEPTED"), "Invalid acceptance scope")
        round_count += 1
        for name, expected in item["sha256"].items():
            require(digest(stage_file(name).read_bytes()) == expected, "Changed scientific evidence: " + name)
            evidence_count += 1
    start = json.loads((STAGE / "_shared/goal_start_receipt.json").read_bytes())
    require(data["goal_at_start"] == "active", "Unexpected start state")
    require(start.get("status_at_start") == "active" and start.get("first_round") == 1364,
            "Start receipt does not identify the active 1364 stage")
    return dict(status="PASS", previous_stage=1363, accepted_through=data["accepted_through"],
                root_navigation=len(ROOT_NAV), round_acceptances=round_count, bound_evidence=evidence_count,
                stage_files=len(data["stage_files_sha256"]), navigation_manifest_sha256=digest(raw))


if __name__ == "__main__":
    print(json.dumps(verify(), ensure_ascii=False))

