"""Locate the active stage and recover the exact navigation bytes accepted at 1815."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research_cognition_physics"
NAMES = {"README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md"}
BEGIN = b"<!-- RESEARCH_STAGE_1818_BEGIN -->"
END = b"<!-- RESEARCH_STAGE_1818_END -->"
RUNNER = ROOT / "scripts/research_navigation_1364.py"
EDITS = (
    (b"from research_sm_1353_closure import original_1353_bytes",
     b"from research_sm_1353_closure import original_1353_bytes\nfrom research_navigation_1818 import before_stage_1818"),
    (b"    raw = original_1353_bytes(path)", b"    raw = original_1353_bytes(path, before_stage_1818(path))"),
    (b"require(digest(path.read_bytes()) == products", b"require(digest(before_stage_1818(path)) == products"),
    (b"and digest(file.read_bytes()) == expected", b"and digest(before_stage_1818(file)) == expected"),
)


def before_stage_1818(path):
    raw = path.read_bytes()
    if path == RUNNER:
        for old, new in EDITS:
            if raw.count(new) != 1:
                raise ValueError("Unrecognized edit in the preceding navigation adapter")
            raw = raw.replace(new, old)
    if path.parent != BASE or path.name not in NAMES:
        return raw
    title_end = raw.index(b"\n") + 1
    nl = b"\r\n" if raw[:title_end].endswith(b"\r\n") else b"\n"
    if not raw[title_end:].startswith(nl + BEGIN + nl) or raw.count(BEGIN) != 1 or raw.count(END) != 1:
        raise ValueError("Unrecognized stage 1818 navigation block in " + path.name)
    after = raw.index(END + nl, title_end) + len(END + nl)
    return raw[:title_end] + raw[after:]
