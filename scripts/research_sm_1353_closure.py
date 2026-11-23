"""Resolve 1353–1363 evidence under its closed archive name and verify its preimages."""
import base64
import hashlib
import json
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research_cognition_physics"
OLD = BASE / "archive_1353_"
CLOSED = BASE / "archive_1353_1363"
MANIFEST = BASE / "_migration/closure_1353_1363/manifest.json"
ROOT_NAV = {"README.md", "research_direction.md", "RESEARCH_STATE.md", "ROADMAP.md"}


@lru_cache(maxsize=1)
def _record():
    return json.loads(MANIFEST.read_bytes())


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def physical_1353(path):
    """Resolve an absolute or project-relative evidence name without following it outside the project."""
    actual = (ROOT / path).resolve()
    if not actual.is_relative_to(ROOT):
        raise ValueError(f"Path outside the project: {path}")
    if actual.is_relative_to(OLD):
        actual = (CLOSED / actual.relative_to(OLD)).resolve()
        if not actual.is_relative_to(CLOSED):
            raise ValueError(f"Path outside the closed archive: {path}")
    return actual


def original_1353_bytes(path, raw=None):
    """Undo precisely this closure layer, leaving earlier layers for their original validators."""
    actual = physical_1353(path)
    supplied = raw is not None
    if raw is None:
        raw = actual.read_bytes()
    if actual.parent == BASE and actual.name in ROOT_NAV:
        return raw.replace(b"archive_1353_1363/", b"archive_1353_/")
    name = actual.relative_to(ROOT).as_posix()
    change = _record()["edits"].get(name)
    if change is None:
        return raw
    digest = _sha(raw)
    if supplied and digest == change["before_sha256"]:
        return raw
    if digest != change["after_sha256"]:
        raise ValueError(f"Unexpected on-disk revision of {name}")
    if change["operation"] == "path":
        raw = raw.replace(b"archive_1353_1363", b"archive_1353_")
    else:
        for at, size, previous in reversed(change["reverse_hunks"]):
            raw = raw[:at] + base64.b64decode(previous, validate=True) + raw[at + size:]
    if _sha(raw) != change["before_sha256"]:
        raise ValueError(f"Original SHA-256 mismatch for {name}")
    return raw


def verify():
    data = _record()
    for name, expected in data["archive_before_sha256"].items():
        if _sha(original_1353_bytes(name)) != expected:
            raise ValueError(f"Changed original evidence: {name}")
    for name, change in data["edits"].items():
        if _sha(original_1353_bytes(name)) != change["before_sha256"]:
            raise ValueError(f"Changed migration evidence: {name}")
    expected = {physical_1353(name) for name in data["archive_before_sha256"]}
    actual = {p for p in CLOSED.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    if expected != actual or OLD.exists():
        raise ValueError("Closed archive inventory differs from the original")
    return {"status": "PASS", "original_files": len(expected), "reversible_files": len(data["edits"]),
            "manifest_sha256": _sha(MANIFEST.read_bytes())}


if __name__ == "__main__":
    print(json.dumps(verify()))
