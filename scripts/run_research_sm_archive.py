"""Replay the archived checks for rounds 1364–1815 in a disposable historical layout.

Examples: python scripts/run_research_sm_archive.py 1364 1404 1809-1811
          python scripts/run_research_sm_archive.py --all --verify-only
          python scripts/run_research_sm_archive.py --all --jobs 4 --report receipts.json
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tempfile

from research_layout import relocation_original_bytes


ROOT = Path(__file__).resolve().parents[1]
BASE = "research_cognition_physics"
OLD = BASE + "/archive_1364_"
NAV = OLD + "/_shared/navigation_manifest.json"
MANIFEST = ROOT / BASE / "_migration/closure_1364_1815/manifest.json"
FOLLOWUP_MANIFEST = ROOT / BASE / "_migration/closure_1353_1363/manifest.json"
FOLLOWUP_NOTE = OLD + "/_admission/round1478_admission.md"
# Earlier archives referenced by the frozen programs or their local input dictionaries.
PREVIOUS = {
    BASE + "/archive_1046_1062/1061/public_inputs.json": "63b9e1e5e8cba7f81526c6b8f384aefd14f1d77f9de2b8f21bc92c90a3df3a8a",
    BASE + "/archive_1323_1352/1351/check.py": "c21afd00799b1773d29a20895edf4a6a48f82b183a6cfbe1b0347085da6e7853",
    BASE + "/archive_1323_1352/1351/proof.md": "465c6b2669b277cb6066c6f9cfbff37075fa7289bf63c94727d8775a5e3e04d9",
    BASE + "/archive_653_701/research_note_666.md": "4b4994593c8cca012a67c8a6030297f963452b0e95aa5513ac4d9e5ef9170e82",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def local_path(root, name):
    parts = PurePosixPath(name)
    require(bool(parts.parts) and parts.as_posix() == name and not parts.is_absolute()
            and ".." not in parts.parts and ":" not in name and "\\" not in name,
            "Invalid project-relative name: " + name)
    path = root / parts
    require(path.resolve().is_relative_to(root), "Path escapes the project: " + name)
    return path


def history(manifest):
    if manifest is None and not MANIFEST.exists():
        live = ROOT / OLD
        require(live.is_dir(), "Closure manifest not found: " + str(MANIFEST))
        files = {file.relative_to(ROOT).as_posix(): file.read_bytes()
                 for file in live.rglob("*") if file.is_file()
                 and "__pycache__" not in file.parts and file.suffix != ".pyc"}
        return files, "pre-closure"
    path = manifest or MANIFEST
    raw = path.read_bytes()
    migration = json.loads(raw)
    require(migration.get("schema") == "research_sm_archive_closure_v1"
            and (migration.get("first_round"), migration.get("accepted_through")) == (1364, 1815),
            "Unsupported Standard Model closure boundary")
    files, destinations = {}, set()
    for item in migration["entries"]:
        name, new = item["original"], item["destination"]
        if not name.startswith(OLD + "/"):
            continue
        local_path(ROOT, name)
        destination = local_path(ROOT, new)
        require(name not in files and destination not in destinations, "Duplicate migration: " + name)
        destinations.add(destination)
        if name == FOLLOWUP_NOTE and FOLLOWUP_MANIFEST.is_file():
            # Closing 1353–1363 subsequently changed one link in this bound document.
            # Reverse that newer edit before checking the first migration's hashes.
            from research_sm_1353_closure import original_1353_bytes
            current = original_1353_bytes(destination)
        else:
            current = destination.read_bytes()
        if item["original_sha256"] == item["current_sha256"]:
            require(not item.get("edits") and len(current) == item["current_bytes"] == item["original_bytes"]
                    and digest(current) == item["current_sha256"], "Changed migration input: " + new)
            files[name] = current
        else:
            # This shared decoder checks both lengths and both digests, preserving BOM and newlines.
            files[name] = relocation_original_bytes(current, item)
        if name == NAV:
            require(new == migration["frozen_navigation"], "Incorrect frozen navigation destination")
    require(NAV in files and json.loads(files[NAV])["accepted_through"] == migration["accepted_through"],
            "Frozen navigation does not cover the migration boundary")
    return files, digest(raw)


def verify(files):
    nav = json.loads(files[NAV])
    require(nav["schema"] == "research_navigation_1364_v1" and nav["first_round"] == 1364,
            "Unrecognized frozen stage")
    last = nav["accepted_through"]
    require(1364 <= last <= 1815, "Unexpected stage boundary")
    hashes = {name: digest(raw) for name, raw in files.items()}
    registrations = {}
    for name, expected in nav["stage_files_sha256"].items():
        require(hashes.get(name) == expected, "Stage navigation hash mismatch: " + name)
    for n in range(1364, last+1):
        folder = f"{OLD}/{n}/"
        name = folder + "acceptance.json"
        require(name in nav["stage_files_sha256"], "Missing frozen acceptance: " + name)
        row = json.loads(files[name])
        require(row["round"] == n and row["status"].startswith("ACCEPTED"), "Unexpected status: " + name)
        for evidence, expected in row["sha256"].items():
            require(hashes.get(evidence) == expected, "Scientific evidence mismatch: " + evidence)
        checks = row["executed_checks"]
        require(checks and len({c["script"] for c in checks}) == len(checks),
                "Missing or duplicate check: " + name)
        for check in checks:
            require(re.fullmatch(r"[\w-]+\.py", check["script"])
                    and folder+check["script"] in row["sha256"]
                    and check["arguments"] in ([], ["--check"])
                    and check["exit_code"] == 0
                    and re.fullmatch(r"[0-9a-f]{6}|[0-9a-f]{12}", check["receipt"]),
                    "Unexpected check registration: " + name)
        registrations[n] = checks
    for name, expected in PREVIOUS.items():
        raw = local_path(ROOT, name).read_bytes()
        require(digest(raw) == expected, "Earlier input hash mismatch: " + name)
        files[name] = raw
    return registrations


def materialize(files, root):
    for name, raw in files.items():
        # Names are authenticated stage files or validated migration paths.
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)


def run_check(root, n, check, result_receipts):
    folder = root / OLD / str(n)
    token = check["receipt"]
    argv = [sys.executable, "-B", str(folder/check["script"]), *check["arguments"]]
    env = {**os.environ, "PYTHONOPTIMIZE": "0", "PYTHONIOENCODING": "utf-8",
           "PYTHONDONTWRITEBYTECODE": "1"}
    env.pop("PYTHONPATH", None)
    # A fresh interpreter per registered entry point keeps all original assertions active.
    run = subprocess.run(argv, cwd=folder, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = run.stdout.decode("utf8", errors="replace")
    receipt_kind = None
    if run.returncode == 0:
        if re.search(r"(?<![0-9a-fA-F])" + token + r"(?![0-9a-fA-F])", output):
            receipt_kind = "printed_by_check"
        elif len(token) == 12:
            competing = re.search(r"(?<![0-9a-fA-F])[0-9a-fA-F]{12}(?![0-9a-fA-F])", output)
            if not competing and token in result_receipts:
                receipt_kind = "verified_result_digest"
        else:
            # Many historical six-hex strings are host execution IDs. Fresh processes
            # cannot regenerate those IDs; the original acceptance has been hash-verified.
            receipt_kind = "historical_execution_reference"
    result = dict(round=n, script=check["script"], arguments=check["arguments"],
                  registered_receipt=token, receipt_kind=receipt_kind, exit_code=run.returncode,
                  output_sha256=digest(run.stdout), passed=bool(receipt_kind))
    if not result["passed"]:
        result["diagnostic"] = output[-4500:] or "No program output"
    return result


def select_rounds(tokens, all_rounds, registered):
    if all_rounds:
        require(not tokens, "Use either round selections or --all")
        return sorted(registered)
    selected = set()
    for token in tokens:
        require(re.fullmatch(r"[0-9]{4}(?:-[0-9]{4})?", token), "Invalid round selection: " + token)
        pieces = token.split("-")
        first, last = int(pieces[0]), int(pieces[-1])
        require(first <= last and first in registered and last in registered,
                "Outside accepted rounds: " + token)
        selected.update(range(first, last+1))
    require(selected, "Supply a round, inclusive range, or --all")
    return sorted(selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rounds", nargs="*", metavar="ROUND|FIRST-LAST")
    parser.add_argument("--all", action="store_true", help="Select all registered rounds")
    parser.add_argument("--verify-only", action="store_true", help="Verify frozen bytes without executing checks")
    parser.add_argument("--jobs", type=int, default=1, help="Number of concurrent Python checks")
    parser.add_argument("--manifest", type=Path, help="Alternative closure manifest for migration verification")
    parser.add_argument("--report", type=Path, help="Write per-check receipts to a new JSON file")
    args = parser.parse_args()
    require(args.jobs >= 1, "--jobs must be positive")
    if args.report:
        require(not args.report.exists(), "The report already exists: " + str(args.report))
    files, boundary = history(args.manifest)
    registrations = verify(files)
    selected = select_rounds(args.rounds, args.all, registrations)
    work = [(n, c) for n in selected for c in registrations[n]]
    results = []
    if not args.verify_only:
        fingerprints = {}
        for n in selected:
            prefix = f"{OLD}/{n}/"
            fingerprints[n] = {digest(data)[:12]
                               for name, raw in files.items()
                               if name.startswith(prefix) and name.endswith("results.json")
                               for data in (raw, raw.replace(b"\r\n", b"\n"))}
        with tempfile.TemporaryDirectory(prefix="sm-replay-") as name:
            view = Path(name)
            materialize(files, view)
            def replay(item):
                n, check = item
                return run_check(view, n, check, fingerprints[n])
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                for result in pool.map(replay, work):
                    results.append(result)
                    if not result["passed"]:
                        print(json.dumps(result, ensure_ascii=False), file=sys.stderr, flush=True)
            for logical, raw in files.items():
                require((view / logical).read_bytes() == raw,
                        "Registered check changed its historical evidence: " + logical)
    failures = sum(not row["passed"] for row in results)
    summary = dict(status="FAIL" if failures else "PASS", mode="verify" if args.verify_only else "execute",
                   navigation_sha256=digest(files[NAV]), closure_manifest=boundary,
                   selected_rounds=len(selected), registered_checks=len(work), executed_checks=len(results),
                   receipt_kinds=dict(Counter(r["receipt_kind"] for r in results if r["passed"])),
                   failures=failures)
    if args.report:
        with args.report.open("x", encoding="utf-8") as out:
            json.dump(dict(summary=summary, checks=results), out, indent=2, ensure_ascii=False)
            out.write("\n")
    print(json.dumps(summary, ensure_ascii=False))
    return int(bool(failures))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, ValueError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
