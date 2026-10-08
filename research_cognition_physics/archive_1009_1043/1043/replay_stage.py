"""Read-only author-verifier replay for rounds 1009--1042.

Only --record writes, exclusively creating this round's aggregate JSON.  No
write option is passed to any historical script.  The default repeats the
read-only regression and compares stable evidence with the saved aggregate.
Passing regression is not a semantic claim that the research goal is complete.
"""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
PROJECT = STAGE.parent.parent
OUT = HERE / "stage_replay_results.json"
ROUNDS = tuple(range(1009, 1043))
WORKERS = 3
TIMEOUT_SECONDS = 600


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(PROJECT).as_posix()


def science_paths(verifier: Path) -> list[Path]:
    tree = ast.parse(verifier.read_text("utf-8-sig"))
    paths = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                candidate = verifier.parent / (alias.name + ".py")
                if candidate.is_file():
                    paths.add(candidate)
    subprocess_science = {
        1040: "nonlinear_common_cone.py",
        1041: "common_relay_time.py",
        1042: "complementary_entropy.py",
    }
    n = int(verifier.parent.name)
    if n in subprocess_science:
        paths.add(verifier.parent / subprocess_science[n])
    return sorted(paths)


def write_preflight(path: Path) -> dict:
    """Inventory known filesystem writes and verify their CLI write guards.

    This is an audited source check, not a general Python sandbox or a proof
    about arbitrary dynamically imported code.  Post-run snapshots separately
    check that all files in these historical round directories stayed intact.
    """
    tree = ast.parse(path.read_text("utf-8-sig"))
    parents = {child: node for node in ast.walk(tree)
               for child in ast.iter_child_nodes(node)}
    calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = ast.unparse(node.func)
        is_write = name.endswith((".write", ".write_text", ".write_bytes"))
        is_write = is_write or name == "json.dump"
        if name.endswith(".open") or name == "open":
            # All relevant open calls use an explicit literal file mode.
            modes = [a.value for a in node.args if isinstance(a, ast.Constant)
                     and isinstance(a.value, str)]
            modes += [k.value.value for k in node.keywords if k.arg == "mode"
                      and isinstance(k.value, ast.Constant)]
            is_write = is_write or any(m in ("x", "w", "a", "xb", "wb", "ab")
                                       for m in modes)
        if not is_write:
            continue
        guards = []
        ancestor = node
        while ancestor in parents:
            ancestor = parents[ancestor]
            if isinstance(ancestor, ast.If):
                guards.append(ast.unparse(ancestor.test))
        write_guards = ("args.write", "args.write_results", "args.write_receipt")
        guarded = any(g in write_guards for g in guards)
        assert guarded, (relative(path), node.lineno, name, guards)
        calls.append({"line": node.lineno, "call": name,
                      "write_flag_guard": next(g for g in guards if g in write_guards)})
    return {"path": relative(path), "sha256": sha(path),
            "known_filesystem_writes": calls,
            "default_write_flags_false": True}


def historical_snapshot() -> dict[str, str]:
    paths = []
    for n in ROUNDS:
        paths.extend(p for p in (STAGE / str(n)).rglob("*") if p.is_file()
                     and "__pycache__" not in p.parts)
        paths.append(STAGE / f"research_note_{n}.md")
    return {relative(p): sha(p) for p in sorted(set(paths))}


def manifest_digest(snapshot: dict[str, str]) -> str:
    return hashlib.sha256(json.dumps(snapshot, sort_keys=True,
                                    ensure_ascii=False).encode("utf8")).hexdigest()


def receipt_path(n: int) -> Path:
    path = STAGE / str(n) / f"research_round_{n}_checks.json"
    assert path.is_file(), path
    return path


def replay_one(n: int, env: dict[str, str]) -> dict:
    verifier = STAGE / str(n) / f"verify_round{n}.py"
    receipt = receipt_path(n)
    before_script, before_receipt = sha(verifier), sha(receipt)
    command = [sys.executable, "-B", "-X", "utf8", str(verifier)]
    started = time.monotonic()
    try:
        completed = subprocess.run(command, cwd=str(PROJECT), env=env,
                                   capture_output=True, text=True, encoding="utf8",
                                   errors="replace", timeout=TIMEOUT_SECONDS)
        exit_code, stdout, stderr = completed.returncode, completed.stdout, completed.stderr
        timeout = False
    except subprocess.TimeoutExpired as exc:
        exit_code, timeout = None, True
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf8", "replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf8", "replace")
    after_script, after_receipt = sha(verifier), sha(receipt)
    unchanged = before_script == after_script and before_receipt == after_receipt
    parsed = None
    try:
        parsed = json.loads(stdout)
    except (ValueError, TypeError):
        pass
    concise = None
    if isinstance(parsed, dict):
        concise = {k: v for k, v in parsed.items()
                   if isinstance(v, (bool, int, float, str)) and
                   any(s in k.lower() for s in
                       ("round", "pass", "reproduc", "mode", "files", "links", "check"))}
    passed = exit_code == 0 and unchanged and not timeout
    return {
        "round": n, "command": command,
        "verifier": relative(verifier), "verifier_sha256_before": before_script,
        "verifier_sha256_after": after_script,
        "receipt": relative(receipt), "receipt_sha256_before": before_receipt,
        "receipt_sha256_after": after_receipt,
        "exit_code": exit_code, "timeout": timeout,
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "verifier_and_receipt_unchanged": unchanged, "regression_passed": passed,
        "brief_result": concise,
        "stdout": stdout, "stderr": stderr,
        "failure_classification": None if passed else "unclassified_requires_read_only_review",
        "failure_does_not_automatically_mean_scientific_refutation": not passed,
    }


def replay() -> dict:
    preflight = []
    for n in ROUNDS:
        verifier = STAGE / str(n) / f"verify_round{n}.py"
        assert verifier.is_file(), verifier
        receipt_path(n)
        for path in [verifier] + science_paths(verifier):
            preflight.append(write_preflight(path))
    before = historical_snapshot()
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1")
    start = datetime.now(timezone.utc).isoformat()
    started = time.monotonic()
    entries = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(replay_one, n, env): n for n in ROUNDS}
        for future in as_completed(futures):
            entry = future.result()
            entries.append(entry)
            print(json.dumps({"round": entry["round"],
                              "regression_passed": entry["regression_passed"],
                              "exit_code": entry["exit_code"],
                              "elapsed_seconds": entry["elapsed_seconds"]}), flush=True)
    entries.sort(key=lambda x: x["round"])
    after = historical_snapshot()
    changed = [{"path": path, "before": before.get(path), "after": after.get(path)}
               for path in sorted(before.keys() | after.keys())
               if before.get(path) != after.get(path)]
    failures = [entry["round"] for entry in entries if not entry["regression_passed"]]
    return {
        "round": 1043, "kind": "read_only_stage_regression_not_new_science",
        "recorded_at_utc": start, "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "rounds": list(ROUNDS), "verifier_count": len(entries),
        "new_scientific_calibration_groups": 0,
        "cumulative_scientific_calibration_groups_unchanged": 3819,
        "semantic_goal_completion_certified_by_this_replay": False,
        "runtime": sys.executable, "runtime_sha256": sha(Path(sys.executable)),
        "replay_script_sha256": sha(Path(__file__)),
        "max_concurrent_verifiers": WORKERS, "timeout_seconds_per_verifier": TIMEOUT_SECONDS,
        "environment_overrides": {key: env[key] for key in
                                  ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                   "PYTHONDONTWRITEBYTECODE", "PYTHONUTF8")},
        "historical_write_flags_passed": [],
        "default_read_only_preflight": preflight,
        "preflight_scope": "Known filesystem writes in 34 author verifiers and 34 direct science modules are CLI-guarded; manually checked helper entry guards. Not a general Python sandbox.",
        "historical_snapshot_files": len(before),
        "historical_snapshot_sha256_before": manifest_digest(before),
        "historical_snapshot_sha256_after": manifest_digest(after),
        "historical_files_changed": changed,
        "historical_snapshot_scope": "Round 1009--1042 directories and numbered notes; excludes __pycache__ and mutable root navigation.",
        "passed_rounds": [entry["round"] for entry in entries if entry["regression_passed"]],
        "failed_rounds": failures,
        "all_regressions_passed": not failures and not changed,
        "live_navigation_incompatibilities": [],
        "live_navigation_incompatibilities_note": "No historic asset is repaired; any failed verifier requires separate read-only classification before attributing its failure to mutable navigation.",
        "results": entries,
    }


def stable_evidence(result: dict) -> dict:
    return {
        "replay_script_sha256": result["replay_script_sha256"],
        "preflight": result["default_read_only_preflight"],
        "all_regressions_passed": result["all_regressions_passed"],
        "historical_files_changed": result["historical_files_changed"],
        "entries": [{key: row[key] for key in
                     ("round", "verifier_sha256_before", "verifier_sha256_after",
                      "receipt_sha256_before", "receipt_sha256_after", "exit_code",
                      "timeout", "regression_passed")}
                    for row in result["results"]],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true",
                        help="Exclusively create only 1043/stage_replay_results.json")
    args = parser.parse_args()
    if args.record and OUT.exists():
        raise FileExistsError(OUT)
    result = replay()
    if args.record:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    elif OUT.exists():
        saved = json.loads(OUT.read_text("utf8"))
        assert stable_evidence(result) == stable_evidence(saved), "Stable replay evidence differs"
    print(json.dumps({"round": 1043, "all_regressions_passed": result["all_regressions_passed"],
                      "passed": len(result["passed_rounds"]), "failed": result["failed_rounds"],
                      "changed_historical_files": len(result["historical_files_changed"]),
                      "new_scientific_calibration_groups": 0,
                      "semantic_goal_completion_certified": False}, ensure_ascii=False))
    if not result["all_regressions_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
