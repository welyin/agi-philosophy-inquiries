"""Verify the closed 1353-1363 evidence. Add --replay to recompute all 23 results."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from research_sm_1353_closure import original_1353_bytes, physical_1353


def checked_file(name, expected):
    path = physical_1353(name)
    if not path.is_relative_to(ROOT):
        raise ValueError(f"Invalid path: {name}")
    if hashlib.sha256(original_1353_bytes(path)).hexdigest() != expected:
        raise ValueError(f"Changed evidence: {name}")
    return path


def verify(replay):
    manifest = json.loads((STAGE / "_shared" / "stage_acceptance.json").read_bytes())
    if set(manifest["rounds"]) != set(map(str,range(1353,1364))):
        raise ValueError("Unexpected round coverage")
    for name, digest in manifest["stage_documents"].items():
        checked_file(name, digest)
    files = commands = 0
    for number, frozen in manifest["rounds"].items():
        folder = STAGE / number
        record = checked_file((folder / "acceptance.json").relative_to(ROOT), frozen)
        acceptance = json.loads(record.read_bytes())
        if acceptance["round"] != int(number) or not acceptance["status"].startswith("ACCEPTED"):
            raise ValueError(f"Unaccepted round {number}")
        for name,digest in acceptance["sha256"].items():
            checked_file(name, digest)
            files += 1
        for check in acceptance["executed_checks"]:
            command, *args = check["script"].split()
            if Path(command).name != command or not command.endswith(".py") or args not in ([], ["--check"]):
                raise ValueError(f"Unexpected replay entry in {number}")
            if check["exit_code"] != 0:
                raise ValueError(f"Unsuccessful receipt in {number}")
            if replay:
                subprocess.run([sys.executable, "-B", "-X", "utf8", str(folder / command), *args],
                               check=True, cwd=ROOT, stdout=subprocess.DEVNULL)
            commands += 1
    if commands != 23 or files != 101:
        raise ValueError(f"Unexpected evidence totals {commands}/{files}")
    return {"verified_rounds":11,"bound_scientific_files":files,"recorded_checks":commands,
            "replayed_checks":commands if replay else 0}


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--replay",action="store_true")
    args=cli.parse_args()
    print(json.dumps(verify(args.replay)))
