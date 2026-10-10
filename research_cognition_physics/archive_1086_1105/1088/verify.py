"""Read-only check of the declared round 1088 assets and finite witness."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def main():
    record = json.loads((HERE / "acceptance.json").read_text(encoding="utf8"))
    assert record["status"] == "PASS"
    assert record["exact_FUCP_plus_CO_model_certified"] is False
    assert record["joint_SR_goal_proved_or_disproved"] is False
    assert record["independent_review"] == "PASS_WITH_DECLARED_SCOPE"
    links = 0
    for relative, expected in record["assets_sha256"].items():
        path = STAGE / relative
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == expected, relative
        if path.suffix == ".md":
            for target in re.findall(r"\]\(([^)]+)\)", raw.decode("utf8")):
                if target.startswith(("https:", "http:", "#")):
                    continue
                assert (path.parent / target.split("#", 1)[0]).exists(), (relative, target)
                links += 1
    if "--recompute" in sys.argv:
        run = subprocess.run([sys.executable, "-B", "-X", "utf8", str(HERE / "check.py")],
                             capture_output=True, text=True, encoding="utf8")
        assert run.returncode == 0, run.stdout + run.stderr
        assert json.loads(run.stdout)["passed"]
    print(json.dumps({"round": 1088, "passed": True,
                      "scope": record["scope"],
                      "assets": len(record["assets_sha256"]), "local_links": links,
                      "recomputed": "--recompute" in sys.argv,
                      "exact_FUCP_plus_CO_model_certified": False,
                      "joint_SR_goal_proved_or_disproved": False,
                      "goal_status": "active"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
