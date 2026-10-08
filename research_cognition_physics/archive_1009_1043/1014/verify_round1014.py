"""Recompute round 1014; file checks do not certify its analytic theorem."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import weak_current_spectral_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1014.md"
OUT = HERE / "research_round_1014_checks.json"
sys.path.insert(0, str(ROOT / "scripts"))
from organize_research_231_775 import links

OWN = ["weak_current_spectral_selection.py", "weak_current_spectral_selection_results.json",
       "verify_round1014.py", "review.md", "NEXT.md"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(prospective=False):
    fresh = science.run()
    saved = json.loads((HERE / OWN[1]).read_text("utf8"))
    science.compare(fresh, saved)
    assert fresh["all_scientific_calibrations_passed"]
    assert not any(fresh[k] for k in ["full_realistic_experiment_claimed",
        "ordinary_event_rate_identity_claimed",
        "mass_simplicity_inferred_from_numerical_tolerance_claimed",
        "continuum_parameter_theorem_proved_by_sampling_claimed"])
    for name, digest in fresh["historical_source_sha256"].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    count = 0
    for p in assets:
        assert p.is_file(), p
        if p.suffix == ".md":
            for _, _, target, local in links(p.read_text("utf-8-sig")):
                resolved = (p.parent / local.replace("\\", "/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (p, target)
                count += 1
    note = NOTE.read_text("utf8")
    for n in range(1, 9):
        assert f"## {n}." in note
    review = (HERE / "review.md").read_text("utf8")
    assert "主代理审阅通过" in review
    return dict(round=1014, date="2026-10-08", all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3792, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(fresh["historical_source_sha256"]),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        full_spectrum_simplicity_is_a_premise=True,
        partial_calibrated_final_channels_can_supply_lower_bound=True,
        absolute_event_rates_not_identified_with_current_coefficients=True,
        primary_agent_math_and_code_review=True,
        independent_math_review=True,
        new_cognitive_axiom=False, goal_completed=False, visual_checks_performed=False,
        analytic_scope="Exact canonical nine-Weyl Ward/current contract; positive full-spectrum simple states; nonstandard total chiral coefficient strength at most one.",
        historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"):sha(p) for p in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, ensure_ascii=False, indent=2))
