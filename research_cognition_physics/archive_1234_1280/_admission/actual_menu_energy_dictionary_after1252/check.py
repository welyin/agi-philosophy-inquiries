"""Exact, read-only replay of the non-numbered energy dictionary witness.

No ODE, sampling, new apparatus, or independent reproof of the old derivative
bounds. The analytic 1248 bounds and frozen 1247 endpoint are distinct inputs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if p.name == "research_cognition_physics")
STAGE = ROOT / "archive_1234_"
OUTPUT = HERE / "results.json"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build() -> dict:
    sources = {
        "1247/results.json": STAGE / "1247/results.json",
        "1248/proof.md": STAGE / "1248/proof.md",
        "1249/canonical.md": STAGE / "1249/canonical.md",
    }
    old = json.loads(sources["1247/results.json"].read_text(encoding="utf-8"))
    point = old["group_finite_memory"]["alpha1"]["checkpoints"]["2"]
    ulo, uhi = map(F, point["state"]["u"])
    hlo, hhi = map(F, point["state"]["H"])
    plo, phi = map(F, point["energy"]["pressure"])
    dt, volume = F(1, 1000), F(1000)
    # These two all-window analytic bounds are inherited from 1248 proof §2.
    abs_Hprime, abs_Pprime = F(9, 40), F(1)
    assert 0 < ulo <= uhi < F(3, 8)
    assert hlo > F(1, 4) + abs_Hprime * dt
    assert phi < -F(4, 625)
    actual_H_lower = hlo - abs_Hprime * dt
    actual_P_upper = phi + abs_Pprime * dt
    assert actual_H_lower > F(1, 4)
    assert actual_P_upper < -F(27, 5000)
    # H>0 and u'=-Hu imply u(t)<=u(2), independently of sampling.
    a3_lower = F(8, 3) ** 3
    energy_rate_lower = 3 * volume * a3_lower * F(1, 4) * F(27, 5000)
    energy_increase_lower = dt * energy_rate_lower
    assert energy_rate_lower == F(384, 5)
    assert energy_increase_lower == F(48, 625)
    assert energy_increase_lower > 0
    return {
        "status": "PASS_EXACT_NON_NUMBERED_ENERGY_DICTIONARY",
        "new_science_groups": 0,
        "fixed_certificate": {
            "alpha": 1, "kappa": "3", "coordinate_volume": str(volume),
            "start": "2", "end": str(F(2) + dt), "duration": str(dt),
            "window_and_coarsenings_are_fixed_certificate_choices": True,
        },
        "frozen_endpoint": {
            "u": [str(ulo), str(uhi)], "H": [str(hlo), str(hhi)],
            "P": [str(plo), str(phi)],
        },
        "analytic_inputs_from_1248_section2": {
            "abs_Hprime_lt": str(abs_Hprime),
            "abs_Pprime_lt": str(abs_Pprime),
            "u_equation": "u'=-H*u",
            "energy_equation": "E_mat'=-3*V*a^3*H*P",
            "these_bounds_are_inherited_proofs_not_numerical_sampling": True,
        },
        "strict_window_bounds": {
            "derived_H_lower": str(actual_H_lower),
            "derived_P_upper": str(actual_P_upper),
            "coarsened_H_gt": "1/4", "coarsened_u_lt": "3/8",
            "coarsened_P_lt": "-27/5000",
            "E_mat_prime_gt": str(energy_rate_lower),
            "E_mat_increase_gt": str(energy_increase_lower),
        },
        "analytic_contrast": {
            "hold_pulses": "w=f=0",
            "J_prime": "0", "H_red": "-V*J",
            "H_red_change": "0",
            "fixed_affine_function_of_H_red_change": "0",
            "rules_out_only_fixed_affine_identification_with_material_energy": True,
            "does_not_rule_out_relational_or_gravitational_parent": True,
        },
        "source_sha256": {name: sha(path) for name, path in sources.items()},
        "code_sha256": sha(Path(__file__).resolve()),
        "scope": "Frozen classical same-history identity check; no physical readout or TR realization certified.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="First save only; exclusive create")
    args = parser.parse_args()
    result = build()
    if args.write:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("SAVED_ONCE_EXACT_NON_NUMBERED_ENERGY_DICTIONARY")
    else:
        saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
        assert result == saved, "Frozen result or source identity mismatch"
        print("READ_ONLY_MATCH_EXACT_NON_NUMBERED_ENERGY_DICTIONARY")


if __name__ == "__main__":
    main()
