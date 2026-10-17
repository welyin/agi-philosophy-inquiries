"""Exact initial-slice screen; no ODE/PDE and no new numbered science group.

Default: recompute in memory and compare results.json without writing.
--write: create results.json exclusively on its first save.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent.parent
RESULTS = HERE / "results.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


# Sparse commutative polynomials over Q. Monomials are sorted variable tuples.
def constant(x: F | int):
    x = F(x)
    return {(): x} if x else {}


def variable(name: str):
    return {(name,): F(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, F(0)) + coefficient
    return {m: c for m, c in out.items() if c}


def scale(c, poly):
    return {m: F(c) * v for m, v in poly.items() if F(c) * v}


def multiply(*polys):
    out = constant(1)
    for poly in polys:
        new = {}
        for m1, c1 in out.items():
            for m2, c2 in poly.items():
                monomial = tuple(sorted(m1 + m2))
                new[monomial] = new.get(monomial, F(0)) + c1 * c2
        out = {m: c for m, c in new.items() if c}
    return out


def equal(left, right):
    return not add(left, scale(-1, right))


def calculate():
    old_path = STAGE / "1247" / "results.json"
    old = json.loads(old_path.read_text(encoding="utf-8"))
    branch = old["group_finite_memory"]["alpha1"]
    tube = branch["all_time_tube"]
    positivity = branch["positive_all_time"]
    u_lo = F(tube["u"][0])
    J_lo = F(tube["J"][0])
    X_abs = max(abs(F(x)) for x in tube["X"])
    assert u_lo > F(1, 3)
    assert J_lo > F(49, 500)
    assert X_abs < F(3, 2)
    assert F(positivity["M_c"][0]) > 0

    t_star, g, y = F(3, 2), F(1, 10), F(1, 1000)
    assert F(tube["theta"][0]) <= t_star <= F(tube["theta"][1])
    f_star = 64 * (t_star - 1) ** 3 * (2 - t_star) ** 3
    fp_star = 192 * ((t_star - 1) ** 2 * (2 - t_star) ** 3
                     - (t_star - 1) ** 3 * (2 - t_star) ** 2)
    assert f_star == 1 and fp_star == 0
    assert t_star > 1  # Original w has zero extension beyond [0,1].

    J_change_abs_upper = g * F(3, 2) * y
    J_strict_lower = F(49, 500) - J_change_abs_upper
    assert J_strict_lower == F(1957, 20000) > 0
    assert J_lo - g * X_abs * y > J_strict_lower

    nu, X, V, Y, J, d, rest, kappa = (
        variable(name) for name in ("nu", "X", "V", "Y", "J", "d", "rest", "kappa")
    )
    gf = g * f_star
    M0 = add(rest, J, scale(gf, multiply(X, Y)))
    source0 = scale(gf, multiply(nu, Y))
    source = {}
    k3 = {}
    identities = {}
    for sign in (-1, 1):
        Ys = add(Y, scale(sign, d))
        Js = add(J, scale(-sign * gf, multiply(X, d)))
        Ms = add(rest, Js, scale(gf, multiply(X, Ys)))
        dM = add(Ms, scale(-1, M0))
        assert not dM
        identities[f"mass_compensation_sign_{sign}"] = True
        # Unchanged U, field data and metric: delta T_c = nu delta M U U.
        assert not multiply(nu, dM)
        assert not scale(F(-1, 3), multiply(kappa, nu, dM))
        identities[f"initial_T_and_constraints_sign_{sign}"] = True

        source[sign] = scale(gf, multiply(nu, Ys))
        dX2 = scale(-sign * gf, multiply(nu, d))
        assert equal(dX2, scale(-1, add(source[sign], scale(-1, source0))))
        k3[sign] = scale(-gf, dX2)
        assert equal(k3[sign], scale(sign * gf * gf, multiply(nu, d)))
        identities[f"X_second_K_third_jet_sign_{sign}"] = True

        dust_energy_prime = scale(sign * gf, multiply(nu, V, d))
        scalar_energy_prime = multiply(V, dX2)
        assert not add(dust_energy_prime, scalar_energy_prime)
        identities[f"complete_initial_energy_exchange_sign_{sign}"] = True

    assert equal(add(source[1], scale(-1, source[-1])),
                 scale(2 * gf, multiply(nu, d)))
    assert equal(add(k3[1], scale(-1, k3[-1])),
                 scale(2 * gf * gf, multiply(nu, d)))
    identities["two_branch_source_gap"] = True
    identities["two_branch_K_third_gap"] = True

    coarse_nu_lower = F(1, 27)
    source_gap = 2 * coarse_nu_lower * gf * y
    k3_gap = 2 * coarse_nu_lower * gf * gf * y
    assert source_gap == F(2, 270000)
    assert k3_gap == F(2, 2700000)
    assert 2 * u_lo ** 3 * gf * y > source_gap
    assert 2 * u_lo ** 3 * gf * gf * y > k3_gap
    # Triangle inequality: max(|estimate-a|, |estimate-b|) >= |a-b|/2.
    assert source_gap / 2 == F(1, 270000)
    assert k3_gap / 2 == F(1, 2700000)

    sources = [
        STAGE / "1245" / "proof.md",
        STAGE / "1247" / "proof.md",
        old_path,
        Path(__file__).resolve(),
    ]
    return {
        "status": "PASS",
        "kind": "non_numbered_exact_initial_slice_screen",
        "new_science_groups": 0,
        "arithmetic": "fractions.Fraction and sparse polynomials over Q",
        "old_process": "1247 alpha1 homogeneous slice at theta=3/2",
        "fixed_inputs": {
            "theta": rational(t_star), "f": rational(f_star),
            "f_prime": rational(fp_star), "w": "0", "g": rational(g),
            "y": rational(y), "chi_abs_upper": "1",
            "chi_label_integral": "0", "witness_point_chi": "1",
            "branches": "Y_plus/minus = Y_star +/- y chi; J_plus/minus = J_star -/+ g X_star y chi",
        },
        "old_tube_checks": {
            "u_enclosure_lower": str(tube["u"][0]),
            "X_enclosure_abs_upper": rational(X_abs),
            "J_enclosure_lower": str(tube["J"][0]),
            "original_Mc_positive_lower": str(positivity["M_c"][0]),
            "strict_coarse_u_lower": "1/3", "strict_coarse_X_abs_upper": "3/2",
            "strict_coarse_J_lower": "49/500",
            "J_change_abs_upper": rational(J_change_abs_upper),
            "new_strict_J_lower": rational(J_strict_lower),
            "Mc_unchanged_exactly": True,
        },
        "polynomial_identities": identities,
        "strict_minimax_bounds": {
            "source_gap_strictly_greater_than": rational(source_gap),
            "source_max_branch_error_strictly_greater_than": rational(source_gap / 2),
            "K_third_gap_strictly_greater_than": rational(k3_gap),
            "K_third_max_branch_error_strictly_greater_than": rational(k3_gap / 2),
            "same_input_to_estimator": "retained zero modes and full initial T, not complete relative data",
        },
        "limits": [
            "Initial-slice source/jet nonclosure, not a finite terminal K separation.",
            "No PDE/ODE integration or explicit future Taylor remainder.",
            "Different declared classical preparations, not an implemented quantum preparation permission.",
            "Does not refute the original known homogeneous input family or any whole cognition programme.",
            "Initial constraints follow unchanged full T and geometry in the original Lagrangian-label dictionary.",
            "Smooth zero-mean chi and pointwise source/jet interpretation are supplied by the analytical screen.",
        ],
        "source_sha256": {p.relative_to(STAGE).as_posix(): sha(p) for p in sources},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusive first save only")
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with RESULTS.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        mode = "FIRST_SAVE_EXCLUSIVE"
    else:
        saved = json.loads(RESULTS.read_text(encoding="utf-8"))
        assert result == saved, "Saved initial-slice certificate differs from read-only recomputation"
        mode = "READ_ONLY_MATCH"
    print(json.dumps({"status": "PASS", "mode": mode,
                      "exact_identities": len(result["polynomial_identities"]),
                      "new_science_groups": 0}, ensure_ascii=False))


if __name__ == "__main__":
    main()
