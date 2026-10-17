"""Finite central-energy-window identities. Non-numbered diagnostic, not a GR certificate."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics").is_dir() and (p / "scripts").is_dir())
RESULT = HERE.with_suffix(".json")
MEMO = HERE.with_suffix(".md")
BASE = ROOT / "research_cognition_physics/archive_1234_/_admission/positive_observer_entropy_after1270/finite_observer_entropy_candidate_after1270.md"
BETA = math.log(2.0)
ENERGIES = (F(0), F(1))
Z = F(3, 2)


def rational(x):
    return {"numerator": x.numerator, "denominator": x.denominator}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tr_corner(n):
    return sum((F(1, 2 ** int(e)) - F(1, 2 ** int(n)))
               for e in ENERGIES if e < n) / Z


def ev2(a, b, c):
    tr = a + b
    rad = math.sqrt((a - b) ** 2 + 4 * abs(c) ** 2)
    return (max(0.0, (tr + rad) / 2), max(0.0, (tr - rad) / 2))


def spectral_integrals(rho, width, upper):
    """Exact piece intervals; only logarithms/eigenvalues evaluated in float."""
    points = sorted(set([F(0), upper] + [
        p for e in ENERGIES for p in (e, e + width) if 0 < p < upper
    ]))
    p = F(0)
    ey = F(0)
    entropy_d = 0.0
    h_y = 0.0
    pieces = []
    for lo, hi in zip(points[:-1], points[1:]):
        mid = (lo + hi) / 2
        active = [e <= mid <= e + width for e in ENERGIES]
        a = rho[0] / width if active[0] else F(0)
        b = rho[1] / width if active[1] else F(0)
        c = rho[2] / width if all(active) else F(0)
        eigenvalues = ev2(float(a), float(b), float(c))
        density = a + b
        dy = hi - lo
        p += dy * density
        ey += (hi * hi - lo * lo) * density / 2
        entropy_d -= float(dy) * sum(x * math.log(x) for x in eigenvalues if x > 1e-15)
        if density:
            h_y -= float(dy * density) * math.log(float(density))
        pieces.append((float(lo), float(hi), eigenvalues))
    return p, ey, entropy_d, h_y, pieces


def entropy_in_corner(rho, width, n):
    p, ey, ent, h_y, pieces = spectral_integrals(rho, width, n)
    t = tr_corner(n)
    c = float(t * Z / p) / BETA
    direct = 0.0
    for lo, hi, eigs in pieces:
        for eigen in eigs:
            if eigen > 1e-15:
                direct -= eigen / float(p) * (
                    (hi - lo) * math.log(c * eigen)
                    + BETA * (hi * hi - lo * lo) / 2
                )
    integral_old = ent - float(p) * math.log(float(Z) / BETA) - BETA * float(ey)
    rhs = integral_old / float(p) + math.log(float(p)) - math.log(float(t))
    assert abs(direct - rhs) < 2e-13
    return p, t, direct, abs(direct - rhs), h_y


def output():
    t3 = tr_corner(F(3))
    assert t3 == F(5, 6)
    s0, s1 = F(2, 3), F(1, 3)
    c0, c1 = F(7, 8), F(3, 4)
    new_reference = (s0 * c0 / t3, s1 * c1 / t3)
    assert new_reference == (F(7, 10), F(3, 10))
    joint_j0_q_gt_2 = s0 * (F(1, 4) - F(1, 8)) / t3
    assert joint_j0_q_gt_2 == F(1, 10)
    correlation = joint_j0_q_gt_2 - new_reference[0] * joint_j0_q_gt_2
    assert correlation == F(3, 100)
    pure = (F(1, 2), F(1, 2), F(1, 2))
    mixed = (F(3, 5), F(2, 5), F(3, 20))
    diagonal = (F(2, 5), F(3, 5), F(0))
    cases = []
    max_error = 0.0
    for name, rho in (("pure_coherent", pure), ("mixed_noncommuting", mixed), ("diagonal", diagonal)):
        for n in (F(2), F(3), F(4)):
            p, t, entropy, err, _ = entropy_in_corner(rho, F(2), n)
            cases.append({"state": name, "N": int(n), "p": rational(p), "trace": rational(t), "entropy": entropy, "identity_residual": err})
            max_error = max(max_error, err)
            if n >= 3:
                assert p == 1
    p3, _, sc3, _, _ = entropy_in_corner(pure, F(2), F(3))
    s_old = math.log(2 * BETA / 3)
    assert abs(sc3 - math.log(4 * BETA / 5)) < 2e-13
    assert abs(sc3 - s_old + math.log(float(t3))) < 2e-13
    p2, _, _, _, _ = entropy_in_corner(pure, F(2), F(2))
    assert p2 == F(3, 4)
    # Bare input after the central projector: min of the two allowed Q windows.
    gram = ((F(1), F(1, 2)), (F(1, 2), F(1, 2)))
    clipped_material = (pure[0] * gram[0][0] / p2, pure[1] * gram[1][1] / p2, pure[2] * gram[0][1] / p2)
    assert clipped_material == (F(2, 3), F(1, 3), F(1, 3))
    # Width below energy gap gives disjoint classical Y supports for a pure coherent input.
    pd, _, _, hy, _ = spectral_integrals(pure, F(1, 2), F(3, 2))
    assert pd == 1
    exact_gap_diagnostic = hy - math.log(0.5)
    assert abs(exact_gap_diagnostic - math.log(2.0)) < 2e-13
    # On every energy pair, N >= B+Emax makes the full MR/O projector identity.
    support_gram = [[min(F(2), F(3)-ej, F(3)-ek) / 2 for ek in ENERGIES] for ej in ENERGIES]
    assert all(x == 1 for row in support_gram for x in row)
    return {
        "status": "PASS_UNNUMBERED_CENTRAL_WINDOW_ONLY",
        "scope": "Finite algebra/window diagnostic; analytic proof carries all-state claims. No geometry, actual preparation, or Gamma recomputation.",
        "source_bindings": {"code": sha(HERE), "memo": sha(MEMO), "published_observer_lemma": sha(BASE)},
        "exact": {
            "trace_N3": rational(t3),
            "conditioned_material_reference": [rational(x) for x in new_reference],
            "reference_energy_correlation": rational(correlation),
            "clipped_work_probability_N2": rational(p2),
            "clipped_material_00_11_01": [rational(x) for x in clipped_material],
            "all_state_support_gram_is_identity_multiplier": True,
        },
        "diagnostics": {
            "cases": cases,
            "max_entropy_identity_residual": max_error,
            "supported_pure_old_entropy": s_old,
            "supported_pure_corner_entropy": sc3,
            "disjoint_support_gap": exact_gap_diagnostic,
            "normalization_shift": -math.log(float(t3)),
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = output()
    if args.write:
        with RESULT.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(actual, out, ensure_ascii=False, indent=2, sort_keys=True)
            out.write("\n")
        print("PASS_FIRST_WRITE", RESULT)
    else:
        expected = json.loads(RESULT.read_text(encoding="utf-8"))
        assert expected == actual, "Saved diagnostic or source binding mismatch"
        print("READ_ONLY_MATCH", len(actual["diagnostics"]["cases"]), "finite cases")

