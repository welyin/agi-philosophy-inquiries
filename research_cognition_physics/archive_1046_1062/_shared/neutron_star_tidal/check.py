"""Exact convention checks only; no EOS solver, waveform fit, or data analysis."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def weighted_coefficients(m1, m2):
    total = m1 + m2
    return (F(16, 13) * (m1 + 12*m2) * m1**4 / total**5,
            F(16, 13) * (m2 + 12*m1) * m2**4 / total**5)


def calculate():
    # Rational substitutions test identities; these are not stellar solutions.
    m1, m2 = F(7, 5), F(6, 5)
    total, mu = m1 + m2, m1*m2/(m1+m2)
    eta = mu/total
    weights = weighted_coefficients(m1, m2)
    old_dimensional = ((m1 + 12*m2)*m1**4/26,
                       (m2 + 12*m1)*m2**4/26)
    assert weights == tuple(32*x/total**5 for x in old_dimensional)
    direct_phase = tuple(-F(9, 16)*26*x/(mu*total**4)
                         for x in old_dimensional)
    new_phase = tuple(-F(117, 256)*x/eta for x in weights)
    assert direct_phase == new_phase
    assert -F(117, 256) == F(3, 128)*(-F(39, 2))
    assert weights == weighted_coefficients(3*m1, 3*m2)
    assert weights == weighted_coefficients(m2, m1)[::-1]
    assert weighted_coefficients(F(1), F(1)) == (F(1, 2), F(1, 2))
    # Dimension vectors have order (length, mass, time).
    grav, light, mass, radius = (3, -1, -2), (1, 0, -1), (0, 1, 0), (1, 0, 0)
    lam = tuple(5*r-g for r, g in zip(radius, grav))
    geometric = tuple(g+l for g, l in zip(grav, lam))
    dimensionless = tuple(10*c+l-4*g-5*m for c, l, g, m
                         in zip(light, lam, grav, mass))
    assert lam == (2, 1, 2)
    assert geometric == (5, 0, 0)
    assert dimensionless == (0, 0, 0)
    return {
        "scope": "count0 exact conventions; no EOS, waveform, or empirical recomputation",
        "rational_test_masses": [str(m1), str(m2)],
        "Lambda_weights": [str(x) for x in weights],
        "phase_coefficients_of_Lambda_times_v5": [str(x) for x in new_phase],
        "lambda_tilde_to_Lambda_tilde_factor": "32/M^5",
        "equal_mass_weights": ["1/2", "1/2"],
        "phase_prefactor_identity": "-117/256 = (3/128)*(-39/2)",
        "SI_dimensions_length_mass_time": {
            "lambda_SI": list(lam), "G_lambda_SI": list(geometric),
            "Lambda": list(dimensionless)},
        "all_assertions_passed": True,
    }


if __name__ == "__main__":
    result = calculate()
    target = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments for read-only checking, or --save-exclusive once.")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print("PASS: exact units, tidal weights and phase conventions; no EOS or data fit.")
