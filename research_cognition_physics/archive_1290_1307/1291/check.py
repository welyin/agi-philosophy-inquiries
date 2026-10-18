"""Round 1291: fixed-temperature heavy KG source baseline.

One author group. All certificate calculations use exact fractions.
Default: read-only recomputation against results.json.
--write: exclusive first creation only. No historical checker is imported.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
ETA = (-1, 1, 1, 1)


def fraction(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def stress_vertices(m, p, ep, q, eq):
    m, ep, eq = map(F, (m, ep, eq))
    p, q = tuple(map(F, p)), tuple(map(F, q))
    assert ep > 0 and eq > 0 and m > 0
    assert ep * ep - sum(x * x for x in p) == m * m
    assert eq * eq - sum(x * x for x in q) == m * m
    pc, qc = (-ep,) + p, (-eq,) + q
    dot = sum(ETA[i] * pc[i] * qc[i] for i in range(4))
    pair, scatter = [], []
    for i in range(4):
        tr, sr = [], []
        for j in range(4):
            sym = (pc[i] * qc[j] + pc[j] * qc[i]) / 2
            metric = F(ETA[i]) if i == j else F(0)
            tr.append(-sym + metric * (dot - m * m) / 2)
            sr.append(sym - metric * (dot + m * m) / 2)
        pair.append(tr)
        scatter.append(sr)
    return pair, scatter


def ward(momentum, vertex):
    return [sum(momentum[i] * vertex[i][j] for i in range(4))
            for j in range(4)]


def vertex_checks():
    # Different non-parallel and reverse rational shells, not one rest-frame row.
    cases = [
        (1, (F(3, 4), 0, 0), F(5, 4), (0, F(4, 3), 0), F(5, 3)),
        (3, (4, 0, 0), 5, (0, 0, 4), 5),
        (2, (0, F(3, 2), 0), F(5, 2), (0, 0, 0), 2),
        (5, (12, 0, 0), 13, (-12, 0, 0), 13),
    ]
    rows = []
    for m, p, ep, q, eq in cases:
        t, s = stress_vertices(m, p, ep, q, eq)
        pp = [F(ep) + F(eq)] + [F(p[i]) + F(q[i]) for i in range(3)]
        pm = [F(ep) - F(eq)] + [F(p[i]) - F(q[i]) for i in range(3)]
        wt, ws = ward(pp, t), ward(pm, s)
        assert wt == ws == [0] * 4
        bound = 3 * F(ep) * F(eq)
        assert all(abs(x) <= bound for v in (t, s) for row in v for x in row)
        opposite, _ = stress_vertices(m, p, ep, tuple(-F(x) for x in p), ep)
        assert opposite[0][0] == 0
        _, diagonal = stress_vertices(m, p, ep, p, ep)
        assert diagonal[0][0] == F(ep) ** 2
        # Delete the explicit mass term from each complete stress vertex.
        badt = [[t[i][j] + (ETA[i] * F(m) ** 2 / 2 if i == j else 0)
                 for j in range(4)] for i in range(4)]
        bads = [[s[i][j] + (ETA[i] * F(m) ** 2 / 2 if i == j else 0)
                 for j in range(4)] for i in range(4)]
        badwt, badws = ward(pp, badt), ward(pm, bads)
        assert any(x != 0 for x in badwt)
        assert any(x != 0 for x in badws)
        rows.append({
            "mass": fraction(m),
            "p": [fraction(x) for x in p],
            "q": [fraction(x) for x in q],
            "energies": [fraction(ep), fraction(eq)],
            "pair_ward": [fraction(x) for x in wt],
            "scatter_ward": [fraction(x) for x in ws],
            "omitted_mass_pair_ward": [fraction(x) for x in badwt],
            "omitted_mass_scatter_ward": [fraction(x) for x in badws],
            "opposite_pair_T00": fraction(opposite[0][0]),
            "zero_transfer_scatter_T00": fraction(diagonal[0][0]),
        })
    return {"rows": rows, "finite_samples_do_not_replace_continuum_Ward_proof": True}


def wick_checks():
    rows = []
    amp, b = F(3, 5), F(7, 4)
    for q in (F(0), F(1, 5), F(1, 3), F(3, 5)):
        # The geometric distribution has factorial moments k! [q/(1-q)]^k.
        n = q / (1 - q)
        first = n
        second = factorial(2) * n * n + first
        # Q = amp*(a^2+a^dagger^2)+b*N. Compute <Q^2> by its
        # transitions out of |j>: amp^2*j(j-1), amp^2*(j+1)(j+2), b^2*j^2.
        second_from_number_transitions = (
            amp * amp * ((second - first) + (second + 3 * first + 2))
            + b * b * second)
        direct_variance = second_from_number_transitions - (b * first) ** 2
        pair = 2 * amp * amp * ((1 + n) ** 2 + n * n)
        scatter = b * b * n * (1 + n)
        assert direct_variance == pair + scatter
        missing_vacuum = 2 * amp * amp * (2 * n + 2 * n * n) + scatter
        missing_bose_factor = pair / 2 + scatter
        assert direct_variance - missing_vacuum == 2 * amp * amp
        assert direct_variance > missing_bose_factor
        if q > 0:
            assert scatter > 0 and direct_variance > pair
        rows.append({
            "geometric_q": fraction(q), "occupation": fraction(n),
            "variance_from_number_transitions": fraction(direct_variance),
            "pair": fraction(pair), "scatter": fraction(scatter),
            "error_omitting_vacuum": fraction(direct_variance - missing_vacuum),
            "error_omitting_pair_factor_two":
                fraction(direct_variance - missing_bose_factor),
            "error_omitting_scatter": fraction(scatter),
        })
    # N evolves at zero frequency but its thermal variance is nonzero.
    n = F(1, 3)
    zero_frequency_variance = n * (1 + n)
    assert zero_frequency_variance == F(4, 9)
    return {
        "amplitude": fraction(amp), "number_coefficient": fraction(b),
        "rows": rows,
        "number_preserving_zero_frequency": "0/1",
        "zero_frequency_variance_at_n_one_third": fraction(zero_frequency_variance),
        "hard_2M_gap_for_scattering_rejected": True,
    }


def gamma_half_coefficient(k):
    """Gamma(k+1/2) divided by sqrt(pi), exactly."""
    return F(factorial(2 * k), 4 ** k * factorial(k))


def p_coefficients(n, beta):
    # Ascending polynomial coefficients, obtained recursively by IBP:
    # P_n(M)=M^n/beta + n P_(n-1)(M)/beta.
    beta = F(beta)
    out = [1 / beta]
    for degree in range(1, n + 1):
        out = [degree * x / beta for x in out] + [1 / beta]
    return out


def evaluate(poly, x):
    total = F(0)
    for c in reversed(poly):
        total = total * x + c
    return total


def integral_checks():
    # int_0^inf x^2/(1+x^2)^4 dx = (1/2) B(3/2,5/2).
    radial_j0_over_pi = (F(1, 2) * gamma_half_coefficient(1)
                        * gamma_half_coefficient(2) / factorial(3))
    assert radial_j0_over_pi == F(1, 32)
    # The angular/measure factor is 1/(2*pi^2).
    j0_times_pi_l3 = radial_j0_over_pi / 2
    radial_j1 = F(1, 2) * F(factorial(1) * factorial(1), factorial(3))
    j1_times_pi2_l4 = radial_j1 / 2
    assert j0_times_pi_l3 == F(1, 64)
    assert j1_times_pi2_l4 == F(1, 24)
    beta44 = sum(F((-1) ** j * comb(3, j), 4 + j) for j in range(4))
    assert beta44 == F(1, 140)
    rows = []
    for beta in (F(1), F(2, 3), F(3, 2)):
        for n in (3, 4):
            coeffs = p_coefficients(n, beta)
            # d[-exp(-beta M)P_n]/dM = exp(-beta M)*M^n.
            defect = [beta * coeffs[j]
                      - ((j + 1) * coeffs[j + 1] if j < n else 0)
                      for j in range(n + 1)]
            assert defect == [0] * n + [1]
            rows.append({"beta": fraction(beta), "degree": n,
                         "coefficients_ascending": [fraction(x) for x in coeffs],
                         "beta_P_minus_P_prime": [fraction(x) for x in defect]})
    assert p_coefficients(3, 1) == [6, 6, 3, 1]
    assert p_coefficients(4, 1) == [24, 24, 12, 4, 1]
    return {"J0_times_pi_l3": fraction(j0_times_pi_l3),
            "J1_times_pi2_l4": fraction(j1_times_pi2_l4),
            "vacuum_beta_4_4": fraction(beta44), "IBP_rows": rows}


def compact_checks():
    width = F(1, 8)
    derivative_rows = []
    for k in range(7):
        atoms = [(width * F(2 * j - k, 2),
                  F((-1) ** j * comb(k, j)) / width ** k)
                 for j in range(k + 1)]
        tv = sum(abs(c) for _, c in atoms)
        assert tv == 16 ** k
        moments = [sum(c * x ** degree for x, c in atoms)
                   for degree in range(k + 1)]
        assert moments[:k] == [0] * k
        assert moments[k] == (-1) ** k * factorial(k)
        derivative_rows.append({
            "order": k, "signed_measure_total_variation": fraction(tv),
            "moments_0_to_order": [fraction(x) for x in moments],
        })
    # Expand (1-Delta)^2 independently as a sparse differential polynomial.
    d = {(0, 0, 0): 1, (2, 0, 0): -1, (0, 2, 0): -1, (0, 0, 2): -1}
    squared = {}
    for a, ca in d.items():
        for b, cb in d.items():
            key = tuple(x + y for x, y in zip(a, b))
            squared[key] = squared.get(key, 0) + ca * cb
    c_bound = sum(abs(c) * 16 ** sum(alpha) for alpha, c in squared.items())
    assert c_bound == 1 + 6 * 16 ** 2 + 9 * 16 ** 4 == 769 ** 2
    time_radius = 6 * width / 2 + F(1, 8)
    space_radius = 4 * width / 2 + F(1, 8)
    assert time_radius == F(1, 2) and space_radius == F(3, 8)

    mass = 2 ** 20
    q_upper = F(1, 2 ** 256)  # e > 2, mass >= 256: never compute exp(-mass).
    assert mass >= 256
    n_upper = 2 * q_upper   # d_M > 1/2
    pair_multiplier = 1 + 2 * n_upper + 2 * n_upper ** 2
    assert pair_multiplier < 2
    p3 = evaluate(p_coefficients(3, 1), mass)
    p4 = evaluate(p_coefficients(4, 1), mass)
    mean_upper = q_upper * p3 / 9  # pi>3 and d_M>1/2
    scattering_upper = 2 * c_bound ** 2 * q_upper * (p4 / 192 + p3 / 216)
    vacuum_upper = F(9, 1120 * 3 ** 4 * 4) * 16 ** 12 / (2 * mass) ** 4
    assert vacuum_upper == F(1, 40320 * 2 ** 36)
    pair_upper = 2 * vacuum_upper
    total_upper = pair_upper + scattering_upper
    assert mean_upper < F(1, 10 ** 59)
    assert scattering_upper < F(1, 10 ** 43)
    assert total_upper < F(1, 10 ** 15)
    assert total_upper + mean_upper ** 2 < F(1, 10 ** 15)
    return {
        "tensor_test": "f00=g(t)*v(x1)*v(x2)*v(x3); all others zero",
        "fixed_beta": "1/1", "fixed_spatial_ell": "1/1",
        "mass": mass, "box_width": fraction(width),
        "time_box_count": 6, "space_boxes_each_axis": 4,
        "fixed_smooth_positive_unit_mass_mollifier_radius": "1/8",
        "time_support_radius": fraction(time_radius),
        "space_support_radius_each_axis": fraction(space_radius),
        "B0": "1/1", "B6_upper": str(16 ** 6),
        "C1_upper": str(c_bound),
        "derivative_measure_checks": derivative_rows,
        "spatial_operator_terms": [
            {"multiindex": list(a), "coefficient": c} for a, c in sorted(squared.items())
        ],
        "exponential_upper": fraction(q_upper),
        "exponential_bound_reason": "M>=256 and e>2; no floating underflow",
        "d_M_strict_lower": "1/2",
        "n_M_strict_upper": fraction(n_upper),
        "pair_multiplier_bound_below_two": fraction(pair_multiplier),
        "P3_at_mass_beta_one": fraction(p3),
        "P4_at_mass_beta_one": fraction(p4),
        "mean_strict_upper": fraction(mean_upper),
        "mean_target": fraction(F(1, 10 ** 59)),
        "scattering_strict_upper": fraction(scattering_upper),
        "scattering_target_power10": -43,
        "vacuum_strict_upper": fraction(vacuum_upper),
        "pair_strict_upper": fraction(pair_upper),
        "total_variance_strict_upper": fraction(total_upper),
        "variance_target": fraction(F(1, 10 ** 15)),
        "variance_margin": fraction(F(1, 10 ** 15) - total_upper),
        "uncentered_second_moment_strict_upper": fraction(total_upper + mean_upper ** 2),
        "requires_mass_within_any_declared_fixed_EFT_upper_bound": True,
    }


def calculate():
    files = {
        "proof.md": HERE / "proof.md",
        "check.py": HERE / "check.py",
        "../_shared/round1291_admission.md": HERE.parent / "_shared/round1291_admission.md",
        "../_admission/heavy_thermal_screen/selection.md":
            HERE.parent / "_admission/heavy_thermal_screen/selection.md",
        "../1290/proof.md": HERE.parent / "1290/proof.md",
        "../1290/results.json": HERE.parent / "1290/results.json",
    }
    return {
        "schema": "round1291-heavy-thermal-fixed-smearing-v1",
        "status": "PASS_AUTHOR_FINITE_THERMAL_SOURCE_BOUND",
        "round": 1291, "author_check_groups": 1,
        "source_bindings": {name: digest(p) for name, p in files.items()},
        "primary_source": {
            "url": "https://arxiv.org/pdf/1209.6068",
            "read_scope": "v2 Theorem 6.2(ii),(iv) and local normality proof",
            "explicit_source_moment_constants_derived_here": True,
        },
        "vertices": vertex_checks(),
        "wick": wick_checks(),
        "integrals": integral_checks(),
        "compact": compact_checks(),
        "scope": {
            "fixed_positive_beta_inertial_physical_time": True,
            "fixed_finite_mass_independent_smooth_spacetime_menu": True,
            "minimal_free_KG_complete_mass_term": True,
            "vacuum_Wick_prescription_not_thermal_normal_ordering": True,
            "nonzero_thermal_mean_retained": True,
            "pair_and_number_preserving_scattering_retained": True,
            "locally_normal_Hadamard_KMS_reference": True,
            "global_finite_energy_thermal_preparation_claimed": False,
            "global_vacuum_Fock_density_matrix_claimed": False,
            "actual_readout_or_preparation_constructed": False,
            "nonzero_portal_history_certified": False,
            "curved_HG_or_self_consistent_Einstein_solution_claimed": False,
            "all_higher_moments_or_all_resolutions_claimed": False,
            "finite_samples_replace_analytic_continuum_proof": False,
            "finite_covariance_implies_unique_species_claimed": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="Exclusively create results.json once.")
    args = parser.parse_args()
    actual = calculate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(actual, stream, ensure_ascii=False, sort_keys=True,
                      indent=2, allow_nan=False)
            stream.write("\n")
        receipt = "CREATED_EXCLUSIVELY"
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        if actual != saved:
            raise AssertionError("Read-only recomputation differs from saved result")
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({
        "status": actual["status"], "receipt": receipt,
        "strict_bounds": {
            "mean": "<1e-59", "thermal_scattering": "<1e-43",
            "full_variance": "<1e-15",
        },
        "result_sha256": digest(path),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
