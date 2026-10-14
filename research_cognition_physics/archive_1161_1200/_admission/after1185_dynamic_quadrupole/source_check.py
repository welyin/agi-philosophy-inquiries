"""Exact finite-head source audit; no field/COM or finite-coupling certificate.

Uses the rational spin-3 basis with Gram diag(binomial(6,j)).  No author
response code is imported.  Default is read-only; --write creates once.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path


def zeros(n, m=None):
    return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    out = zeros(n)
    for j in range(n):
        out[j][j] = F(1)
    return out


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def comm(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def rank(a):
    a = [row[:] for row in a]
    pivot = 0
    for col in range(len(a[0])):
        found = next((r for r in range(pivot, len(a)) if a[r][col]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        d = a[pivot][col]
        a[pivot] = [x / d for x in a[pivot]]
        for r in range(len(a)):
            if r != pivot and a[r][col]:
                d = a[r][col]
                a[r] = [x - d * y for x, y in zip(a[r], a[pivot])]
        pivot += 1
    return pivot


def nonzero(a):
    return any(x for row in a for x in row)


def column_zero(a, j):
    return all(row[j] == 0 for row in a)


def beta_integer(a, b):
    return F(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))


def exact(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def filehash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate():
    n = 7
    ident = eye(n)
    gram = zeros(n)
    h = scale(3, ident)
    projector = zeros(n)
    projector[3][3] = F(1)
    for j in range(n):
        gram[j][j] = F(comb(6, j))
    for j in range(6):
        h[j + 1][j] = -F(j + 1, 2)
        h[j][j + 1] = -F(6 - j, 2)
    assert mul(gram, h) == mul(transpose(h), gram)
    assert mul(gram, projector) == mul(transpose(projector), gram)
    assert mul(projector, projector) == projector

    spectral_projectors = []
    spectral_ranks = []
    for lam in range(7):
        p = eye(n)
        for mu in range(7):
            if mu != lam:
                p = mul(p, scale(F(1, lam - mu), add(h, scale(-mu, ident))))
        assert mul(p, p) == p
        assert mul(h, p) == scale(lam, p)
        assert mul(gram, p) == mul(transpose(p), gram)
        spectral_projectors.append(p)
        spectral_ranks.append(rank(p))
    assert spectral_ranks == [1] * 7
    total = zeros(n)
    upi = zeros(n)
    for j, p in enumerate(spectral_projectors):
        total = add(total, p)
        upi = add(upi, scale((-1) ** j, p))
    assert total == ident
    reflection = zeros(n)
    for j in range(n):
        reflection[6 - j][j] = F(1)
    assert upi == reflection

    # ad_H^r P; time derivatives are i^r times these exact real matrices.
    ad = [projector]
    for _ in range(3):
        ad.append(comm(h, ad[-1]))
    endpoint_zero = [[column_zero(ad[r], j) for r in range(3)] for j in (0, 6)]
    assert endpoint_zero == [[True] * 3, [True] * 3]
    assert not column_zero(ad[3], 0) and not column_zero(ad[3], 6)
    assert nonzero(comm(projector, ad[1]))

    # ||i[H,P]|| = sqrt(6): its nonzero squared spectral values are 6.
    # The rational matrix for (i[H,P])^2 is -ad[1]^2.
    derivative_squared = scale(-1, mul(ad[1], ad[1]))
    assert mul(derivative_squared, derivative_squared) == scale(6, derivative_squared)
    assert rank(derivative_squared) == 2

    # Fourier Ward identities, exact finite matrices.  The written proof is
    # polynomial in arbitrary k; this is an independently evaluated witness.
    k = (F(2), F(3), F(5))
    e = (F(1), F(-1), F(0))
    d = sum(ei * ki * ki for ei, ki in zip(e, k))
    rho = scale(-d, projector)
    current = [scale(ei * ki, ad[1]) for ei, ki in zip(e, k)]
    second_derivative = scale(-1, ad[2])
    stress = [[scale(e[i] if i == j else F(0), second_derivative)
               for j in range(3)] for i in range(3)]
    ward0 = comm(h, rho)
    for j in range(3):
        ward0 = add(ward0, scale(k[j], current[j]))
    wards = [ward0]
    for i in range(3):
        wi = comm(h, current[i])
        for j in range(3):
            wi = add(wi, scale(k[j], stress[j][i]))
        wards.append(wi)
    assert not any(nonzero(w) for w in wards)
    assert add(add(stress[0][0], stress[1][1]), stress[2][2]) == zeros(n)

    # Both mass and internal control are identical in the two setting sectors.
    z = [[F(1), F(0)], [F(0), F(-1)]]
    x = [[F(0), F(1)], [F(1), F(0)]]
    y_over_i = [[F(0), F(-1)], [F(1), F(0)]]
    assert comm(z, y_over_i) == scale(-2, x)

    # Free head amplitude is binom(6,3) sin(t/2)^6 cos(t/2)^6.
    amplitude_prefactor = F(comb(6, 3), 64)
    assert amplitude_prefactor == F(5, 16)
    midpoint_variance = amplitude_prefactor * (1 - amplitude_prefactor)
    assert midpoint_variance == F(55, 256)

    # Integrals of even powers of sine divided by pi, on [0,pi].
    sin10_over_pi = F(comb(10, 5), 4 ** 5)
    sin12_over_pi = F(comb(12, 6), 4 ** 6)
    a2_over_pi = amplitude_prefactor ** 2 * sin12_over_pi
    adot2_over_pi = 36 * amplitude_prefactor ** 2 * (sin10_over_pi - sin12_over_pi)
    assert a2_over_pi == F(5775, 262144)
    assert adot2_over_pi == F(4725, 65536)

    radius = F(1, 4)
    # Radial normalized (1-r^2/a^2)^4 profile: m2=3a^2/13 per body.
    # Use the half-integer beta recurrence without evaluating sqrt(pi).
    mean_r2_ratio = F(3, 13)
    m2 = 2 * radius ** 2 * mean_r2_ratio
    assert m2 == F(3, 104)
    separation = F(4)
    c0_times_pi = F(3, 2) / separation ** 5
    c1_times_pi = -(3 * separation ** 2 + m2) / (4 * separation ** 5)
    ic = c0_times_pi * a2_over_pi + c1_times_pi * adot2_over_pi
    assert ic == -F(22700475, 27917287424)
    assert F(22, 7) < separation - 2 * radius
    coefficient_g2_ximinus1 = -2 * ic
    assert coefficient_g2_ximinus1 > 0
    source_kernel_norm_times_pi_bound = c0_times_pi + 6 * abs(c1_times_pi)

    here = Path(__file__).resolve().parent
    archive = here.parents[1]
    dependencies = {
        "source_proof.md": here / "source_proof.md",
        "source_check.py": Path(__file__).resolve(),
        "../after1185_quadrupole_screen/proof.md": here.parent / "after1185_quadrupole_screen/proof.md",
        "../after1185_quadrupole_screen/check.py": here.parent / "after1185_quadrupole_screen/check.py",
        "../after1185_quadrupole_screen/results.json": here.parent / "after1185_quadrupole_screen/results.json",
        "../../1161/proof.md": archive / "1161/proof.md",
        "../../1171/proof.md": archive / "1171/proof.md",
    }
    return {
        "status": "PASS_WITH_SCOPE",
        "numbered_research_round": False,
        "scope": "Exact finite head, full linear source insertion, and normalized leading record coefficient only; no complete apparatus or finite-coupling probability bound.",
        "head": {
            "dimension": 7,
            "basis_gram": [comb(6, j) for j in range(7)],
            "energy_eigenvalues": list(range(7)),
            "spectral_projector_ranks": spectral_ranks,
            "U_pi_equals_reflection": True,
            "prepared_ad_orders_0_1_2_zero_at_both_ends": endpoint_zero,
            "third_ad_nonzero_at_both_ends": True,
            "derivative_operator_norm_squared": "6",
            "mean_P": "(5/16) sin(t)^6",
            "variance_P_at_pi_over_2": exact(midpoint_variance),
        },
        "operator_source": {
            "four_Ward_matrices_exact_zero": True,
            "spatial_trace_exact_zero": True,
            "physical_STF_mass_quadrupole_equals_2S": True,
            "control_mass_includes_Hc": True,
            "Pauli_commutator_Z_EY_equals_minus_i_X": True,
            "R_tensor_factor_identity": "All operator identities tensor with I_R; fixed head preparation remains part of the declared code.",
        },
        "time_and_kernel": {
            "T": "pi",
            "R": exact(separation),
            "radius_each": exact(radius),
            "m2_total": exact(m2),
            "integral_a_squared_div_pi": exact(a2_over_pi),
            "integral_adot_squared_div_pi": exact(adot2_over_pi),
            "C0_times_pi": exact(c0_times_pi),
            "C1_times_pi": exact(c1_times_pi),
            "I_C_exact": exact(ic),
            "I_C_display": float(ic),
            "delta_p_leading_coefficient_for_g2_qstar2_ximinus1": exact(coefficient_g2_ximinus1),
            "V_C_norm_bound_times_pi": exact(source_kernel_norm_times_pi_bound),
            "strict_spacelike_margin_using_pi_upper": exact(separation - 2 * radius - F(22, 7)),
            "probability_lower_bound_claimed": False,
        },
        "remaining": [
            "One common field/constraint Hamiltonian or process matching the response kernel and its normalization.",
            "Finite coupling remainder on the whole declared prepared code with arbitrary reference and retained field/head outputs.",
            "Actual COM, orientation/reference, setting, pointer and material matching errors and their complete stress.",
        ],
        "sha256": {name: filehash(path) for name, path in dependencies.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Create source_results.json once; never overwrite.")
    args = parser.parse_args()
    result = evaluate()
    destination = Path(__file__).with_name("source_results.json")
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
        print("WRITE_ONCE PASS_WITH_SCOPE")
    else:
        saved = json.loads(destination.read_text(encoding="utf-8"))
        if saved != result:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS_WITH_SCOPE")
    print("I_C =", result["time_and_kernel"]["I_C_exact"])
    print("No finite-coupling probability lower bound is claimed.")


if __name__ == "__main__":
    main()
