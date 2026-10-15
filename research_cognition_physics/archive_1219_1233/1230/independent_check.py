"""1230 independent rational and finite-CAR audit.

Does not import the author's program or receipt. The continuum Schwinger
inequality and domain statements are proved analytically; the finite CAR
fixture verifies its two-mode algebra, not a cutoff approximation to QFT.
Default and --verify are read-only. --write-results creates a missing receipt.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = HERE / "independent_results.json"
SOURCES = (
    "archive_1120_1158/1126/proof.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1228/results.json",
    "archive_1219_/_admission/after1229_geometry_roles_and_area_priority.md",
    "archive_1219_/1230/independent_check.py",
)
GAMMA = F(1, 10**8)
SIGMA = F(1, 10**120)
DELTA = F(1, 10**24)
M = 10**260
T = 2
# Deliberately retain the wider old guard, independent of author's C=3000.
C = 16 * 10**6
L = 16 / DELTA**2
PI_LO = F(157, 50)
PI_HI = F(22, 7)


def display(x):
    return {"exact": str(x), "display": float(x)}


def integrate_even_power(n):
    # Integral of (1-u^2)^n on [-1,1], exact polynomial evaluation.
    return sum((F(2 * (-1)**j * math.comb(n, j), 2*j+1)
                for j in range(n+1)), F(0))


def atan_partial(x, n):
    return sum(((-1)**j * x**(2*j+1) / (2*j+1)
                for j in range(n)), F(0))


def exact_integrals():
    j = integrate_even_power(8)
    i0 = integrate_even_power(16)
    assert j == F(2**17 * math.factorial(8)**2, math.factorial(17))
    assert i0 == F(2**33 * math.factorial(16)**2, math.factorial(33))
    derivative_integral = 256 * sum(
        (F(2 * (-1)**k * math.comb(14, k), 2*k+3) for k in range(15)), F(0))
    assert derivative_integral / i0 == F(44, 5)
    p_lo = 16*atan_partial(F(1, 5), 16) - 4*atan_partial(F(1, 239), 5)
    p_hi = 16*atan_partial(F(1, 5), 15) - 4*atan_partial(F(1, 239), 6)
    assert PI_LO < p_lo < p_hi < PI_HI
    # heat kernel: pi^-1/2 * (4*pi)^-3/2 = 1/(8*pi^2).
    # derivative and projector factors give 1/(32*pi^2).
    assert F(1, 8) * F(1, 2) * F(1, 2) == F(1, 32)
    return dict(status="PASS", J=display(j), I0=display(i0),
                derivative_ratio="44/5", L1_product_coefficient=display(j**6/i0**3),
                pi_Machin_lower=display(p_lo), pi_Machin_upper=display(p_hi),
                selected_pi_lower=str(PI_LO), selected_pi_upper=str(PI_HI),
                Schwinger_prefactor_without_pi="1/32"), j**6/i0**3


def exp_negative_lower(x):
    # Odd Taylor degree: e^-x >= sum_{j=0}^9 (-x)^j/j!.
    # Lagrange remainder is positive for all x>=0.
    assert 0 <= x <= 1
    result = sum(((-x)**j / math.factorial(j) for j in range(10)), F(0))
    assert result > 0
    return result


def kernel_geometry(q3):
    rows = []
    for a0 in (F(1, 10), F(1, 5)):
        am, ap = a0-GAMMA, a0+GAMMA
        bm, bp = F(1, 100)-GAMMA, F(1, 100)+GAMMA
        dm, dp = 1-GAMMA-2*bp, 1+GAMMA+2*bp
        r2max = dp**2+8*ap**2
        assert 0 < dm < dp and r2max < 4
        # I(r) = integral u exp(-alpha*u-1/u)du.
        i_lower = 16/r2max**2 - 4/r2max
        alpha = dm**2/4
        # e^-1/u <= u/(u+1); lower bound integral u/(1+u)e^-alpha*u
        # on [1,2],[2,4],[4,infinity] gives this independently tighter upper.
        j_lower = (F(1, 2)*exp_negative_lower(alpha)
                   + F(1, 6)*exp_negative_lower(2*alpha)
                   + F(2, 15)*exp_negative_lower(4*alpha)) / alpha
        i_upper = 1/alpha**2-j_lower
        assert 0 < i_lower < i_upper < 1/alpha**2
        cl = am**2*bm*q3*dm*i_lower/(32*PI_HI**2)
        cu = ap**2*bp*q3*dp*i_upper/(32*PI_LO**2)
        assert 0 < cl < cu < F(1, 2)
        h2 = 1+F(44, 5)*(2/am**2+1/bm**2)
        assert h2 < 400**2
        rows.append(dict(a0=a0, c_lower=cl, c_upper=cu, area_lower=4*am**2,
                         area_upper=4*ap**2, h2_upper=h2,
                         d_min=dm, d_max=dp, r2_max=r2max))
    output = [{k: display(v) for k, v in row.items()} for row in rows]
    return dict(status="PASS", independent_upper="piecewise u/(1+u), Taylor-exp lower",
                all_actual_geometry_parameters_transported=["a_x", "a_y", "b", "D"],
                shared_actual_geometry_in_each_menu=True,
                continuum_cutoff_used=False, rows=output), rows


def annihilator(k):
    result = np.zeros((4, 4), dtype=np.complex128)
    for state in range(4):
        if (state >> k) & 1:
            sign = (-1)**((state & ((1 << k)-1)).bit_count())
            result[state ^ (1 << k), state] = sign
    return result


def finite_car():
    a, b = annihilator(0), annihilator(1)
    ident = np.eye(4)
    nf, ng = a.conj().T@a, b.conj().T@b
    assert np.array_equal(a@b+b@a, np.zeros((4, 4)))
    assert np.array_equal(a@a.conj().T+a.conj().T@a, ident)
    assert np.array_equal(nf@ng, ng@nf)
    # Four field states times the two full four-dimensional P/E registers.
    iso = np.zeros((64, 4))
    for state in range(4):
        x, y = state & 1, (state >> 1) & 1
        iso[(state*4+3*x)*4+3*y, state] = 1
    assert np.array_equal(iso.T@iso, ident)
    cases = []
    for c in (F(1, 8), F(1, 4)):
        # Dyadic entries: these finite products and sums are exact in binary64.
        cf = float(c)
        diagonal = np.array([.25-cf*cf, .25+cf*cf, .25+cf*cf, .25-cf*cf])
        rho = np.diag(diagonal).astype(complex)
        rho[1, 2], rho[2, 1] = -1j*cf, 1j*cf
        assert np.trace(rho@nf) == .5 and np.trace(rho@ng) == .5
        assert np.trace(rho@(a.conj().T@b)) == 1j*cf
        assert np.trace(rho@(nf@ng)) == .25-cf*cf
        assert F(1, 4)-c*c > 0 and (F(1, 2)-c)**2 > 0
        output = iso@rho@iso.T
        joint = [0., 0., 0., 0.]
        for field in range(4):
            for pe1 in range(4):
                for pe2 in range(4):
                    index = (field*4+pe1)*4+pe2
                    x, y = pe1 & 1, pe2 & 1
                    joint[2*x+y] += float(output[index, index].real)
        expected = [.25-cf*cf, .25+cf*cf, .25+cf*cf, .25-cf*cf]
        assert joint == expected
        cases.append(dict(c=str(c), joint=joint,
                          exact_covariance=str(-c*c)))
    return dict(status="PASS", fixture="exact-dyadic two-mode CAR restriction identity",
                full_field_spectra_in_continuum_proof=True,
                finite_Fock_dimension=4, full_PE_dimension=16,
                isometry_shape=[64, 4], all_input_isometry=True,
                cases=cases, numerical_continuum_vacuum_replacement=False)


def binary_information_identity():
    cases = []
    for r, s, kappa in ((F(1, 3), F(2, 5), F(1, 100)),
                        (F(1, 2), F(1, 2), -F(1, 64)),
                        (F(2, 3), F(3, 7), -F(1, 50))):
        product = [(1-r)*(1-s), (1-r)*s, r*(1-s), r*s]
        p = [product[0]+kappa, product[1]-kappa,
             product[2]-kappa, product[3]+kappa]
        assert min(p) > 0 and sum(p) == 1
        tv = sum(abs(x-y) for x, y in zip(p, product))/2
        chi2 = sum((x-y)**2/y for x, y in zip(p, product))
        assert tv == 2*abs(kappa)
        assert chi2 == kappa**2/(r*(1-r)*s*(1-s))
        # This display is a check only; the continuum guards use exact
        # Pinsker and ln x <= x-1, not floating evaluation of tiny MI.
        mi = sum(float(x)*math.log(float(x/y)) for x, y in zip(p, product))
        assert float(2*tv**2) <= mi <= float(chi2)
        cases.append(dict(marginals=[str(r), str(s)], covariance=str(kappa),
                          tv=str(tv), chi2=str(chi2), MI_display=mi))
    return dict(status="PASS", cases=cases, lower="8 covariance^2",
                upper="covariance^2/[r(1-r)s(1-s)]",
                ideal_regional_channel=True, actual_channel_exactly_regional=False,
                actual_region_comparison_requires_covariance_error=True)


def mi_bounds(row, error):
    # Deliberately coarser than author's 2e+e^2.
    xi = 3*error+error**2
    assert row["c_lower"]**2 > xi
    lower = 8*(row["c_lower"]**2-xi)**2
    upper = (row["c_upper"]**2+xi)**2/(F(1, 4)-error**2)**2
    return lower, upper, xi


def probability_separation(rows):
    b_t = SIGMA+F(T, 2)/M/SIGMA
    ec = L*T*b_t
    ep = 2*C*DELTA
    e = 2*(ec+ep)
    bounds = [mi_bounds(row, e) for row in rows]
    gap = bounds[1][0]-bounds[0][1]
    assert gap > F(1, 10**21)
    s = F(1, 10**14)
    empirical = [mi_bounds(row, e+2*s) for row in rows]
    empirical_gap = empirical[1][0]-empirical[0][1]
    assert empirical_gap > F(1, 10**21)
    n = 4*10**28
    assert 2*n*s*s == 8
    exp8_lower = sum((F(8)**j/math.factorial(j) for j in range(20)), F(0))
    assert exp8_lower > 1600  # 16 exp(-8) < .01.
    return dict(status="PASS", independent_C=C, covariance_guard="3e+e^2",
                probability_error=display(e), actual_intervals=[
                    [display(x), display(y)] for x, y, _ in bounds],
                actual_large_minus_small_lower=display(gap),
                optional_empirical_intervals=[
                    [display(x), display(y)] for x, y, _ in empirical],
                empirical_large_minus_small_lower=display(empirical_gap),
                per_menu_independent_trials=n, per_letter_error=str(s),
                two_menu_failure_bound="16 exp(-8) < .01",
                physical_repetitions_performed=False), ec, ep, n


def complete_resources(rows, ec, ep, n):
    assert 2*PI_HI*(400+14) < C
    assert max(row["h2_upper"] for row in rows) < 400**2
    a0 = 2**28
    rh, rq = 2*C, a0+2*C
    rha, rqa = 2*T*C/DELTA, a0+2*T*C/DELTA
    eh = (rha+rh)*ec+2*rh*ep
    eq = (rqa+rq)*ec+2*rq*ep
    assert eh < F(21, 10**10) and eq < F(2, 10**8)
    rend = L*(3*SIGMA+L*T*T/(2*M))
    nbar = (SIGMA*L*T+L*T*T/(4*M*SIGMA))**2
    assert rend < F(1, 10**60) and 0 < nbar < F(1, 10**140)
    # g(n) <= n(1+n-ln n), increasing in this range. ln(10)<3.
    n_cap = F(1, 10**140)
    entropy = 2*n_cap*(1+n_cap+420)
    assert entropy < F(1, 10**136)
    per_run_free_energy = M+1/(4*M*SIGMA**2)+6
    return dict(status="PASS", h2_guard=400**2,
                C_is_deliberately_wider_than_author=True,
                field_energy_error=display(eh),
                branch_free_extension_TQ_error=display(eq),
                endpoint_interaction_energy_budget=display(rend),
                per_clock_occupation=display(nbar),
                two_clock_entropy_upper=display(entropy),
                free_initial_energy=display(per_run_free_energy),
                exact_initial_energy_also_includes_nonzero_V0=True,
                PE_single_registers_per_run=4, clocks_per_run=2,
                retained_clocks_two_menu_sampling=4*n,
                retained_PE_single_registers_two_menu_sampling=8*n,
                sampling_free_baseline=display(2*n*per_run_free_energy),
                same_terminal_field_state=True,
                future_source_is_free_reference_extension=True,
                actual_continued_future_source_history=False,
                force_and_clock_backreaction_included=True,
                microscopic_local_device_or_full_boost_conservation=False)


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write-results", action="store_true")
    modes.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    g1, q3 = exact_integrals()
    g2, rows = kernel_geometry(q3)
    g3 = finite_car()
    g4 = binary_information_identity()
    g5, ec, ep, n = probability_separation(rows)
    g6 = complete_resources(rows, ec, ep, n)
    result = dict(round=1230, status="PASS", check_groups=6,
                  groups=dict(exact_integrals_and_kernel_prefactor=g1,
                              transported_geometry_and_independent_kernel_bound=g2,
                              complete_CAR_PE_instrument=g3,
                              classical_information_identity_and_scope=g4,
                              actual_information_and_finite_sampling=g5,
                              complete_source_clock_and_record_budget=g6),
                  limitations=dict(no_universal_area_law=True,
                                   no_full_region_entropy_numerics=True,
                                   no_UV_limit_or_Newton_coefficient=True,
                                   no_dynamical_area_backreaction=True,
                                   no_empirical_experiment=True,
                                   not_stage_completion=True),
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                 for p in SOURCES})
    if args.write_results:
        with TARGET.open("x", encoding="utf8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
    elif json.loads(TARGET.read_text(encoding="utf8")) != result:
        raise SystemExit("Saved independent receipt differs; no files changed.")
    print(json.dumps(dict(round=1230, status="PASS", check_groups=6,
                          actual_MI_gap_lower=g5["actual_large_minus_small_lower"]["display"],
                          empirical_MI_gap_lower=g5["empirical_large_minus_small_lower"]["display"],
                          read_only=not args.write_results), ensure_ascii=False))


if __name__ == "__main__":
    main()
