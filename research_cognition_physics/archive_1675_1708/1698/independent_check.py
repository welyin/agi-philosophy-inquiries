"""Independent neutral supertrace and charged gauge thresholds for R1698."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def root(func, lo, hi):
    assert func(lo) * func(hi) < 0
    ascending = func(hi) > func(lo)
    for _ in range(75):
        mid = (lo + hi) / 2
        if (func(mid) < 0) == ascending:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def logmass(x, q2):
    return x * (math.log(x / q2) - 1) if x else 0.


def coefficient(a, ell, n, q2):
    # In units y^2/(16 pi^2); ell and n are soft squared masses.
    b = a + n
    divided = (logmass(b, q2) - logmass(ell, q2)) / (b - ell) if b != ell else math.log(b / q2)
    return logmass(b, q2) + logmass(ell, q2) + a * divided - 2 * logmass(a, q2)


def scalar_potential(a, ell, n, q2, h2):
    # In units 1/(32 pi^2); one complex scalar and one Weyl component per eigenvalue.
    off = math.sqrt(a * h2)
    hess = np.array([[h2, off], [off, a + h2]])
    fermions = np.linalg.eigvalsh(hess)
    scalars = np.linalg.eigvalsh(hess + np.diag([ell, n]))

    def tr(values):
        return sum(x*x*(math.log(x / q2) - 1.5) for x in values if x > 0)
    return tr(scalars) - tr(fermions)


def direct_gauge(ms, charged, low):
    """For a trial weak mass, separately cross couplings 1,2, then impose 3."""
    mz = low["mu_ref_GeV"]
    beta = np.array([33/5, 1, -3.])
    beta_sm = np.array([41/10, -19/6, -7.])
    color = np.array([2/5, 0, 1.])
    weak = 1 - color
    inv0 = 4*math.pi/np.array([low[k][0] for k in ("g1_GUT", "g2", "g3")])**2
    conv = np.array([0, 2, 3.]) / (12*math.pi)

    def solve_two(mw):
        # a_i(U)=intercept_i-(b_i+2)*ln(U/MZ)/(2 pi).
        intercept = inv0 - conv + (beta-beta_sm)*math.log(ms/mz)/(2*math.pi)
        for mf, rep in [(m, color) for m in charged] + [(mw, weak)]*2:
            weighted_log = (2*math.log(mf/mz) + math.log(math.hypot(mf, ms)/mz))/3
            intercept += rep*weighted_log/(2*math.pi)
        length = 2*math.pi*(intercept[0]-intercept[1])/(beta[0]-beta[1])
        final = intercept - (beta+2)*length/(2*math.pi)
        return mz*math.exp(length), final

    mw = root(lambda mw: np.diff(solve_two(mw)[1])[1], mz, 1e5)
    unif, inverse = solve_two(mw)
    assert np.ptp(inverse) < 1e-11
    return dict(equal_weak_fermions_GeV=mw, weak_scalars_GeV=math.hypot(mw, ms),
                meeting_GeV=unif, common_inverse_alpha=float(inverse[0]))


def compute():
    inp = read("inputs.json")
    physical = read("results.json")["shared_witness"]
    gauge_inputs = read(inp["gauge_inputs"])
    m = np.array(physical["neutral_poles_GeV"])
    y2 = np.array(physical["neutral_pole_Yukawa_squared"])
    a = m*m
    ref = float(y2 @ a / (8*math.pi**2))
    mh2 = read(inp["higgs"])["atlas"]["mass_GeV"][0]**2
    assert abs(ref-mh2) < 1e-5
    ms = inp["ordinary_soft_scale_GeV"]
    uv = inp["common_flavor_cutoff_GeV"]
    charged = physical["colored_fermions_GeV"]

    def neutral(ell, n, q):
        return sum(share*coefficient(ai, ell, n, q*q) for ai,share in zip(a, y2))/(16*math.pi**2)

    def univ(scale, q):
        s = scale*scale
        return float(y2 @ ((a+s)*np.log1p(s/a) - s + s*np.log(a/q**2)) / (8*math.pi**2))

    # Exact univariate second-order coefficient, including coinciding scalar levels.
    for l, n in [(0.6, 0.3), (2.3, 0.3), (3., 0.), (0., 1.), (1., 1.)]:
        eps = 0.000002
        f = lambda u: scalar_potential(2., l, n, 9., u)
        derivative = (-3*f(0) + 4*f(eps) - f(2*eps))/(2*eps)
        np.testing.assert_allclose(derivative/2, coefficient(2., l, n, 9.), atol=1e-6)
    for q in [ms, uv]:
        np.testing.assert_allclose(neutral(ms**2, ms**2, q), univ(ms, q), rtol=1e-12)
    # Scale dependence is canceled by a soft counterterm with the opposite slope.
    test_l, test_n = ms**2, 0.17*ms**2
    shift = -sum(y2)*(test_l+test_n)*math.log(uv**2 / ms**2)/(16*math.pi**2)
    np.testing.assert_allclose(neutral(test_l, test_n, uv) - neutral(test_l, test_n, ms), shift)
    baseline = {}
    for name, q in [("soft_scale", ms), ("flavor_cutoff", uv)]:
        c = univ(ms, q)
        baseline[name] = dict(Q_GeV=q, neutral_CW_coefficient_GeV2=c,
                              required_other_and_tree_coefficient_GeV2=-mh2/2-c,
                              neutral_CW_at_same_gauge_but_zero_singlet_soft_GeV2=neutral(ms**2, 0, q))
    gauge0 = direct_gauge(ms, charged, gauge_inputs)
    np.testing.assert_allclose(gauge0["equal_weak_fermions_GeV"], physical["weak_fermions_GeV"][0],rtol=2e-9)
    points = {}
    # Boundary C_tree(Q)=0 for each distinct construction; no other sector in this curvature.
    for label, choose_q in [("boundary_at_common_soft_scale", lambda scale: scale),
                            ("boundary_at_fixed_flavor_cutoff", lambda _: uv)]:
        scale = root(lambda v: univ(v, choose_q(v)) + mh2/2, 100, 3e4)
        q = choose_q(scale)
        assert max(np.hypot(m, scale)) < uv
        points[label] = dict(soft_GeV=scale, Q_GeV=q, mass_ratio_soft_over_largest_neutral=scale/max(m),
                             neutral_CW_coefficient_GeV2=univ(scale, q), gauge=direct_gauge(scale,charged,gauge_inputs))
    zstar = root(lambda z: (1+z)*math.log1p(z)-z*(math.log(z)+1)+0.5, 1, 5)
    return dict(fixed_neutral_fermion_positive_moment_GeV2=ref, target_potential_coefficient_GeV2=-mh2/2,
                h_equals_minus_half_squared_scale_ratio=zstar,
                frozen_30TeV_universal_soft=baseline, actual_zero_boundary_solutions=points)


def tidy(data):
    if isinstance(data, dict):
        return {key:tidy(value) for key,value in data.items()}
    return float(f"{data:.9g}") if isinstance(data, (float, np.floating)) else data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    data = tidy(compute())
    if parser.parse_args().check:
        assert data == read("independent_results.json")
        print("1698 independent: neutral soft matrices, scheme transfer and three-coupling crossings verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(data, indent=2)+"\n",encoding="utf-8")
