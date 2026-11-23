"""Invariant relative weak phase from two mass-ratio slopes and a colored determinant."""
import argparse
import cmath
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(filename):
    return json.loads((HERE / filename).read_bytes())


def diagonal_projection(source, target):
    left, masses, vh = np.linalg.svd(source)
    return (np.diag(left.conj().T @ target @ vh.conj().T) / masses)[::-1]


def ratio_slope(source, derivative):
    step = 2e-5

    def ratio_at(x):
        masses = np.linalg.svd(source + x*derivative, compute_uv=False)
        return math.log(masses[1]/masses[0])

    return (ratio_at(step)-ratio_at(-step))/(2*step)


def slopes(d, e, z0):
    edge_e = z0*diagonal_projection(e, d)
    edge_d = diagonal_projection(d, e)/z0
    a = -.6*(edge_e[1]-edge_e[2])
    b = -.4*(edge_d[1]-edge_d[2])
    return np.array([[a.real, -a.imag], [b.real, b.imag]])


def pair(z):
    return [z.real, z.imag]


def run():
    previous = runpy.run_path(str(HERE/"../1715/check.py"))
    _, report = previous["charged_background"](.1)
    d = previous["unpack"](report["down_texture_10PeV"])
    e = previous["unpack"](report["charged_texture_10PeV"]).T
    z0, zt = e[0,1]/d[0,1], e[2,2]/d[2,2]
    t = zt/z0
    m = abs(t)
    matrix = slopes(d, e, z0)
    # A common lepton phase changes the high-scale texture, both pencil roots
    # and raw z together, but leaves the relative port in this matrix unchanged.
    chi = cmath.exp(.731j)
    np.testing.assert_allclose(slopes(d, chi*e, chi*z0), matrix, atol=2e-13)
    eps = []
    for shifts in itertools.product((-.01, .01), repeat=2):
        eps.append(np.linalg.norm(np.linalg.solve(matrix, np.array(shifts)*[1,m*m])))
    theta_work = 1.7475728155339805e-10
    f = read("../1726/results.json")["canonical_f_GeV_at_crossing"]
    data = {}
    for name, z in (("real_root", z0), ("complex_root", zt),
                    ("reflected_relative_phase", z0*t.conjugate())):
        u = z/z0
        dd = -.4*(e/z-d)
        ee = .6*(e-z*d)
        truth = np.array([u.real, u.imag])
        output = matrix @ truth/np.array([1, abs(u)**2])
        np.testing.assert_allclose([ratio_slope(e, ee), ratio_slope(d, dd)], output, atol=1e-8)
        # Both rates ignore arbitrary independent aligned complex drifts.
        np.testing.assert_allclose([ratio_slope(e, ee+(.37+.44j)*e),
                                    ratio_slope(d, dd+(.23-.27j)*d)], output, atol=1e-8)
        np.testing.assert_allclose(np.linalg.solve(matrix, output*[1, abs(u)**2]), truth, atol=2e-15)
        deriv = .4*(3-(2+t)/u)
        np.testing.assert_allclose(np.trace(np.linalg.solve(d,dd)), deriv)
        modulus_checks = []
        for x in (-.01,.01):
            actual = np.linalg.det(d+x*dd)/np.linalg.det(d)
            expect = (1+.4*x*(1-1/u))**2*(1+.4*x*(1-t/u))
            np.testing.assert_allclose(actual, expect, atol=1e-12)
            modulus_checks.append(pair(actual))
        data[name] = {"relative_zeta": pair(u), "angle_deg": math.degrees(cmath.phase(u)),
                      "slopes_mu_over_tau_s_over_b": output.tolist(),
                      "diagonal_coupling_contrasts_at_v_over_f_unit_projection":
                          (output*246./f).tolist(),
                      "complex_derivative_log_det_down": pair(deriv),
                      "det_ratio_at_x_minus_plus_001": modulus_checks,
                      "small_displacement_for_single_source_EDM_budget_linear":
                          theta_work/abs(deriv.imag)}
    return {"real_readout_matrix": matrix.tolist(),
            "normalized_determinant_abs": abs(np.linalg.det(matrix))/np.prod(np.linalg.norm(matrix, axis=1)),
            "relative_modulus_for_phase_test": m,
            "worst_relative_complex_coordinate_error_per_0p01_both_absolute_slopes": max(eps),
            "max_angle_error_deg_per_0p01_both_absolute_slopes_at_fixed_modulus":
                math.degrees(math.asin(max(eps)/m)),
            "vacuum_reference_f_Sigma_GeV": f, "single_source_theta_budget": theta_work,
            "cases": data}


if __name__ == "__main__":
    arg = argparse.ArgumentParser()
    arg.add_argument("--check", action="store_true")
    args = arg.parse_args()
    result = json.loads(json.dumps(run()), parse_float=lambda z: float(f"{float(z):.11g}"))
    path = HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
        print("1728 independent phase and strong-angle check passed")
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
