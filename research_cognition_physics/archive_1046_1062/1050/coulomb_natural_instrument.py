"""1050: finite calibration of a genuine Coulomb-parent instrument bound.

No spatial grid or surrogate six-dimensional Hamiltonian evolves Coulomb states.
The continuous propagation/instrument estimate is proved in proof.md.
Default is read-only comparison; --write creates results exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULT = HERE / "coulomb_natural_instrument_results.json"
HISTORY = [
    "archive_956_989/research_note_956.md",
    "archive_956_989/research_note_960.md",
    "archive_956_989/research_note_965.md",
    "archive_956_989/research_note_981.md",
    "archive_1046_/research_note_1048.md",
    "archive_1046_/1048/proof.md",
    "archive_1046_/1048/material_matching_results.json",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(a):
    return float(np.linalg.norm(a, 2))


def analytic_raw(R):
    e = math.exp(-R)
    S = e * (1 + R + R * R / 3)
    kappa = e * R * (R + 1) / 3
    xx = e * (1 + R + 9 * R**2 / 20 + 7 * R**3 / 60 + R**4 / 60)
    pp = e * (1 + R - R * R) / 3
    grad = e * (1 + R - R * R / 3)
    return {
        "S": S,
        "G": np.array([[1, S], [S, 1.]]),
        "x": np.diag([-R / 2, R / 2]),
        "x2": np.array([[R * R / 4 + 1, xx], [xx, R * R / 4 + 1]]),
        "p": np.array([[0, -1j * kappa], [1j * kappa, 0]]),
        "p2": np.array([[1 / 3, pp], [pp, 1 / 3]]),
        "grad2": np.array([[1., grad], [grad, 1.]]),
    }


def prolate_cross(R, order=40):
    # u = R(mu-1): Gauss-Laguerre; nu: Gauss-Legendre.
    u, wu = np.polynomial.laguerre.laggauss(order)
    v, wv = np.polynomial.legendre.leggauss(order)
    mu = 1 + u[:, None] / R
    nu = v[None, :]
    weight = R * R * math.exp(-R) / 4 * wu[:, None] * wv[None, :]
    jac = mu * mu - nu * nu
    px_left = (mu * nu + 1) / (mu + nu)
    px_right = (mu * nu - 1) / (mu - nu)
    x = R * mu * nu / 2
    return {
        "S": float(np.sum(weight * jac)),
        "x": float(np.sum(weight * jac * x)),
        "x2": float(np.sum(weight * jac * x * x)),
        "p": 1j * float(np.sum(weight * jac * px_right)),
        "p2": float(np.sum(weight * jac * px_left * px_right)),
        "grad2": float(np.sum(weight * (mu * mu + nu * nu - 2))),
    }


def creation_operators():
    creators = []
    for j in range(4):
        c = np.zeros((16, 16), complex)
        for bits in range(16):
            if not (bits >> j) & 1:
                c[bits | (1 << j), bits] = (-1) ** ((bits & ((1 << j) - 1)).bit_count())
        creators.append(c)
    sector = np.array([i for i in range(16) if i.bit_count() == 2])
    return creators, sector


def second_quantize(a, creators, sector):
    full = sum(a[i, j] * creators[i] @ creators[j].conj().T
               for i in range(4) for j in range(4))
    return full[np.ix_(sector, sector)]


def geometry(R, creators, sector):
    raw = analytic_raw(R)
    val, vec = np.linalg.eigh(raw["G"])
    L = (vec * val**(-.5)) @ vec.T
    mats = {name: L @ raw[name] @ L for name in ["x", "x2", "p", "p2", "grad2"]}
    d = R / (2 * math.sqrt(1 - raw["S"]**2))
    Q = second_quantize(np.diag([-.5, -.5, .5, .5]), creators, sector)
    X = second_quantize(np.kron(mats["x"], np.eye(2)), creators, sector)
    J = second_quantize(np.kron(mats["p"], np.eye(2)), creators, sector)
    cx = mats["x2"] - mats["x"] @ mats["x"]
    cp = mats["p2"] - mats["p"] @ mats["p"]
    # Compression of full one-body sums squared includes outside-orbital variance.
    leak_gram = second_quantize(np.kron(cx, np.eye(2)), creators, sector)
    current_second = J @ J + second_quantize(np.kron(cp, np.eye(2)), creators, sector)
    b_exact = math.sqrt(max(0., float(np.linalg.eigvalsh(leak_gram)[-1])))
    j_initial = math.sqrt(max(0., float(np.linalg.eigvalsh(current_second)[-1])))
    r = math.sqrt(2 / (1 - raw["S"]))
    errors = [norm(L @ raw["G"] @ L - np.eye(2)), norm(X - 2 * d * Q)]
    kappa = math.exp(-R) * R * (R + 1) / 3
    expected_p = kappa / math.sqrt(1 - raw["S"]**2) * np.array([[0, -1j], [1j, 0]])
    errors.append(norm(mats["p"] - expected_p))
    assert np.linalg.eigvalsh(cx)[0] > -2e-12
    assert np.linalg.eigvalsh(cp)[0] > -2e-12
    assert b_exact <= 2 * r + 1e-12
    assert j_initial <= math.sqrt(8 * (32 + 7 * r))
    # A one-body microscopic current never directly moves a doubly occupied pair.
    dl, dr = list(sector).index(3), list(sector).index(12)
    assert abs(J[dl, dr]) < 1e-14
    assert max(errors) < 2e-12
    return {
        "R": R, "overlap": raw["S"], "ell": 2 * d,
        "normalization_dipole_current_max_residual": max(errors),
        "actual_position_leakage_norm": b_exact,
        "coarse_position_leakage_bound": 2 * r,
        "actual_initial_current_second_moment_norm": j_initial,
        "energy_norm_upper_bound": 7 * r,
        "microscopic_pair_current_matrix_element_abs": float(abs(J[dl, dr])),
    }, Q


def instrument(Q, sector):
    D = Q @ Q
    phase = np.diag(np.exp(-1j * np.pi * np.diag(Q)))
    K = [(np.eye(6) + phase) / 2, (np.eye(6) - phase) / 2]
    ideal = [np.eye(6) - D, D]
    errors = [norm(K[i] - ideal[i]) for i in range(2)]
    errors += [norm(sum(k.conj().T @ k for k in K) - np.eye(6))]
    # Passive six-dimensional reference, maximum entanglement; compare full branches.
    ent = np.eye(6, dtype=complex).reshape(36) / math.sqrt(6)
    branch = [np.kron(k, np.eye(6)) @ ent for k in K]
    ideal_branch = [np.kron(k, np.eye(6)) @ ent for k in ideal]
    errors += [float(np.linalg.norm(a - b)) for a, b in zip(branch, ideal_branch)]
    dl, single = list(sector).index(3), list(sector).index(5)
    psi = (np.eye(6)[:, dl] + np.eye(6)[:, single]) / math.sqrt(2)
    rho = np.outer(psi, psi)
    dephased = sum(k @ rho @ k.conj().T for k in K)
    disturbance = float(np.sum(np.abs(np.linalg.eigvalsh(rho - dephased))) / 2)
    assert abs(disturbance - .5) < 1e-13
    assert max(errors) < 1e-12
    return {
        "Q_diagonal_in_binary_two_particle_basis": np.diag(Q).real.tolist(),
        "kraus_ranks_plus_minus": [int(np.linalg.matrix_rank(k, tol=1e-10)) for k in K],
        "full_reference_max_residual": max(errors),
        "ideal_double_vs_single_read_contrast": 1,
        "ideal_coherent_input_nonselective_disturbance_half_trace": disturbance,
        "natural_Coulomb_output_isometry": "W_T = exp(-i H T) W; not numerically replaced",
    }


def rational_window():
    # exp(10) > its positive degree-7 Taylor polynomial already exceeds 100(1+10+100/3).
    e10lower = sum(F(10**j, math.factorial(j)) for j in range(8))
    assert e10lower > 100 * (1 + 10 + F(100, 3))
    assert F(200, 99) < F(100, 49)  # sqrt(2/(1-S)) < 10/7
    assert 8 * 42 < F(37, 2)**2
    R, T = F(1000), F(1000, 137)
    eps = F(22, 7) * (F(20, 7) / R + F(37, 2) * T / (2 * R))
    contrast = 1 - 2 * eps
    assert eps == F(74239, 335650) < F(9, 40)
    assert contrast == F(93586, 167825) > F(11, 20)
    # Independent quadrature checks pulse area and first moment, not propagation error.
    z, w = np.polynomial.legendre.leggauss(64)
    t = (z + 1) * float(T) / 2
    beta = 2 * np.pi / (float(R) * float(T)) * np.sin(np.pi * t / float(T))**2
    area = float(w @ beta * float(T) / 2)
    moment = float(w @ (t * beta) * float(T) / 2)
    pulse_error = max(abs(area - np.pi / float(R)), abs(moment - float(T) * area / 2))
    assert pulse_error < 1e-14
    return {
        "R_atomic_units": 1000, "T_atomic_units_exact": str(T),
        "S_less_than": "1/100", "b_X_less_than": "20/7",
        "electronic_H_W_norm_less_than": 10, "E_star_less_than": 42,
        "J0_less_than": "37/2", "pi_upper_bound": "22/7",
        "instrument_error_strict_upper_bound_exact": str(eps),
        "instrument_error_strict_upper_bound_decimal": float(eps),
        "actual_read_contrast_strict_lower_bound_exact": str(contrast),
        "actual_read_contrast_strict_lower_bound_decimal": float(contrast),
        "pulse_area_first_moment_quadrature_residual": pulse_error,
        "microscopic_current_expectation_error_upper_bound":
            float(eps * (37 + 4 * F(22, 7000))),
        "comparison_time_not_less_than_R_over_c_if_c_ge_137": True,
        "time_check_is_not_a_Maxwell_or_long_wavelength_realization": True,
    }


def run():
    creators, sector = creation_operators()
    car = max(norm(c.conj().T @ d + d @ c.conj().T - (np.eye(16) if i == j else 0))
              for i, c in enumerate(creators) for j, d in enumerate(creators))
    assert car < 1e-14
    rows, quad_rows = [], []
    for R in [1., 2., 4., 8., 16.]:
        raw, q = analytic_raw(R), prolate_cross(R)
        discrepancies = [abs(q["S"] - raw["S"]), abs(q["x"])]
        for name in ["x2", "p", "p2", "grad2"]:
            discrepancies.append(abs(q[name] - raw[name][0, 1]))
        assert max(discrepancies) < 2e-11, (R, discrepancies)
        quad_rows.append({"R": R, "Slater_prolate_cross_max_residual": max(discrepancies)})
        row, Q = geometry(R, creators, sector)
        rows.append(row)
    # Single normalized Slater: radial u=2r; angular v=cos(theta).
    u, wu = np.polynomial.laguerre.laggauss(40)
    v, wv = np.polynomial.legendre.leggauss(40)
    weight = wu[:, None] * wv[None, :] * u[:, None]**2 / 4
    one = {"norm": float(np.sum(weight)),
           "relative_x2": float(np.sum(weight * (u[:, None] * v[None, :] / 2)**2)),
           "px2": float(np.sum(weight * v[None, :]**2))}
    assert max(abs(one["norm"] - 1), abs(one["relative_x2"] - 1),
               abs(one["px2"] - 1 / 3)) < 2e-12
    return {
        "round": 1050, "date": "2026-10-08", "all_passed": True,
        "new_scientific_groups": 1, "new_cognitive_axioms": 0,
        "full_roadmap_completed": False,
        "historical_sha256": {rel: digest(ROOT / rel) for rel in HISTORY},
        "single_Slater_quadrature": one,
        "prolate_quadrature_calibration": quad_rows,
        "Lowdin_and_full_six_dimensional_geometry": rows,
        "CAR_max_residual": car,
        "instrument": instrument(Q, sector),
        "finite_window": rational_window(),
        "scope": {
            "parent": "fixed-nuclei two-electron 3D electronic Coulomb H; total BO energy H+1/R",
            "preparation": "normal Slater/Lowdin six-dimensional two-electron isometry W",
            "P_spectral_or_invariant_required": False,
            "continuous_Coulomb_propagation_established_by_analysis_not_grid": True,
            "actual_propagation_numerically_simulated": False,
            "task": "initial charge parity instrument with complete natural output W_T",
            "arbitrary_unknown_input_and_passive_reference": True,
            "source": "actual microscopic current first moment with independent second-moment bound",
            "old_956_U_v_long_time_matching": False,
            "assumed_control": "finite-time homogeneous controlled dipole; sin-squared pulse",
            "Maxwell_controller_or_autonomous_energy_closure": False,
            "full_stress_or_GR_or_P981_matching": False,
            "new_dimension_or_Coulomb_law_generation": False,
            "R1000_overlap_underflow_used_as_exact_zero": False,
        },
    }


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path + "." + key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for j, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{j}]")
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-10, abs_tol=3e-12), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1050, "all_passed": True,
                      "mode": "exclusive_write" if args.write else "read_only",
                      "error_bound": result["finite_window"]["instrument_error_strict_upper_bound_exact"],
                      "new_scientific_groups": 1}, ensure_ascii=False))


if __name__ == "__main__":
    main()
