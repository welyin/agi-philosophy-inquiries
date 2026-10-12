"""1122 algebra checks. Default recomputes and compares results.json read-only.

The finite positive-energy Fock model is an algebra calibration, NOT a
discretization/truncation of the continuum Dirac field. Continuum compact
support, the covariance identity and the heat-kernel sign are proved in proof.md.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TOL = 2e-11
N = 12  # 2 flavours x (3 particle modes + 3 hole modes)
VAC = {0: 1.0 + 0j}
X = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                 [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
H1 = np.array([[2, -.2, math.sqrt(.8)],
               [-.2, 2, math.sqrt(.8)],
               [math.sqrt(.8), math.sqrt(.8), 4]], dtype=complex)
H = np.kron(np.eye(4), H1)
EVAL, EVEC = np.linalg.eigh(H1)


def combine(*terms):
    out = {}
    for coefficient, vector in terms:
        for key, value in vector.items():
            out[key] = out.get(key, 0j) + coefficient * value
    return {key: value for key, value in out.items() if abs(value) > 1e-15}


def inner(left, right):
    return sum(value.conjugate() * right.get(key, 0j)
               for key, value in left.items())


def distance(left, right):
    delta = combine((1, left), (-1, right))
    return float(math.sqrt(max(0, inner(delta, delta).real)))


def fermi(vector, mode, creation):
    result = {}
    for mask, value in vector.items():
        occupied = bool(mask & (1 << mode))
        if occupied == creation:
            continue
        sign = (-1) ** ((mask & ((1 << mode) - 1)).bit_count())
        result[mask ^ (1 << mode)] = sign * value
    return result


def c(vector, flavour, position, dagger=False):
    base = 6 * flavour
    return combine((1 / math.sqrt(2), fermi(vector, base + position, dagger)),
                   (1 / math.sqrt(2), fermi(vector, base + 3 + position, not dagger)))


def empty_projection(vector, flavour, position):
    return combine((1, vector), (-1, c(c(vector, flavour, position),
                                     flavour, position, True)))


def empty_state(positions):
    state = VAC
    for position in positions:
        for flavour in (0, 1):
            state = empty_projection(state, flavour, position)
    probability = float(inner(state, state).real)
    return combine((1 / math.sqrt(probability), state)), probability


def code(positions):
    zero, probability = empty_state(positions)
    columns = []
    for flavours in itertools.product((0, 1), repeat=len(positions)):
        vector = zero
        for flavour, position in reversed(list(zip(flavours, positions))):
            vector = c(vector, flavour, position, True)
        columns.append(vector)
    return columns, probability


def q(vector, position):
    return combine((1, c(c(vector, 1, position), 0, position, True)),
                   (1, c(c(vector, 0, position), 1, position, True)))


def local_gate(vector, position, theta):
    first = q(vector, position)
    return combine((1, vector), (math.cos(theta) - 1, q(first, position)),
                   (-1j * math.sin(theta), first))


def logical_vector(columns, amplitudes):
    return combine(*zip(amplitudes, columns))


def gamma_apply(vector, matrix):
    """Exterior-power action from creation operators, independent of minors."""
    result = {}
    for mask, amplitude in vector.items():
        output = VAC
        for mode in reversed([j for j in range(N) if mask & (1 << j)]):
            output = combine(*[(matrix[row, mode], fermi(output, row, True))
                               for row in range(N) if abs(matrix[row, mode]) > 1e-15])
        result = combine((1, result), (amplitude, output))
    return result


def exterior_element(left, matrix, right):
    """Full Fock group matrix element using occupied-mode determinants."""
    answer = 0j
    for lm, lv in left.items():
        rows = [j for j in range(N) if lm & (1 << j)]
        for rm, rv in right.items():
            if lm.bit_count() != rm.bit_count():
                continue
            cols = [j for j in range(N) if rm & (1 << j)]
            determinant = 1.0 if not rows else np.linalg.det(matrix[np.ix_(rows, cols)])
            answer += lv.conjugate() * rv * determinant
    return answer


def compressed(columns, matrix):
    return np.array([[exterior_element(left, matrix, right) for right in columns]
                     for left in columns])


def group(t):
    u = (EVEC * np.exp(1j * t * EVAL)) @ EVEC.conj().T
    return u, np.kron(np.eye(4), u)


def global_rotation(theta):
    e = math.cos(theta) * I2 - 1j * math.sin(theta) * X
    result = np.eye(N, dtype=complex)
    for mode in range(3):
        result[np.ix_([mode, mode + 6], [mode, mode + 6])] = e
        result[np.ix_([mode + 3, mode + 9], [mode + 3, mode + 9])] = e.conj()
    return e, result


def expected_kernel(u):
    p0, pL = u[0, 0], u[0, 1]
    determinant = p0 * p0 - pL * pL
    ff = p0 * (1 + determinant) / 2
    fg = pL * (1 - determinant) / 2
    empty = (1 + 2 * p0 * p0 + 2 * pL * pL + determinant * determinant) / 4
    a, b = ff * ff, -fg * fg
    return a * np.eye(4) + b * SWAP, (p0, pL, determinant, ff, fg, empty, a, b)


def close_scalar(left, right, path="root"):
    if isinstance(left, dict):
        assert left.keys() == right.keys(), path
        for key in left:
            close_scalar(left[key], right[key], path + "." + key)
    elif isinstance(left, list):
        assert len(left) == len(right), path
        for index, value in enumerate(left):
            close_scalar(value, right[index], path + f"[{index}]")
    elif isinstance(left, (float, int)) and not isinstance(left, bool):
        assert abs(left - right) <= 1e-10 * max(1, abs(left), abs(right)), path
    else:
        assert left == right, path


def run():
    one, probability_one = code((0,))
    two, probability_two = code((0, 1))
    theta = .71
    e, g = global_rotation(theta)
    v = np.kron(e, I2)
    gram1 = np.array([[inner(a, b) for b in one] for a in one])
    gram2 = np.array([[inner(a, b) for b in two] for a in two])
    code_residual = max(np.linalg.norm(gram1 - I2), np.linalg.norm(gram2 - np.eye(4)))
    local_residual = 0.0
    global_residual = 0.0
    for columns, logical_g in [(one, e), (two, np.kron(e, e))]:
        for index, column in enumerate(columns):
            global_residual = max(global_residual,
                                  distance(gamma_apply(column, g),
                                           logical_vector(columns, logical_g[:, index])))
    for position, logical_u in [(0, v), (1, np.kron(I2, e))]:
        for index, column in enumerate(two):
            local_residual = max(local_residual,
                                 distance(local_gate(column, position, theta),
                                          logical_vector(two, logical_u[:, index])))
    # Explicit entangled three-dimensional passive reference, in addition to
    # checking the entire four-dimensional operator identity above.
    rng = np.random.default_rng(1122)
    amplitudes = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
    amplitudes /= np.linalg.norm(amplitudes)
    reference_residual = max(distance(local_gate(logical_vector(two, amplitudes[:, k]), 0, theta),
                                      logical_vector(two, (v @ amplitudes)[:, k])) for k in range(3))
    parity_ok = all(mask.bit_count() % 2 == 0 for col in two for mask in col)
    assert max(code_residual, local_residual, global_residual, reference_residual) < TOL
    assert abs(probability_one - .25) < TOL and abs(probability_two - .0625) < TOL
    assert parity_ok

    times = [.13, .37, .81]
    formula_residual = 0.0
    determinant_vs_creation_residual = 0.0
    single_scalar_residual = 0.0
    global_kernel_commutator = 0.0
    selective_commutator_residual = 0.0
    nonunitarity = []
    selective_defects = []
    wrong_sign_defects = []
    wrong_modulus_defects = []
    same_flavour_residual = 0.0
    for t in times:
        u, full = group(t)
        actual = compressed(two, full)
        expected, (_, _, d, ff, fg, empty, a, b) = expected_kernel(u)
        formula_residual = max(formula_residual, float(np.linalg.norm(actual - expected)))
        same_flavour_residual = max(same_flavour_residual, abs(a + b - d * empty))
        direct_column = gamma_apply(two[1], full)
        determinant_vs_creation_residual = max(determinant_vs_creation_residual,
            max(abs(inner(left, direct_column) - actual[k, 1]) for k, left in enumerate(two)))
        single = compressed(one, full)
        single_scalar_residual = max(single_scalar_residual,
                                    float(np.linalg.norm(single - np.trace(single) / 2 * I2)))
        total_logical = np.kron(e, e)
        global_kernel_commutator = max(global_kernel_commutator,
                                      float(np.linalg.norm(actual @ total_logical - total_logical @ actual)))
        defect = float(np.linalg.norm(actual @ v - v @ actual, 2))
        selective_defects.append(defect)
        selective_commutator_residual = max(selective_commutator_residual,
                                           abs(defect - 2 * abs(b) * abs(math.sin(theta))))
        compressed_u = u[:2, :2]
        nonunitarity.append(float(np.linalg.norm(compressed_u.conj().T @ compressed_u - I2, 2)))
        wrong_sign_defects.append(float(np.linalg.norm(actual - (a * np.eye(4) - b * SWAP))))
        wrong_modulus_defects.append(float(np.linalg.norm(actual -
                                    (abs(ff) ** 2 * np.eye(4) - abs(fg) ** 2 * SWAP))))
    assert max(formula_residual, same_flavour_residual, determinant_vs_creation_residual,
               single_scalar_residual, global_kernel_commutator, selective_commutator_residual) < TOL
    assert min(selective_defects) > 1e-5 and min(nonunitarity) > .01
    assert min(wrong_sign_defects) > 1e-5 and min(wrong_modulus_defects) > .1

    expected_quartic = -(H1[0, 0].real * H1[0, 1].real) ** 2
    asymptotic = []
    for t in [.02, .01, .005, .0025]:
        u, full = group(t)
        actual = compressed(two, full)[1, 2]
        ratio = actual / t ** 4
        asymptotic.append({"t": t, "B_over_t4_real": float(ratio.real),
                           "B_over_t4_imag": float(ratio.imag),
                           "error_from_leading_coefficient": float(abs(ratio - expected_quartic))})
    errors = [row["error_from_leading_coefficient"] for row in asymptotic]
    assert all(b < .6 * a for a, b in zip(errors, errors[1:]))
    assert errors[-1] < .004
    assert abs((H1 @ H1)[0, 1]) < TOL
    full_conservation = float(np.linalg.norm(H @ g - g @ H))
    assert full_conservation < TOL and min(EVAL) > 0
    u, full = group(.13)
    k = compressed(two, full)
    trivial = np.kron(-I2, I2)
    trivial_residual = float(np.linalg.norm(k @ trivial - trivial @ k))
    assert trivial_residual < TOL

    return {
        "round": 1122,
        "scope": "Finite positive-energy Fock algebra calibration; not a continuum Dirac truncation",
        "groups_passed": 4,
        "continuum_obligations": "Compact support, Dirac covariance I/2, and e_L<0 are analytic in proof.md",
        "code_and_addressability": {
            "single_empty_probability": probability_one,
            "double_empty_probability": probability_two,
            "max_isometry_residual": float(code_residual),
            "max_local_gate_residual": local_residual,
            "max_global_gate_residual": global_residual,
            "entangled_reference_residual": reference_residual,
            "double_code_fixed_even_parity": parity_ok,
        },
        "full_fock_propagation": {
            "times": times,
            "max_formula_residual": formula_residual,
            "max_same_flavour_identity_residual": float(same_flavour_residual),
            "max_determinant_vs_creation_residual": float(determinant_vs_creation_residual),
            "compressed_one_particle_nonunitarity": nonunitarity,
            "wrong_exchange_sign_defects": wrong_sign_defects,
            "wrong_absolute_square_defects": wrong_modulus_defects,
        },
        "conservation_and_spectrum": {
            "one_particle_energy_eigenvalues": [float(x) for x in EVAL],
            "global_energy_commutator_residual": full_conservation,
            "max_single_code_scalar_residual": single_scalar_residual,
            "max_double_global_rotation_commutator": global_kernel_commutator,
            "selective_rotation_commutators": selective_defects,
            "max_exact_commutator_norm_formula_residual": selective_commutator_residual,
            "scalar_theta_pi_exception_residual": trivial_residual,
            "trace_kernel_at_first_time_abs": float(abs(np.trace(k))),
        },
        "quartic_exchange": {
            "finite_model_e0": float(H1[0, 0].real),
            "finite_model_eL": float(H1[0, 1].real),
            "H1_squared_fg_abs": float(abs((H1 @ H1)[0, 1])),
            "expected_B_over_t4": float(expected_quartic),
            "samples": asymptotic,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write initial results artifact")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        status = "written"
    else:
        close_scalar(result, json.loads(path.read_text(encoding="utf-8")))
        status = "matched read-only"
    print(json.dumps({"round": 1122, "groups_passed": result["groups_passed"],
                      "saved_results": status,
                      "max_formula_residual": result["full_fock_propagation"]["max_formula_residual"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
