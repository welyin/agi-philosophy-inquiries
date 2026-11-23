"""Magnetic translations on an actual symmetric bilinear and heavy-current holdouts."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OMEGA = np.exp(2j * np.pi / 3)


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def spectrum_and_perturbation(table, count):
    central, error = np.array(table, dtype=float).T
    lo = central - count * error
    hi = central + count * error
    assert lo[1] > hi[0] and lo[2] > hi[1]
    lower = min(lo[1] - hi[0], lo[2] - hi[1]) / 2
    return {
        "central_spectral_distance": float(np.min(np.diff(central)) / 2),
        "full_box_spectral_distance_lower": float(lower),
        "full_box_relative_to_largest_lower": float(lower / hi[2]),
        "full_box_high_gap_lower": float(lo[2] - hi[1]),
    }


def response_box(quoted, count, response):
    value, *errors = quoted
    return [(value - count * sum(errors)) / response[1], (value + count * sum(errors)) / response[0]]


def translations():
    x = np.roll(np.eye(3, dtype=complex), 1, axis=0)
    z = np.diag([1, OMEGA, OMEGA ** 2])
    tensor = np.eye(9, dtype=complex).reshape(3, 3, 9)
    # Each copy has the same charge 2i and is shifted by i -> i+1.
    bases = [tensor[i, i] for i in range(3)]
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        bases.append((tensor[j, k] + tensor[k, j]) / math.sqrt(2))
    columns = np.stack(bases, axis=1)
    for gate, copy in [(x, x), (z, z ** 2)]:
        output = columns.conj().T @ np.kron(gate, gate) @ columns
        np.testing.assert_allclose(output, np.kron(np.eye(2), copy), atol=2e-15)
    directions = [z, x, x @ z, x @ z @ z]
    bases = []
    for gate in directions:
        cube_root_phases = OMEGA ** np.arange(3)
        eigenvalues, vectors = np.linalg.eig(gate)
        indices = [np.argmin(abs(eigenvalues - phase)) for phase in cube_root_phases]
        basis = vectors[:, indices]
        np.testing.assert_allclose(basis.conj().T @ basis, np.eye(3), atol=2e-15)
        bases.append(basis)
        for character in range(3):
            allowed = np.fromfunction(lambda i, j: (i + j - character) % 3 == 0, (3, 3))
            free = np.where(allowed, np.full((3, 3), 1 + 2j), 0)
            np.fill_diagonal(free, np.diag(free) * (3 + 4j))
            physical = basis.conj() @ free @ basis.conj().T
            np.testing.assert_allclose(gate.T @ physical @ gate, OMEGA ** character * physical, atol=2e-14)
            np.testing.assert_allclose(np.linalg.svd(physical, compute_uv=False), np.sqrt(5) * np.array([5, 1, 1]))
    for left, right in itertools.combinations(bases, 2):
        overlap = left.conj().T @ right
        np.testing.assert_allclose(abs(overlap) ** 2, np.full((3, 3), 1/3), atol=2e-15)
        quartet = overlap[0, 0] * overlap[1, 1] * (overlap[0, 1] * overlap[1, 0]).conj()
        np.testing.assert_allclose(abs(quartet.imag), 1 / (6 * math.sqrt(3)))
    return {"distinct_projective_lines": len(bases), "distinct_line_pairs_tested": 6,
            "flat_modulus": 1 / math.sqrt(3), "weak_gram_only_abs_J": 1 / (6 * math.sqrt(3))}


def two_characters():
    up = np.array([[0, 1, 0], [1, 2, 0], [0, 0, 4]], dtype=complex)
    down = np.array([[0, 1, 1], [1, 3, 0], [1, 0, 5j]], dtype=complex)
    grams = [y @ y.conj().T for y in (up, down)]
    commutator = grams[0] @ grams[1] - grams[1] @ grams[0]
    eigensystems = [np.linalg.eigh(h) for h in grams]
    overlap = eigensystems[0][1].conj().T @ eigensystems[1][1]
    quartet = overlap[0, 0] * overlap[1, 1] * (overlap[0, 1] * overlap[1, 0]).conj()
    assert abs(quartet.imag) > 1e-7
    return {"positive_Yu_singular_values": np.sqrt(eigensystems[0][0]).tolist(),
            "positive_Yd_singular_values": np.sqrt(eigensystems[1][0]).tolist(),
            "imaginary_determinant_of_gram_commutator": float(np.linalg.det(commutator).imag),
            "absolute_quartet_imaginary_part": float(abs(quartet.imag))}


def compute():
    cfg = read("inputs.json")
    count = cfg["individual_half_widths"]
    response = cfg["extracted_over_model_response"]
    up = read(cfg["one_TeV_up"])["one_TeV_DR"]["up"]
    down = read(cfg["one_TeV_down"])["down_yukawa_dr_1TeV"]
    cb = response_box(read(cfg["tree_BelleII"])["BelleII_Vcb_BGL"], count, response)
    ratio = response_box(read(cfg["tree_heavy_baryon"])["ub_over_cb"], count, response)
    sin_angles = [c * math.hypot(1, r) for c, r in zip(cb, ratio)]
    angles = [math.asin(s) for s in sin_angles]
    target = {}
    for tag, start in (("same_line_same_heavy", 0.),
                       ("different_lines", math.acos(1 / math.sqrt(3))),
                       ("same_line_other_heavy", math.pi / 2)):
        change = max(0., start - angles[1], angles[0] - start)
        sine = math.sin(change / 2)
        target[tag] = {"sum_of_heavy_line_rotations_deg_at_least": math.degrees(change),
                       "max_operator_error_over_physical_high_gap_at_least": sine / (1 + sine)}
    return {"translation_algebra": translations(),
            "up": spectrum_and_perturbation(up, count),
            "down": spectrum_and_perturbation(down, count),
            "heavy_tree_inputs": {"cb_box": cb, "ub_over_cb_box": ratio,
                                  "tb_from_column_unitarity": [math.cos(angles[1]), math.cos(angles[0])],
                                  "physical_heavy_angle_deg": list(map(math.degrees, angles)),
                                  "absolute_J_upper_from_two_tree_inputs": cb[1] ** 2 * ratio[1] / 2},
            "residual_breaking": target, "two_character_competitor": two_characters()}


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    else:
        np.testing.assert_allclose(a, b, rtol=2e-10, atol=1e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    results = compute()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf8") as target:
            json.dump(results, target, indent=2)
            target.write("\n")
    else:
        compare(results, read("results.json"))
    print("Verified the magnetic translations, spectrum, heavy-line bounds and CP competitor.")
