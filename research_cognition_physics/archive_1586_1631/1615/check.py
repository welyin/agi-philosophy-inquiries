"""Check discrete flavour lifts on SM zero modes and the ordered duql operator."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def group_elements(order):
    return [(k, reflected) for reflected in (0, 1) for k in range(order)]


def planar(order, harmonic, k, reflected):
    angle = 2 * np.pi * harmonic * k / order
    cosine, sine = np.cos(angle), np.sin(angle)
    rotation = np.array([[cosine, -sine], [sine, cosine]])
    return rotation @ np.diag([1, (-1) ** reflected])


def family(order, harmonic, k, reflected, rotation_sign=1, reflection_sign=1):
    result = np.zeros((3, 3))
    result[0, 0] = rotation_sign ** k * reflection_sign ** reflected
    result[1:, 1:] = planar(order, harmonic, k, reflected)
    return result


def sign_det(matrix):
    return round(np.linalg.det(matrix))


def reynolds(quarks, leptons):
    average = np.zeros((81, 81))
    for quark, lepton in zip(quarks, leptons):
        pair = np.kron(quark, quark)
        average += np.kron(pair, np.kron(quark, lepton))
    average /= len(quarks)
    assert np.max(np.abs(average @ average - average)) < 1e-12
    assert np.max(np.abs(average - average.T)) < 1e-12
    return average


def calc(config):
    order = config["polygon_order"]
    elements = group_elements(order)
    original = [family(order, config["quark_planar_harmonic"], k, f) for k, f in elements]
    odd = [(-1)**f * a for (_, f), a in zip(elements, original)]
    inert = [np.eye(3) for _ in elements]
    projectors = [reynolds(q, inert) for q in (original, odd)]
    all_singlets = np.zeros((3, 3, 3, 3))
    all_singlets[0, 0, 0, 0] = 1
    singlet = all_singlets.flatten()
    area = np.zeros_like(all_singlets)
    area[1, 2, 0, 0], area[2, 1, 0, 0] = 1, -1
    dot = np.zeros_like(area)
    dot[1, 1, 0, 0] = dot[2, 2, 0, 0] = 1
    tensors = (singlet, dot.flatten(), area.flatten())
    survival = [[round(np.linalg.norm(p @ v)**2 / np.linalg.norm(v)**2)
                 for v in tensors] for p in projectors]
    assert survival == [[1, 1, 0], [0, 0, 1]]
    ranks = [round(np.trace(p)) for p in projectors]
    assert ranks == [12, 9]
    y = np.array([[1j, 2, 0], [1, -3, 0.1j], [2j, 0.5, 4]])
    for a, b in zip(original, odd):
        np.testing.assert_allclose(a.T @ y @ a, b.T @ y @ b)
    copies = config["weak_colours"]
    phases = []
    for q in (original, odd):
        weak = [sign_det(a)**copies for a in q]
        strong = [sign_det(a)**4 for a in q]
        phases.append({"weak_negative_elements": weak.count(-1),
                       "strong_negative_elements": strong.count(-1)})
    assert phases[0]["weak_negative_elements"] == order
    assert phases[1]["weak_negative_elements"] == 0
    assert all(item["strong_negative_elements"] == 0 for item in phases)

    common = []
    for rot in config["lepton_singlet_rotation"]:
        for ref in config["lepton_singlet_reflection"]:
            for j in config["lepton_doublet_harmonics"]:
                les = [family(order, j, k, f, rot, ref) for k, f in elements]
                quarks = [rot**k * ref**f * a for (k, f), a in zip(elements, original)]
                for q, ell in zip(quarks, les):
                    assert sign_det(q)**copies * sign_det(ell) == 1
                    assert sign_det(q)**4 == 1
                dim = round(np.trace(reynolds(quarks, les)))
                harmonic = j if rot == 1 else order//2 - j
                predicted = 4 + {1: 6, 2: 3, 3: 1}.get(harmonic, 0)
                assert dim == predicted
                common.append({"lepton_singlet_r": rot, "lepton_singlet_s": ref,
                               "lepton_planar_j": j, "invariant_duql_dimension": dim})
    # Keep only the planar lepton of the independently tested muon assignment.
    # Two identical scalar insertions transform in the symmetric square.
    muon_counts = np.zeros(3)
    for (k, f), quark in zip(elements, odd):
        lepton = planar(order, 5, k, f)
        scalar = planar(order, 1, k, f)
        sym2 = (np.trace(scalar)**2 + np.trace(scalar @ scalar))/2
        muon_counts += np.trace(quark)**3 * np.trace(lepton) * np.array([1, np.trace(scalar), sym2])
    muon_counts = np.rint(muon_counts / len(elements)).astype(int).tolist()
    assert muon_counts == [0, 0, 1]
    return {"muon_planar_duql_with_zero_one_two_fundamental_insertions": muon_counts,
            "inert_leptons": {"duql_dimensions": ranks, "instantons": phases,
                              "all_singlets_dot_area_projector_survival": survival},
            "planar_leptons": common,
            "plane_two_slots_squared_neighbour_ratio_area_over_orthogonal":
                round(float(np.sin(np.pi/order)**2), 12)}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = calc(json.loads((HERE / "inputs.json").read_text(encoding="utf-8")))
    output = HERE / "results.json"
    if args.write:
        output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    else:
        assert json.loads(output.read_text(encoding="utf-8")) == data
    print("1615: flavour lifts, both instanton phases and ordered duql projectors verified")


if __name__ == "__main__":
    main()
