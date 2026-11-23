"""Fixed-pole, fixed-light-line test of a three-step charged transfer."""
import argparse
import itertools
import json
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent


def prior(round_number, name="results.json"):
    return json.loads((ARCHIVE / str(round_number) / name).read_text(encoding="utf-8"))


def complete_frame(columns):
    """Orthonormal frame adapted to ordered, linearly independent seed columns."""
    result = []
    for vector in list(columns) + list(np.eye(5)):
        orth = np.array(vector, dtype=float, copy=True)
        for fixed in result:
            orth -= fixed.dot(orth) * fixed
        norm = np.linalg.norm(orth)
        if norm > 1e-10:
            result.append(orth / norm)
    assert len(result) == 5
    return np.column_stack(result)


def norm_budget(null, source, exit_):
    p = null[source] * null[exit_]
    return p, np.sqrt(1-null[source]**2), np.sqrt(1-null[exit_]**2)


def align_left(null, source, exit_, residues):
    """Explicit equal-Gram identification in the fixed, tagged five-doublet basis."""
    p, norm_a, norm_b = norm_budget(null, source, exit_)
    total = np.linalg.norm(residues, 1)
    active = np.flatnonzero(abs(residues) > 1e-12)
    spare = next(i for i in range(4) if i not in active)
    spectral_a, spectral_b = np.zeros(5), np.zeros(5)
    spectral_a[0], spectral_b[0] = null[source], null[exit_]
    spectral_a[active+1] = norm_a*np.sqrt(abs(residues[active])/total)
    spectral_b[active+1] = residues[active]/spectral_a[active+1]
    spectral_b[spare+1] = np.sqrt(norm_b**2-(total/norm_a)**2)
    axes = np.eye(5)
    n1 = complete_frame([null, axes[source], axes[exit_]])
    n2 = complete_frame([axes[0], spectral_a, spectral_b])
    left = n1 @ n2.T
    assert np.allclose(left[:, 0], null)
    assert np.allclose(left[source, :] * left[exit_, :], np.r_[p, residues])
    return left


def old_mass(overlap):
    """R1714 equal bare scales, with the H_d row appended before the four old rows."""
    unit_s = np.array([1., 0., 0.])
    unit_u = np.array([overlap, np.sqrt(1.-overlap**2), 0.])
    mass = np.zeros((5, 5))
    mass[1:4, 0], mass[1:4, 1:4] = -unit_s, np.eye(3)
    mass[4, 1:4], mass[4, 4] = -unit_u, 1.
    return mass


def number(x):
    return 0. if abs(x) < 8.e-14 else float(f"{x:.11g}")


def list_numbers(xs):
    return [number(x) for x in xs]


def compute():
    old = prior(1714)
    charged = prior(1732)["exploratory_k3"]["registered_with_lower_us"]
    unit_gev = prior(1714, "inputs.json")["example_heavy_doublet_mass_units_GeV"][0]
    c = charged["Vud_10PeV"]
    null = np.array([1., c, c**2, c**3, 0.])
    null /= np.linalg.norm(null)
    assert abs(null[3]/null[0] - charged["positive_x"]) < 2.e-10
    source, exit_ = 0, 3
    p, norm_a, norm_b = norm_budget(null, source, exit_)
    budget = norm_a * norm_b / p
    mass = old_mass(old["Vus_static"])
    eig, right = np.linalg.eigh(mass.T @ mass)
    assert abs(eig[0]) < 1.e-14
    eig[0] = 0.
    poles = eig[1:]
    archived_gev = np.array(old["same_1712_point"]["heavy_masses_GeV"])
    assert np.allclose(np.sqrt(poles) * unit_gev, archived_gev, atol=.002, rtol=0.)

    choices = []
    for active in itertools.combinations(range(4), 3):
        ds = poles[list(active)]
        weights = np.linalg.solve(np.array([np.ones(3), ds, ds**2]), [-1., 0., 0.])
        residues = np.zeros(4)
        residues[list(active)] = weights
        choices.append((np.linalg.norm(weights, 1), active, residues))
    best, used, three = min(choices, key=lambda entry: entry[0])
    assert sum(weight < budget for weight, _, _ in choices) == 1
    assert best < budget
    direct = np.zeros(4)
    direct[-1] = -1.
    ts = [(1e6 / unit_gev)**2, (1.e7 / unit_gev)**2, 1.]
    cases = {}
    for name, weights, active in [
            ("three_moments", three, used), ("direct_largest_pole", direct, (3,))]:
        left = align_left(null, source, exit_, p*weights)
        extended = left @ np.diag(np.sqrt(eig)) @ right.T
        gram = extended @ extended.T
        assert np.allclose(extended.T @ extended, mass.T @ mass, atol=2e-14)
        assert np.allclose(extended.T @ null, 0., atol=2e-14)
        assert np.allclose(np.linalg.eigvalsh(gram), eig, atol=2e-14)
        moments = np.array([np.dot(weights, poles**n) for n in range(1, 4)])
        expected = len(active)
        assert np.max(np.abs(moments[:expected-1]), initial=0) < 2.e-14
        assert abs(moments[expected-1]) > .01
        shape = []
        slope = np.dot(weights, 1/poles)
        for t in ts:
            form = t * np.linalg.inv(t*np.eye(5) + gram)[exit_, source] / p
            target = float(np.prod(poles[list(active)]/(poles[list(active)]+t)))
            assert abs(form - target) < 1.e-11
            shape.append(number(target))
        cases[name] = {
            "normalized_massive_cross_residues": list_numbers(weights),
            "unit_mass_squared_times_dlogF_dQ_squared_at_zero": number(slope),
            "signed_cross_moments_1_2_3_in_unit_mass_squared": list_numbers(moments),
            "F_over_F0_at_1_10_30_PeV": shape,
            "left_gram_in_tagged_basis_in_bare_mass_squared": [
                list_numbers(row) for row in gram],
        }
    amplitudes = {key: value["F_over_F0_at_1_10_30_PeV"] for key,value in cases.items()}
    ratio = np.array(amplitudes["three_moments"])/amplitudes["direct_largest_pole"]
    joint_cfg = prior(1732, "inputs.json")
    offsets = dict(joint_cfg["exploratory_weak_offsets"]["registered_with_lower_us"])
    offsets.update(zip(("strange", "bottom"), charged["down_strange_bottom_marginal_widths"][1:]))
    engine = runpy.run_path(str(ARCHIVE/"1731"/"check.py"))
    joint = engine["endpoint"]("explicit_open", joint_cfg["log_step"]/2, explicit_offsets=offsets)
    mixing = np.linalg.svd(joint["D"])[0][:, ::-1]
    joint_us = abs(mixing[0, 1])
    assert abs(abs(mixing[0, 0])-c) < 1e-8
    joint_mass = old_mass(joint_us)
    joint_poles = np.linalg.eigvalsh(joint_mass.T @ joint_mass)[1:]
    jp = joint_poles[list(used)]
    joint_w = np.linalg.solve(np.array([np.ones(3), jp, jp**2]), [-1., 0., 0.])
    assert np.linalg.norm(joint_w, 1) < budget
    t10 = (1.e7/unit_gev)**2
    joint_rate = np.prod(jp/(t10+jp))**2/(joint_poles[-1]/(t10+joint_poles[-1]))**2
    return {
        "charged_point": "1732 exploratory_k3.registered_with_lower_us",
        "positive_side": "1714 equal bare doublets; unbroken SUSY quadratic tree level",
        "charged_light_line_site_amplitudes": list_numbers(null),
        "c": c, "exit_to_source_light_ratio": number(null[exit_]/null[source]),
        "pole_unit_GeV": unit_gev, "weak_masses_GeV": archived_gev.tolist(),
        "light_cross_pole_residue": number(p),
        "available_relative_total_variation": number(budget),
        "required_relative_total_variation_by_triple": {
            ",".join(str(i+1) for i in triple): number(value)
            for value,triple,_ in choices
        },
        "minimum_required": number(best),
        "replacing_old_overlap_with_actual_1732_10PeV_charged_directions": {
            "abs_Vus": number(joint_us), "masses_at_same_bare_scale_GeV": list_numbers(np.sqrt(joint_poles)*unit_gev),
            "triple_134_variation": number(np.linalg.norm(joint_w, 1)), "rate_ratio_10PeV": number(joint_rate)},
        "models": cases, "three_over_direct_reduced_rates_1_10_30_PeV": list_numbers(ratio**2),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    dest = HERE / "independent_results.json"
    results = compute()
    if args.check:
        if json.loads(dest.read_text(encoding="utf-8")) != results:
            raise SystemExit("1733 independent results differ")
        print("1733 independent: fixed neutral Gram, poles, light residues and two charged shapes verified")
    else:
        dest.write_text(json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
