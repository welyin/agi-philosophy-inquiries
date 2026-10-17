"""Independent finite diagnostics for round 1234; no actual permissions added.

Uses the older 1111 construction plus its explicit work-mass offset, checks it
against the 1112 material, and constructs separate energy-shell comparisons.
All-n reference stability and fixed-domain conservativity are analytic claims,
not inferred from these finite diagnostics. Default is a read-only comparison.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SOURCES = {
    "archive_1106_1119/1111/check.py": "704ae4ca016c77e0df83fc9c2d39e90d79e9fbf3d6469c8c5b9bf1dad8a5bd44",
    "archive_1106_1119/1112/check.py": "6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f",
    "archive_1161_1200/1166/proof.md": "02335c08c46865b438b075a0d8475e143fc113e01eb711f725b96e4290bb0812",
    "archive_1201_1218/1209/proof.md": "9c2eecc8d9a372996f3c317f65613ecbecfad54486fc4a5f23bcf72fa7897e47",
    "archive_1201_1218/_shared/thermal_roles_adoption_after1218.md": "0bfd3f7494d731874ffbad260f83ba66885dbd0e4841efbc2a430b03c649553c",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, BASE / relative)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def small(value, label, tolerance=2e-10):
    value = float(np.linalg.norm(value))
    assert value < tolerance, (label, value)
    return value


def entropy(rho):
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    assert vals.min() > -2e-12
    vals = vals[vals > 1e-14]
    return float(-vals @ np.log(vals))


def tdist(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)


def same_material():
    earlier = module("independent1234_1111", "archive_1106_1119/1111/check.py").model()
    later = module("independent1234_1112", "archive_1106_1119/1112/check.py").Material(1)
    # 1166 uses the complete work-extended material, with rest offset +2.
    h = earlier["hout"] + 2 * np.eye(288)
    expected_h = later.baseline * np.eye(288) + later.hfk
    material_error = small(h - expected_h, "same full 1166 physical H")
    eig, vec = np.linalg.eigh(h)
    expected_eig = np.repeat(later.baseline + np.arange(4), 72)
    spectrum_error = small(eig - expected_eig, "four full energy layers")
    rows = []
    for beta in (math.log(2), 0.0, -math.log(3)):
        weights = np.exp(-beta * (eig - eig.min()))
        weights /= weights.sum()
        rho = (vec * weights) @ vec.conj().T
        residual = small(h @ rho - rho @ h, "same H stationary Gibbs family")
        assert np.linalg.eigvalsh(rho).min() > 0
        rows.append(dict(beta=beta, trace=float(np.trace(rho).real),
                         minimum_eigenvalue=float(weights.min()),
                         commutator_residual=residual,
                         complete_mean_energy=float(np.trace(h @ rho).real)))
    beta = math.log(2)
    weights = np.exp(-beta * (eig - eig.min())); weights /= weights.sum()
    rho = (vec * weights) @ vec.conj().T
    expected_probs = [Fraction(8, 15), Fraction(4, 15), Fraction(2, 15), Fraction(1, 15)]
    exact_mean = sum(Fraction(k) * p for k, p in enumerate(expected_probs))
    assert exact_mean == Fraction(11, 15)
    mean_error = small(np.trace(h @ rho).real - later.baseline - float(exact_mean), "complete mean")
    # Independently assemble the old rank-three captured ready code.
    gamma = earlier["gamma"]
    hc = earlier["hc"]
    ev, vv = np.linalg.eigh(hc)
    cap = vv @ (np.exp(-1j * (math.pi / 4) * ev) * (vv.conj().T @ np.eye(4)[:, 0]))
    aux = np.kron(np.array([1, 1j]) / math.sqrt(2), np.kron([1, 0], [1, 0]))
    columns = np.column_stack([gamma @ np.kron(cap, np.kron(aux, np.kron(np.eye(3)[:, q], [1, 0, 0])))
                               for q in range(3)])
    small(columns.conj().T @ columns - np.eye(3), "rank-three original prepared code")
    p = float(np.trace(columns.conj().T @ rho @ columns).real)
    support_error = small(p - float(Fraction(3, 320)), "thermal mass on original code")
    exact_clock_spectral_weights = [Fraction(1, 8), Fraction(3, 8), Fraction(3, 8), Fraction(1, 8)]
    exact_p = Fraction(3, 72) * sum(a * b for a, b in zip(expected_probs, exact_clock_spectral_weights))
    assert exact_p == Fraction(3, 320)
    return dict(full_dimension=288, complete_H_difference=material_error,
                spectrum_difference=spectrum_error, family_checks=rows,
                beta_log2_layer_probabilities=[str(p) for p in expected_probs],
                full_mean_energy_above_independent_baseline=str(exact_mean), mean_error=mean_error,
                original_code_rank=3, Gibbs_probability_on_old_code=str(exact_p),
                support_probability_error=support_error,
                all_old_code_states_half_trace_distance_lower_bound=str(1 - exact_p),
                this_does_not_assert_no_other_actual_preparation_exists=True)


def mathematical_comparisons():
    rng = np.random.default_rng(1234071)
    rows = []
    for n in (1, 2):
        energies = np.zeros(1, dtype=int)
        for _ in range(n):
            energies = (energies[:, None] + np.arange(4)[None, :]).ravel()
        ds, da = len(energies), 4
        aux_e = np.arange(da)
        total_e = (energies[:, None] + aux_e).ravel()
        unitary = np.zeros((ds * da, ds * da), complex)
        for level in np.unique(total_e):
            ids = np.where(total_e == level)[0]
            z = rng.normal(size=(len(ids), len(ids))) + 1j * rng.normal(size=(len(ids), len(ids)))
            q, _ = np.linalg.qr(z)
            unitary[np.ix_(ids, ids)] = q
        small(unitary.conj().T @ unitary - np.eye(ds * da), "full finite shell unitary")
        comm = small(unitary * total_e[None, :] - total_e[:, None] * unitary, "strict total energy")
        z = rng.normal(size=(da, da)) + 1j * rng.normal(size=(da, da))
        b = z @ z.conj().T + np.eye(da)
        b /= np.trace(b)
        for beta in (math.log(2), 0.0, -math.log(3)):
            p = np.exp(-beta * energies); p /= p.sum()
            reference = np.diag(p)
            original = np.kron(reference, b)
            omega = unitary @ original @ unitary.conj().T
            cube = omega.reshape(ds, da, ds, da)
            omega_s = np.trace(cube, axis1=1, axis2=3)
            omega_a = np.trace(cube, axis1=0, axis2=2)
            de_s = float(energies @ (omega_s.diagonal().real - p))
            de_a = float(aux_e @ (omega_a.diagonal().real - b.diagonal().real))
            ds_s = entropy(omega_s) - entropy(reference)
            ds_a = entropy(omega_a) - entropy(b)
            mutual = entropy(omega_s) + entropy(omega_a) - entropy(omega)
            relative = -entropy(omega_s) - float(omega_s.diagonal().real @ np.log(p))
            energy_error = small(de_s + de_a, "full energy account")
            identity_error = small(relative - (-beta * de_a + ds_a - mutual), "TR exact entropy account")
            distance = tdist(omega_s, reference)
            bound = math.sqrt(max(0.0, (abs(beta) * abs(de_a) + abs(ds_a)) / 2))
            assert distance <= bound + 1e-10
            rows.append(dict(n=n, beta=beta, factored_complete_material_dimension=288 ** n,
                             explicitly_computed_clock_aux_dimension=ds * da,
                             degeneracy_identity_dimension=72 ** n,
                             strict_energy_commutator=comm, energy_account_residual=energy_error,
                             entropy_identity_residual=identity_error, relative_entropy=relative,
                             retained_mutual_information=mutual,
                             auxiliary_state_change=tdist(omega_a, b),
                             material_half_trace_distance=distance, TR_budget_bound=bound,
                             actual_operation_permission_claimed=False))
    return rows


def deletion_and_budget():
    p = np.array([8., 4., 2., 1.]) / 15
    wrong_energies = np.array([0., 1., 2., 4.])
    design = np.column_stack([np.ones(4), wrong_energies])
    residual = float(np.linalg.norm(design @ np.linalg.lstsq(design, -np.log(p), rcond=None)[0] + np.log(p)))
    assert residual > .2
    # A mere bounded reference label does not mean uniform apparatus resources.
    budget_rows = []
    for power in (2, 6, 12):
        delta = Fraction(1, 10 ** power)
        kappa = Fraction(2, 10 ** power)
        b = math.sqrt((math.log(2) * float(delta) + float(kappa)) / 2)
        budget_rows.append(dict(delta=str(delta), kappa=str(kappa), stability_bound=b))
    assert budget_rows[-1]["stability_bound"] < 2e-6
    return dict(wrong_physical_H_Gibbs_fit_residual=residual,
                independently_fixed_H_must_be_retained=True,
                all_n_analytic_bound_budget_examples=budget_rows,
                actual_state_preparation_or_motion_not_inferred=True)


def calculate():
    for rel, expected in SOURCES.items():
        assert sha(BASE / rel) == expected, rel
    result = dict(schema="round1234_independent_diagnostics_v1", status="PASS_WITH_SCOPE",
                  independent_check_groups=3, code_sha256=sha(Path(__file__)), source_sha256=SOURCES,
                  complete_material_and_nonempty_ideal_reference=same_material(),
                  ideal_complete_energy_shell_comparisons=mathematical_comparisons(),
                  deletion_and_budget=deletion_and_budget(),
                  scope=dict(no_new_actual_permission=True, old_actual_domain_not_simulated_anew=True,
                             universal_semantics_requires_analytic_review=True,
                             new_adopted_axioms=0, GR_deliverables_completed=0))
    for rel, expected in SOURCES.items():
        assert sha(BASE / rel) == expected, rel
    return result


def compare(actual, saved, path="root"):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys(), path
        for key in actual:
            compare(actual[key], saved[key], path + "." + key)
    elif isinstance(actual, list):
        assert len(actual) == len(saved), path
        for index, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, f"{path}[{index}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, abs_tol=2e-9, rel_tol=2e-9), path
    else:
        assert actual == saved, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-once", action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = HERE / "independent_results.json"
    if args.write_once:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=result["status"], groups=3,
                          mode="WRITE_ONCE" if args.write_once else "READ_ONLY_MATCH",
                          same_H_dimension=288, old_code_distance_lower_bound="317/320",
                          new_actual_permissions=0, GR_deliverables_completed=0)))


if __name__ == "__main__":
    main()
