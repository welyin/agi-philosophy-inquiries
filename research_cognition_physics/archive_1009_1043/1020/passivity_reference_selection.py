"""1020: exact finite-copy passivity and reference-representation boundaries.

Occupation triples represent every product-basis population including its
degeneracy; equal-energy reordering is not a work-extraction witness. GNS
generators are constructed as actual matrices on the supported HS space.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE/"passivity_reference_selection_results.json"


def exact(value):
    return str(F(value))


def occupations(n, eta):
    z = F(3, 2)+eta/4
    rows = []
    for n2 in range(n+1):
        for n1 in range(n-n2+1):
            n0 = n-n1-n2
            multiplicity = math.comb(n, n2)*math.comb(n-n2, n1)
            weight = F(1, 2)**n1*(eta/4)**n2
            rows.append(((n0, n1, n2), n1+2*n2, weight, multiplicity))
    assert sum(item[3] for item in rows) == 3**n
    assert sum((weight*mult for _, _, weight, mult in rows), F(0)) == z**n
    return rows, z


def population_passivity(n, eta):
    rows, z = occupations(n, eta)
    # Compare every pair of occupied product-basis population classes with
    # different energies. No 3**n matrix and no sorting inside degeneracies.
    compared = 0
    violations = 0
    worst_gain = F(0)
    worst_pair = None
    for low in rows:
        for high in rows:
            if low[1] >= high[1]:
                continue
            compared += 1
            if low[2] < high[2]:
                violations += 1
                gain = (high[1]-low[1])*(high[2]-low[2])/(z**n)
                if gain > worst_gain:
                    worst_gain = gain
                    worst_pair = dict(low_occupation=list(low[0]), low_energy=low[1],
                        high_occupation=list(high[0]), high_energy=high[1])
    # Separate adjacent-energy population-envelope check.
    groups = {}
    for _, energy, weight, _ in rows:
        groups.setdefault(energy, []).append(weight)
    envelope_passive = all(min(groups[e]) >= max(groups[e+1]) for e in range(2*n))
    passive = violations == 0
    assert envelope_passive == passive
    criterion = eta**((n-1)//2) >= F(1, 2)
    assert passive == criterion
    return dict(copies=n, occupation_classes=len(rows), represented_basis_states=3**n,
                distinct_energy_comparisons=compared, inversion_pairs=violations,
                passive=passive, adjacent_energy_envelopes_passive=envelope_passive,
                exact_criterion=criterion,
                eta_to_floor_n_minus1_over2=exact(eta**((n-1)//2)),
                largest_single_configuration_swap_work=exact(worst_gain),
                worst_pair=worst_pair,
                equal_energy_pairs_excluded=True)


def activation_case(eta, guaranteed_window=None):
    k = 1
    while eta**k >= F(1, 2):
        k += 1
    first = 2*k+1
    scans = [population_passivity(n, eta) for n in range(1, first+2)]
    assert all(row["passive"] for row in scans[:first-1])
    assert not scans[first-1]["passive"] and not scans[first]["passive"]
    z = F(3, 2)+eta/4
    high_p = F(1, 2)**first/z**first
    low_p = (eta/4)**k/z**first
    work = high_p-low_p
    assert work > 0 and eta**(k-1) >= F(1, 2) > eta**k
    if guaranteed_window is not None:
        assert all(row["passive"] for row in scans[:guaranteed_window])
    return dict(eta=exact(eta), beta="log(2)", delta="-log(eta)",
                guaranteed_finite_window=guaranteed_window, k_first_strict_crossing=k,
                first_active_copy_count=first,
                exact_last_nonstrict_power=exact(eta**(k-1)),
                exact_first_strict_power=exact(eta**k), scans=scans,
                normalized_single_copy_population=[exact(F(1)/z), exact(F(1, 2)/z), exact(eta/(4*z))],
                exact_swap_witness=dict(high_occupation=[0, first, 0], high_energy=first,
                    low_occupation=[k+1, 0, k], low_energy=2*k,
                    high_individual_population=exact(high_p), low_individual_population=exact(low_p),
                    extracted_work=exact(work), extracted_work_decimal=float(work),
                    degeneracy_factor_not_multiplied_into_individual_population=True,
                    physical_battery_or_cyclic_controller_implemented=False))


def gns_matrix(energies, support, probabilities):
    h = np.diag(np.asarray(energies, float))
    d = len(energies)
    basis = []
    for j in support:
        for i in range(d):
            b = np.zeros((d, d), complex); b[i, j] = 1
            basis.append(b)
    L = np.zeros((len(basis), len(basis)), complex)
    for col, b in enumerate(basis):
        action = h @ b-b @ h
        for row, a in enumerate(basis):
            L[row, col] = np.vdot(a, action)
    rho = np.diag(np.asarray(probabilities, float))
    omega_matrix = np.diag(np.sqrt(np.asarray(probabilities, float)))
    omega = np.array([np.vdot(b, omega_matrix) for b in basis])
    expected = sorted(float(energies[i]-energies[j]) for j in support for i in range(d))
    actual = np.linalg.eigvalsh(L)
    assert np.max(np.abs(actual-np.array(expected))) < 2e-14
    assert np.linalg.norm(L @ omega) < 2e-14
    assert abs(np.vdot(omega, omega)-1) < 2e-14
    assert np.linalg.norm(h @ rho-rho @ h) < 2e-14
    positive = bool(min(actual) >= -2e-14)
    ground_support = all(energies[j] == min(energies) for j in support)
    assert positive == ground_support
    return dict(energies=list(energies), support=list(support),
                probabilities=list(map(float, probabilities)), GNS_dimension=len(basis),
                matrix_eigenvalues=actual.tolist(), expected_energy_differences=expected,
                GNS_generator_positive=positive, support_in_ground_eigenspace=ground_support,
                stationary=True, physical_H_positive=min(energies) >= 0,
                state_purity=float(np.trace(rho @ rho).real)), L, omega


def reference_spectrum_cases():
    cases = []
    for energies in ([0, 1], [0, 0, 2], [0, 1, 1, 3], [0, 0, 0]):
        d = len(energies)
        for mask in range(1, 2**d):
            support = [j for j in range(d) if mask & (1 << j)]
            total = sum(j+1 for j in support)
            probabilities = [F(j+1, total) if j in support else F(0) for j in range(d)]
            case, _, _ = gns_matrix(energies, support, probabilities)
            cases.append(case)
    assert len(cases) == 32
    mixed_ground = [row for row in cases if len(row["support"]) > 1 and row["GNS_generator_positive"]]
    assert mixed_ground and all(row["state_purity"] < 1 for row in mixed_ground)
    thermal, L, omega = gns_matrix([0, 1], [0, 1], [F(2, 3), F(1, 3)])
    assert thermal["matrix_eigenvalues"] == [-1., 0., 0., 1.]
    thermal_copy_samples = []
    for n in range(1, 7):
        p = [F(2, 3)**(n-k)*F(1, 3)**k for k in range(n+1)]
        assert all(p[k] >= p[k+1] for k in range(n))
        assert sum((math.comb(n,k)*p[k] for k in range(n+1)), F(0)) == 1
        thermal_copy_samples.append(dict(copies=n, energy_levels=list(range(n+1)),
                                         individual_populations=list(map(exact, p)), passive=True))
    return dict(support_enumeration=cases, support_cases=len(cases),
                mixed_ground_positive_examples=len(mixed_ground),
                equivalence_proof_scope="finite full matrix algebra, stationary density matrix, supported Hilbert-Schmidt GNS representation",
                thermal_qubit=thermal, finite_thermal_copy_calibrations=thermal_copy_samples,
                thermal_all_copy_passivity_is_analytic_Gibbs_identity=True,
                canonical_GNS_generator_not_same_as_physical_H=True,
                full_QFT_forward_cone_derived=False)


def purification_case():
    h = np.diag([0., 1.]); eye = np.eye(2)
    minus = np.kron(h, eye)-np.kron(eye, h)
    plus = np.kron(h, eye)+np.kron(eye, h)
    omega = np.array([math.sqrt(2/3), 0., 0., math.sqrt(1/3)])
    pure = np.outer(omega, omega)
    reduced = np.einsum("ijkj->ik", pure.reshape(2, 2, 2, 2))
    mean = float(omega @ plus @ omega)
    variance = float(omega @ plus @ plus @ omega)-mean**2
    assert np.linalg.norm(minus @ omega) == 0
    assert abs(variance-8/9) < 2e-14
    assert abs(float(np.trace(pure @ pure))-1) < 2e-14
    assert abs(float(np.trace(reduced @ reduced))-5/9) < 2e-14
    assert np.max(np.abs(reduced-np.diag([2/3, 1/3]))) < 2e-14
    return dict(omega=omega.tolist(), purification_purity=float(np.trace(pure @ pure)),
                physical_reduced_state=reduced.tolist(), physical_state_purity=float(np.trace(reduced @ reduced)),
                minus_generator_eigenvalues=np.linalg.eigvalsh(minus).tolist(),
                minus_generator_fixes_omega=True,
                plus_generator_eigenvalues=np.linalg.eigvalsh(plus).tolist(),
                plus_mean_energy=mean, plus_energy_variance=variance,
                plus_generator_fixes_omega=False,
                vector_purity_does_not_imply_algebraic_state_purity=True,
                auxiliary_copy_not_treated_as_free_physical_resource=True)


def remaining_gap_freedom():
    rows = []
    for lam in (F(2), F(3), F(5, 2)):
        energies = [F(0), F(1), lam]
        weights = np.exp(-math.log(2)*np.array(energies, float))
        probabilities = weights/weights.sum()
        state = np.array([1., 0., 1.], complex)/math.sqrt(2)
        evolved = np.exp(-1j*math.pi*np.array(energies, float))*state
        plus_probability = float(abs(np.vdot(state, evolved))**2)
        assert abs(plus_probability-(1+math.cos(math.pi*float(lam)))/2) < 2e-14
        rows.append(dict(lambda_gap_ratio=exact(lam), physical_energies=list(map(exact, energies)),
            beta="log(2)", Gibbs_reference=probabilities.tolist(),
            fixed_clock_time="pi", common_preparation="(|0>+|2>)/sqrt(2)",
            common_plus_read_probability=plus_probability,
            Gibbs_all_copy_passivity_by_analytic_identity=True,
            all_FUCP_and_all_physical_recovery_claimed=False))
    assert np.max(np.abs(np.array([r["common_plus_read_probability"] for r in rows])-[1., 0., .5])) < 2e-14
    return rows


def compare(fresh, saved, path="root"):
    if isinstance(fresh, float):
        assert isinstance(saved, (int, float)) and math.isclose(fresh, saved, rel_tol=3e-11, abs_tol=2e-13), path
    elif isinstance(fresh, dict):
        assert isinstance(saved, dict) and fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert isinstance(saved, list) and len(fresh) == len(saved), path
        for i, (a, b) in enumerate(zip(fresh, saved)):
            compare(a, b, path+f"[{i}]")
    else:
        assert fresh == saved, (path, fresh, saved)


def run():
    families = [activation_case(F(2*N, 2*N+1), N) for N in (1, 2, 3, 6)]
    assert [row["first_active_copy_count"] for row in families] == [5, 9, 11, 19]
    equality = activation_case(F(1, 2))
    assert equality["scans"][2]["passive"] and equality["scans"][3]["passive"]
    assert not equality["scans"][4]["passive"]
    sources = [BASE/"archive_1009_/research_note_1019.md",
        BASE/"archive_1009_/1019/spin_statistics_locality_results.json",
        BASE/"archive_301_341/research_note_305.md",
        BASE/"archive_370_428/research_note_422.md",
        BASE/"archive_923_934/research_note_929.md",
        BASE/"archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
        BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return dict(round=1020, new_calibration_groups=1, cumulative_test_groups=3798,
        new_cognitive_axioms=0, all_scientific_calibrations_passed=True,
        scope="fixed finite Hamiltonians, full cyclic-unitary operation menu, exact finite-copy diagonal passivity and canonical GNS generator on the actual reference support",
        complete_passivity_generated_from_internal_resource_accounting=False,
        mature_complete_passivity_theorem_claimed_as_new=False,
        finite_copy_samples_claimed_to_prove_complete_passivity=False,
        exact_family_general_criterion="delta*floor((n-1)/2)<=beta; first active n=2*floor(beta/delta)+3 for beta,delta>0",
        exact_finite_copy_families=families, equality_boundary=equality,
        GNS_reference_spectrum=reference_spectrum_cases(),
        purification_boundary=purification_case(), remaining_gap_freedom=remaining_gap_freedom(),
        physical_H_positive_in_all_matrix_examples=True,
        unrestricted_operation_menu_is_an_extra_input=True,
        all_independent_copies_is_stronger_than_fixed_finite_task_budget=True,
        positive_reference_generator_requires_ground_support_in_stated_matrix_scope=True,
        no_unique_gap_ratio_selected=True, full_FUCP_nonimplication_proved=False,
        analytic_obligations=[
            "occupation energy-population ordering gives the exact all-n condition for the declared qutrit family",
            "same energy degeneracies cannot yield work by a population swap",
            "Gibbs tensor powers are passive for every finite copy count by the analytic Gibbs identity",
            "on supported HS space L=H_left-H_right is positive iff reference support lies in the ground eigenspace",
            "positive physical H, positive canonical GNS L, and pure purification vector are different conditions",
            "complete passivity relative to independently chosen H does not fix its measurable gap ratios"],
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\", "/"):
                                  hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); start = time.perf_counter(); result = run()
    if args.write:
        serialized = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
        with OUT.open("x", encoding="utf8") as stream:
            stream.write(serialized)
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1020, passed=result["all_scientific_calibrations_passed"],
        first_active_copies=[row["first_active_copy_count"] for row in result["exact_finite_copy_families"]],
        work=[row["exact_swap_witness"]["extracted_work_decimal"] for row in result["exact_finite_copy_families"]],
        equality_first_active=result["equality_boundary"]["first_active_copy_count"],
        GNS_support_cases=result["GNS_reference_spectrum"]["support_cases"],
        elapsed_seconds=time.perf_counter()-start), indent=2))
