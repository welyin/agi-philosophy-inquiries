"""1201 independent candidate audit; native-basis 16D QQ dynamics.
Default: recompute, compare saved JSON, print JSON; no writes.
--write: create the initial results file, never overwrite it.
Floating diagnostics only: not the main certified/resource-complete admission.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "independent_results_1201.json"


def op(a):
    return float(np.linalg.norm(a, 2))


def entropy(rho):
    values = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    assert min(values) > -1e-12
    values = values[values > 1e-14]
    return float(-np.sum(values * np.log(values)))


def marginals(rho):
    t = rho.reshape(4, 4, 4, 4)
    return np.einsum("abcb->ac", t), np.einsum("abad->bd", t)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-8, abs_tol=2e-11), (a, b)
    else:
        assert a == b, (a, b)


def run():
    # Reconstruct native matrices explicitly; do not import 980 or main 1201 code.
    # Native basis agrees with 968/material() and 980's exact h_eta.
    g, time, eta = 0.1, 600.0, 0.05
    h = np.array([[0, -.2, 0, -eta * math.sqrt(3) / 4],
                  [-.2, 1, 0, 0], [0, 0, 1, 0],
                  [-eta * math.sqrt(3) / 4, 0, 0, -eta / 2]], float)
    q = np.zeros((4, 4))
    q[1, 2] = q[2, 1] = 1
    eye = np.eye(4)
    ha, hb = np.kron(h, eye), np.kron(eye, h)
    interaction = np.kron(q, q)
    H = ha + hb + g * interaction
    energies, R = np.linalg.eigh(h)
    values, W = np.linalg.eigh(H)
    U = (W * np.exp(-1j * time * values)) @ W.T
    center = np.eye(4)[:, 2]
    assert abs(energies[2] - 1) < 1e-13
    assert op(h @ center - center) == 0
    assert abs(abs(R[2, 2]) - 1) < 1e-13

    # Input is native |c><c| on A. B population basis is h's ordered energy basis.
    # Output A projectors use that same independently specified h.
    amplitudes = (np.kron(R.T, eye) @ U[:, 8:12] @ R).reshape(4, 4, 4)
    effects = [amplitudes[j].conj().T @ amplitudes[j] for j in range(4)]
    M = np.array([np.diag(effect).real for effect in effects])
    inverse = np.linalg.inv(M)
    inverse8 = np.round(inverse, 8)
    residual8 = float(np.linalg.norm(np.eye(4) - inverse8 @ M, np.inf))
    inverse_norm = float(np.linalg.norm(inverse, np.inf))
    neumann = float(np.linalg.norm(inverse8, np.inf) / (1 - residual8))
    assert max(abs(M.sum(axis=0) - 1)) < 1e-12
    assert min(M.ravel()) > 0
    assert inverse_norm < 3.14 and residual8 < 2e-8
    # This is an explicitly conditional budget, not a certified numerical enclosure.
    # If every true M entry is within 1e-7 of this matrix, Neumann gives < 3.14.
    assumed_entry_error = 1e-7
    q_bound = residual8 + np.linalg.norm(inverse8, np.inf) * 4 * assumed_entry_error
    conditional_inverse_bound = float(np.linalg.norm(inverse8, np.inf) / (1 - q_bound))
    assert conditional_inverse_bound < 3.14

    columns = np.array([[1, -1, 0, 0], [1, 1, -2, 0],
                        [1, 1, 1, -3]], float).T
    columns /= np.linalg.norm(columns, axis=0)
    singular = np.linalg.svd(M @ columns, compute_uv=False)
    states = []
    r_test = np.array([3/8, 1/8, 1/8, 3/8])
    beta = 1.3
    calibration = np.exp(-beta * (energies - energies[0]))
    calibration /= sum(calibration)
    for name, r in [("unknown_stationary_test", r_test),
                    ("beta_1_3_calibration_only", calibration)]:
        rho_b = (R * r) @ R.T
        rho_initial = np.kron(np.outer(center, center), rho_b)
        rho_final = U @ rho_initial @ U.conj().T
        a, b = marginals(rho_final)
        output = np.diag(R.T @ a @ R).real
        da = float(np.trace(h @ a).real - 1)
        db = float(np.trace(h @ (b - rho_b)).real)
        initial_interaction = float(g * np.trace(interaction @ rho_initial).real)
        final_interaction = float(g * np.trace(interaction @ rho_final).real)
        work_off = -final_interaction
        sa, sb = entropy(a), entropy(b)
        sj0, sj1 = entropy(rho_initial), entropy(rho_final)
        energy_residual = da + db + final_interaction - initial_interaction
        assert max(abs(output - M @ r)) < 1e-12
        assert max(abs(inverse @ output - r)) < 1e-11
        assert abs(energy_residual) < 1e-11
        assert abs(sj0 - sj1) < 1e-11
        assert abs(initial_interaction) < 1e-12
        states.append(dict(
            label=name, source_populations=r.tolist(), output=output.tolist(),
            inverse_recovery=(inverse @ output).tolist(),
            population_beta_triplet=[
                float(math.log(r[i] / r[2]) / (energies[2] - energies[i]))
                for i in (0, 1, 3)],
            energy=dict(A_initial=1.0, A_final=1.0 + da,
                B_initial=float(energies @ r), B_final=float(energies @ r) + db,
                delta_A=da, delta_B=db, initial_interaction=initial_interaction,
                final_interaction=final_interaction, switch_on_mean_work=0.0,
                switch_off_mean_work=work_off, balance_residual=energy_residual),
            entropy_nats=dict(joint_initial=sj0, joint_final=sj1,
                A_final=sa, B_initial=sj0, B_final=sb,
                B_change=sb-sj0, mutual_information=sa+sb-sj1),
            offdiagonal_norm_in_energy_basis=dict(
                A=op(R.T @ a @ R - np.diag(output)),
                B=op(R.T @ b @ R - np.diag(np.diag(R.T @ b @ R))))
        ))

    # Exact symmetry of the actual menu, not an added physical assumption.
    # S h S = h, S q S = -q, hence [H, S_A S_B] = 0.
    parity = np.diag([1., 1., -1., 1.])
    parity_joint = np.kron(parity, parity)
    assert op(parity @ h - h @ parity) == 0
    assert op(parity @ q @ parity + q) == 0
    assert op(H @ parity_joint - parity_joint @ H) == 0
    energy_parity = R.T @ parity @ R
    assert max(op(N @ energy_parity - energy_parity @ N) for N in effects) < 1e-12
    X = np.zeros((4, 4))
    X[0, 2] = X[2, 0] = 1
    coherent = np.diag(calibration) + .01 * X
    coherent_output = np.array([np.trace(N @ coherent).real for N in effects])
    kernel_outputs = [float(np.trace(N @ X).real) for N in effects]
    nonstationarity = op(np.diag(energies) @ coherent - coherent @ np.diag(energies))
    assert min(np.linalg.eigvalsh(coherent)) > 0
    assert max(abs(coherent_output - M @ calibration)) < 1e-12
    assert nonstationarity > .01
    # One predeclared holdout, not a search: same preparation and full parent at t=300.
    hold_time = 300.0
    U_hold = (W * np.exp(-1j * hold_time * values)) @ W.T
    amp_hold = (np.kron(R.T, eye) @ U_hold[:, 8:12] @ R).reshape(4, 4, 4)
    M_hold = np.sum(abs(amp_hold)**2, axis=1)
    energy_coefficients = {"delta_A": [], "delta_B": [], "switch_off_mean_work": []}
    for j in range(4):
        psi0 = np.kron(center, R[:, j])
        psi1 = U_hold @ psi0
        da = float(np.vdot(psi1, ha @ psi1).real - 1)
        db = float(np.vdot(psi1, hb @ psi1).real - energies[j])
        wo = float(-g * np.vdot(psi1, interaction @ psi1).real)
        assert abs(da + db - wo) < 1e-11
        energy_coefficients["delta_A"].append(da)
        energy_coefficients["delta_B"].append(db)
        energy_coefficients["switch_off_mean_work"].append(wo)
    holdout_states = []
    for name, r in [("unknown_stationary_test", r_test),
                    ("beta_1_3_calibration_only", calibration)]:
        rho_b = (R * r) @ R.T
        joint0 = np.kron(np.outer(center, center), rho_b)
        joint1 = U_hold @ joint0 @ U_hold.conj().T
        a, b = marginals(joint1)
        direct_p = np.diag(R.T @ a @ R).real
        prediction = {key: float(np.array(coef) @ r)
                      for key, coef in energy_coefficients.items()}
        direct = dict(
            delta_A=float(np.trace(h @ a).real - 1),
            delta_B=float(np.trace(h @ (b-rho_b)).real),
            switch_off_mean_work=float(-g * np.trace(interaction @ joint1).real))
        assert max(abs(direct_p-M_hold @ r)) < 1e-12
        assert max(abs(prediction[key]-direct[key]) for key in direct) < 1e-12
        holdout_states.append(dict(label=name, predicted_probability=(M_hold @ r).tolist(),
                                   direct_probability=direct_p.tolist(),
                                   predicted_energy=prediction, direct_energy=direct))
    assert max(abs(M_hold.sum(axis=0)-1)) < 1e-12
    holdout = dict(time=hold_time, matrix=M_hold.tolist(),
                   energy_coefficients_by_initial_B_population=energy_coefficients,
                   checked_states=holdout_states,
                   control_timing_and_switching_remain_conditional=True)
    old = ROOT / "archive_956_989"
    sources = [old / "968/internal_relay.py", old / "980/finite_thermal_records.py",
               old / "research_note_980.md"]
    assert all(p.is_file() for p in sources)
    return dict(
        status="1201 in-progress candidate; floating audit, not final admission",
        date="2026-10-14",
        all_candidate_assertions_passed=True,
        fixed_parameters=dict(g=g, time=time, eta=eta, A_native_preparation_index=2,
                              joint_dimension=16),
        native_material=dict(h=h.tolist(), Q=q.tolist(), energies=energies.tolist()),
        numerical_diagnostics=dict(
            H_eigen_residual=op(H @ W - W * values),
            H_eigenvector_gram_error=op(W.T @ W - np.eye(16)),
            h_eigen_residual=op(h @ R - R * energies),
            unitary_gram_error=op(U.conj().T @ U - np.eye(16)),
            effect_sum_error=op(sum(effects) - np.eye(4))),
        population_readout=dict(
            matrix_rows_A_output_columns_B_population=M.tolist(),
            full_inverse=inverse.tolist(), inverse_infinity_norm=inverse_norm,
            inverse_rounded_8_decimals=inverse8.tolist(),
            rounded_inverse_residual_infinity=residual8,
            floating_neumann_inverse_upper=neumann,
            zero_sum_singular_values=singular.tolist(),
            conditional_main_certificate_budget=dict(
                assumed_M_entry_error=assumed_entry_error,
                neumann_residual_upper=float(q_bound),
                inverse_infinity_upper=conditional_inverse_bound,
                selected_conservative_bound=3.14,
                assumption_certified_here=False)),
        checked_states=states,
        predeclared_holdout=holdout,
        test_state_all_real_beta_gibbs_distance=dict(
            exact_infinity_distance=0.125,
            argument="beta>=0 requires r[2]>=r[3]; beta<=0 requires r[1]>=r[0]; both violated by 1/4; beta=0 attains distance 1/8"),
        remove_stationarity=dict(
            actual_joint_parity_commutator=op(H @ parity_joint - parity_joint @ H),
            X="|e0><c|+|c><e0|", offdiagonal_coefficient=.01,
            kernel_outputs=kernel_outputs,
            coherent_state_min_eigenvalue=float(min(np.linalg.eigvalsh(coherent))),
            coherent_state_commutator_norm=nonstationarity,
            same_output_max_error=float(max(abs(coherent_output-M @ calibration)))),
        scope=dict(
            original_full_QQ_parent_used=True,
            rotating_wave_or_weak_coupling_replacement=False,
            parameter_scan_in_this_script=False,
            unknown_reference_assumed_gibbs=False,
            quantum_full_marginal_entropies_not_diagonal_shannon=True,
            rigorous_main_numerical_certificate_completed=False,
            original_clock_contact_timing_implemented=False,
            energy_measurement_instrument_implemented=False,
            controller_preparation_cost_completed=False,
            fresh_copy_preparation_cost_completed=False,
            geometry_region_identity_completed=False,
            final_1201_admission=False),
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        if OUT.exists():
            raise FileExistsError("Refusing to overwrite independent results")
        OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
