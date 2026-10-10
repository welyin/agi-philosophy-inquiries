"""1077 independent diagnostic from archived 480 and a rebuilt fixed reader.

Default stdout JSON; no evidence writes. Spectral/Frechet differentiation and
finite differences are numerical diagnostics, not the author's strict integer
certificate. The reader is rebuilt directly from the original Hamiltonian,
fixed ancillary preparation and fixed source-Z measurement.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts"))
from research_layout import Layout, ResearchRuntime

GROUPS = []


def checked(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    GROUPS.append(name)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tensor(values):
    result = np.ones((1, 1), complex)
    for value in values:
        result = np.kron(result, value)
    return result


def run_independent():
    with ResearchRuntime(Layout(), "late").installed():
        import correlated_reference_local_access as source
        import relational_direction_distance_response as archived_reader
        h = source.frame.system()[1]
        k = source.prior.s12()
        columns = source.pure_columns()
        eigenvalues, eigenvectors = np.linalg.eigh(h)
        p0 = np.array([.2, .4, math.pi/2])
        tau = 1e-10
        qgraph = source.prior.qgraph()
        identity = np.eye(2, dtype=complex)
        x, y, z = archived_reader.P
        pauli = (identity, x, y, z)
        raw = np.asarray([
            tensor([identity, s, identity, identity+x, identity+y, identity+z])
            for s in pauli])
        m = np.kron(tensor([identity, (identity+z)/2,
                            identity, identity, identity, identity]), np.eye(6))

        def uh(t):
            return (eigenvectors*np.exp(-1j*t*eigenvalues)) @ eigenvectors.conj().T

        def state(parameters, pulse_time=tau, derivatives=False):
            t, u, theta = parameters
            right, left = uh(t), uh(u)
            if pulse_time == 0:
                pulse = math.cos(theta)*np.eye(384)-1j*math.sin(theta)*k
                derivative_pulse = -1j*k @ pulse
            else:
                values, vectors = np.linalg.eigh(pulse_time*h+theta*k)
                pulse = (vectors*np.exp(-1j*values)) @ vectors.conj().T
                midpoint = (values[:, None]+values[None, :])/2
                difference = values[:, None]-values[None, :]
                kernel = -1j*np.exp(-1j*midpoint)*np.sinc(difference/(2*math.pi))
                derivative_pulse = vectors @ (
                    kernel*(vectors.conj().T @ k @ vectors)) @ vectors.conj().T
            w = left @ pulse @ right
            output = (w @ columns).reshape(64, 6, 28)
            graph = np.einsum("dga,dha->gh", output, output.conj())/384
            if not derivatives:
                return graph
            dw = [w @ (-1j*h), (-1j*h) @ w, left @ derivative_pulse @ right]
            derivatives_g = []
            for derivative in dw:
                dy = (derivative @ columns).reshape(64, 6, 28)
                term = np.einsum("dga,dha->gh", dy, output.conj())/384
                derivatives_g.append(term+term.conj().T)
            return graph, np.asarray(derivatives_g), output, w

        checked("archived source normalization and original generators",
                columns.shape == (384, 28)
                and abs(np.sum(abs(columns)**2)-384) < 1e-12
                and np.array_equal(h, h.conj().T)
                and np.allclose(k @ k, np.eye(384), atol=1e-14)
                and np.max(abs(h-archived_reader.system()[1])) == 0)

        # Independent fixed reader, graph-order convention inherited from 474.
        reader_t = .7
        reader_u = uh(reader_t)
        reader_o = reader_u.conj().T @ m @ reader_u
        effect_coefficients = np.einsum(
            "diej,ked->kij", reader_o.reshape(64, 6, 64, 6), raw)/64
        joint_effect = sum(np.kron(s, b) for s, b in zip(pauli, effect_coefficients))
        joint_spectrum = np.linalg.eigvalsh(joint_effect)
        checked("rebuilt original fixed reader is a physical joint effect",
                np.max(abs(joint_effect-joint_effect.conj().T)) < 2e-12
                and joint_spectrum[0] > -2e-12
                and joint_spectrum[-1] < 1+2e-12)

        graph, derivatives_g, output, w = state(p0, derivatives=True)
        graph_eigenvalues = np.linalg.eigvalsh(graph)
        graph_coherence = float(np.linalg.norm(graph-np.diag(np.diag(graph))))
        plus, minus = [5, 4, 2], [0, 1, 3]
        pair_total = np.diag(graph).real[plus]+np.diag(graph).real[minus]
        derivative_pair_total = np.asarray([
            np.diag(d).real[plus]+np.diag(d).real[minus] for d in derivatives_g])
        checked("full coherent graph state and matching-reference statistics",
                graph_eigenvalues[0] > .08
                and abs(np.trace(graph)-1) < 2e-12
                and graph_coherence > .38
                and np.max(abs(pair_total-1/3)) < 2e-12
                and np.max(abs(derivative_pair_total)) < 2e-12
                and max(abs(np.trace(d)) for d in derivatives_g) < 2e-12)

        # The two retained classical preparation records remain distinct.
        branch_graphs = [
            np.einsum("dga,dha->gh", output[:, :, sl], output[:, :, sl].conj())/192
            for sl in (slice(0, 4), slice(4, 28))]
        branch_residual = max(
            float(np.max(abs(np.diag(g).real[plus]+np.diag(g).real[minus]-1/3)))
            for g in branch_graphs)
        checked("both original correlated-source branches remain reference-safe",
                branch_residual < 2e-12
                and np.max(abs((branch_graphs[0]+branch_graphs[1])/2-graph)) < 2e-12)

        jq = np.asarray([qgraph @ np.diag(d).real for d in derivatives_g]).T
        je = np.einsum("mgh,khg->mk", effect_coefficients, derivatives_g).real
        effect_center = np.einsum("mgh,hg->m", effect_coefficients, graph).real
        diagonal_derivatives = np.asarray([np.diag(np.diag(d)) for d in derivatives_g])
        diagonal_je = np.einsum(
            "mgh,khg->mk", effect_coefficients, diagonal_derivatives).real
        coherence_je = je-diagonal_je
        old_w, old_jq, old_graph = source.numerical_control()
        coordinate_x = np.diag(graph).real[plus]-np.diag(graph).real[minus]
        current_q = qgraph @ np.diag(graph).real
        checked("existing original distance chart is recovered",
                np.max(abs(jq-old_jq)) < 2e-12
                and np.max(abs(graph-old_graph)) < 2e-12
                and np.max(abs(current_q-(np.ones((3, 3))-np.eye(3)) @ coordinate_x)) < 2e-12)
        checked("coherence is retained rather than replacing the source by its diagonal",
                np.linalg.norm(coherence_je[1:]) > .03)

        finite_differences = []
        for delta in (1e-3, 1e-4, 1e-5):
            numerical_dg = []
            for axis in range(3):
                offset = np.eye(3)[axis]*delta
                numerical_dg.append((state(p0+offset)-state(p0-offset))/(2*delta))
            fd_j = np.einsum(
                "mgh,khg->mk", effect_coefficients, np.asarray(numerical_dg)).real
            finite_differences.append({
                "step": delta, "effect_jacobian": fd_j.tolist(),
                "max_difference_from_Frechet": float(np.max(abs(fd_j-je)))})
        fd_errors = [row["max_difference_from_Frechet"] for row in finite_differences]
        checked("three-scale finite differences agree with exact Frechet formula",
                fd_errors[0] < 2e-7 and fd_errors[1] < 2e-9
                and fd_errors[2] < 3e-11
                and fd_errors[1] < fd_errors[0]/80
                and fd_errors[2] < fd_errors[1]/40)

        ideal_g, ideal_dg, _, _ = state(p0, pulse_time=0, derivatives=True)
        ideal_je = np.einsum("mgh,khg->mk", effect_coefficients, ideal_dg).real
        ideal_jq = np.asarray([qgraph @ np.diag(d).real for d in ideal_dg]).T
        finite_pulse_difference = float(np.max(abs(je-ideal_je)))
        checked("finite H-on pulse agrees with the old Duhamel error scale",
                finite_pulse_difference < 324*tau
                and np.max(abs(jq-ideal_jq)) < 324*tau)

        # The active old data are stored, not destroyed. For fresh independent Q,
        # contraction of the retained source columns equals the graph effect.
        rng = np.random.default_rng(1077)
        retained_probability_residual = 0.
        original_reader_probability_residual = 0.
        for _ in range(6):
            r = rng.normal(size=3)
            r *= rng.uniform(.05, .95)/np.linalg.norm(r)
            weights = np.concatenate(([1.], r))
            graph_effect = sum(a*b for a, b in zip(weights, effect_coefficients))
            # Full old data plus retained 28-column purification label.
            stored_probability = np.einsum(
                "dga,gh,dha->", output.conj(), graph_effect, output).real/384
            graph_probability = np.trace(graph_effect @ graph).real
            retained_probability_residual = max(retained_probability_residual,
                                                 abs(stored_probability-graph_probability))
            fresh_data = sum(a*c for a, c in zip(weights, raw))/64
            direct_reader_probability = np.trace(
                reader_o @ np.kron(fresh_data, graph)).real
            effect_probability = float(weights @ effect_center)
            original_reader_probability_residual = max(
                original_reader_probability_residual,
                abs(direct_reader_probability-effect_probability))
        checked("old-data storage and independent new-probe reader contractions",
                retained_probability_residual < 2e-12
                and original_reader_probability_residual < 2e-12)

        # Diagnostics only. Strict invertibility is the separate author's result.
        q_singular = np.linalg.svd(jq, compute_uv=False)
        b_singular = np.linalg.svd(je[1:], compute_uv=False)
        checked("well-resolved nondegenerate numerical candidate",
                q_singular[-1] > 1e-4 and b_singular[-1] > .0012)

    return {
        "round": 1077, "status": "passed_diagnostic",
        "assertion_groups": len(GROUPS), "groups": GROUPS,
        "parameters": {"p0": p0.tolist(), "pulse_time": tau,
                       "reader_time": reader_t, "reader": "source 1 fixed Z+"},
        "graph_eigenvalues": graph_eigenvalues.tolist(),
        "graph_coherence_frobenius": graph_coherence,
        "matching_pair_totals": pair_total.tolist(),
        "matching_pair_derivative_max": float(np.max(abs(derivative_pair_total))),
        "branch_matching_residual": branch_residual,
        "old_Q": current_q.tolist(),
        "old_Q_Jacobian": jq.tolist(),
        "old_Q_Jacobian_reproduction_residual": float(np.max(abs(jq-old_jq))),
        "old_Q_singular_values_diagnostic": q_singular.tolist(),
        "old_Q_determinant_diagnostic": float(np.linalg.det(jq)),
        "effect_center": effect_center.tolist(),
        "effect_Jacobian": je.tolist(),
        "Bloch_singular_values_diagnostic": b_singular.tolist(),
        "Bloch_determinant_diagnostic": float(np.linalg.det(je[1:])),
        "Bloch_vs_Q_local_derivative_diagnostic": (je[1:] @ np.linalg.inv(jq)).tolist(),
        "population_only_effect_Jacobian": diagonal_je.tolist(),
        "coherence_effect_Jacobian": coherence_je.tolist(),
        "finite_difference_diagnostics": finite_differences,
        "ideal_effect_Jacobian": ideal_je.tolist(),
        "finite_pulse_effect_Jacobian_difference": finite_pulse_difference,
        "old_Duhamel_entry_bound": 324*tau,
        "retained_source_probability_residual": float(retained_probability_residual),
        "fixed_reader_probability_residual": float(original_reader_probability_residual),
        "source_hashes": {
            "source480": digest(ROOT/"research_cognition_physics/archive_467_530/480/correlated_reference_local_access.py"),
            "reader474": digest(ROOT/"research_cognition_physics/archive_467_530/474/relational_direction_distance_response.py"),
            "frozen1076proof": digest(ROOT/"research_cognition_physics/archive_1063_/1076/proof.md"),
            "independent_check": digest(Path(__file__)),
        },
        "scope": {
            "strict_rank_proof_from_float": False,
            "all_graph_coherences_retained": True,
            "coordinate_to_Pauli_feedback": False,
            "new_dephasing_device": False,
            "old_data_and_records_destroyed": False,
            "complete_location_encoded_in_one_qubit": False,
            "same_Q_arbitrary_history_reader_equivalence": False,
            "fixed_actual_source_control_context": True,
            "spatial_dimension_derived": False,
            "independent_probe_and_storage_permissions_are_inputs": True,
        },
    }


if __name__ == "__main__":
    result = run_independent()
    # Only after all independent calculations: compare against saved exact
    # source/reader contraction and its finite-pulse enclosure.
    author_results_path = Path(__file__).with_name("results.json")
    author_code_path = Path(__file__).with_name("check.py")
    author_result = json.loads(author_results_path.read_text(encoding="utf-8"))
    author_ideal = np.asarray(
        author_result["ideal_pulse_effect_jacobian_4x3_diagnostic"], float)
    bounds = np.asarray(
        author_result["finite_pulse_jacobian_entry_error_by_column_float"], float)
    independent_ideal = np.asarray(result["ideal_effect_Jacobian"])
    actual = np.asarray(result["effect_Jacobian"])
    ideal_difference = float(np.max(abs(independent_ideal-author_ideal)))
    actual_column_differences = np.max(abs(actual-author_ideal), axis=0)
    checked("final comparison to author exact ideal source and finite-pulse enclosure",
            author_result["status"] == "strict_certificate_passed"
            and ideal_difference < 2e-12
            and np.all(actual_column_differences < bounds))
    result["assertion_groups"] = len(GROUPS)
    result["author_comparison"] = {
        "performed_only_after_independent_calculation": True,
        "ideal_intermediate_max_difference": ideal_difference,
        "actual_finite_pulse_column_max_difference":
            actual_column_differences.tolist(),
        "author_finite_pulse_entry_bounds": bounds.tolist(),
        "actual_inside_author_strict_enclosure": True,
        "author_code_sha256": digest(author_code_path),
        "author_results_sha256": digest(author_results_path),
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
