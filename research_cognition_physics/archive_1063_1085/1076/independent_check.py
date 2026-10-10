"""Independent round 1076 check; default stdout JSON, no evidence writes.

Uses the archived 474/442 model through the existing migration runtime.
The author's exact_reader is imported only after every independent calculation
and certificate, solely to compare values after explicit graph-order mapping.
No image, control synthesis, external data, or new dependency is used.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import math
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from research_layout import Layout, ResearchRuntime

GROUPS = []
TOL = 2e-11


def checked(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    GROUPS.append(name)


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tensor(matrices):
    answer = np.ones((1, 1), complex)
    for a in matrices:
        answer = np.kron(answer, a)
    return answer


def maxrow(a):
    return float(np.max(np.sum(np.abs(a), axis=1)))


def rational_certificate(a, entry_error=F(1, 100000)):
    """An exact rational Neumann test, NOT a numerical rank threshold."""
    n = len(a)
    center = [[F(round(float(x) * 10**12), 10**12) for x in row] for row in a]
    approximate_inverse = [
        [F(round(float(x) * 10**6), 10**6) for x in row]
        for row in np.linalg.inv(a)]
    residual = [[
        sum(approximate_inverse[i][k] * center[k][j] for k in range(n))
        - int(i == j) for j in range(n)] for i in range(n)]
    residue = max(sum(abs(x) for x in row) for row in residual)
    inverse_norm = max(sum(abs(x) for x in row) for row in approximate_inverse)
    upper = residue + inverse_norm * n * entry_error
    if upper >= 1:
        raise AssertionError("rational inverse does not certify nonsingularity")
    return {
        "center": [[str(x) for x in row] for row in center],
        "approximate_inverse": [[str(x) for x in row] for row in approximate_inverse],
        "entry_error": str(entry_error),
        "exact_residual": str(residue),
        "inverse_infinity_norm": str(inverse_norm),
        "neumann_bound": str(upper),
        "neumann_bound_float": float(upper),
        "strictly_below_one": True,
    }


def det3(a):
    return (a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
            - a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
            + a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]))


def main():
    with ResearchRuntime(Layout(), "late").installed():
        import relational_direction_distance_response as archived
        trees, h, _, _, graph_q, eigenvalues, eigenvectors = archived.system()
        h = h.astype(np.int64)
        identity = np.eye(2, dtype=complex)
        pauli = tuple(np.asarray(p, complex) for p in archived.P)
        x, y, z = pauli
        raw = np.asarray([
            tensor([identity, source, identity, identity+x, identity+y, identity+z])
            for source in (identity, x, y, z)])
        m = np.kron(tensor([identity, (identity+z)/2,
                            identity, identity, identity, identity]), np.eye(6))
        choices = [
            tuple(j for j in (0, 3, 4, 5)
                  if archived.old.edge(1, j) in tree)
            for tree in trees]
        canonical = list(itertools.combinations((0, 3, 4, 5), 2))
        canonical_to_old = [choices.index(c) for c in canonical]
        plus = [choices.index((0, j)) for j in (3, 4, 5)]
        minus = [choices.index(tuple(k for k in (0, 3, 4, 5) if k not in (0, j)))
                 for j in (3, 4, 5)]
        checked("archived Hamiltonian and graph-coordinate convention",
                np.array_equal(h, h.T)
                and np.all(h.sum(axis=1) == 9)
                and canonical_to_old == [5, 4, 2, 3, 1, 0]
                and plus == [5, 4, 2] and minus == [0, 1, 3])

        def coefficients(o):
            # O indexed (data,input graph,data,output graph).
            # Pauli-source coefficient is raw_mu/64, including source sigma_mu/2.
            return np.einsum("diej,ked->kij", o.reshape(64, 6, 64, 6), raw)/64

        def diagonal_response(b):
            return ((b[:, plus, plus]-b[:, minus, minus])/2).real

        def evolve(t):
            u = (eigenvectors*np.exp(-1j*t*eigenvalues)) @ eigenvectors.conj().T
            return u.conj().T @ m @ u

        t = 0.7
        o = evolve(t)
        b = coefficients(o)
        n = sum(np.kron(s, c) for s, c in zip((identity, x, y, z), b))
        spectrum = np.linalg.eigvalsh(n)
        checked("same-process qubit-graph effect is Hermitian and physical",
                np.max(np.abs(b-b.conj().transpose(0, 2, 1))) < TOL
                and spectrum[0] >= -TOL and spectrum[-1] <= 1+TOL)

        # Direct probabilities on arbitrary independently injected qubit and graph.
        rng = np.random.default_rng(1076)
        direct_probability_residual = 0.
        for _ in range(8):
            v = rng.normal(size=3)
            v *= rng.uniform(0, .9)/np.linalg.norm(v)
            source = (identity+sum(a*p for a, p in zip(v, pauli)))/2
            g = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
            sigma = g @ g.conj().T
            sigma /= np.trace(sigma)
            data = (raw[0]+sum(a*c for a, c in zip(v, raw[1:])))/64
            direct = np.trace(o @ np.kron(data, sigma))
            factored = np.trace(n @ np.kron(source, sigma))
            direct_probability_residual = max(direct_probability_residual,
                                             abs(direct-factored))
        checked("independent-input probability identity",
                direct_probability_residual < TOL)

        # Independent Taylor algorithm: full 384x384 Heisenberg recursion,
        # distinct from the author's 384x48 exact-integer Schrödinger columns.
        degree = 80
        term = m.astype(complex)
        polynomial = term.copy()
        gamma = 1e-12
        length = 12.6
        current_error = accumulated_error = addition_error = 0.
        for k in range(1, degree+1):
            previous_norm = maxrow(term)
            total_norm = maxrow(polynomial)
            current_error = (length/k)*current_error + gamma*(length/k)*previous_norm
            term = 1j*(t/k)*(h @ term-term @ h)
            accumulated_error += current_error
            addition_error += gamma*(total_norm+maxrow(term))
            polynomial += term
        # Conservative factor 2 covers floating evaluation of the error ledger,
        # including row sums. Explicit rational tail uses exp(L)<3^13.
        tail = F(3**13) * F(63, 5)**81 / math.factorial(81)
        doubled_error = 2*(accumulated_error+addition_error)+float(tail)
        series_b = coefficients(polynomial)
        series_response = diagonal_response(series_b)
        raw_entry_l1 = [float(np.sum(np.abs(a))/64) for a in raw]
        # Each contraction is bounded by 4 times max-row effect error.
        # 1e-9 covers coefficient contraction/decimal-centering arithmetic;
        # unused margin is > three times the full bound.
        entry_error_budget = 4*doubled_error + 1e-9
        checked("independent series enclosure fits chosen coefficient budget",
                entry_error_budget < 1e-5
                and maxrow(polynomial-o) < doubled_error
                and raw_entry_l1 == [4., 4., 4., 4.])
        rank3 = rational_certificate(series_response[1:])
        checked("rank three certified by rational approximate inverse",
                rank3["strictly_below_one"])

        # Canonical author's W01 is old W54; never silently reuse graph labels.
        a, c = canonical_to_old[0], canonical_to_old[1]
        w = np.zeros((6, 6), complex)
        w[a, c] = w[c, a] = .5
        canonical_coh = series_b[:, a, c].real
        combined = np.column_stack((series_response, canonical_coh))
        combined_certificate = rational_certificate(combined)
        checked("three diagonal directions plus one coherence certify rank four",
                combined_certificate["strictly_below_one"])

        # All-zero population coherence subspace alone spans four coefficients.
        labels, columns, perturbations = [], [], []
        for g in range(6):
            for j in range(g+1, 6):
                wr = np.zeros((6, 6), complex)
                wi = wr.copy()
                wr[g, j] = wr[j, g] = .5
                wi[g, j], wi[j, g] = -.5j, .5j
                for label, perturbation in (("Re", wr), ("Im", wi)):
                    labels.append((g, j, label))
                    perturbations.append(perturbation)
                    columns.append(np.array([np.trace(bb @ perturbation).real
                                             for bb in series_b]))
        coherence_matrix = np.asarray(columns).T
        chosen_labels = [(0, 1, "Im"), (0, 3, "Re"),
                         (1, 3, "Im"), (2, 4, "Re")]
        selected = [labels.index(item) for item in chosen_labels]
        selected_matrix = coherence_matrix[:, selected]
        coherence_certificate = rational_certificate(selected_matrix)
        checked("zero-population coherences alone certify all four coefficients",
                coherence_certificate["strictly_below_one"])

        # Pure trace direction exists exactly because the preceding matrix is
        # invertible. The following floating combination is an illustration,
        # not a claim of exact preparation from decimal coefficients.
        weights = np.linalg.solve(selected_matrix, np.array([1., 0, 0, 0]))
        weights /= max(abs(weights))
        pure_trace_w = sum(value*perturbations[index]
                           for value, index in zip(weights, selected))
        pure_trace_response = selected_matrix @ weights
        trace_pair_gap = pure_trace_response[0]/6
        checked("trace-only coherence illustration and positive graph states",
                np.linalg.norm(pure_trace_response[1:]) < TOL
                and pure_trace_response[0] > .03
                and np.allclose(np.diag(pure_trace_w), 0)
                and np.linalg.eigvalsh(np.eye(6)/6+pure_trace_w/12)[0] > 0
                and np.linalg.eigvalsh(np.eye(6)/6-pure_trace_w/12)[0] > 0)

        sigma_plus, sigma_minus = np.eye(6)/6+.1*w, np.eye(6)/6-.1*w
        old_q_difference = graph_q @ (
            np.diag(sigma_plus).real-np.diag(sigma_minus).real)
        effect_difference = np.array(
            [np.trace(bb @ (sigma_plus-sigma_minus)).real for bb in series_b])
        probability_gap_bound = .2*(abs(canonical_coh[3])-1e-5)
        fixed_plus_z_gap_bound = .2*(canonical_coh[0]+canonical_coh[3]-2e-5)
        checked("same original position means give distinguishable fixed-reader effects",
                np.linalg.norm(old_q_difference) == 0
                and np.linalg.eigvalsh(sigma_plus)[0] > 0
                and np.linalg.eigvalsh(sigma_minus)[0] > 0
                and probability_gap_bound > .0048 and fixed_plus_z_gap_bound > .0074)

        # Finite exact-symbolic coefficient audit: at these orders Gaussian
        # integers remain exactly representable, including each contraction.
        nested = m.astype(complex)
        taylor = []
        maximum_nested_entry = 0.
        for order in range(9):
            checked_integer = (np.array_equal(nested.real, np.rint(nested.real))
                               and np.array_equal(nested.imag, np.rint(nested.imag)))
            if not checked_integer:
                raise AssertionError("nested commutator lost integer exactness")
            maxentry = float(np.max(np.abs(nested.real))
                             + np.max(np.abs(nested.imag)))
            maximum_nested_entry = max(maximum_nested_entry, maxentry)
            # H row sum 9; raw contraction has <=4096 entries of magnitude <=2.
            if 2*9*maxentry >= 2**53 or 8192*maxentry >= 2**53:
                raise AssertionError("integer representability certificate failed")
            derivative = diagonal_response(coefficients(nested))[1:]
            if any(F(float(v)).denominator > 1024 for v in derivative.flat):
                raise AssertionError("dyadic exact coefficient denominator exceeded")
            taylor.append([[F(float(v)).limit_denominator(1024)/math.factorial(order)
                            for v in row] for row in derivative])
            nested = 1j*(h @ nested-nested @ h)
        determinant_coefficients = []
        for order in range(9):
            value = F(0)
            for i in range(order+1):
                for j in range(order-i+1):
                    k = order-i-j
                    value += det3([taylor[i][0], taylor[j][1], taylor[k][2]])
            determinant_coefficients.append(value)
        checked("nonzero eighth-order determinant coefficient",
                determinant_coefficients[:8] == [F(0)]*8
                and determinant_coefficients[8] == F(-5, 1536))

        diagnostics = []
        for ti in (.001, .01, .1, .3, .7, 1., 2.):
            bt = coefficients(evolve(ti))
            ar = diagonal_response(bt)
            diagnostics.append({
                "t": ti,
                "trace_response": ar[0].tolist(),
                "bloch_singular_values_diagnostic":
                    np.linalg.svd(ar[1:], compute_uv=False).tolist(),
            })

    # Last step only: compare with the independently authored exact algorithm.
    author_path = Path(__file__).with_name("check.py")
    spec = importlib.util.spec_from_file_location("round1076_author_comparison", author_path)
    author = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(author)
    _, _, author_h, author_response, author_coh, _, author_error, _ = author.exact_reader()
    index = [6*d+g for d in range(64) for g in canonical_to_old]
    author_response_residual = float(np.max(np.abs(
        np.asarray(author_response, float)-series_response)))
    author_coherence_residual = float(np.max(np.abs(
        np.asarray(author_coh, float)-canonical_coh)))
    checked("final comparison after explicit graph-order permutation",
            np.array_equal(h[np.ix_(index, index)], author_h)
            and author_response_residual < TOL
            and author_coherence_residual < TOL)

    output = {
        "round": 1076, "status": "passed",
        "assertion_groups": len(GROUPS), "groups": GROUPS,
        "source": "archived 474/442 model; independent Heisenberg contraction",
        "graph_order": [list(c) for c in choices],
        "canonical_to_archived_permutation": canonical_to_old,
        "effect_convention": "N=sum sigma_mu tensor B_mu; E_sigma=sum Tr(B_mu sigma) sigma_mu",
        "t": "7/10", "Hilbert_dimension": 384,
        "effect_extreme_eigenvalues": [float(spectrum[0]), float(spectrum[-1])],
        "direct_probability_residual": float(direct_probability_residual),
        "response_4x3": series_response.tolist(),
        "canonical_W01_response": canonical_coh.tolist(),
        "old_W01_response": series_b[:, 0, 1].real.tolist(),
        "series": {
            "degree": degree, "H_infinity_norm": 9, "gamma": gamma,
            "arithmetic_assumptions": "IEEE binary64; no overflow/underflow; 384-term real-integer times complex dot products; gamma deliberately exceeds their operation-count bound",
            "recurrence_error": accumulated_error,
            "addition_error": addition_error,
            "double_margin_error": doubled_error,
            "rational_tail_bound": str(tail),
            "coefficient_error_budget": entry_error_budget,
            "certificate_entry_radius": "1/100000",
            "series_vs_spectral_maxrow": maxrow(polynomial-o),
        },
        "rank3_certificate": rank3,
        "diagonal_plus_coherence_rank4_certificate": combined_certificate,
        "coherence_only_rank4_certificate": coherence_certificate,
        "coherence_only_columns": chosen_labels,
        "coherence_only_response": selected_matrix.tolist(),
        "trace_only_weights_illustrative": weights.tolist(),
        "trace_only_response_illustrative": pure_trace_response.tolist(),
        "trace_only_pair_probability_gap_illustrative": float(trace_pair_gap),
        "same_population_effect_difference": effect_difference.tolist(),
        "same_population_probability_gap_lower": probability_gap_bound,
        "same_population_fixed_plus_Z_gap_lower": fixed_plus_z_gap_bound,
        "determinant_taylor_coefficients_0_to_8":
            [str(v) for v in determinant_coefficients],
        "integer_nested_maximum_entry": maximum_nested_entry,
        "time_diagnostics": diagnostics,
        "author_comparison": {
            "used_only_after_independent_calculations": True,
            "response_residual": author_response_residual,
            "coherence_residual": author_coherence_residual,
            "author_strict_effect_error": float(author_error),
            "author_code_sha256": hash_file(author_path),
        },
        "source_hashes": {
            "archived474": hash_file(ROOT/"research_cognition_physics/archive_467_530/474/relational_direction_distance_response.py"),
            "archived442": hash_file(ROOT/"research_cognition_physics/archive_429_466/442/branching_tree_distance_audit.py"),
            "independent_check": hash_file(Path(__file__)),
        },
        "scope": {
            "fixed_readout_bloch_response_on_declared_diagonal_section": 3,
            "effect_response_with_same_population_coherences": 4,
            "descends_to_entire_original_distance_mean_quotient": False,
            "full_unknown_location_encoded_as_qubit": False,
            "actual_spatial_dimension_derived": False,
            "new_control_device_built": False,
            "pure_trace_exact_existence_from_invertibility": True,
            "decimal_weights_exact_preparation_claimed": False,
            "same_population_coherence_obstruction_previously_known": True,
        },
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
