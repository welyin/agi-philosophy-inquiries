"""1078 independent spectral-invariant audit. Default stdout JSON only.

First recomputes the center by the frozen 1077 independent spectral/Frechet
algorithm. Exact propagation is reused only at the final interval stage through
author original_data(); interval extrema are independently computed by grouped
multiaffine corner enumeration, not by the author's interval implementation.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import numpy as np

BASE = Path(__file__).resolve().parent
GROUPS = []


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    GROUPS.append(name)


def endpoints(center, radius):
    return center-radius, center+radius


def outward(pair, denominator=10**12):
    low, high = pair
    left = F((low*denominator).numerator//(low*denominator).denominator, denominator)
    x = high*denominator
    right = F(-((-x.numerator)//x.denominator), denominator)
    assert left <= low <= high <= right
    return [str(left), str(right)]


def grouped_corner_range(values, derivatives):
    """Exact extrema on the full independent box of the u/theta minor.

For fixed c_u,c_theta the three terms
2*b_k*(c_u*j_k_theta-c_theta*j_k_u) depend on disjoint triples.
Enumerate their eight corners separately, sum extrema, then enumerate the
four shared c corners. This is equivalent to all 2^11 full-box corners.
"""
    ranges = []
    evaluations = 0
    for cu, ct in itertools.product(derivatives[0][1], derivatives[0][2]):
        low, high = F(0), F(0)
        for k in (1, 2, 3):
            candidates = [
                2*b*(cu*jt-ct*ju)
                for b, ju, jt in itertools.product(
                    values[k], derivatives[k][1], derivatives[k][2])]
            evaluations += len(candidates)
            low += min(candidates)
            high += max(candidates)
        ranges.append((low, high))
    return (min(pair[0] for pair in ranges), max(pair[1] for pair in ranges)), evaluations


def positive_control():
    """Mature abstract legal effect family; NOT a new physical reader."""
    identity = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    pauli = (x, y, z)
    strength = F(1, 5)
    e0 = identity/2+float(strength)*z
    # Tangent sphere directions at north pole give two independent effects.
    tangent = np.column_stack((np.array([float(strength), 0, 0]),
                               np.array([0, float(strength), 0])))
    ux = (identity-1j*x)/np.sqrt(2)
    uy = (identity-1j*y)/np.sqrt(2)
    qab = ux @ uy @ e0 @ uy.conj().T @ ux.conj().T
    qba = uy @ ux @ e0 @ ux.conj().T @ uy.conj().T
    plus_x = (identity+x)/2
    gap = float(abs(np.trace(plus_x @ (qab-qba)).real))
    spectra = [np.linalg.eigvalsh(a) for a in (e0, qab, qba)]
    checked("abstract radial fixed-trace control admits angular rank two",
            np.linalg.matrix_rank(tangent) == 2
            and max(np.max(abs(v-np.array([.3, .7]))) for v in spectra) < 2e-14
            and abs(gap-.2) < 2e-14)
    return {
        "effect_family": "I/2 + (1/5) n.sigma, |n|=1",
        "fixed_spectral_invariants": ["1/2", "1/25"],
        "tangent_rank_diagnostic": 2,
        "fixed_plus_X_order_gap": gap,
        "physical_source_or_spatial_dimension_claimed": False,
    }


def main():
    previous = load("frozen1077_independent_center", BASE.parent/"1077/independent_check.py")
    old = previous.run_independent()
    center = np.asarray(old["effect_center"], float)
    jacobian = np.asarray(old["effect_Jacobian"], float)
    invariant_jacobian = np.vstack((jacobian[0], 2*center[1:] @ jacobian[1:]))
    minors = []
    for cols in itertools.combinations(range(3), 2):
        minors.append({
            "columns": [("t", "u", "theta")[i] for i in cols],
            "value": float(np.linalg.det(invariant_jacobian[:, cols]))})
    checked("independently recomputed center has well-resolved spectral minor",
            abs(minors[2]["value"]+1.4559584855481872e-5) < 2e-15
            and all(abs(item["value"]) > 8e-6 for item in minors))
    checked("reused propagation retained source coherence and fixed-reader context",
            old["graph_coherence_frobenius"] > .38
            and old["scope"]["all_graph_coherences_retained"]
            and not old["scope"]["coordinate_to_Pauli_feedback"])
    control = positive_control()

    # Only now reuse the exact propagation data. Do not call author's run()
    # or any of the author's interval routines.
    author = load("round1078_exact_data_only", BASE/"check.py")
    exact_values, exact_jacobian, source_error, reader_error = author.original_data()
    tau = F(1, 10**10)
    radius = F(1, 246240000)
    weights = (9, 9, 1)
    value_error = (reader_error
                   +(1+reader_error)*source_error*(2+source_error)+18*tau)
    derivative_error = [
        2*m*reader_error+(1+reader_error)*2*m*source_error*(2+source_error)+36*m*tau
        for m in weights]
    value_boxes = [endpoints(c, value_error+38*radius) for c in exact_values]
    derivative_boxes = [[endpoints(c, derivative_error[a]+76*weights[a]*radius)
                         for a, c in enumerate(row)] for row in exact_jacobian]
    exact_range, corner_evaluations = grouped_corner_range(value_boxes, derivative_boxes)
    checked("independent exact corner enclosure is strictly negative on the full old cube",
            F(-1, 50000) < exact_range[0] <= exact_range[1] < F(-1, 100000))
    # Confirm the fresh numerical center is within the center enclosure,
    # separately from the uniform-cube enclosure.
    numerical_allowance = 2e-14
    checked("independent finite-pulse center lies in the reused exact error boxes",
            max(abs(center[k]-float(exact_values[k]))-float(value_error)
                for k in range(4)) < numerical_allowance
            and max(abs(jacobian[k, a]-float(exact_jacobian[k][a]))
                    -float(derivative_error[a])
                    for k in range(4) for a in range(3)) < numerical_allowance)

    author_result = json.loads((BASE/"results.json").read_text(encoding="utf-8"))
    author_outer = tuple(F(v) for v in author_result["uniform_minor_enclosure"])
    checked("independent corner result fits author's outward interval",
            author_outer[0] <= exact_range[0] <= exact_range[1] <= author_outer[1])
    checked("old local Bloch embedding budget is reused exactly",
            F(1, 2000)-2052*radius == F(59, 120000))

    output = {
        "round": 1078, "status": "passed",
        "assertion_groups": len(GROUPS), "groups": GROUPS,
        "center_effect": center.tolist(),
        "center_effect_jacobian": jacobian.tolist(),
        "center_spectral_invariant_jacobian": invariant_jacobian.tolist(),
        "center_spectral_singular_values_diagnostic":
            np.linalg.svd(invariant_jacobian, compute_uv=False).tolist(),
        "center_minors_diagnostic": minors,
        "chosen_minor_columns": ["u", "theta"],
        "domain": "unchanged entire closed 480/1077 control cube",
        "radius": str(radius),
        "independent_uniform_minor_enclosure": outward(exact_range),
        "independent_uniform_minor_float": [float(v) for v in exact_range],
        "author_outward_enclosure": [str(v) for v in author_outer],
        "independent_interval_algorithm":
            "exact grouped multiaffine corners; four shared derivative corners and three independent eight-corner terms",
        "corner_term_evaluations": corner_evaluations,
        "equivalent_full_box_corner_count": 2**11,
        "strict_minor_upper": "-1/100000",
        "center_error_budgets": {
            "effect_float": float(value_error),
            "jacobian_by_column_float": [float(v) for v in derivative_error]},
        "abstract_positive_control": control,
        "source_reuse": {
            "independent_spectral_Frechet_source": "1077/independent_check.py",
            "exact_propagation_source_reused_at_final_stage":
                "1078.original_data -> frozen 1077 Gaussian source and 1076 reader",
            "author_interval_algorithm_reused": False,
            "new_propagation_model_or_parameters": False,
        },
        "hashes": {
            "independent_check": sha(Path(__file__)),
            "frozen1077_independent_check": sha(BASE.parent/"1077/independent_check.py"),
            "author_check": sha(BASE/"check.py"),
            "author_results": sha(BASE/"results.json"),
            "frozen1077_proof": sha(BASE.parent/"1077/proof.md"),
        },
        "scope": {
            "fixed_context_spectral_level_is_locally_one_dimensional": True,
            "local_angular_rank_two_conjugacy_patch_in_this_context": False,
            "argument_relies_only_on_small_cube_not_fitting_a_full_sphere": False,
            "later_departure_from_cube_rescues_a_local_patch_through_it": False,
            "all_readers_or_actual_context_transport_excluded": False,
            "spatial_dimension_derived_or_refuted": False,
            "numerical_rank_is_universal_proof": False,
            "positive_control_is_an_abstract_legality_check": True,
        },
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
