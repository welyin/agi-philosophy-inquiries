"""Independent exact audit; standard library only, no author-code import."""
from pathlib import Path
from itertools import permutations
import argparse
import json


def zeros(n=3):
    return [[0 for _ in range(n)] for _ in range(n)]


def add(*matrices):
    n = len(matrices[0])
    return [[sum(a[i][j] for a in matrices) for j in range(n)]
            for i in range(n)]


def scale(c, a):
    return [[c * entry for entry in row] for row in a]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def commutator(a, b):
    return add(multiply(a, b), scale(-1, multiply(b, a)))


def determinant(a):
    n = len(a)
    total = 0
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = -1 if inversions % 2 else 1
        for i in range(n):
            term *= a[i][p[i]]
        total += term
    return total


def frobenius_squared(a):
    return sum(x * x for row in a for x in row)


def linear_combination(coefficients, matrices):
    return add(*(scale(c, a) for c, a in zip(coefficients, matrices)))


def compute():
    # Enter the displayed qutrit matrices directly, then independently use
    # Sylvester minors and a nonzero diagonal minor, not the author's RREF.
    terms = [
        [[10, 1, 0], [1, 11, 1], [0, 1, 14]],
        [[10, 1, 0], [1, 9, 1], [0, 1, 6]],
        [[10, -2, 0], [-2, 11, -2], [0, -2, 12]],
    ]
    identity = [[int(i == j) for j in range(3)] for i in range(3)]
    d = [[i if i == j else 0 for j in range(3)] for i in range(3)]
    a = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    d2 = multiply(d, d)
    assert terms == [
        add(scale(10, identity), a, d2),
        add(scale(10, identity), a, scale(-1, d2)),
        add(scale(10, identity), d, scale(-2, a)),
    ]
    assert all(h[i][j] == h[j][i]
               for h in terms for i in range(3) for j in range(3))
    leading_minors = [
        [determinant([row[:k] for row in h[:k]]) for k in range(1, 4)]
        for h in terms
    ]
    assert leading_minors == [[10, 109, 1516], [10, 89, 524], [10, 106, 1232]]
    assert all(value > 0 for group in leading_minors for value in group)
    lower_bounds = [
        min(h[i][i] - sum(abs(h[i][j]) for j in range(3) if j != i)
            for i in range(3)) for h in terms
    ]
    assert lower_bounds == [9, 5, 7]
    diagonal_column_minor = [[h[i][i] for h in terms] for i in range(3)]
    independence_minor = determinant(diagonal_column_minor)
    assert independence_minor == 40

    h_total = add(*terms)
    weights = [3, 1, 2]
    source = linear_combination(weights, terms)
    assert h_total == add(scale(30, identity), d)
    assert source == add(scale(60, identity), scale(2, d), scale(2, d2))
    pairs = [(0, 1), (0, 2), (1, 2)]
    pair_commutators = [commutator(terms[i], terms[j]) for i, j in pairs]
    flows = [commutator(h_total, term) for term in terms]
    pair_norms = [frobenius_squared(c) for c in pair_commutators]
    flow_norms = [frobenius_squared(c) for c in flows]
    assert pair_norms == [80, 116, 52]
    assert flow_norms == [4, 4, 16]
    assert all(n > 0 for n in pair_norms + flow_norms)
    # Explicit dependence is stronger than merely reporting a rank.
    assert flows[0] == flows[1]
    assert flows[2] == scale(-2, flows[0])
    assert flows[0][0][1] == -1
    assert linear_combination([-2, -1, 1], pair_commutators) == zeros()
    assert determinant([[pair_commutators[0][0][1], pair_commutators[1][0][1]],
                        [pair_commutators[0][1][2], pair_commutators[1][1][2]]]) == 4
    null_basis = [[1, 1, 1], [1, -1, 0]]
    assert all(linear_combination(v, flows) == zeros() for v in null_basis)
    assert determinant([[1, 1], [1, -1]]) != 0
    assert weights == [2 * x + y for x, y in zip(*null_basis)]
    assert linear_combination([weights[j] - weights[i] for i, j in pairs],
                              pair_commutators) == commutator(h_total, source)
    assert commutator(h_total, source) == zeros()
    source_eigenvalues = [source[i][i] for i in range(3)]
    second_difference = source_eigenvalues[2] - 2 * source_eigenvalues[1] + source_eigenvalues[0]
    assert second_difference == 4 and min(source_eigenvalues) > 0
    assert source == add(scale(2, h_total), scale(2, d2))
    # In the normalized state (|0> + i|1>)/sqrt(2),
    # <i[H,h_a]> = -[H,h_a]_(0,1). All values remain exact integers.
    witness_currents = [-c[0][1] for c in flows]
    assert witness_currents == [1, 1, -2]
    assert sum(witness_currents) == sum(g * v for g, v in zip(weights, witness_currents)) == 0

    # Compare only after independently deriving all quantities above.
    author = json.loads(Path(__file__).with_name("results.json").read_text(encoding="utf-8"))
    comparisons = {
        "hamiltonian_terms": terms,
        "strict_positive_gershgorin_bounds": lower_bounds,
        "term_span_rank": 3,
        "H": h_total,
        "source_weights": weights,
        "J": source,
        "pair_commutator_frobenius_squared": pair_norms,
        "actual_exchange_frobenius_squared": flow_norms,
        "weight_map_rank": 1,
        "weight_kernel_dimension": 2,
        "kernel_equation": "g1 + g2 - 2*g3 = 0",
        "source_commutator_squared": 0,
        "source_spectral_second_difference": second_difference,
    }
    assert all(author[key] == value for key, value in comparisons.items())
    return {
        "status": "PASS_WITH_SCOPE",
        "audit_type": "Unnumbered independent exact matrix audit; not a new scientific round",
        "implementation": "Python standard library; no numpy and no import or execution of author check.py",
        "hamiltonian_terms": terms,
        "strict_positivity_sylvester_leading_minors": leading_minors,
        "strict_positive_gershgorin_bounds": lower_bounds,
        "linear_independence_diagonal_minor": independence_minor,
        "H": h_total,
        "source_weights": weights,
        "J": source,
        "pair_commutators": pair_commutators,
        "pair_commutator_frobenius_squared": pair_norms,
        "pair_commutator_rank": 2,
        "pair_commutator_zero_combination": [-2, -1, 1],
        "actual_exchange_commutators": flows,
        "actual_exchange_frobenius_squared": flow_norms,
        "normal_qutrit_witness_current_expectations": witness_currents,
        "weight_map_rank": 1,
        "weight_kernel_dimension": 2,
        "weight_kernel_basis": null_basis,
        "kernel_equation": "g1 + g2 - 2*g3 = 0",
        "source_eigenvalues": source_eigenvalues,
        "source_spectral_second_difference": second_difference,
        "source_commutator_squared": 0,
        "nonaffine_identity": "J = 2*H + 2*(H - 30*I)^2",
        "author_result_fields_independently_matched": sorted(comparisons),
        "scope": {
            "same_qutrit_not_three_spatial_subsystems": True,
            "positive_linear_source_class_counterexample_only": True,
            "physical_monopole_not_constructed": True,
            "1047_TT_is_not_total_monopole_counterexample": True,
            "Lloyd_assessment_is_literature_admission_only": True,
            "new_axioms": 0,
            "new_numbered_round": False,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = compute()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    print(json.dumps(result, ensure_ascii=False))