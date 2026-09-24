"""Round 415: compact bounded-control budgets need not be continuous costs.

The infinite-dimensional compactness theorem is analytic. These finite checks
audit Dyson tails, control limits, and an exact injective recurrent example.
"""
import argparse
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

TARGET = Path(__file__).with_name("bounded_control_budget_audit_results.json")
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
GENERATORS = (X, Z)


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def endpoint(protocol):
    out = np.eye(2, dtype=complex)
    for duration, control in protocol:
        h = sum(x * g for x, g in zip(control, GENERATORS))
        out = unitary(h, duration) @ out
    return out


def graded_dyson(protocol, order):
    """Grade by total generator degree, not separately truncated pulse order."""
    grades = [np.eye(2, dtype=complex)] + [np.zeros((2, 2), complex) for _ in range(order)]
    for duration, control in protocol:
        a = -1j * duration * sum(x * g for x, g in zip(control, GENERATORS))
        powers = [np.eye(2, dtype=complex)]
        for k in range(1, order + 1):
            powers.append(powers[-1] @ a / k)
        grades = [sum((powers[j] @ grades[k-j] for j in range(k+1)), np.zeros((2, 2), complex))
                  for k in range(order+1)]
    return sum(grades), grades


def tail_bound(a, order):
    # Positive terms, avoiding cancellation in exp(a) - truncated sum.
    assert order + 2 > a >= 0
    return a ** (order + 1) / math.factorial(order + 1) / (1 - a / (order + 2))


def word_coefficients(protocol, order):
    """Real iterated-integral coefficients; later operator words prepend."""
    coeff = {(): 1.}
    for duration, control in protocol:
        segment = {(): 1.}
        for k in range(1, order+1):
            for word in product(range(2), repeat=k):
                segment[word] = duration**k / math.factorial(k) * math.prod(control[j] for j in word)
        new = {}
        for a, x in segment.items():
            for b, y in coeff.items():
                if len(a) + len(b) <= order:
                    new[a+b] = new.get(a+b, 0.) + x*y
        coeff = new
    return coeff


def word_matrix(word):
    out = np.eye(2, dtype=complex)
    for j in word:
        out = out @ GENERATORS[j]
    return (-1j)**len(word) * out


def from_words(coeff):
    return sum((x * word_matrix(word) for word, x in coeff.items()), np.zeros((2, 2), complex))


def trace_distance(a, b):
    return float(np.sum(abs(np.linalg.eigvalsh(a-b))) / 2)


def projector(v):
    return np.outer(v, v.conj())


def choi(u):
    d = len(u)
    return projector(u.reshape(-1) / np.sqrt(d))


def pell_pairs(count):
    p, q = 1, 1
    out = []
    for _ in range(count):
        out.append((p, q))
        p, q = 3*p + 4*q, 2*p + 3*q
    return out


def recurrent_endpoint(p, q, fraction=1.):
    # sqrt(2)*q - p = 1/(sqrt(2)*q+p), exactly for p^2-2q^2=-1.
    delta = 2 * np.pi / (p + np.sqrt(2)*q)
    if fraction == 1.:
        return np.diag([1., 1., 1., np.exp(-1j*delta)]), delta
    assert fraction == .5 and p % 2 == q % 2 == 1
    return np.diag([1., -1., 1., -np.exp(-.5j*delta)]), delta


@lru_cache(None)
def report():
    rng = np.random.default_rng(415)
    radius, order = 1.7, 8
    max_error = max_unitarity = max_grade_ratio = 0.
    for _ in range(18):
        durations = rng.random(6)
        durations *= radius / durations.sum()
        controls = rng.normal(size=(6, 2))
        controls /= np.sum(abs(controls), axis=1, keepdims=True)
        protocol = list(zip(durations, controls))
        exact = endpoint(protocol)
        truncated, grades = graded_dyson(protocol, order)
        max_error = max(max_error, opnorm(exact-truncated))
        max_unitarity = max(max_unitarity, opnorm(exact.conj().T@exact-np.eye(2)))
        max_grade_ratio = max(max_grade_ratio, max(opnorm(grades[k])/(radius**k/math.factorial(k))
                                                    for k in range(1, order+1)))

    protocol = [(.3, (.7, .3)), (.25, (-.4, .6)), (.35, (.2, -.8))]
    word_order, mesh = 6, 1e-4
    coeff = word_coefficients(protocol, word_order)
    reconstructed = from_words(coeff)
    matrix_polynomial, _ = graded_dyson(protocol, word_order)
    rounded = {w: (v if not w else mesh*round(v/mesh)) for w, v in coeff.items()}
    quantized = from_words(rounded)
    word_count = len(coeff)-1
    radius_words = sum(t for t, _ in protocol)
    degree_bounds = [sum(abs(v) for w, v in coeff.items() if len(w)==k) <=
                     radius_words**k/math.factorial(k)+1e-14 for k in range(1, word_order+1)]
    inverse_protocol = [(t, tuple(-v for v in u)) for t, u in reversed(protocol)]
    inverse_error = opnorm(endpoint(inverse_protocol)@endpoint(protocol)-np.eye(2))

    horizon = 1.3
    mixed = unitary((X+Z)/2, horizon)
    switch_rows = []
    comm_norm = opnorm(X@Z-Z@X)
    for n in (4, 16, 64):
        switched = np.linalg.matrix_power(unitary(Z, horizon/(2*n))@unitary(X, horizon/(2*n)), n)
        switch_rows.append(dict(cycles=n, error=opnorm(switched-mixed),
                                analytic_bound=horizon*horizon*comm_norm/(8*n),
                                integrated_l1_cost=horizon))

    rows = []
    choi_formula_error = witness_formula_error = direct_formula_error = 0.
    direct_count = 0
    full_return_probe = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    half_return_probe = np.array([1, 1, 0, 0], complex)/np.sqrt(2)
    for p, q in pell_pairs(8):
        u, delta = recurrent_endpoint(p, q)
        half, _ = recurrent_endpoint(p, q, .5)
        measured = trace_distance(choi(u), choi(np.eye(4)))
        choi_expected = np.sqrt(3)/2 * abs(np.sin(delta/2))
        witness = trace_distance(projector(u@full_return_probe),projector(full_return_probe))
        diamond_half = abs(np.sin(delta/2))
        half_witness = trace_distance(projector(half@half_return_probe),projector(half_return_probe))
        choi_formula_error = max(choi_formula_error, abs(measured-choi_expected))
        witness_formula_error = max(witness_formula_error, abs(witness-diamond_half))
        if q <= 1000:
            direct = np.diag(np.exp(-1j*2*np.pi*q*np.array([0.,1.,0.,np.sqrt(2)])))
            direct_formula_error = max(direct_formula_error, opnorm(direct-u))
            direct_count += 1
        rows.append(dict(p=p, q=q, pell_residual=p*p-2*q*q,
                         exact_minimum_cost=float(2*np.pi*q), phase_remainder=float(delta),
                         full_process_half_diamond=float(diamond_half),
                         choi_trace_distance=measured, half_scaled_witness_distance=half_witness))
    parameters = [(.2,-.4),(-.3,.7),(.8,.1)]
    h1 = np.diag([0.,1.,0.,np.sqrt(2)])
    h2 = np.diag([0.,0.,1.,np.sqrt(3)])
    composition_error = 0.
    for s,t in parameters:
        for a,b in parameters:
            uu = unitary(h1,s)@unitary(h2,t)
            vv = unitary(h1,a)@unitary(h2,b)
            together = unitary(h1,s+a)@unitary(h2,t+b)
            composition_error = max(composition_error,opnorm(uu@vv-together))
    return dict(round=415,scientific_baseline=414,
        scope="Finite bounded control generators give compact budget sublevels under an explicit measurable convex-control contract; compactness does not make exact endpoint cost or scaling continuous in process distinguishability",
        dyson=dict(protocol_count=18,horizon=radius,order=order,maximum_error=max_error,
                   analytic_tail_bound=tail_bound(radius,order),maximum_unitarity_error=max_unitarity,
                   largest_grade_to_majorant_ratio=max_grade_ratio),
        finite_word_cover=dict(order=word_order,nonempty_word_count=word_count,mesh=mesh,
             coefficient_degree_bounds_passed=all(degree_bounds),
             independent_polynomial_error=opnorm(reconstructed-matrix_polynomial),
             quantized_endpoint_error=opnorm(endpoint(protocol)-quantized),
             analytic_cover_error=tail_bound(radius_words,word_order)+word_count*mesh/2),
        inverse_protocol=dict(endpoint_error=inverse_error,cost_forward=radius_words,
                              cost_inverse=sum(t for t,_ in inverse_protocol)),
        measurable_control_limit=switch_rows,
        recurrent_injective_example=dict(generator_count=2,hilbert_dimension=4,
             parameter_cost="abs(s)+abs(t)",channel_endpoint_injectivity_proved_analytically=True,
             sequence=rows,choi_formula_error=choi_formula_error,witness_formula_error=witness_formula_error,
             independent_direct_exponential_count=direct_count,independent_direct_exponential_error=direct_formula_error,
             composition_error=composition_error),
        finite_bounded_generators_are_extra_input=True,measurable_convex_control_access_is_extra_input=True,
        compactness_requires_finite_hilbert_dimension=False,
        optimal_cost_attained_in_declared_measurable_control_class=True,
        finite_program_exact_optimality_claimed=False,
        endpoint_resource_cost_inferred_continuous=False,
        scaling_continuity_inferred_from_amplitude_adjustment=False,
        single_device_autonomous_implementation_proved=False,
        budget_identified_with_total_internal_resource_cost=False,
        actual_position_or_dimension_generated=False,new_cognitive_axiom_adopted=False,
        phase_closure_triggered=False)


class Audit(unittest.TestCase):
    def test_01_noncommuting_dyson_tail(self):
        r=report()["dyson"]
        self.assertLess(r["maximum_error"],r["analytic_tail_bound"])
        self.assertLess(r["maximum_unitarity_error"],1e-13)
        self.assertLessEqual(r["largest_grade_to_majorant_ratio"],1+1e-12)

    def test_02_finite_word_representation_and_cover(self):
        r=report()["finite_word_cover"]
        self.assertEqual(r["nonempty_word_count"],126)
        self.assertTrue(r["coefficient_degree_bounds_passed"])
        self.assertLess(r["independent_polynomial_error"],1e-13)
        self.assertLess(r["quantized_endpoint_error"],r["analytic_cover_error"])

    def test_03_actual_reverse_control(self):
        r=report()["inverse_protocol"]
        self.assertLess(r["endpoint_error"],1e-13)
        self.assertAlmostEqual(r["cost_forward"],r["cost_inverse"])

    def test_04_closed_control_class_limit(self):
        rows=report()["measurable_control_limit"]
        self.assertTrue(all(r["error"]<=r["analytic_bound"] for r in rows))
        self.assertTrue(all(a["error"]>b["error"] for a,b in zip(rows,rows[1:])))

    def test_05_exact_integer_pell_certificate(self):
        rows=report()["recurrent_injective_example"]["sequence"]
        self.assertTrue(all(r["pell_residual"]==-1 and r["p"]%2==r["q"]%2==1 for r in rows))
        self.assertTrue(all(a["exact_minimum_cost"]<b["exact_minimum_cost"] for a,b in zip(rows,rows[1:])))

    def test_06_independent_choi_and_optimal_witness(self):
        r=report()["recurrent_injective_example"]
        self.assertLess(r["choi_formula_error"],1e-13)
        self.assertLess(r["witness_formula_error"],1e-13)
        self.assertLess(r["independent_direct_exponential_error"],1e-10)

    def test_07_injective_but_discontinuous_half_scaling(self):
        rows=report()["recurrent_injective_example"]["sequence"]
        self.assertTrue(all(abs(r["half_scaled_witness_distance"]-1)<1e-13 for r in rows))
        self.assertTrue(all(a["full_process_half_diamond"]>b["full_process_half_diamond"] for a,b in zip(rows,rows[1:])))
        self.assertLess(rows[-1]["full_process_half_diamond"],1e-5)
        self.assertGreater(rows[-1]["exact_minimum_cost"],1e6)

    def test_08_independent_endpoint_composition(self):
        self.assertLess(report()["recurrent_injective_example"]["composition_error"],1e-13)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    tests=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result=dict(report(),checks=dict(run=tests.testsRun,failures=len(tests.failures),errors=len(tests.errors)),
                runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
