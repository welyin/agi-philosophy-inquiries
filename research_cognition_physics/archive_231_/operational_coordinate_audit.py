"""Round 412: budget-calibrated coordinates under the explicit round-385 contract.

The reconstruction receives only endpoint tokens and measured pair budgets.
Hidden vectors occur solely in the model oracle and independent truth checks.
No claim of generating the endpoint contract from cognitive axioms is made.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name("operational_coordinate_audit_results.json")


def rotation(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    k = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    return np.eye(3) + np.sin(angle)*k + (1-np.cos(angle))*(k@k)


class EndpointModel:
    """An explicit input model, not a discovered physical space."""
    def __init__(self):
        self.linear = np.array([[1.2, .3, -.2], [0., .8, .25], [0., 0., 1.5]])
        self.beta = 1.7
        self.twist_axis = np.array([1., 2., -1.])
        self.twist_rate = .83
        self._vectors = []
        self.calls = 0
        self.zero = self.store(np.zeros(3))

    def store(self, v):
        self._vectors.append(np.asarray(v, dtype=float).copy())
        return len(self._vectors)-1

    def from_truth(self, x):
        return self.store(np.linalg.solve(self.linear, np.asarray(x, dtype=float)))

    def truth(self, token):
        return self.linear@self._vectors[token]

    def combine(self, a, b):
        return self.store(self._vectors[a]+self._vectors[b])

    def inverse(self, a):
        return self.store(-self._vectors[a])

    def redirect(self, a, r):
        return self.from_truth(r@self.truth(a))

    def dilate(self, a, r):
        rot = rotation(self.twist_axis, self.twist_rate*np.log(r))
        return self.from_truth(r**(1/self.beta)*rot@self.truth(a))

    def budget(self, a, b):
        self.calls += 1
        return float(np.linalg.norm(self.linear@(self._vectors[b]-self._vectors[a]))**self.beta)

    def direction(self, a):
        x = self.truth(a)
        return x/np.linalg.norm(x)


def calibrate_exponent(budget, origin, unit, doubled):
    if not np.isclose(budget(origin, unit), 1., atol=1e-12, rtol=0):
        raise ValueError("calibration anchor must have unit budget")
    return float(np.log(budget(origin, doubled))/np.log(2.))


class CoordinateChart:
    """Only a pair-budget callable, tokens, and an observed exponent enter."""
    def __init__(self, budget, beta, origin, anchors):
        if beta <= 0 or len(anchors) != 3:
            raise ValueError("positive calibrated exponent and three anchors required")
        self.budget, self.beta = budget, beta
        self.origin, self.anchors = origin, tuple(anchors)
        self.radii2 = np.array([self.distance2(origin, a) for a in anchors])
        self.gram = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                self.gram[i, j] = (self.radii2[i]+self.radii2[j]-self.distance2(anchors[i], anchors[j]))/2
        self.cholesky = np.linalg.cholesky(self.gram)

    def distance2(self, a, b):
        return self.budget(a, b)**(2/self.beta)

    def covector(self, token):
        radius2 = self.distance2(self.origin, token)
        return (self.radii2+radius2-np.array([self.distance2(a, token) for a in self.anchors]))/2

    def locate(self, token):
        return np.linalg.solve(self.cholesky, self.covector(token))

    def coefficients(self, token):
        return np.linalg.solve(self.gram, self.covector(token))


def transition(a, b):
    t = b.locate(a.origin)
    columns_a = np.column_stack([a.locate(p) for p in a.anchors])
    columns_b = np.column_stack([b.locate(p)-t for p in a.anchors])
    return columns_b@np.linalg.inv(columns_a), t


@lru_cache(None)
def report():
    m = EndpointModel()
    unit = m.from_truth([1., 0., 0.])
    beta = calibrate_exponent(m.budget, m.zero, unit, m.combine(unit, unit))
    anchors = [m.from_truth(v) for v in ([1.,0.,0.],[.3,1.,0.],[-.2,.4,1.2])]
    chart = CoordinateChart(m.budget, beta, m.zero, anchors)
    points = [m.from_truth(v) for v in ([0.,0.,0.], [.7,-.4,.9], [-1.3,.2,.4], [.01,.02,-.03], [2.,-1.,.3])]
    anchor_matrix = np.column_stack([m.truth(a) for a in anchors])
    true_frame = np.linalg.solve(chart.cholesky, anchor_matrix.T)
    reconstruction = max(np.linalg.norm(chart.locate(p)-true_frame@m.truth(p)) for p in points)
    distance_error = max(abs(np.linalg.norm(chart.locate(p)-chart.locate(q))-np.linalg.norm(m.truth(p)-m.truth(q))) for p in points for q in points)
    add_error = max(np.linalg.norm(chart.locate(m.combine(p,q))-chart.locate(p)-chart.locate(q)) for p in points for q in points)
    inverse_error = max(np.linalg.norm(chart.locate(m.inverse(p))+chart.locate(p)) for p in points)
    calls = m.calls
    chart.locate(points[1])
    calls_per_point = m.calls-calls

    contract_error = 0.
    scaling_add_error = 0.
    for r in (.2, .7, 2.3):
        for p in points[1:]:
            contract_error = max(contract_error, abs(m.budget(m.zero,m.dilate(p,r))-r*m.budget(m.zero,p)))
        p, q = points[1:3]
        scaling_add_error = max(scaling_add_error, np.linalg.norm(m.truth(m.dilate(m.combine(p,q),r))-m.truth(m.combine(m.dilate(p,r),m.dilate(q,r)))))
    rot = rotation([2.,-1.,1.], .7)
    covariance_error = max(np.linalg.norm(m.direction(m.redirect(p,rot))-rot@m.direction(p)) for p in points[1:])

    charts = [chart]
    for origin_x, frame in (([.4,-.2,.5], rotation([1.,2.,3.],.5)), ([-.6,.3,.2], rotation([3.,1.,-1.],-.8))):
        origin = m.from_truth(origin_x)
        anchor_tokens = [m.from_truth(np.asarray(origin_x)+frame@anchor_matrix[:,i]) for i in range(3)]
        charts.append(CoordinateChart(m.budget,beta,origin,anchor_tokens))
    ab, tab = transition(charts[0], charts[1])
    bc, tbc = transition(charts[1], charts[2])
    ac, tac = transition(charts[0], charts[2])
    transition_error = max(np.linalg.norm(charts[1].locate(p)-ab@charts[0].locate(p)-tab) for p in points)

    # Error guarantee in a fixed actual anchor basis, not falsely aligned noisy axes.
    p = points[1]
    q = chart.covector(p)
    coefficient = np.linalg.solve(chart.gram,q)
    dg = 1e-4*np.array([[1.,.2,-.3],[.2,-.8,.4],[-.3,.4,.6]])
    dq = 1e-4*np.array([.7,-.9,.3])
    tau = float(np.linalg.norm(dg,2))
    nu = float(np.linalg.norm(dq))
    eigenvalues = np.linalg.eigvalsh(chart.gram)
    bound = (nu+tau*np.linalg.norm(coefficient))/(eigenvalues[0]-tau)
    noisy = np.linalg.solve(chart.gram+dg,q+dq)
    coefficient_error = float(np.linalg.norm(noisy-coefficient))

    # Budget normalization can spiral; even exponent-corrected polar labels fail.
    def naive(p):
        b = m.budget(m.zero,p)
        if b == 0:
            return np.zeros(3)
        return b**(1/beta)*m.direction(m.dilate(p,1/b))
    p,q = points[1:3]
    naive_error = np.linalg.norm(naive(m.combine(p,q))-naive(p)-naive(q))
    return dict(round=412,
        scope="Conditional operational coordinates from round-385 endpoint, dilation, redirection and qubit contracts; no unconditional spatial generation",
        calibration=dict(true_beta=m.beta,measured_beta=beta,unit_doubling_budget=m.budget(m.zero,m.combine(unit,unit))),
        contract=dict(homogeneity_error=float(contract_error),dilation_additivity_error=float(scaling_add_error),direction_covariance_error=float(covariance_error)),
        coordinates=dict(reconstruction_error=float(reconstruction),distance_error=float(distance_error),additivity_error=float(add_error),inverse_error=float(inverse_error),budget_queries_per_new_position=calls_per_point,anchor_gram=chart.gram.tolist(),gram_minimum_eigenvalue=float(eigenvalues[0])),
        overlap=dict(orthogonality_error=float(np.linalg.norm(ab.T@ab-np.eye(3))),point_transition_error=float(transition_error),three_chart_matrix_error=float(np.linalg.norm(bc@ab-ac)),three_chart_origin_error=float(np.linalg.norm(bc@tab+tbc-tac))),
        uncertainty=dict(gram_error_norm=tau,covector_error_norm=nu,coefficient_error=coefficient_error,coefficient_upper_bound=float(bound),physical_position_error=float(np.linalg.norm(anchor_matrix@(noisy-coefficient))),physical_position_upper_bound=float(np.sqrt(eigenvalues[-1])*bound)),
        twisted_scaling=dict(rate=m.twist_rate,naive_polar_additivity_error=float(naive_error)),
        actual_endpoint_contract_still_input=True,conditional_coordinate_construction=True,
        direction_injectivity_added_as_axiom=False,actual_so3_action_added_as_axiom=False,
        real_group_roots_required_by_readout=False,finite_statistics_certify_global_contract=False,
        globally_flat_universe_derived=False,spatial_dimension_generated_from_cognitive_axioms=False,
        new_cognitive_axiom_adopted=False,phase_closure_triggered=False)


class Audit(unittest.TestCase):
    def test_01_explicit_contract_with_twisting_dilation(self):
        for v in report()["contract"].values():
            self.assertLess(v, 2e-13)

    def test_02_exponent_obtained_by_actual_doubling(self):
        r=report()["calibration"]
        self.assertAlmostEqual(r["measured_beta"],r["true_beta"],places=13)
        self.assertGreater(r["unit_doubling_budget"],1)

    def test_03_four_budget_readout_reconstructs_oblique_coordinates(self):
        r=report()["coordinates"]
        self.assertLess(r["reconstruction_error"],1e-13)
        self.assertEqual(r["budget_queries_per_new_position"],4)
        self.assertGreater(r["gram_minimum_eigenvalue"],.4)

    def test_04_addition_inverse_and_distance(self):
        r=report()["coordinates"]
        for k in ("additivity_error","inverse_error","distance_error"):
            self.assertLess(r[k],1e-13)

    def test_05_actual_common_anchor_chart_transition(self):
        for v in report()["overlap"].values():
            self.assertLess(v,1e-13)

    def test_06_calibrated_error_bound(self):
        r=report()["uncertainty"]
        self.assertLessEqual(r["coefficient_error"],r["coefficient_upper_bound"])
        self.assertLessEqual(r["physical_position_error"],r["physical_position_upper_bound"])

    def test_07_spiral_normalization_is_not_an_additive_coordinate(self):
        self.assertGreater(report()["twisted_scaling"]["naive_polar_additivity_error"],.1)

    def test_08_dependent_anchors_are_rejected(self):
        m=EndpointModel()
        anchors=[m.from_truth([1.,0.,0.]),m.from_truth([0.,1.,0.]),m.from_truth([1.,0.,0.])]
        with self.assertRaises(np.linalg.LinAlgError):
            CoordinateChart(m.budget,m.beta,m.zero,anchors)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    tests=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result=dict(report(),checks=dict(run=tests.testsRun,failures=len(tests.failures),errors=len(tests.errors)),runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
