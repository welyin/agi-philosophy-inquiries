"""Round 242: coherent schedules versus duplicate labels of one history.

All circuits and event interfaces are supplied examples. The quantum switch is
an existing higher-order construction, not a derived spacetime growth rule.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np

from quantum_history_record_audit import I, X, Z, H, P, PLUS, density, channel


def choi(ks):
    d = ks[0].shape[1]
    return sum(density(k.reshape(-1)) for k in ks)/d


def controlled(u):
    return np.kron(P[0], I)+np.kron(P[1], u)


def switch_kraus(a, b):
    return [np.kron(P[0], bj@ai)+np.kron(P[1], ai@bj)
            for ai in a for bj in b]


def control_plus_probability(ks):
    rho = np.kron(density(PLUS), I/2)
    return float(np.trace(np.kron(density(PLUS), I)@channel(ks, rho)).real)


def mix_kraus(ks, mixing):
    return [sum(mixing[i,j]*k for j,k in enumerate(ks)) for i in range(len(ks))]


def linear_extensions():
    # Events 0<1 on A, 2<3 on B; cross-chain events are incomparable.
    return [p for p in itertools.permutations(range(4))
            if p.index(0)<p.index(1) and p.index(2)<p.index(3)]


def event_operators(outcomes=None):
    theta = .37
    rotation = np.array([[math.cos(theta), -math.sin(theta)],
                         [math.sin(theta), math.cos(theta)]], dtype=complex)
    end_a, end_b = Z, X
    if outcomes is not None:
        a, b = outcomes
        end_a = P[a]
        end_b = math.sqrt([.7,.3][b])*[I,X][b]
    return [np.kron(H,I),np.kron(end_a,I),
            np.kron(I,rotation),np.kron(I,end_b)]


def product(operators, schedule):
    result = np.eye(len(operators[0]),dtype=complex)
    for e in schedule:
        result = operators[e]@result
    return result


def schedule_products(outcomes=None):
    return [product(event_operators(outcomes),s) for s in linear_extensions()]


def coherent_schedule_output():
    products = schedule_products()
    m, d = len(products), len(products[0])
    basis = np.eye(m)
    w = sum(np.kron(density(basis[:,j]),u) for j,u in enumerate(products))
    uniform = np.ones(m)/math.sqrt(m)
    # Entangled A,B input; this is a known controlled finite circuit.
    psi = np.array([1,0,0,1],dtype=complex)/math.sqrt(2)
    target = density(products[0]@psi)
    output = w@np.kron(density(uniform),density(psi))@w.conj().T
    expected = np.kron(density(uniform),target)
    return output,expected,float(np.trace(np.kron(density(uniform),np.eye(d))@output).real)


def report():
    phases = [0.,math.pi/2,math.pi]
    phase_rows = [{'phase':p,
                   'ordinary_channel_max_difference':float(np.max(np.abs(choi([np.exp(1j*p)*I])-choi([I])))),
                   'controlled_plus_probability':control_plus_probability([controlled(np.exp(1j*p)*I)])}
                  for p in phases]
    exact = schedule_products()
    selective_error = max(float(np.max(np.abs(u-schedule_products(outcome)[0])))
                          for outcome in itertools.product(range(2),repeat=2)
                          for u in schedule_products(outcome))
    a = [math.sqrt(.7)*I,math.sqrt(.3)*X]
    b = [math.sqrt(.6)*I,math.sqrt(.4)*Z]
    mixing = np.array([[1,1j],[1j,1]],dtype=complex)/math.sqrt(2)
    mixed_error = np.max(np.abs(choi(switch_kraus(a,b))-choi(switch_kraus(mix_kraus(a,H),mix_kraus(b,mixing)))))
    output,expected,prob = coherent_schedule_output()
    return {'round':242,'date':'2026-09-22',
            'controlled_phase_rows':phase_rows,
            'pauli_order':{
                'ordinary_composite_channel_difference':float(np.max(np.abs(choi([Z@X])-choi([X@Z])))),
                'switch_plus_probability':control_plus_probability(switch_kraus([X],[Z])),
                'commuting_comparison_plus_probability':control_plus_probability(switch_kraus([X],[X])),
                'switch_kraus_mixing_difference':float(mixed_error)},
            'two_independent_chains':{
                'natural_schedules':[list(s) for s in linear_extensions()],
                'unitary_product_max_difference':float(max(np.max(np.abs(u-exact[0])) for u in exact)),
                'selective_product_max_difference':selective_error,
                'coherent_schedule_factorization_error':float(np.max(np.abs(output-expected))),
                'uniform_control_readout_probability':prob},
            'duplicate_label_counterexample':{
                'physical_probabilities':[.5,.5], 'label_multiplicities':[1,2],
                'renormalized_probability_counting':[1/3,2/3],
                'renormalized_equal_phase_amplitude_counting':[.2,.8]},
            'extra_inputs':['Specified finite event poset and operators',
                            'Known coherent control implementation or supplied quantum switch',
                            'Exact operator interchange for incomparable events, including outcome operators'],
            'conclusion':'Operator interchange makes finite natural schedules equivalent; ordinary channel equality alone does not identify histories in coherent contexts.',
            'not_claimed':['A quantum switch has been derived from fixed-order black-box calls',
                           'Natural-label multiplicity is physical interference',
                           'The event order, amplitudes, or a causal-set quantum growth law has been selected'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class LabelTests(unittest.TestCase):
    def test_same_channel_different_controlled_phase(self):
        for p in [0.,.2,math.pi/2,math.pi]:
            u = np.exp(1j*p)*I
            np.testing.assert_allclose(choi([u]),choi([I]),atol=1e-14)
            self.assertAlmostEqual(control_plus_probability([controlled(u)]),(1+math.cos(p))/2)

    def test_commuting_channels_need_not_have_same_coherent_order_readout(self):
        np.testing.assert_allclose(choi([Z@X]),choi([X@Z]),atol=1e-14)
        self.assertAlmostEqual(control_plus_probability(switch_kraus([X],[Z])),0)
        self.assertAlmostEqual(control_plus_probability(switch_kraus([X],[X])),1)

    def test_switch_is_invariant_under_phases_of_each_called_unitary(self):
        original = choi(switch_kraus([X],[Z]))
        phased = choi(switch_kraus([np.exp(.31j)*X],[np.exp(-.72j)*Z]))
        np.testing.assert_allclose(original,phased,atol=1e-14)

    def test_switch_channel_is_independent_of_kraus_basis(self):
        a = [math.sqrt(.7)*I,math.sqrt(.3)*X]
        b = [math.sqrt(.6)*I,math.sqrt(.4)*Z]
        mixing = np.array([[1,1j],[1j,1]],dtype=complex)/math.sqrt(2)
        ks = switch_kraus(a,b)
        np.testing.assert_allclose(sum(k.conj().T@k for k in ks),np.eye(4),atol=1e-14)
        np.testing.assert_allclose(choi(ks),choi(switch_kraus(mix_kraus(a,H),mix_kraus(b,mixing))),atol=1e-14)

    def test_all_natural_schedules_agree_with_operator_interchange(self):
        schedules = linear_extensions()
        self.assertEqual(len(schedules),6)
        ops = event_operators()
        self.assertGreater(np.max(np.abs(ops[0]@ops[1]-ops[1]@ops[0])),1)
        for a in [0,1]:
            for b in [2,3]:
                np.testing.assert_allclose(ops[a]@ops[b],ops[b]@ops[a],atol=1e-14)
        for u in schedule_products():
            np.testing.assert_allclose(u,schedule_products()[0],atol=1e-14)

    def test_interchange_includes_selective_outcomes_and_complete_instruments(self):
        ks = []
        for outcome in itertools.product(range(2),repeat=2):
            products = schedule_products(outcome)
            for k in products:
                np.testing.assert_allclose(k,products[0],atol=1e-14)
            ks.append(products[0])
        np.testing.assert_allclose(sum(k.conj().T@k for k in ks),np.eye(4),atol=1e-14)

    def test_normalized_coherent_schedule_register_does_not_multiply_weight(self):
        output,expected,prob = coherent_schedule_output()
        np.testing.assert_allclose(output,expected,atol=1e-14)
        self.assertAlmostEqual(prob,1)
        self.assertAlmostEqual(np.trace(output).real,1)

    def test_counting_duplicate_labels_biases_physical_distribution(self):
        probabilities = np.array([.5,.5])
        multiplicity = np.array([1,2])
        naive_prob = probabilities*multiplicity
        naive_amplitude = probabilities*multiplicity**2
        np.testing.assert_allclose(naive_prob/naive_prob.sum(),[1/3,2/3])
        np.testing.assert_allclose(naive_amplitude/naive_amplitude.sum(),[.2,.8])
        self.assertFalse(np.allclose(naive_prob/naive_prob.sum(),probabilities))
        self.assertFalse(np.allclose(naive_amplitude/naive_amplitude.sum(),probabilities))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LabelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('coherent_label_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
