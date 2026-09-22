"""Round 241: bridge finite operations to a histories decoherence functional.

The interferometer is a supplied finite interface example, not a spacetime
growth law. Ordinary channel erasure, conditioning and coherent feedforward
are kept distinct. D[h,k] = Tr(C_h^dagger C_k rho).
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Z = np.diag([1,-1]).astype(complex)
H = (X+Z)/math.sqrt(2)
P = [np.diag([1,0]).astype(complex), np.diag([0,1]).astype(complex)]
PLUS = np.array([1,1],dtype=complex)/math.sqrt(2)


def density(v):
    return np.outer(v, v.conj())


def channel(ks, rho):
    return sum(k@rho@k.conj().T for k in ks)


def trace_distance(a,b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def gram(operators, rho):
    return np.array([[np.trace(a.conj().T@b@rho) for b in operators] for a in operators])


def measure(d, subset):
    if not subset:
        return 0.
    return float(d[np.ix_(subset,subset)].sum().real)


def history_model(record_overlap=1.):
    histories = [(a,x) for x in range(2) for a in range(2)]
    operators = [P[x]@H@P[a] for a,x in histories]
    d = gram(operators, density(PLUS))
    # Record states |r0>=|0>, |r1>=gamma|0>+sqrt(1-|gamma|^2)|1>.
    gamma = complex(record_overlap)
    if abs(gamma)>1+1e-14:
        raise ValueError('Record overlap must have modulus <=1')
    record_gram = np.array([[1,gamma],[gamma.conjugate(),1]],dtype=complex)
    factors = np.array([[record_gram[a,b] for b,y in histories] for a,x in histories])
    return histories,operators,d*factors


def copy_marker(rho):
    # rho acts on target A and arbitrary reference R; output axes are A,R,M.
    d = len(rho)
    ref = d//2
    v = sum(np.kron(np.kron(P[a],np.eye(ref)),np.eye(2)[:,a:a+1]) for a in range(2))
    return v@rho@v.conj().T


def trace_last(rho, d):
    return np.trace(rho.reshape(d,2,d,2),axis1=1,axis2=3)


def marker_channel(joint, d, ks):
    return channel([np.kron(np.eye(d),k) for k in ks],joint)


def erase_marker(joint, d, feedback=False):
    ref = d//2
    output = []
    for sign in [1,-1]:
        bra = np.array([[1,sign]],dtype=complex)/math.sqrt(2)
        k = np.kron(np.eye(d),bra)
        block = k@joint@k.conj().T
        if feedback and sign == -1:
            correction = np.kron(Z,np.eye(ref))
            block = correction@block@correction.conj().T
        output.append(block)
    return output


def all_i3_error(d):
    errors = []
    for assignment in itertools.product(range(4),repeat=len(d)):
        a,b,c = ([i for i,label in enumerate(assignment) if label == j] for j in range(3))
        value = (measure(d,a+b+c)-measure(d,a+b)-measure(d,b+c)-measure(d,a+c)
                 +measure(d,a)+measure(d,b)+measure(d,c))
        errors.append(abs(value))
    return max(errors)


def report():
    histories,ops,d = history_model()
    _,_,recorded = history_model(0.)
    plus = density(PLUS)
    tagged = copy_marker(plus)
    reset = [np.array([[1,0],[0,0]],dtype=complex),np.array([[0,1],[0,0]],dtype=complex)]
    reset_state = trace_last(marker_channel(tagged,2,reset),2)
    erased = erase_marker(tagged,2)
    feedback = sum(erase_marker(tagged,2,True))
    bell = density(np.array([1,0,0,1],dtype=complex)/math.sqrt(2))
    bell_tagged = copy_marker(bell)
    copied_elsewhere = marker_channel(bell_tagged,4,P)
    failed_recovery = sum(erase_marker(copied_elsewhere,4,True))
    invalid = gram([I/math.sqrt(2)]*2,plus)
    return {'round':241,'date':'2026-09-22',
            'convention':'D[h,k]=Tr(C_h^dagger C_k rho); histories ordered as (path,final port)',
            'histories':[list(h) for h in histories],
            'unrecorded_D_real':d.real.tolist(),
            'recorded_D_real':recorded.real.tolist(),
            'history_hilbert_ranks':{'unrecorded':int(np.linalg.matrix_rank(d)),
                                     'orthogonal_marker':int(np.linalg.matrix_rank(recorded))},
            'normalization':{'unrecorded':float(d.sum().real),'recorded':float(recorded.sum().real),
                             'bare_Kraus_Gram_counterexample':float(invalid.sum().real)},
            'unrecorded_nonadditivity':{'port0_union':measure(d,[0,1]),
                                        'port0_diagonal_sum':float(d[0,0].real+d[1,1].real),
                                        'three_history_subset_measure':measure(d,[0,1,2]),
                                        'all_256_disjoint_triples_I3_error':all_i3_error(d)},
            'final_port0_probabilities':{
                'no_marker':float(np.trace(P[0]@H@plus@H).real),
                'orthogonal_marker_ignored':float(np.trace(P[0]@H@trace_last(tagged,2)@H).real),
                'marker_reset_only':float(np.trace(P[0]@H@reset_state@H).real),
                'eraser_conditioned_plus':float(np.trace(P[0]@H@erased[0]@H).real/np.trace(erased[0]).real),
                'eraser_plus_probability':float(np.trace(erased[0]).real),
                'eraser_with_feedback':float(np.trace(P[0]@H@feedback@H).real)},
            'unknown_state_reference_recovery':{
                'coherent_marker_feedback_trace_distance':trace_distance(sum(erase_marker(bell_tagged,4,True)),bell),
                'unavailable_extra_record_trace_distance':trace_distance(failed_recovery,bell)},
            'extra_inputs':['Given finite interferometer and projective history decomposition',
                            'Specified marker coupling, accessibility and feedback'],
            'conclusion':'Strongly positive histories are available conditionally; CP instruments alone do not fix coherent history weights or record treatment.',
            'not_claimed':['A natural event-growth measure has been selected',
                           'Ignoring or resetting a record locally restores interference',
                           'All coherent erasure requires postselection: feedback can be deterministic when records remain accessible'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class HistoryTests(unittest.TestCase):
    def test_gram_is_strongly_positive_normalized_and_has_quantum_sum_rule(self):
        for overlap in [0.,.4,1.,.5j]:
            _,ops,d = history_model(overlap)
            np.testing.assert_allclose(d,d.conj().T,atol=1e-14)
            self.assertGreater(np.linalg.eigvalsh(d).min(),-1e-14)
            self.assertAlmostEqual(d.sum().real,1)
            self.assertLess(all_i3_error(d),1e-14)
        np.testing.assert_allclose(sum(ops),H,atol=1e-14)

    def test_history_space_rank_and_measure_need_not_be_probability(self):
        _,_,d = history_model()
        _,_,r = history_model(0.)
        self.assertEqual(np.linalg.matrix_rank(d),2)
        self.assertEqual(np.linalg.matrix_rank(r),4)
        self.assertAlmostEqual(measure(d,[0,1,2]),1.25)
        self.assertAlmostEqual(measure(d,[0,1]),1)
        self.assertAlmostEqual(d[0,0].real+d[1,1].real,.5)

    def test_instrument_completeness_does_not_normalize_bare_history_gram(self):
        ks = [I/math.sqrt(2)]*2
        np.testing.assert_allclose(sum(k.conj().T@k for k in ks),I,atol=1e-14)
        d = gram(ks,density(PLUS))
        self.assertAlmostEqual(d.sum().real,2)
        self.assertAlmostEqual(np.trace(d).real,1)
        self.assertAlmostEqual(np.diag(np.diag(d)).sum().real,1)

    def test_record_overlap_controls_interference(self):
        for gamma in [0.,.25,1.,.5j,-.7]:
            _,_,d = history_model(gamma)
            self.assertAlmostEqual(measure(d,[0,1]),(1+complex(gamma).real)/2)
            self.assertAlmostEqual(measure(d,[2,3]),(1-complex(gamma).real)/2)

    def test_marker_only_trace_preserving_operations_do_not_change_marginal(self):
        bell = density(np.array([1,0,0,1],dtype=complex)/math.sqrt(2))
        tagged = copy_marker(bell)
        marginal = trace_last(tagged,4)
        damp = [np.diag([1,math.sqrt(.3)]),np.array([[0,math.sqrt(.7)],[0,0]])]
        reset = [np.array([[1,0],[0,0]]),np.array([[0,1],[0,0]])]
        for ks in [[H],P,damp,reset]:
            np.testing.assert_allclose(trace_last(marker_channel(tagged,4,ks),4),marginal,atol=1e-14)

    def test_eraser_subensembles_keep_success_weights(self):
        blocks = erase_marker(copy_marker(density(PLUS)),2)
        for block in blocks:
            self.assertAlmostEqual(np.trace(block).real,.5)
        self.assertAlmostEqual(np.trace(P[0]@H@blocks[0]@H).real,.5)
        self.assertAlmostEqual(np.trace(P[0]@H@blocks[1]@H).real,0)
        np.testing.assert_allclose(sum(blocks),I/2,atol=1e-14)

    def test_feedback_recovers_unknown_state_including_reference(self):
        corrected = [I/math.sqrt(2),Z@Z/math.sqrt(2)]
        # Choi input checks the complete channel, not just one uncorrelated input.
        bell = density(np.array([1,0,0,1],dtype=complex)/math.sqrt(2))
        np.testing.assert_allclose(sum(erase_marker(copy_marker(bell),4,True)),bell,atol=1e-14)
        np.testing.assert_allclose(channel(corrected,density(PLUS)),density(PLUS),atol=1e-14)

    def test_an_unavailable_copy_of_the_record_blocks_that_recovery(self):
        bell = density(np.array([1,0,0,1],dtype=complex)/math.sqrt(2))
        tagged = copy_marker(bell)
        decohered = marker_channel(tagged,4,P)
        recovered = sum(erase_marker(decohered,4,True))
        target = channel([np.kron(a,I) for a in P],bell)
        np.testing.assert_allclose(recovered,target,atol=1e-14)
        self.assertAlmostEqual(trace_distance(recovered,bell),.5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(HistoryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('quantum_history_record_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
