"""Round 245: rank limits for generalized records in a histories Gram space.

Projector records and their physical persistence/accessibility are different
requirements. The no-record statement preserves the original decoherence
functional; it does not forbid a changed model or another interpretation.
"""
import argparse
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import complex_growth_extension_audit as ext
from quantum_history_record_audit import I, X, H, P, density


def gram(vectors):
    return vectors.conj().T@vectors


def correlation(d):
    scale = np.sqrt(np.outer(d.diagonal().real,d.diagonal().real))
    return np.divide(np.abs(d),scale,out=np.zeros_like(scale),where=scale>1e-30)


def record_error(vectors, projectors):
    return max(float(np.linalg.norm(r@vectors[:,j]-(vectors[:,j] if i==j else 0)))
               for i,r in enumerate(projectors) for j in range(vectors.shape[1]))


def weak_example():
    a = np.array([z for z in ext.level(3,2).values() if abs(z)>0])
    return a.reshape(1,-1)


def weak_capacity_example(rank):
    return np.column_stack([phase*np.eye(rank)[:,j]/math.sqrt(2*rank)
                            for j in range(rank) for phase in [1,1j]])


def trine_vectors():
    return np.array([[1,1,1],[1,np.exp(2j*math.pi/3),np.exp(4j*math.pi/3)]])/math.sqrt(2)


def report():
    weak = gram(weak_example())
    a = np.array(list(ext.level(4,2).values()))
    d = np.outer(a.conj(),a)
    active = np.flatnonzero(abs(a)>1e-14)
    corr = correlation(d[np.ix_(active,active)])
    off = corr[~np.eye(len(active),dtype=bool)]
    r0,r1 = np.kron(P[0],I),np.kron(P[1],I)
    vectors = np.column_stack([np.array([1,0,0,0]),np.array([0,0,1,0])])/math.sqrt(2)
    stable = np.kron(P[0],H)+np.kron(P[1],X)
    mixing = np.kron(H,I)
    return {'round':245,'date':'2026-09-22',
            'scope':'Generalized exact projector records for specified history vectors; original Gram functional fixed',
            'weak_consistency_counterexample':{
                'D_real':weak.real.tolist(),'D_imag':weak.imag.tolist(),
                'rank':int(np.linalg.matrix_rank(weak)),
                'sum_diagonal':float(np.trace(weak).real),
                'normalized_branch_density_distance':0.,
                'equal_prior_best_discrimination_success':.5},
            'rank_one_growth_N4':{
                'nonzero_cylinders':len(active),
                'normalized_offdiagonal_min':float(off.min()),
                'normalized_offdiagonal_max':float(off.max()),
                'original_total_weight':float(d.sum().real),
                'weight_after_deleting_offdiagonals':float(np.trace(d).real)},
            'capacity_bounds':{
                'exact_complex_records':'m <= rank(D)',
                'weak_real_consistency':'m <= 2*rank(D), does not imply records',
                'uniform_normalized_overlap':'delta^2 >= (m-r)/(r*(m-1)) when m>r',
                'trine_r2_m3_overlap':float(abs(gram(trine_vectors())[0,1]))},
            'record_persistence_example':{
                'block_preserving_future_error':record_error(stable@vectors,[r0,r1]),
                'record_mixing_future_error':record_error(mixing@vectors,[r0,r1])},
            'extra_inputs':['Gram representation of the supplied decoherence functional',
                            'Exact orthogonal projector criterion for generalized records',
                            'Specified future dynamics when testing persistence'],
            'not_claimed':['Generalized record projectors are automatically locally accessible',
                           'No model can have records after changing its functional',
                           'All histories interpretations must adopt this record criterion'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class RankTests(unittest.TestCase):
    def test_exact_orthogonal_records_saturate_rank_capacity(self):
        probabilities = np.array([.2,.3,.5])
        vectors = np.diag(np.sqrt(probabilities)).astype(complex)
        records = [np.diag(np.eye(3)[j]) for j in range(3)]
        self.assertEqual(np.linalg.matrix_rank(gram(vectors)),3)
        self.assertLess(record_error(vectors,records),1e-14)
        self.assertAlmostEqual(gram(vectors).sum().real,1)

    def test_weak_consistency_can_have_identical_conditional_states(self):
        v = weak_example()
        d = gram(v)
        self.assertAlmostEqual(d[0,1].real,0)
        self.assertAlmostEqual(abs(d[0,1]),.5)
        states = [density(v[:,j]/np.linalg.norm(v[:,j])) for j in range(2)]
        np.testing.assert_allclose(states[0],states[1],atol=1e-14)
        self.assertEqual(np.linalg.matrix_rank(d),1)

    def test_all_nonzero_rank_one_growth_branches_have_unit_correlation(self):
        a = np.array([z for z in ext.level(5,2).values() if abs(z)>0])
        np.testing.assert_allclose(correlation(np.outer(a.conj(),a)),np.ones((len(a),len(a))),atol=1e-13)

    def test_blank_ancilla_isometry_preserves_the_obstruction(self):
        v = weak_example()
        embedding = np.array([[1],[1j],[2],[-1j]],dtype=complex)/math.sqrt(7)
        np.testing.assert_allclose(gram(embedding@v),gram(v),atol=1e-14)
        out = embedding@v
        rho = [density(out[:,j]/np.linalg.norm(out[:,j])) for j in range(2)]
        np.testing.assert_allclose(rho[0],rho[1],atol=1e-14)

    def test_weak_consistency_capacity_is_twice_complex_rank(self):
        for r in [1,2,3]:
            d = gram(weak_capacity_example(r))
            np.testing.assert_allclose(d.real,np.eye(2*r)/(2*r),atol=1e-14)
            self.assertEqual(np.linalg.matrix_rank(d),r)
            self.assertEqual(np.linalg.matrix_rank(d.real),2*r)
            self.assertAlmostEqual(d.sum().real,1)

    def test_frame_overlap_bound_and_saturating_example(self):
        v = trine_vectors()
        d = gram(v)
        off = np.abs(d)[~np.eye(3,dtype=bool)]
        np.testing.assert_allclose(off,.5,atol=1e-14)
        self.assertAlmostEqual(off.max()**2,(3-2)/(2*(3-1)))
        self.assertAlmostEqual(np.sum(abs(d)**2),3**2/2)

    def test_deleting_interference_is_a_change_of_normalized_model(self):
        row = report()['rank_one_growth_N4']
        self.assertAlmostEqual(row['original_total_weight'],1)
        self.assertAlmostEqual(row['weight_after_deleting_offdiagonals'],.5)

    def test_future_evolution_needs_a_record_preservation_condition(self):
        row = report()['record_persistence_example']
        self.assertLess(row['block_preserving_future_error'],1e-14)
        self.assertGreater(row['record_mixing_future_error'],.4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RankTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('history_record_rank_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
