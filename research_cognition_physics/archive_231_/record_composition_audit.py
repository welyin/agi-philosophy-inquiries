"""Round232: a known single-use consistent process and its composition failure.

The Boolean process is Baumeler--Wolf (2016), arXiv:1507.01714, Eq.(98).
We reproduce it as a restricted consistency diagnostic, not as a full F+U+C+P
model. Copy aggregation is an explicit higher-order interface assumption.
"""
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import platform
import unittest

import numpy as np

BITS=list(product(range(2),repeat=3))
FUNCS=list(product(range(2),repeat=2))
STATES2=np.array(list(product(range(4),repeat=3)),dtype=int)


def process(out):
    a,b,c=map(int,out)
    return ((1-b)*c,(1-c)*a,(1-a)*b)


def fixed_points(local):
    return [i for i in BITS if i==process(tuple(f[x] for f,x in zip(local,i)))]


def joint_distribution(kernels):
    """Joint over each local input and output; kernels[j][o,i]."""
    joint=np.zeros((2,)*6)
    for outputs in BITS:
        inputs=process(outputs)
        joint[inputs+outputs]=np.prod([k[o,i] for k,o,i in zip(kernels,outputs,inputs)])
    return joint


def process_two(out):
    low=process(tuple(int(x)%2 for x in out))
    high=process(tuple(int(x)//2 for x in out))
    return tuple(a+2*b for a,b in zip(low,high))


TABLE2=np.array([process_two(o) for o in STATES2])


def fixed_points_two(local):
    # The table is indexed lexicographically by the three base-four outputs.
    outputs=np.stack([np.asarray(f)[STATES2[:,j]] for j,f in enumerate(local)],axis=1)
    indices=16*outputs[:,0]+4*outputs[:,1]+outputs[:,2]
    return STATES2[np.all(STATES2==TABLE2[indices],axis=1)]


def reversible_two_bit_maps():
    # Independent flips, optionally composed with exchange of the two ports.
    return [tuple((((i%2)*2+i//2) if swap else i)^mask for i in range(4))
            for swap in (False,True) for mask in range(4)]


def product_maps():
    return [tuple(f[i%2]+2*g[i//2] for i in range(4)) for f,g in product(FUNCS,repeat=2)]


def composition_certificate():
    controls=reversible_two_bit_maps()
    hist=Counter(len(fixed_points_two(fs)) for fs in product(controls,repeat=3))
    witness=((1,0,3,2),(0,2,1,3),(0,2,1,3))
    return {'reversible_local_settings_checked':len(controls)**3,
            'fixed_point_count_histogram':{str(k):hist[k] for k in sorted(hist)},
            'zero_weight_witness':{'local_truth_tables':witness,
                                   'fixed_points':fixed_points_two(witness).tolist(),
                                   'normalization_weight':len(fixed_points_two(witness))}}


def single_certificate():
    counts=Counter(len(fixed_points(fs)) for fs in product(FUNCS,repeat=3))
    return {'deterministic_local_settings_checked':64,
            'fixed_point_count_histogram':{str(k):counts[k] for k in sorted(counts)},
            'truth_table':[{'outputs':o,'inputs':process(o)} for o in BITS]}


def record_certificate(seed=232):
    rng=np.random.default_rng(seed)
    normalization_error=record_marginal_error=0.
    for _ in range(128):
        kernels=[]
        for _party in range(3):
            t=rng.random((2,2));t/=t.sum(axis=0,keepdims=True);kernels.append(t)
        joint=joint_distribution(kernels)
        # Three explicit classical record bits r_j=i_j.  Copies are readouts,
        # not extra information-bearing systems sent around the feedback.
        with_records=np.zeros((2,)*9)
        for inputs in BITS:
            for outputs in BITS: with_records[inputs+outputs+inputs]=joint[inputs+outputs]
        normalization_error=max(normalization_error,abs(float(joint.sum())-1))
        record_marginal_error=max(record_marginal_error,
                                  float(np.max(np.abs(with_records.sum(axis=(6,7,8))-joint))))
    return {'random_stochastic_settings':128,'explicit_record_bits':3,
            'maximum_normalization_error':normalization_error,
            'maximum_record_discard_error':record_marginal_error,
            'memory_scope':'Three local classical input records per use; no global time order or free infinite storage.'}


def quantum_diagonal_certificate(seed=233):
    rng=np.random.default_rng(seed);worst=0.
    for _ in range(64):
        kernels=[]
        for _party in range(3):
            z=rng.normal(size=(6,2))+1j*rng.normal(size=(6,2))
            q,_=np.linalg.qr(z)
            kraus=q.reshape(3,2,2)
            # <o| Phi(|i><i|) |o>, valid for every local CPTP channel.
            kernels.append(np.sum(np.abs(kraus)**2,axis=0))
        worst=max(worst,abs(float(joint_distribution(kernels).sum())-1))
    return {'random_local_CPTP_settings':64,'maximum_normalization_error':worst}


class RecordCompositionTests(unittest.TestCase):
    def test_single_use_all_deterministic_operations_have_one_fixed_point(self):
        self.assertEqual(single_certificate()['fixed_point_count_histogram'],{'1':64})

    def test_stochastic_operations_and_record_marginals(self):
        c=record_certificate()
        self.assertLess(c['maximum_normalization_error'],2e-15)
        self.assertEqual(c['maximum_record_discard_error'],0.)

    def test_record_queries_are_compatible_without_a_time_label(self):
        kernels=[np.array([[.2,.7],[.8,.3]])]*3
        joint=joint_distribution(kernels)
        records=joint.sum(axis=(3,4,5))
        # Two query orders are ordinary marginals of the same joint record.
        np.testing.assert_allclose(records.sum(axis=2).sum(axis=1),
                                   records.sum(axis=1).sum(axis=1))
        self.assertAlmostEqual(float(records.sum()),1.)

    def test_no_party_can_be_first_for_all_independent_constant_settings(self):
        answers=np.array([process(s) for s in BITS])
        for party in range(3):
            for own in (0,1):
                values=answers[[s[party]==own for s in BITS],party]
                self.assertEqual(sorted(values.tolist()),[0,0,0,1])
        # First-party correct-guess probability bounds joint success.
        self.assertEqual(3/4,.75)

    def test_diagonal_quantum_channels_reduce_to_stochastic_kernels(self):
        self.assertLess(quantum_diagonal_certificate()['maximum_normalization_error'],4e-15)

    def test_two_uses_with_unmixed_local_ports_remain_consistent(self):
        fs=product_maps()
        self.assertEqual(len(fs),16)
        for local in product(fs,repeat=3):
            self.assertEqual(len(fixed_points_two(local)),1)

    def test_two_uses_can_fail_with_only_reversible_local_controls(self):
        c=composition_certificate()
        self.assertEqual(c['fixed_point_count_histogram'],{'0':66,'1':384,'2':60,'4':2})
        self.assertEqual(c['zero_weight_witness']['normalization_weight'],0)
        for f in c['zero_weight_witness']['local_truth_tables']:
            self.assertEqual(sorted(f),list(range(4)))

    def test_boolean_contradiction_for_explicit_witness(self):
        # If all six consistency equations held, elimination would force
        # y_B=x_C=y_A=0, y_C=x_B=x_A and x_B=1-x_B.
        witness=((1,0,3,2),(0,2,1,3),(0,2,1,3))
        direct=[]
        for x,y in product(BITS,repeat=2):
            o=(1-x[0],y[1],y[2]);p=(y[0],x[1],x[2])
            if x==process(o) and y==process(p): direct.append((x,y))
        self.assertEqual(direct,[])
        self.assertEqual(fixed_points_two(witness).tolist(),[])
        self.assertTrue(all(b!=1-b for b in (0,1)))


def report():
    return {'round':232,'date':'2026-09-22',
            'question':'Does compatible finite recording imply event order, and does single-use consistency survive independent composition with local regrouping?',
            'source':'Baumeler and Wolf, NJP18 013036 (2016), arXiv:1507.01714 Eq.(98).',
            'single_use':single_certificate(),'records':record_certificate(),
            'quantum_diagonal_check':quantum_diagonal_certificate(),
            'two_use_product_operations':{'settings_checked':4096,'all_fixed_point_counts':1},
            'two_use_regrouped_operations':composition_certificate(),
            'analytic_results':['Finite compatible records alone do not force a fixed global event order.',
                                'This single-use process fails independent two-copy composition with arbitrary local regrouping.',
                                'A local X on the first A port and SWAP at B and C give normalization zero.',
                                'Dense local permissions already exclude a normalization discontinuity escape.'],
            'scope':'Restricted classical/diagonal process diagnostic, not a complete F+U+C+P countermodel or a derivation of spacetime.',
            'extra_input':'Treating a whole process as an independently repeatable resource with jointly accessible ports is an explicit higher-order interface contract.',
            'next':'Specify physically legitimate port regrouping and test that it does not reject ordinary multi-event or adaptive circuits by misidentifying a laboratory with one event.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RecordCompositionTests))
    if not result.wasSuccessful(): raise SystemExit(1)
    data=report();data['checks']={'run':result.testsRun,'failures':0,'errors':0}
    if args.write_results:
        target=Path(__file__).with_name('record_composition_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8'))!=json.loads(json.dumps(data)):
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
