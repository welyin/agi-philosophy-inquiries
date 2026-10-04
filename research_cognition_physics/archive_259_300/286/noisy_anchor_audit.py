"""Round 286: finite-confidence anchor mapping, including correlated noise.

Noise acts only on marker readout. Port transport, route memory and the graph
remain reliable. Exact binomial risks and first-error bounds carry the proof;
seeded Monte Carlo is only an end-to-end implementation check.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from anchored_position_audit import PortProbe, map_ball, certify_map, bounds
from observable_cover_audit import cube, complete_four


def ball_upper(radius):
    if radius < 0:
        raise ValueError('nonnegative radius required')
    return 3*2**radius-2


def query_budget(radius):
    return 3*ball_upper(radius)**2


def majority_error(p, repetitions):
    """Binomial tail, evaluated in log space; repetitions must be odd."""
    if not 0 <= p < .5 or repetitions < 1 or repetitions % 2 != 1:
        raise ValueError('require 0 <= p < 1/2 and positive odd repetitions')
    if p == 0:
        return 0.
    m = repetitions
    return math.fsum(math.exp(math.lgamma(m+1)-math.lgamma(k+1)-math.lgamma(m-k+1)
                              +k*math.log(p)+(m-k)*math.log1p(-p))
                     for k in range((m+1)//2,m+1))


def exponential_bound(p, repetitions):
    if p == 0:
        return 0.
    return math.exp(repetitions*math.log(2*math.sqrt(p*(1-p))))


def choose_repetitions(p, radius, delta, exact=False):
    if not 0 <= p < .5 or not 0 < delta < 1:
        raise ValueError('require p < 1/2 and a probability error budget')
    if p == 0:
        return 1
    qmax = query_budget(radius)
    if exact:
        m = 1
        while qmax*majority_error(p,m) > delta:
            m += 2
        return m
    m = max(1,math.ceil(math.log(qmax/delta)/(-math.log(2*math.sqrt(p*(1-p))))))
    return m if m % 2 else m+1


class MajorityDetector:
    def __init__(self, probe, p, repetitions, seed, correlated=False):
        majority_error(p,repetitions)  # Validate the parameter domain.
        self.probe, self.p, self.m = probe,p,repetitions
        self.rng = np.random.default_rng(seed)
        self.correlated = correlated
        self.blocks = self.raw_reads = self.errors = 0

    def __call__(self):
        truth = self.probe.at_anchor()
        # Every raw read is represented in the probe's physical read ledger.
        observations = [truth] + [self.probe.at_anchor() for _ in range(self.m-1)]
        assert all(value == truth for value in observations)
        if self.correlated:
            shared = bool(self.rng.random() < self.p)
            flips = np.full(self.m,shared,dtype=bool)
        else:
            flips = self.rng.random(self.m) < self.p
        votes = np.logical_xor(np.asarray(observations,dtype=bool),flips)
        answer = bool(np.count_nonzero(votes) > self.m//2)
        self.blocks += 1
        self.raw_reads += self.m
        self.errors += answer != truth
        return answer


def risk_case(radius,p,delta=.01):
    size = ball_upper(radius)
    qmax = query_budget(radius)
    chosen = choose_repetitions(p,radius,delta)
    exact = choose_repetitions(p,radius,delta,exact=True)
    return {'radius':radius,'ball_size_upper':size,'queries_before_first_error_upper':qmax,
            'single_read_error':p,'target_failure_probability':delta,
            'repetitions_exponential_bound':chosen,'repetitions_exact_union_bound':exact,
            'single_block_exact_error':majority_error(p,chosen),
            'whole_map_union_bound':qmax*majority_error(p,chosen),
            'whole_map_exponential_bound':qmax*exponential_bound(p,chosen),
            'good_event_raw_reads_upper':qmax*chosen,
            'good_event_traversals_upper':bounds(radius,size)['steps_upper'],
            'perfectly_correlated_block_error':p}


def experiment(repetitions,correlated=False,trials=160):
    rows = cube()
    correct = comparisons = raw = errors = 0
    for seed in range(trials):
        probe = PortProbe(rows)
        detector = MajorityDetector(probe,.1,repetitions,10000+seed,correlated)
        value = map_ball(probe,2,detector)
        correct += certify_map(rows,0,2,value)
        comparisons += value['comparisons']
        raw += detector.raw_reads
        errors += detector.errors
        assert probe._position == 0 and probe.reads == detector.raw_reads
    return {'graph':'cube','radius':2,'single_read_error':.1,'repetitions':repetitions,
            'correlation':'one shared flip per read block' if correlated else 'independent fresh flips',
            'trials':trials,'exact_maps':correct,'incorrect_maps':trials-correct,
            'comparison_calls_total':comparisons,'raw_marker_reads_total':raw,
            'wrong_read_blocks_total':errors,'seed_start':10000,
            'purpose':'Finite implementation check; confidence theorem does not use empirical success rate.'}


def report():
    chosen = choose_repetitions(.1,2,.01)
    return {'round':286,
            'scope':'Finite-confidence classical position mapping under explicit independent marker noise; stable ports and route memory remain additional assumptions.',
            'risk_table':[risk_case(r,p) for r in (1,3,6) for p in (.1,.25,.4)],
            'experiments':[experiment(1),experiment(chosen),
                           experiment(chosen,correlated=True)],
            'noise_counterexample':{
                'p':.1,'repetitions':[1,5,25,101],
                'independent_majority_errors':[majority_error(.1,m) for m in (1,5,25,101)],
                'shared_flip_majority_errors':[.1]*4,
                'reason':'All repeated reads inherit one latent inversion; majority cannot remove this error.'},
            'resource_scope':'An error guarantee requires a justified p bound and fresh read noise; one anchor still needs retention and physical reads. No energy or universal clock conversion is supplied.',
            'dimension_scope':'A correct finite graph ball gives finite-scale counts and root hop distances, not an infinite-limit dimension or selection of three.'}


class Audit(unittest.TestCase):
    def test_01_exact_binomial(self):
        from itertools import product
        for p in (.05,.2,.4):
            for m in (1,3,5,7):
                enumeration = math.fsum(p**sum(bits)*(1-p)**(m-sum(bits))
                                        for bits in product((0,1),repeat=m)
                                        if sum(bits) > m//2)
                self.assertAlmostEqual(majority_error(p,m),enumeration,places=13)

    def test_02_exponential_bound(self):
        for p in (.01,.1,.25,.4,.49):
            for m in (1,3,11,51,151):
                self.assertLessEqual(majority_error(p,m),exponential_bound(p,m)*(1+1e-12))

    def test_03_certified_and_minimal_repetitions(self):
        for row in report()['risk_table']:
            self.assertLessEqual(row['whole_map_exponential_bound'],.01*(1+1e-12))
            q = row['queries_before_first_error_upper']
            p = row['single_read_error']
            m = row['repetitions_exact_union_bound']
            self.assertLessEqual(q*majority_error(p,m),.01)
            if m > 1:
                self.assertGreater(q*majority_error(p,m-2),.01)

    def test_04_zero_noise_and_read_ledger(self):
        for rows in (cube(),complete_four()):
            probe = PortProbe(rows)
            detector = MajorityDetector(probe,0.,5,42)
            value = map_ball(probe,2,detector)
            self.assertTrue(certify_map(rows,0,2,value))
            self.assertEqual(probe.reads,detector.blocks*5)
            self.assertEqual(detector.errors,0)
            self.assertLessEqual(detector.blocks,query_budget(2))

    def test_05_end_to_end_independent(self):
        m = choose_repetitions(.1,2,.01)
        row = experiment(m)
        self.assertEqual(row['incorrect_maps'],0)  # Fixed seed regression, not a probability proof.
        self.assertGreater(experiment(1)['incorrect_maps'],0)

    def test_06_correlated_repetition_counterexample(self):
        rows = [experiment(m,correlated=True,trials=40) for m in (1,5,25)]
        self.assertEqual(len({r['incorrect_maps'] for r in rows}),1)
        self.assertGreater(rows[0]['incorrect_maps'],0)
        self.assertEqual(len({r['wrong_read_blocks_total'] for r in rows}),1)
        self.assertEqual(rows[2]['raw_marker_reads_total'],25*rows[0]['raw_marker_reads_total'])

    def test_07_adaptive_first_error_accounting(self):
        # Sum all possible positions of the first error on a tiny correct trace.
        # Probability of any first error is 1-(1-q)^Q, regardless of what a
        # corrupted continuation would do after it; Q*q is the union bound.
        probe = PortProbe([[1],[0]])
        value = map_ball(probe,1)
        blocks = value['marker_reads']
        for p in (.1,.3):
            q = majority_error(p,5)
            first_errors = math.fsum((1-q)**i*q for i in range(blocks))
            self.assertAlmostEqual(first_errors,1-(1-q)**blocks,places=14)
            self.assertLessEqual(first_errors,query_budget(1)*q)

    def test_08_domain_and_budget(self):
        for p,m in ((.5,3),(-.1,3),(.1,2),(.1,0)):
            with self.assertRaises(ValueError):
                majority_error(p,m)
        self.assertEqual(choose_repetitions(0,3,.01),1)
        for r in range(5):
            self.assertEqual(ball_upper(r),1+3*(2**r-1))
            self.assertEqual(query_budget(r),3*ball_upper(r)**2)


if __name__ == '__main__':
    main(__name__,'noisy_anchor_audit',report)
