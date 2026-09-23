"""Round 291: baseline delay, queue observability and probe load.
Classical finite network measurements only; packet type and clock rates fixed.
"""
import itertools
import math
import unittest
import numpy as np
from growing_stream_audit import main


def fcfs(arrivals, service=1., propagation=2., available=0.):
    """One nonpreemptive FIFO server; timestamps taken at queue admission."""
    times = np.asarray(arrivals, dtype=float)
    if service <= 0 or propagation < 0 or np.any(np.diff(times) < 0):
        raise ValueError('positive service and ordered admissions required')
    rows = []
    for arrival in times:
        start = max(float(arrival), available)
        finish = start + service
        rows.append({'admit':float(arrival), 'start':start, 'finish':finish,
                     'receive':finish+propagation, 'queue_wait':start-arrival,
                     'delay':finish+propagation-arrival})
        available = finish
    return rows


def held_echo(hold):
    """Forward delay 1; reverse server busy until 10, service 1, flight 1."""
    if hold < 0:
        raise ValueError('nonnegative hold required')
    row = fcfs([1.+hold], service=1., propagation=1., available=10.)[0]
    return {'hold':hold, 'raw_roundtrip':row['receive'],
            'corrected_roundtrip':row['receive']-hold,
            'unobserved_reverse_queue_wait':row['queue_wait']}


def baseline_interval(observed, queue_threshold, error, drift=0.):
    """Conditional on at least one near-empty sample and uniform drift bound."""
    if min(queue_threshold,error,drift) < 0 or not len(observed):
        raise ValueError('nonempty data and nonnegative bounds required')
    minimum = float(min(observed))
    return [max(0., minimum-queue_threshold-error-drift),
            minimum+error+drift]


def samples_per_link(links, p, failure):
    if links < 1 or not 0 < p <= 1 or not 0 < failure < 1:
        raise ValueError('invalid probability or link count')
    return 1 if p == 1 else math.ceil(math.log(links/failure)/(-math.log1p(-p)))


def budget(links=6, p=.25, failure=.01, service_work=2.,
           capacity=1., fraction=.1, drift_rate=.0001,
           target_error=.1, queue_threshold=.02, error=.005):
    n = samples_per_link(links,p,failure)
    work = links*n*service_work
    lower = work/(fraction*capacity)
    upper = (target_error-queue_threshold-error)/drift_rate
    return {'links':links, 'samples_per_link':n, 'total_echoes':links*n,
            'service_work':work, 'measurement_fraction':fraction,
            'minimum_window_from_work':lower, 'maximum_window_from_drift':upper,
            'certificate_intervals_overlap':bool(lower<=upper),
            'union_failure_bound':links*(1-p)**n,
            'error_bound_at_minimum_window':queue_threshold+error+drift_rate*lower}


def conditional_failure(n, p):
    """Exact binary history enumeration; success chance depends on the past."""
    failure = total = 0.
    for bits in itertools.product((0,1), repeat=n):
        probability = 1.
        for k,bit in enumerate(bits):
            chance = p if sum(bits[:k]) % 2 == 0 else (1+p)/2
            probability *= chance if bit else 1-chance
        total += probability
        if not any(bits):
            failure += probability
    return total,failure


def report():
    measurements = np.array([9.,7.,8.,7.])
    alternatives = [{'baseline':b, 'unseen_queue':(measurements-b).tolist()}
                    for b in (3.,6.,7.)]
    streams=[]
    for gap in (.5,1.,2.):
        rows=fcfs(np.arange(8)*gap)
        streams.append({'admission_gap':gap, 'delays':[r['delay'] for r in rows],
                        'total_service_work':8.,
                        'observation_end':rows[-1]['receive']})
    total,failure=conditional_failure(8,.25)
    return {'round':291,
            'scope':'Conditional baseline-delay inference with finite samples, endogenous probe queues and drift budgets; no spacetime curvature, quantum backaction or energy closure derived.',
            'same_observations':measurements.tolist(),
            'indistinguishable_decompositions':alternatives,
            'self_loaded_probe_streams':streams,
            'held_echoes':[held_echo(w) for w in (0.,4.,9.,20.)],
            'bounded_noise_example':{
                'true_baseline':6., 'observed':[7.95,6.07,6.46],
                'queue_threshold':.1, 'absolute_error_bound':.05,
                'conditional_baseline_interval':baseline_interval([7.95,6.07,6.46],.1,.05)},
            'conditional_probability_example':{
                'samples':8, 'per_history_near_empty_lower_bound':.25,
                'enumerated_total_probability':total, 'all_bad_probability':failure,
                'product_bound':.75**8,
                'common_mode_all_bad_probability_without_conditional_assumption':.75},
            'window_certificates':[budget(),budget(fraction=.5)],
            'resource_limits':['The work bound is necessary for this certificate, not a sufficient scheduling theorem.',
                               'The near-empty probability must hold under the actual probing policy; changing probe load can change it.',
                               'Unobserved queue wait cannot be subtracted from endpoint timestamps alone.',
                               'The baseline includes fixed transmission/processing for a specified probe type; it is not pure spatial flight time.']}


class Audit(unittest.TestCase):
    def test_01_probe_created_queue_closed_form(self):
        for service,gap in itertools.product((.5,1.,2.),(.25,.5,1.,3.)):
            rows=fcfs(np.arange(12)*gap,service=service,propagation=2.)
            expected=service+2.+np.arange(12)*max(service-gap,0.)
            np.testing.assert_allclose([r['delay'] for r in rows],expected)
            self.assertAlmostEqual(sum(r['finish']-r['start'] for r in rows),12*service)

    def test_02_baseline_not_identifiable_even_without_noise(self):
        y=np.array([9.,7.,8.,7.])
        for b in np.linspace(.01,7.,31):
            q=y-b
            self.assertTrue(np.all(q>=0))
            np.testing.assert_allclose(b+q,y)
        self.assertGreater(min(y),6.)

    def test_03_hold_subtraction_changes_experienced_queue(self):
        rows=[held_echo(w) for w in (0.,4.,9.,20.)]
        self.assertEqual([r['corrected_roundtrip'] for r in rows],[12.,8.,3.,3.])
        for r in rows:
            self.assertEqual(r['corrected_roundtrip'],
                             3.+r['unobserved_reverse_queue_wait'])

    def test_04_minimum_interval_with_noise_and_slow_drift(self):
        for q in ((.1,2.,.7),(3.,.05,.2)):
            for errors in itertools.product((-.05,.05),repeat=3):
                for shifts in itertools.product((-.03,.03),repeat=3):
                    y=6.+np.array(q)+np.array(errors)+np.array(shifts)
                    low,high=baseline_interval(y,.1,.05,.03)
                    self.assertLessEqual(low,6.+1e-12)
                    self.assertGreaterEqual(high,6.-1e-12)
                    self.assertGreaterEqual(min(y),6.-.08-1e-12)
                    self.assertLessEqual(min(y),6.+.18+1e-12)

    def test_05_history_dependent_sampling_product_bound(self):
        for p in (.1,.25,.7):
            for n in (1,4,8):
                total,failure=conditional_failure(n,p)
                self.assertAlmostEqual(total,1.)
                self.assertLessEqual(failure,(1-p)**n+1e-14)

    def test_06_common_mode_does_not_average_away(self):
        p=.25
        # One shared Bernoulli draw fixes the queue regime for the entire run.
        for n in (2,8,100):
            all_bad=1-p
            self.assertGreater(all_bad,(1-p)**n)
        self.assertEqual(1-p,.75)

    def test_07_work_and_drift_certificates(self):
        slow,fast=budget(),budget(fraction=.5)
        self.assertEqual(slow['samples_per_link'],23)
        self.assertLessEqual(slow['union_failure_bound'],.01)
        self.assertGreater(6*.75**22,.01)
        self.assertFalse(slow['certificate_intervals_overlap'])
        self.assertTrue(fast['certificate_intervals_overlap'])
        self.assertEqual(slow['minimum_window_from_work'],2760.)

    def test_08_added_queue_telemetry_identifies_baseline(self):
        rows=fcfs(np.arange(8)*.5,available=4.)
        corrected=[r['delay']-r['queue_wait'] for r in rows]
        np.testing.assert_allclose(corrected,np.full(8,3.))
        self.assertGreater(min(r['delay'] for r in rows),3.)


if __name__ == '__main__':
    main(__name__,'probe_load_audit',report)

