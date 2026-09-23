"""Round 297: causal, expiring clock certificates with explicit receipt times."""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from clock_window_certificate_audit import envelope, path_enclosure


def received(records, now):
    # A record is (nominal launch, residual reading, safe availability time).
    if any(s > a for s, y, a in records):
        raise ValueError('availability must follow the nominal launch')
    return [row for row in records if row[2] <= now]


def online_interval(records, now, query, error, lipschitz, domain):
    if not domain[0] <= now <= query <= domain[1]:
        raise ValueError('decision and query must lie inside the prior domain')
    ready = received(records, now)
    if not ready or any(s < domain[0] for s, y, a in ready):
        raise ValueError('no valid received history')
    lo, hi = envelope([r[0] for r in ready], [r[1] for r in ready],
                      error, lipschitz, [query])
    return float(lo[0]), float(hi[0])


def expiry(records, now, tolerance, error, lipschitz, domain):
    lo, hi = online_interval(records, now, now, error, lipschitz, domain)
    worst = max(abs(lo), abs(hi))
    if worst > tolerance + 1e-12:
        return None
    horizon = math.inf if lipschitz == 0 else max(0., (tolerance-worst)/lipschitz)
    return min(domain[1], now+horizon)


def period_budget(channels, work, capacity, fraction, delay, horizon,
                  observed_cap, error, lipschitz, tolerance):
    if min(channels, work, capacity, fraction, lipschitz) <= 0 or fraction > 1:
        raise ValueError('positive service and drift parameters required')
    minimum = channels*work/(fraction*capacity)
    maximum = (tolerance-observed_cap-error)/lipschitz-delay-horizon
    return {'minimum_period_by_work':minimum, 'maximum_period_by_certificate':maximum,
            'this_sufficient_policy_has_feasible_period':minimum <= maximum}


def example_records():
    sites = np.arange(0., 21., 4.)
    return [(float(s), float(.005*np.sin(.2*s)), float(s+1.1)) for s in sites]


def report():
    records = example_records()
    now, horizon, error, L, tolerance = 10., 2., .0002, .001, .02
    lo, hi = online_interval(records, now, now+horizon, error, L, (0., 30.))
    last = max(s for s, y, a in received(records, now))
    budget = period_budget(6, 2, 12, .25, 1.1, horizon, .005, error, L, tolerance)
    return {'round':297,
            'scope':'Causal prediction for a prescribed clock under supplied drift, error, delivery and service bounds; not a derivation of physical time or a universal impossibility test.',
            'inputs':{'lipschitz':L, 'sample_error':error, 'domain':[0.,30.],
                      'reference_shift':1., 'tolerance':tolerance, 'decision_time':now,
                      'prediction_horizon':horizon, 'period':4., 'delivery_bound':1.1},
            'records':[list(r) for r in records],
            'received_at_decision':[list(r) for r in received(records, now)],
            'latest_received_launch':last, 'age_at_decision':now-last,
            'future_residual_interval':[lo,hi], 'actual_future_residual':float(.005*np.sin(.2*(now+horizon))),
            'expiry_without_new_data':expiry(records, now, tolerance, error, L, (0.,30.)),
            'period_budget':budget,
            'uniform_policy_bound':.005+error+L*(4.+1.1+horizon),
            'seven_unit_window_policy_bound':.005+error+L*(4.+1.1+7.),
            'six_hop_certificate':list(path_enclosure([0,1,0,1,0,1,0], 10.,
                 {(0,1):1.,(1,0):1.}, {(0,1):.0173,(1,0):.0173}, (10.,17.))),
            'resource_ledger':{'channels':6, 'work_per_update':2, 'service_per_unit':12,
                               'reserved_fraction':.25, 'update_work_rate':3.,
                               'delivery_bound_must_be_established_separately':True},
            'limits':['Receipt filtering uses fixed candidate coordinates and conservative availability times, not a newly derived absolute time.',
                      'The Lipschitz bound must extend over the future certification domain and the actual probe policy.',
                      'An empty sufficient period interval rejects this certificate/policy combination, not all synchronization schemes.',
                      'More past data cannot remove future drift uncertainty without stronger assumptions.']}


class Audit(unittest.TestCase):
    def test_01_unreceived_records_cannot_change_decisions(self):
        rows = [(0.,0.,1.),(8.,.003,9.)]
        a = online_interval(rows,10.,12.,.001,.01,(0.,30.))
        self.assertEqual(a, online_interval(rows+[(9.,900.,11.)],10.,12.,.001,.01,(0.,30.)))
        with self.assertRaises(ValueError):
            online_interval(rows,0.,2.,.001,.01,(0.,30.))

    def test_02_future_envelopes_cover_bounded_noise_and_drift(self):
        for seed in range(6):
            rng = np.random.default_rng(seed)
            sites = np.arange(21.)
            rows = [(float(s), float(.005*np.sin(.2*s)+rng.uniform(-.0002,.0002)),
                     float(s+rng.uniform(1.1,2.))) for s in sites]
            for now in (3.,8.,14.):
                for q in np.linspace(now,now+5,17):
                    lo,hi = online_interval(rows,now,float(q),.0002,.001,(0.,30.))
                    self.assertLessEqual(lo-1e-12,.005*np.sin(.2*q))
                    self.assertGreaterEqual(hi+1e-12,.005*np.sin(.2*q))

    def test_03_expiry_is_sharp_for_indistinguishable_extensions(self):
        rows = [(0.,0.,1.1),(4.,0.,5.1),(8.,0.,9.1)]
        end = expiry(rows,10.,.02,.0002,.001,(0.,40.))
        self.assertAlmostEqual(end,27.8)
        # r_+/- = +/-[b+L max(0,z-8)] both fit every received reading.
        for q in (10.,end,end+.1):
            lo,hi = online_interval(rows,10.,q,.0002,.001,(0.,40.))
            self.assertAlmostEqual(hi,.0002+.001*(q-8.))
            self.assertAlmostEqual(lo,-hi)
        self.assertGreater(hi,.02)

    def test_04_out_of_order_receipt_does_not_reset_age(self):
        rows = [(8.,0.,9.),(0.,0.,10.)]
        self.assertEqual(max(r[0] for r in received(rows,10.)),8.)
        self.assertEqual(online_interval(rows,10.,12.,0.,.001,(0.,20.)),(-.004,.004))

    def test_05_domain_and_zero_drift_and_expired_cases(self):
        rows = [(0.,0.,1.)]
        self.assertEqual(expiry(rows,2.,.1,.01,0.,(0.,10.)),10.)
        self.assertIsNone(expiry(rows,2.,.001,.01,.01,(0.,10.)))
        with self.assertRaises(ValueError):
            online_interval(rows,2.,11.,.01,.01,(0.,10.))

    def test_06_period_age_bound_including_receipt_delay(self):
        rows = example_records()
        for s,y,a in rows:
            self.assertGreaterEqual(a,s+1.+.005*np.sin(.2*s)+.05)
        for now in np.linspace(1.1,21.099,1001):
            age = now-max(s for s,y,a in received(rows,float(now)))
            self.assertLessEqual(age,5.1+1e-12)
            lo,hi = online_interval(rows,float(now),float(now+2.),.0002,.001,(0.,30.))
            self.assertLessEqual(max(abs(lo),abs(hi)),.0123+1e-12)

    def test_07_service_and_precision_constraints(self):
        good = period_budget(6,2,12,.25,1.1,2,.005,.0002,.001,.02)
        self.assertAlmostEqual(good['minimum_period_by_work'],4.)
        self.assertAlmostEqual(good['maximum_period_by_certificate'],11.7)
        self.assertTrue(good['this_sufficient_policy_has_feasible_period'])
        bad = period_budget(6,2,1,.25,1.1,2,.005,.0002,.001,.02)
        self.assertFalse(bad['this_sufficient_policy_has_feasible_period'])

    def test_08_multihop_use_checks_actual_possible_entry_window(self):
        tau = {(0,1):1.,(1,0):1.}
        cert = {e:.005+.0002+.001*(4.+1.1+7.) for e in tau}
        path = [0,1,0,1,0,1,0]
        interval = path_enclosure(path,10.,tau,cert,(10.,17.))
        arrival = 10.
        for _ in path[1:]:
            arrival += 1.+.005*np.sin(.2*arrival)
        self.assertTrue(interval[0] <= arrival <= interval[1])
        with self.assertRaises(ValueError):
            path_enclosure(path,10.,tau,cert,(10.,12.))

    def test_09_old_zero_history_does_not_determine_next_event(self):
        # Without a common supplied L, an arbitrarily steep nonnegative future
        # ramp keeps T(z)=z+1+r(z) causal and increasing while matching all past.
        for speed in (1.,10.,100.):
            r = lambda z: speed*max(z-8.,0.)
            self.assertTrue(all(r(s)==0 for s in (0.,4.,8.)))
            self.assertAlmostEqual(r(9.),speed)


if __name__ == '__main__':
    main(__name__,'online_clock_certificate_audit',report)
