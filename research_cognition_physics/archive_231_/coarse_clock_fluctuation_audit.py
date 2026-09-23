"""Round 298: when nonuniform tick intervals admit a stable coarse rate.
Event counts index observations; all intervals are relative to a supplied reference.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def mean_variance(n, sigma=.2, rho=0.):
    if n < 1 or not -1 <= rho <= 1:
        raise ValueError('positive block size and valid correlation required')
    lag = np.arange(1,n,dtype=float)
    return float(sigma**2/n*(1+2*np.sum((1-lag/n)*rho**lag)))


def enumerate_markov(n, rho, mu=1., sigma=.2):
    means, probabilities = [], []
    for signs in itertools.product((-1.,1.),repeat=n):
        signs = np.array(signs)
        weights = (1+rho*signs[1:]*signs[:-1])/2
        probabilities.append(float(.5*np.prod(weights)))
        means.append(float(mu+sigma*np.mean(signs)))
    return np.array(means), np.array(probabilities)


def sample_means(n, paths=4096, rho=0., sigma=.2, seed=298):
    rng = np.random.default_rng(seed)
    state = rng.choice([-1.,1.],size=paths)
    total = state.copy()
    for _ in range(1,n):
        state *= np.where(rng.random(paths) < (1+rho)/2,1.,-1.)
        total += state
    return 1.+sigma*total/n


def alternating_blocks(blocks=12, mu=1., sigma=.2):
    count, signed_total = 0, 0
    rows=[]
    for k in range(blocks):
        length = 1 if k == 0 else 10*count
        sign = 1 if k%2 == 0 else -1
        count += length
        signed_total += sign*length
        rows.append({'block':k+1, 'events':count,
                     'block_interval':mu+sigma*sign,
                     'cumulative_mean':mu+sigma*signed_total/count})
    return rows


def allan_variance_iid(n, sigma=.2):
    # Half the variance of two adjacent, disjoint block-mean differences.
    return sigma**2/n


def report():
    rows=[]
    for n in (1,16,256,4096):
        iid = mean_variance(n)
        correlated = mean_variance(n,rho=.9)
        rows.append({'events_per_block':n, 'iid_rate_sd':float(np.sqrt(iid)),
                     'markov_rho_0_9_rate_sd':float(np.sqrt(correlated)),
                     'common_offset_ensemble_rate_sd':.2,
                     'iid_accumulated_time_sd':float(n*np.sqrt(iid)),
                     'two_independent_clocks_rate_difference_sd':float(np.sqrt(2*iid))})
    simulations=[]
    for rho in (0.,.9):
        values = sample_means(256,rho=rho)
        simulations.append({'rho':rho,'events_per_block':256,'paths':len(values),'seed':298,
                            'empirical_variance':float(np.var(values)),
                            'analytic_variance':mean_variance(256,rho=rho)})
    return {'round':298,
            'scope':'Conditional statistical coarse graining of operational clock comparisons; neither fundamental time fluctuations nor a shared physical time have been derived.',
            'inputs':{'mean_interval':1.,'fluctuation_amplitude':.2,
                      'all_tick_intervals_positive':True,
                      'covariance_law_and_reference_readout_are_supplied':True},
            'block_statistics':rows, 'monte_carlo_cross_checks':simulations,
            'bounded_nonstationary_counterexample':alternating_blocks(),
            'stationary_common_offset':{'possible_single_history_rates':[.8,1.2],
                                        'ensemble_variance_about_prescribed_rate':.04,
                                        'allan_variance_within_each_history':0.,
                                        'removable_by_run_specific_rate_calibration':True},
            'iid_next_tick':{'conditional_variance_even_after_arbitrarily_long_past':.04,
                             'possible_intervals':[.8,1.2]},
            'limits':['Averaging a rate is different from bounding absolute accumulated time error or predicting the next event.',
                      'A common offset is a calibration ambiguity, not proof that a single history lacks a uniform rate.',
                      'The growing alternating-block example has no asymptotic constant rate despite bounded positive intervals.',
                      'These are classical statistical comparisons; clock noise cannot by itself identify microscopic spacetime ontology.']}


class Audit(unittest.TestCase):
    def test_01_covariance_formula_matches_exact_markov_enumeration(self):
        for n in range(1,9):
            for rho in (-.7,0.,.6,.9):
                means,prob = enumerate_markov(n,rho)
                self.assertAlmostEqual(float(sum(prob)),1.)
                self.assertAlmostEqual(float(prob@means),1.)
                self.assertAlmostEqual(float(prob@(means-1.)**2),mean_variance(n,rho=rho))

    def test_02_summable_correlations_give_explicit_mean_square_bound(self):
        for rho in (-.8,0.,.5,.9):
            covariance_sum=.04*(1+abs(rho))/(1-abs(rho))
            for n in (1,2,17,256,4096):
                self.assertLessEqual(mean_variance(n,rho=rho),covariance_sum/n+1e-12)
                self.assertGreaterEqual(mean_variance(n,rho=rho),-1e-12)

    def test_03_empirical_covariance_check_is_not_the_proof(self):
        for rho in (0.,.9):
            empirical=np.var(sample_means(256,rho=rho))
            exact=mean_variance(256,rho=rho)
            self.assertLess(abs(empirical/exact-1),.09)

    def test_04_common_mode_does_not_average_to_prescribed_universal_rate(self):
        for n in (1,16,256):
            means,prob=enumerate_markov(min(n,8),1.)
            self.assertAlmostEqual(float(prob@(means-1.)**2),.04)
            self.assertAlmostEqual(mean_variance(n,rho=1.),.04)

    def test_05_bounded_nonstationary_intervals_have_two_subsequence_limits(self):
        rows=alternating_blocks(20)
        self.assertTrue(all(.8 <= r['block_interval'] <= 1.2 for r in rows))
        self.assertAlmostEqual(rows[-2]['cumulative_mean'],1.+1/6,places=12)
        self.assertAlmostEqual(rows[-1]['cumulative_mean'],1.-1/6,places=12)
        self.assertGreater(rows[-2]['cumulative_mean']-rows[-1]['cumulative_mean'],.33)

    def test_06_absolute_error_and_future_event_do_not_vanish_with_rate_error(self):
        for n in (1,16,256,4096):
            self.assertAlmostEqual(n*n*mean_variance(n),.04*n)
        # Exact factorization of independent signs leaves the next sign unbiased.
        histories=list(itertools.product((-1.,1.),repeat=5))
        for past in itertools.product((-1.,1.),repeat=4):
            next_intervals=np.array([1.+.2*h[-1] for h in histories if h[:4]==past])
            self.assertEqual(len(next_intervals),2)
            self.assertAlmostEqual(float(np.var(next_intervals)),.04)

    def test_07_allan_difference_separates_noise_from_constant_offset(self):
        for n in (1,2,3,4):
            differences=[]
            for signs in itertools.product((-1.,1.),repeat=2*n):
                a=.2*np.mean(signs[:n]); b=.2*np.mean(signs[n:])
                differences.append((b-a)**2/2)
            self.assertAlmostEqual(float(np.mean(differences)),allan_variance_iid(n))
        for offset in (-.2,.2):
            self.assertEqual(((1.+offset)-(1.+offset))**2/2,0.)

    def test_08_two_independent_clocks_need_comparison_not_absolute_background(self):
        n=3
        values,prob=enumerate_markov(n,0.)
        variance=float(np.sum(prob[:,None]*prob[None,:]*(values[:,None]-values[None,:])**2))
        self.assertAlmostEqual(variance,2*mean_variance(n))
        # The event index is a count, not itself a supplied physical duration.
        self.assertTrue(np.all(values>0))


if __name__ == '__main__':
    main(__name__,'coarse_clock_fluctuation_audit',report)
