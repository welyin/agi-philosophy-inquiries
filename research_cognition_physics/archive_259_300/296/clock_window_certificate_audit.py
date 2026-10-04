"""Round 296: bounded-data certificates for approximate clock translations.
Uses standard Lipschitz/set-membership envelopes on a stated finite domain.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main
from temporal_distance_audit import path_arrival, echo_matrix, reference_metric


def residual_error_bound(lipschitz,launch_error,read_error):
    return (1+lipschitz)*launch_error+read_error


def cover_radius(sites,domain):
    points=np.sort(np.asarray(sites,dtype=float))
    a,b=domain
    if not len(points) or a>b or points[0]<a or points[-1]>b:
        raise ValueError('sample sites must be inside the domain')
    gaps=[points[0]-a,b-points[-1]]
    if len(points)>1:
        gaps.extend(np.diff(points)/2)
    return float(max(gaps))


def envelope(sites,values,error,lipschitz,queries):
    sites=np.asarray(sites,dtype=float)
    values=np.asarray(values,dtype=float)
    if len(sites)!=len(values) or min(error,lipschitz)<0:
        raise ValueError('invalid samples or bounds')
    discrepancies=abs(values[:,None]-values[None,:])
    admissible=2*error+lipschitz*abs(sites[:,None]-sites[None,:])
    if np.any(discrepancies>admissible+1e-12):
        raise ValueError('samples inconsistent with supplied noise/regularity bounds')
    distances=abs(np.atleast_1d(queries)[:,None]-sites[None,:])
    lower=np.max(values[None,:]-error-lipschitz*distances,axis=1)
    upper=np.min(values[None,:]+error+lipschitz*distances,axis=1)
    if np.any(lower>upper+1e-12):
        raise ValueError('empty feasible envelope')
    return lower,upper


def uniform_bound(sites,values,error,lipschitz,domain):
    envelope(sites,values,error,lipschitz,sites)
    return float(max(abs(np.asarray(values)))+error+lipschitz*cover_radius(sites,domain))


def fixture(samples=17,seed=0):
    rng=np.random.default_rng(seed)
    tau={(0,1):.6,(1,0):.4,(1,2):.8,(2,1):.7,(0,2):1.5,(2,0):1.2}
    amplitude,omega=.005,.2
    lipschitz=amplitude*omega
    eta=read_error=.0001
    domain=(0.,20.)
    sites=np.linspace(*domain,samples)
    arcs,observed,certificates={},{},{}
    error=residual_error_bound(lipschitz,eta,read_error)
    for index,(edge,d) in enumerate(tau.items()):
        phase=float(rng.uniform(0,2*np.pi))
        arcs[edge]=lambda z,d=d,p=phase:z+d+amplitude*np.sin(omega*z+p)
        perturb=(-1.)**(np.arange(samples)+index)
        launches=sites+eta*perturb
        returns=arcs[edge](launches)+read_error*np.roll(perturb,1)
        values=returns-sites-d
        observed[edge]=values
        certificates[edge]=uniform_bound(sites,values,error,lipschitz,domain)
    return tau,arcs,sites,observed,certificates,error,lipschitz,domain


def path_enclosure(path,start,tau,certificates,domain):
    low=high=float(start)
    a,b=domain
    for edge in zip(path,path[1:]):
        if low<a-1e-12 or high>b+1e-12:
            raise ValueError('whole possible entry interval must stay inside verified domain')
        low+=tau[edge]-certificates[edge]
        high+=tau[edge]+certificates[edge]
    return low,high


def report():
    tau,arcs,sites,obs,cert,error,L,domain=fixture()
    epsilon=max(cert[e]/tau[e] for e in tau)
    ref=reference_metric(3,tau)
    actual=echo_matrix(3,arcs,10.)
    path=[0,1,2,0,1,2,0]
    low,high=path_enclosure(path,2.,tau,cert,domain)
    exact=path_arrival(arcs,path,2.)
    sample_max=max(float(max(abs(v))) for v in obs.values())
    return {'round':296,
            'scope':'Conditional finite-window approximate translation certificates from noisy samples plus supplied regularity bounds and candidate clocks; no exact all-time common clock, natural dimension or gravitational dynamics derived.',
            'inputs':{'domain':list(domain),'directed_channels':len(tau),
                      'samples_per_channel':len(sites),'launch_error':.0001,
                      'return_reading_error':.0001,'residual_lipschitz_bound':L,
                      'effective_sample_error':error,'covering_radius':cover_radius(sites,domain),
                      'candidate_clocks_and_units_fixed_before_validation':True},
            'maximum_observed_residual':sample_max,
            'edge_certificates':[{'edge':list(e),'reference_shift':tau[e],
                                  'uniform_residual_bound':cert[e]} for e in tau],
            'relative_delay_bound':epsilon,
            'verified_metric_example':{'departure':10.,'reference_metric':ref.tolist(),
                                       'actual_half_roundtrips':actual.tolist(),
                                       'conservative_roundtrip_horizon':2*(1+epsilon)*2*max(tau.values())},
            'six_hop_path':{'path':path,'start':2.,
                            'predicted_arrival':2.+sum(tau[e] for e in zip(path,path[1:])),
                            'certified_arrival_interval':[low,high],
                            'actual_arrival':exact},
            'probe_ledger':{'one_way_samples':len(tau)*len(sites),
                            'endpoint_timestamp_samples':2*len(tau)*len(sites),
                            'maximum_sample_launch':float(sites[-1]),
                            'all_returns_may_arrive_after_last_launch':True,
                            'data_exchange_storage_and_clock_energy_not_included':True},
            'limits':['Regularity and absolute error bounds must be justified independently; consistency with samples is not proof that they hold.',
                      'The certificate concerns a prescribed candidate clock, not optimality or uniqueness among all monotone clocks.',
                      'Window coverage is retrospective unless stationarity or a predictive bound extends it to future measurements.',
                      'Path errors add per edge; long paths can exceed tolerance or exit the verified domain.',
                      'Uniform rescaling changes absolute errors and reference shifts together, leaving the relative certificate unchanged.']}


class Audit(unittest.TestCase):
    def test_01_set_membership_envelopes_cover_all_noise_corners(self):
        sites=np.linspace(0,2,4)
        queries=np.linspace(0,2,101)
        L=.02
        eta=epsilon=.001
        b=residual_error_bound(L,eta,epsilon)
        for signs in itertools.product((-1.,1.),repeat=8):
            actual_sites=sites+eta*np.array(signs[:4])
            # A known L-Lipschitz residual, not used by the envelope itself.
            y=actual_sites+1.+.02*np.sin(actual_sites)+epsilon*np.array(signs[4:])
            residuals=y-sites-1.
            lower,upper=envelope(sites,residuals,b,L,queries)
            truth=.02*np.sin(queries)
            self.assertTrue(np.all(lower<=truth+1e-12))
            self.assertTrue(np.all(upper>=truth-1e-12))

    def test_02_inconsistent_prior_bounds_are_rejected(self):
        with self.assertRaises(ValueError):
            envelope([0.,1.],[0.,1.],.01,.1,[.5])

    def test_03_cover_radius_counts_boundaries_and_internal_gaps(self):
        self.assertEqual(cover_radius([0.,1.,2.],(0.,2.)),.5)
        self.assertEqual(cover_radius([.2,.4,.8],(0.,1.)),.2)
        self.assertEqual(cover_radius([.5],(0.,1.)),.5)

    def test_04_dense_model_evaluation_is_inside_proved_uniform_bounds(self):
        for seed in range(6):
            tau,arcs,sites,obs,cert,error,L,domain=fixture(seed=seed)
            grid=np.linspace(*domain,2001)
            for edge in tau:
                truth=arcs[edge](grid)-grid-tau[edge]
                self.assertLessEqual(float(max(abs(truth))),cert[edge]+1e-12)
                low,high=envelope(sites,obs[edge],error,L,grid)
                self.assertTrue(np.all(low<=truth+1e-12))
                self.assertTrue(np.all(high>=truth-1e-12))

    def test_05_path_enclosures_contain_actual_compositions(self):
        tau,arcs,_,_,cert,_,_,domain=fixture()
        for path in ([0,1,2],[0,2,1,0],[0,1,2,0,1,2,0]):
            for start in (0.,2.,10.):
                low,high=path_enclosure(path,start,tau,cert,domain)
                end=path_arrival(arcs,path,start)
                self.assertLessEqual(low,end+1e-12)
                self.assertGreaterEqual(high,end-1e-12)

    def test_06_data_certificate_transfers_to_round292_metric_bound(self):
        for seed in range(6):
            tau,arcs,_,_,cert,_,_,domain=fixture(seed=seed)
            epsilon=max(cert[e]/tau[e] for e in tau)
            ref=reference_metric(3,tau)
            self.assertLess(epsilon,1.)
            for start in (0.,2.,10.):
                self.assertLessEqual(start+4*(1+epsilon)*max(tau.values()),domain[1])
                actual=echo_matrix(3,arcs,start)
                self.assertTrue(np.all(actual>=(1-epsilon)*ref-1e-12))
                self.assertTrue(np.all(actual<=(1+epsilon)*ref+1e-12))

    def test_07_per_step_error_can_accumulate_linearly(self):
        tau={(0,1):1.,(1,0):1.}
        certificates={e:.02 for e in tau}
        arcs={e:lambda z:z+1.02 for e in tau}
        for hops in (2,4,8):
            path=[k%2 for k in range(hops+1)]
            low,high=path_enclosure(path,0.,tau,certificates,(0.,20.))
            end=path_arrival(arcs,path,0.)
            self.assertAlmostEqual(end,high)
            self.assertAlmostEqual(end-hops,.02*hops)

    def test_08_out_of_window_path_is_not_certified(self):
        tau,_,_,_,cert,_,_,domain=fixture()
        with self.assertRaises(ValueError):
            path_enclosure([0,1,2,0],19.5,tau,cert,domain)

    def test_09_unit_rescaling_cannot_fake_relative_accuracy(self):
        tau,_,sites,obs,cert,error,L,domain=fixture()
        alpha=.001
        for edge in tau:
            scaled=uniform_bound(alpha*sites,alpha*obs[edge],
                                 alpha*error,L,tuple(alpha*np.array(domain)))
            self.assertAlmostEqual(scaled/(alpha*tau[edge]),cert[edge]/tau[edge],places=12)


if __name__ == '__main__':
    main(__name__,'clock_window_certificate_audit',report)

