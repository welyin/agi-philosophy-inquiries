"""Round 299: three-cornered hats, common modes and clock/link ambiguity.
Classical aligned scalar observations; no microscopic spacetime model inferred.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def incidence(n, edges):
    b = np.zeros((len(edges),n))
    for k,(i,j) in enumerate(edges):
        b[k,i],b[k,j] = -1.,1.
    return b


def pair_variances(covariance):
    s=np.asarray(covariance,dtype=float)
    return np.diag(s)[:,None]+np.diag(s)[None,:]-2*s


def three_cornered_hat(pair):
    v=np.asarray(pair,dtype=float)
    return np.array([(v[0,1]+v[0,2]-v[1,2])/2,
                     (v[0,1]+v[1,2]-v[0,2])/2,
                     (v[0,2]+v[1,2]-v[0,1])/2])


def relative_covariance(b, observed):
    inverse=np.linalg.pinv(b)
    return inverse@observed@inverse.T


def cycle_projector(b):
    return np.eye(len(b))-b@np.linalg.pinv(b)


def report():
    base=np.diag([1.,4.,9.])
    hidden=base+2*np.ones((3,3))
    correlated=np.outer([1.,2.,-1.],[1.,2.,-1.])
    b=incidence(3,[(0,1),(1,2),(0,2)])
    cycle=np.array([1.,1.,-1.])
    link=np.array([.1,.2,.8])
    p=np.eye(3)-np.ones((3,3))/3
    return {'round':299,
            'scope':'Identifiability of aligned classical clock comparisons under explicit covariance and link models; no identification of microscopic time fluctuations.',
            'three_cornered_hat':{'independent_true_variances':np.diag(base).tolist(),
                'pair_variances':pair_variances(base).tolist(),
                'same_observed_pairs_with_true_variances':np.diag(hidden).tolist(),
                'hat_for_both_models':three_cornered_hat(pair_variances(hidden)).tolist(),
                'hidden_common_variance':2.},
            'correlated_counterexample':{'true_covariance':correlated.tolist(),
                'true_variances':np.diag(correlated).tolist(),
                'incorrect_uncorrelated_hat':three_cornered_hat(pair_variances(correlated)).tolist()},
            'observable_covariance':{'incidence':b.tolist(),
                'edge_covariance':(b@hidden@b.T).tolist(),
                'recoverable_centered_covariance':relative_covariance(b,b@hidden@b.T).tolist(),
                'matches_centered_truth':bool(np.allclose(relative_covariance(b,b@hidden@b.T),p@hidden@p))},
            'link_cycle_example':{'link_residuals':link.tolist(),'cycle_coefficients':cycle.tolist(),
                'cycle_value':float(cycle@link),'per_reading_bound':.01,
                'certified_interval':[float(cycle@link-.03),float(cycle@link+.03)],
                'cycle_component':(cycle_projector(b)@link).tolist()},
            'limits':['Simultaneous or correctly aligned scalar data are a supplied interface; transport and alignment costs are not removed.',
                      'Positive three-cornered-hat outputs do not establish independence or identify a common mode.',
                      'Unknown link terms in the incidence image can be exchanged with clock terms without changing any record.',
                      'Cycle residuals reject a node-only model but can arise from links, alignment or measurement errors.']}


class Audit(unittest.TestCase):
    def test_01_hat_recovers_independent_covariances(self):
        for v in ([1.,4.,9.],[.1,.2,.3],[0.,1.,2.]):
            self.assertTrue(np.allclose(three_cornered_hat(pair_variances(np.diag(v))),v))

    def test_02_correlation_can_give_negative_hat_without_negative_true_variance(self):
        s=np.outer([1.,2.,-1.],[1.,2.,-1.])
        self.assertTrue(np.allclose(three_cornered_hat(pair_variances(s)),[-2.,3.,6.]))
        self.assertGreaterEqual(np.linalg.eigvalsh(s).min(),-1e-12)

    def test_03_positive_hat_cannot_certify_absence_of_common_noise(self):
        base=np.diag([1.,4.,9.])
        for q in (.1,2.,100.):
            s=base+q*np.ones((3,3))
            self.assertTrue(np.allclose(pair_variances(base),pair_variances(s)))
            self.assertTrue(np.all(three_cornered_hat(pair_variances(s))>0))
            self.assertFalse(np.allclose(np.diag(s),three_cornered_hat(pair_variances(s))))

    def test_04_connected_graph_recovers_only_centered_covariance(self):
        rng=np.random.default_rng(299)
        for n in range(2,9):
            edges=[(i,i+1) for i in range(n-1)]
            if n>2:
                edges.append((0,n-1))
            b=incidence(n,edges)
            a=rng.normal(size=(n,n)); s=a@a.T
            p=np.eye(n)-np.ones((n,n))/n
            self.assertEqual(np.linalg.matrix_rank(b),n-1)
            self.assertTrue(np.allclose(relative_covariance(b,b@s@b.T),p@s@p))
            a=rng.normal(size=n)
            invisible=np.ones((n,1))*a[None,:]+a[:,None]*np.ones((1,n))
            self.assertTrue(np.allclose(b@invisible@b.T,0.))

    def test_05_clock_link_reassignment_preserves_all_readings(self):
        b=incidence(3,[(0,1),(1,2),(0,2)])
        rng=np.random.default_rng(29)
        x,p,a=[rng.normal(size=(3,57)) for _ in range(3)]
        self.assertTrue(np.allclose(b@x+p,b@(x+a)+(p-b@a)))
        self.assertTrue(np.allclose(cycle_projector(b)@p,cycle_projector(b)@(p-b@a)))

    def test_06_tree_has_no_snapshot_cycle_diagnostic(self):
        b=incidence(5,[(0,1),(1,2),(1,3),(3,4)])
        self.assertTrue(np.allclose(cycle_projector(b),0.))
        link=np.array([.1,.2,-.1,.4])
        self.assertTrue(np.allclose(b@np.linalg.pinv(b)@link,link))

    def test_07_cycle_error_certificate_covers_every_error_corner(self):
        cycle=np.array([1.,1.,-1.]); link=np.array([.1,.2,.8])
        for signs in itertools.product((-1.,1.),repeat=3):
            measured=float(cycle@(link+.01*np.array(signs)))
            self.assertLessEqual(abs(measured-float(cycle@link)),.03+1e-12)
            self.assertLess(measured,0.)

    def test_08_pathwise_common_mode_cancellation_survives_all_block_sizes(self):
        b=incidence(3,[(0,1),(1,2),(0,2)])
        rng=np.random.default_rng(298)
        x=rng.normal(size=(3,128)); common=10*np.sin(np.arange(128)*.17)
        y=b@x; altered=b@(x+common)
        self.assertTrue(np.allclose(y,altered))
        for n in (1,2,8,32,128):
            self.assertTrue(np.allclose(y.reshape(3,-1,n).mean(axis=2),
                                        altered.reshape(3,-1,n).mean(axis=2)))

    def test_09_exact_bounded_sign_ensemble_realizes_the_covariance_ambiguity(self):
        signs=np.array(list(itertools.product((-1.,1.),repeat=4)))
        x=signs[:,:3]*np.array([1.,2.,3.])
        shifted=x+np.sqrt(2)*signs[:,3,None]
        b=incidence(3,[(0,1),(1,2),(0,2)])
        self.assertTrue(np.allclose(x@b.T,shifted@b.T))
        self.assertTrue(np.allclose(x.T@x/len(x),np.diag([1.,4.,9.])))
        self.assertTrue(np.allclose(shifted.T@shifted/len(x),np.diag([1.,4.,9.])+2))
        self.assertTrue(np.all(10+shifted>0))


if __name__ == '__main__':
    main(__name__,'clock_source_identifiability_audit',report)
