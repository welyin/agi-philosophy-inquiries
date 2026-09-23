"""Round 295: finite robust obstructions and a hidden-gap countermodel.
Classical local arrival maps; no all-time certification from finite queries.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main
from local_clock_gauge_audit import loop_f, loop_g
from interval_density_audit import crossing_pair


def contrast_interval(y_first,y_second,k_first,k_second,launch_error,output_error):
    """Interval for second(x)-first(x), with separately perturbed launches."""
    if min(k_first,k_second,launch_error,output_error)<0:
        raise ValueError('nonnegative error and Lipschitz bounds required')
    gap=float(y_second-y_first)
    radius=2*output_error+(k_first+k_second)*launch_error
    return gap-radius,gap+radius


def certified_sign(interval):
    lo,hi=interval
    if lo>1e-12:
        return 1
    if hi<-1e-12:
        return -1
    return 0


def perturbed_observations(first,second,x,eta,epsilon):
    for d1,d2,e1,e2 in itertools.product((-1.,1.),repeat=4):
        yield first(x+d1*eta)+e1*epsilon,second(x+d2*eta)+e2*epsilon


def adaptive_transcript(g,steps=18):
    """A finite policy choosing the next input after seeing previous returns."""
    records=[]
    x=0.
    for k in range(steps):
        name='G' if k%3 else 'F'
        f=g if name=='G' else lambda z:z+1.
        y=float(f(x))
        records.append({'step':k,'map':name,'input':x,'output':y})
        x=float((.37*y+.11*(k%4))%1.)
    return records


def empty_gap(records):
    sites=sorted(set([0.,1.]+[r['input'] for r in records]))
    left,right=max(zip(sites,sites[1:]),key=lambda p:p[1]-p[0])
    center=(left+right)/2
    halfwidth=(right-left)/3
    return center,halfwidth


def bump(x,center,halfwidth):
    z=(np.asarray(x)-center)/halfwidth
    return np.maximum(1-z*z,0.)**3


def bump_prime(x,center,halfwidth):
    z=(np.asarray(x)-center)/halfwidth
    return np.where(abs(z)<1,-6*z*np.maximum(1-z*z,0.)**2/halfwidth,0.)


def bump_second(x,center,halfwidth):
    z=(np.asarray(x)-center)/halfwidth
    p=np.maximum(1-z*z,0.)
    return np.where(abs(z)<1,(-6*p*p+24*z*z*p)/halfwidth**2,0.)


def hidden_model(records):
    center,halfwidth=empty_gap(records)
    amplitude=halfwidth/12
    altered=lambda x:np.asarray(x)+1.+amplitude*bump(x,center,halfwidth)
    return altered,center,halfwidth,amplitude


def report():
    fg=lambda x:loop_g(loop_f(x))
    gf=lambda x:loop_f(loop_g(x))
    ci=contrast_interval(gf(0.),fg(0.),1.2,1.2,.01,.005)
    f,g,_=crossing_pair()
    cross=[]
    for x in (.25,.75):
        interval=contrast_interval(f(x),g(x),1.,1+.2*np.pi,.005,.005)
        cross.append({'launch':x,'observed_gap':float(g(x)-f(x)),
                      'certified_gap_interval':list(interval),
                      'sign':certified_sign(interval)})
    records=adaptive_transcript(lambda x:x+1.)
    altered,c,h,a=hidden_model(records)
    other=adaptive_transcript(altered)
    return {'round':295,
            'scope':'Finite robust rejection witnesses and non-identifiability of exact global constant-delay structure; bounded input/output errors and map regularity are explicit assumptions.',
            'noncommutation_witness':{
                'launch_uncertainty':.01,'output_uncertainty':.005,
                'word_lipschitz_bounds':[1.2,1.2],
                'observed_gap':float(fg(0.)-gf(0.)),
                'certified_gap_interval':list(ci),'sign':certified_sign(ci),
                'actual_edge_transmissions':8},
            'crossing_witness':cross,
            'finite_query_countermodel':{
                'query_count':len(records),'transcript':records,
                'same_adaptive_transcript':records==other,
                'hidden_support':[c-h,c+h], 'center':c,'amplitude':a,
                'global_derivative_conservative_bounds':[.5,1.5],
                'hidden_commutator_gap':float(altered(c+1.)-(altered(c)+1.)),
                'first_derivatives_at_sample_inputs':[
                    float(1+a*bump_prime(r['input'],c,h)) for r in records],
                'perturbation_is_C2_and_strictly_future_directed':True},
            'limits':['A finite obstruction disproves simultaneous constant-delay clocks for the tested reproducible model, not time coordinates in general.',
                      'A finite data fit does not prove an all-time functional identity; even bounded first derivatives leave unsampled countermodels.',
                      'Repeated reads do not reduce an adversarial deterministic error bound without a new statistical assumption.',
                      'Both compared experiments must share a specified nominal local launch event; launch jitter and probe interference must be accounted for.']}


class Audit(unittest.TestCase):
    def test_01_input_and_output_error_interval_covers_true_contrast(self):
        f,g,_=crossing_pair()
        for x in np.linspace(-1,1,21):
            truth=float(g(x)-f(x))
            for a,b in perturbed_observations(f,g,x,.005,.005):
                lo,hi=contrast_interval(a,b,1.,1+.2*np.pi,.005,.005)
                self.assertLessEqual(lo,truth+1e-12)
                self.assertGreaterEqual(hi,truth-1e-12)

    def test_02_noncommutation_witness_survives_all_error_corners(self):
        first=lambda x:loop_f(loop_g(x))
        second=lambda x:loop_g(loop_f(x))
        for a,b in perturbed_observations(first,second,0.,.01,.005):
            self.assertEqual(certified_sign(contrast_interval(a,b,1.2,1.2,.01,.005)),1)

    def test_03_two_start_crossing_is_robust(self):
        f,g,_=crossing_pair()
        for x,sign in ((.25,1),(.75,-1)):
            for a,b in perturbed_observations(f,g,x,.005,.005):
                interval=contrast_interval(a,b,1.,1+.2*np.pi,.005,.005)
                self.assertEqual(certified_sign(interval),sign)

    def test_04_exactly_commuting_model_never_falsely_rejected(self):
        first=lambda x:x+3.
        second=lambda x:x+3.
        for a,b in perturbed_observations(first,second,0.,.01,.005):
            self.assertEqual(certified_sign(contrast_interval(a,b,1.,1.,.01,.005)),0)
        self.assertEqual(certified_sign(contrast_interval(3.,3.1,1.,1.,.1,.1)),0)

    def test_05_hidden_bump_is_smooth_enough_and_causal(self):
        records=adaptive_transcript(lambda x:x+1.)
        altered,c,h,a=hidden_model(records)
        grid=np.linspace(c-2*h,c+2*h,1001)
        derivative=1+a*bump_prime(grid,c,h)
        self.assertTrue(np.all(derivative>=.5))
        self.assertTrue(np.all(derivative<=1.5))
        self.assertTrue(np.all(altered(grid)>grid))
        for x in (c-h,c+h):
            self.assertAlmostEqual(float(bump(x,c,h)),0.,places=12)
            self.assertAlmostEqual(float(bump_prime(x,c,h)),0.,places=9)
            self.assertAlmostEqual(float(bump_second(x,c,h)),0.,places=6)

    def test_06_full_adaptive_transcript_is_identical(self):
        for n in (6,18,40,100):
            records=adaptive_transcript(lambda x:x+1.,steps=n)
            altered,_,_,_=hidden_model(records)
            self.assertEqual(records,adaptive_transcript(altered,steps=n))

    def test_07_hidden_model_has_a_real_translation_obstruction(self):
        records=adaptive_transcript(lambda x:x+1.)
        altered,c,h,a=hidden_model(records)
        self.assertAlmostEqual(float(altered(c+1.)-(altered(c)+1.)),-a,places=12)
        self.assertGreater(a,0.)
        self.assertLess(h,.25)

    def test_08_local_derivative_samples_also_miss_the_hidden_gap(self):
        records=adaptive_transcript(lambda x:x+1.)
        altered,c,h,a=hidden_model(records)
        for row in records:
            x=row['input']
            self.assertGreater(abs(x-c),h)
            self.assertEqual(float(altered(x)),x+1.)
            self.assertEqual(float(1+a*bump_prime(x,c,h)),1.)
            self.assertEqual(float(a*bump_second(x,c,h)),0.)


if __name__ == '__main__':
    main(__name__,'finite_clock_witness_audit',report)

