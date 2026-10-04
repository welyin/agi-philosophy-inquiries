"""Round 300: reuse time-delay interferometry algebra on a two-stream model.
Retarded records are classical data; no LISA model or gravitational signal derived.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def retarded_maps(a=1., b=2., va=0., vb=0.):
    return {'A':lambda t:(1-va)*np.asarray(t)-a,
            'B':lambda t:(1-vb)*np.asarray(t)-b}


def word_time(word,t,maps,domain=None):
    value=np.asarray(t,dtype=float)
    for name in word:
        previous=value
        value=maps[name](previous)
        if domain is not None:
            if np.any(value<domain[0]) or np.any(previous>domain[1]) or np.any(value>previous):
                raise ValueError('retarded record unavailable or map not causal on requested domain')
    return value


def delay(word,f,t,maps):
    # D_A D_B f(t) = f(R_B(R_A(t))); operator order is intentionally explicit.
    return f(word_time(word,t,maps))


def streams(c,sa,sb,maps):
    return (lambda t:delay('A',c,t,maps)-c(t)+sa(t),
            lambda t:delay('B',c,t,maps)-c(t)+sb(t))


def first(y_a,y_b,t,maps):
    return delay('B',y_a,t,maps)-y_a(t)-delay('A',y_b,t,maps)+y_b(t)


def lifted(y_a,y_b,t,maps):
    y_ab=lambda q:y_a(q)+delay('A',y_b,q,maps)
    y_ba=lambda q:y_b(q)+delay('B',y_a,q,maps)
    return delay('BA',y_ab,t,maps)-y_ab(t)-delay('AB',y_ba,t,maps)+y_ba(t)


def word_gaps(t,maps):
    return (word_time('BA',t,maps)-word_time('AB',t,maps),
            word_time('BAAB',t,maps)-word_time('ABBA',t,maps))


def zero(t):
    return np.zeros_like(np.asarray(t),dtype=float)


def report():
    constant=retarded_maps()
    grid=np.linspace(8.,30.,257)
    c=lambda t:10*np.sin(.7*np.asarray(t))+.2*np.asarray(t)**2
    sa=lambda t:.1*np.sin(.43*np.asarray(t))
    ya,yb=streams(c,sa,zero,constant)
    no_signal=streams(c,zero,zero,constant)
    static_noise=first(*no_signal,grid,constant)
    signal_response=first(sa,zero,grid,constant)
    rows=[]
    for scale in (1.,.5,.25,.125):
        va,vb=.02*scale,.03*scale
        maps=retarded_maps(va=va,vb=vb)
        ramp=lambda t:10*np.asarray(t)
        y=streams(ramp,zero,zero,maps)
        gap1,gap2=word_gaps(10.,maps)
        rows.append({'rate_scale':scale,'first_word_time_gap':float(gap1),
                     'lifted_word_time_gap':float(gap2),
                     'first_source_residual':float(first(*y,10.,maps)),
                     'lifted_source_residual':float(lifted(*y,10.,maps)),
                     'predicted_suppression_factor':va+vb-va*vb})
    maps=retarded_maps(va=.02,vb=.03)
    return {'round':300,
            'scope':'Exact delay-operator identities and affine-delay tests of known TDI-style cancellation; no derived gravity, universal clock or separation of arbitrary physical signals from clock noise.',
            'constant_delay_example':{'delays':[1.,2.],
                'maximum_source_only_residual':float(np.max(abs(static_noise))),
                'maximum_surviving_signal_response':float(np.max(abs(signal_response))),
                'maximum_noise_plus_signal_identity_error':float(np.max(abs(first(ya,yb,grid,constant)-signal_response)))},
            'varying_delay_examples':rows,
            'record_ledger':{'raw_measurement_terms_first':4,'raw_measurement_terms_lifted':8,
                'maximum_raw_stream_shift_depth_first':1,'maximum_raw_stream_shift_depth_lifted':3,
                'maximum_source_history_depth_lifted':4,
                'minimum_source_time_on_test_grid':float(min(word_time(w,8.,maps) for w in ('BAAB','ABBA'))),
                'per_term_error_example':.0001,'first_readout_error_bound':.0004,
                'lifted_readout_error_bound':.0008,
                'acquisition_transmission_interpolation_and_buffer_costs_are_extra':True},
            'limits':['Delay maps, aligned records and an additive two-stream response model are supplied inputs.',
                      'The lifted combination is exact algebraically, but first-order source cancellation is not exact for general varying delays.',
                      'Signals with the same response pattern as the removable source may be canceled as well.',
                      'A nonzero residual can be source leakage, interpolation error or channel response; it does not uniquely identify time fluctuations.']}


class Audit(unittest.TestCase):
    def test_01_constant_delays_cancel_arbitrary_source_functions(self):
        maps=retarded_maps(a=1.,b=1.7)
        grid=np.linspace(8,30,301)
        for c in (lambda t:np.sin(t*t),lambda t:t**3,lambda t:np.exp(.05*t)):
            y=streams(c,zero,zero,maps)
            self.assertTrue(np.allclose(first(*y,grid,maps),0,atol=2e-10))
            self.assertTrue(np.allclose(lifted(*y,grid,maps),0,atol=2e-10))

    def test_02_noncommuting_delays_leave_the_exact_commutator(self):
        maps=retarded_maps(va=.02,vb=.03); grid=np.linspace(8,30,151)
        c=lambda t:3*np.sin(.7*t)+.2*t*t
        y=streams(c,zero,zero,maps)
        rhs=delay('BA',c,grid,maps)-delay('AB',c,grid,maps)
        self.assertTrue(np.allclose(first(*y,grid,maps),rhs))
        self.assertGreater(np.max(abs(rhs)),.01)

    def test_03_affine_word_gaps_match_exact_formulas(self):
        for va,vb in ((.02,.03),(.01,.015),(-.01,.02),(0.,0.)):
            maps=retarded_maps(a=1.,b=2.,va=va,vb=vb)
            kappa=2*va-vb; q=(1-va)*(1-vb)
            g1,g2=word_gaps(np.linspace(8,30,71),maps)
            self.assertTrue(np.allclose(g1,kappa,atol=1e-12))
            self.assertTrue(np.allclose(g2,-(1-q)*kappa,atol=1e-12))

    def test_04_lifted_identity_and_second_order_suppression(self):
        c=lambda t:10*np.asarray(t)
        for scale in (1.,.5,.25,.125):
            va,vb=.02*scale,.03*scale
            maps=retarded_maps(va=va,vb=vb); y=streams(c,zero,zero,maps)
            rhs=delay('BAAB',c,10.,maps)-delay('ABBA',c,10.,maps)
            self.assertAlmostEqual(float(lifted(*y,10.,maps)),float(rhs),places=11)
            self.assertAlmostEqual(float(first(*y,10.,maps)),.1*scale,places=11)
            self.assertAlmostEqual(float(lifted(*y,10.,maps)),-.1*scale*(va+vb-va*vb),places=11)

    def test_05_differential_channel_signal_can_survive_cancellation(self):
        maps=retarded_maps(); t=np.linspace(8,30,201)
        c=lambda q:10*np.sin(q); sa=lambda q:.1*np.sin(.43*q)
        y=streams(c,sa,zero,maps)
        signal=first(sa,zero,t,maps)
        self.assertTrue(np.allclose(first(*y,t,maps),signal))
        self.assertGreater(np.max(abs(signal)),.08)

    def test_06_a_signal_in_the_source_response_image_is_also_removed(self):
        maps=retarded_maps(); q=lambda t:.3*np.sin(.61*t)
        sa,sb=streams(q,zero,zero,maps)
        self.assertTrue(np.allclose(first(sa,sb,np.linspace(8,30,211),maps),0.))

    def test_07_source_signal_reassignment_preserves_both_raw_streams(self):
        maps=retarded_maps(va=.02,vb=.03); t=np.linspace(8,30,101)
        c=lambda z:10*np.sin(z); q=lambda z:.2*np.cos(.3*z)
        sa=lambda z:.1*np.sin(.7*z); sb=lambda z:.2*np.cos(.5*z)
        old=streams(c,sa,sb,maps)
        new=streams(lambda z:c(z)+q(z),
                    lambda z:sa(z)-(delay('A',q,z,maps)-q(z)),
                    lambda z:sb(z)-(delay('B',q,z,maps)-q(z)),maps)
        self.assertTrue(all(np.allclose(a(t),b(t)) for a,b in zip(old,new)))

    def test_08_finite_record_error_bounds_include_all_terms(self):
        maps=retarded_maps(va=.02,vb=.03); eps=1e-4
        for function,terms in ((first,4),(lifted,8)):
            # Each invocation is an independently bounded stored-term evaluation.
            for signs in itertools.product((-1.,1.),repeat=terms):
                errors=iter(eps*np.array(signs))
                sample=lambda t:next(errors)
                value=float(function(sample,sample,10.,maps))
                self.assertLessEqual(abs(value),terms*eps+1e-12)
                with self.assertRaises(StopIteration):
                    next(errors)
        c=lambda t:10*np.asarray(t)
        sa=lambda t:.1*np.sin(.43*t)
        ya,yb=streams(c,sa,zero,maps)
        target=float(first(sa,zero,10.,maps))
        # |y_A'| <= .2+.043, |y_B'| = .3 in this affine example.
        bound=10*.01+4*eps+.243*.002+.3*.001
        for signs in itertools.product((-1.,1.),repeat=6):
            da,db=.001*signs[0],.002*signs[1]
            values=[ya(maps['B'](10.)+db),ya(10.),yb(maps['A'](10.)+da),yb(10.)]
            observed=float(np.array([1.,-1.,-1.,1.])@(np.array(values)+eps*np.array(signs[2:])))
            self.assertLessEqual(abs(observed-target),bound+1e-12)

    def test_09_retarded_history_is_available_and_short_buffers_are_rejected(self):
        maps=retarded_maps(va=.02,vb=.03)
        for word in ('BA','AB','BAAB','ABBA'):
            for t in np.linspace(8,30,71):
                self.assertLessEqual(word_time(word,t,maps,(0.,30.)),t)
        with self.assertRaises(ValueError):
            word_time('BAAB',8.,maps,(5.,30.))
        used=[]
        y=lambda t:used.append(float(t)) or 0.
        lifted(y,y,10.,maps)
        self.assertEqual(len(used),8)
        self.assertTrue(all(0 <= t <= 10 for t in used))

    def test_10_order_obstruction_survives_monotone_recharting(self):
        maps=retarded_maps(va=.02,vb=.03)
        phi=lambda t:np.sinh(.1*t); inverse=lambda t:10*np.arcsinh(t)
        changed={k:(lambda z,f=f:phi(f(inverse(z)))) for k,f in maps.items()}
        for t in (8.,10.,20.):
            original=word_gaps(t,maps)[0]
            transformed=word_gaps(phi(t),changed)[0]
            self.assertEqual(np.sign(original),np.sign(transformed))
            self.assertAlmostEqual(float(word_time('BA',phi(t),changed)),float(phi(word_time('BA',t,maps))))


if __name__ == '__main__':
    main(__name__,'delayed_clock_comparison_audit',report)
