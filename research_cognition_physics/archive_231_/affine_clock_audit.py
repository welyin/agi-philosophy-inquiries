"""Round 289: affine-clock calibration and timestamp identifiability.

Applies Freris-Graham-Kumar (2011) and the four-timestamp exchange underlying
RFC 5905. The common affine representation and fixed delays are model inputs,
not a derived universal physical clock.
"""
from collections import deque
from itertools import product
import unittest
import numpy as np
from growing_stream_audit import main


class AffineNetwork:
    """Hidden simulator parameters; inference receives timestamp pairs only."""
    def __init__(self,rates,offsets,delays):
        self.rates = np.asarray(rates,dtype=float)
        self.offsets = np.asarray(offsets,dtype=float)
        self.delays = dict(delays)
        assert len(self.rates)==len(self.offsets) and min(self.rates)>0
        assert all(v>0 for v in self.delays.values())

    def receive(self,i,j,sent):
        departure = (sent-self.offsets[i])/self.rates[i]
        return self.rates[j]*(departure+self.delays[i,j])+self.offsets[j]

    def exchange(self,i,j,sent,wait_ticks):
        received = self.receive(i,j,sent)
        reply = received+wait_ticks
        returned = self.receive(j,i,reply)
        return np.array([sent,received,reply,returned])

    def shifted(self,potentials):
        f = np.asarray(potentials,dtype=float)
        return AffineNetwork(self.rates,self.offsets+self.rates*f,
                             {(i,j):d+f[i]-f[j] for (i,j),d in self.delays.items()})


def timestamp_data(net,span=20.):
    return {(i,j):np.array([[0.,net.receive(i,j,0.)],
                            [span,net.receive(i,j,span)]])
            for i,j in net.delays}


def fit_edge(data):
    (x0,y0),(x1,y1) = np.asarray(data,dtype=float)
    if x1<=x0:
        raise ValueError('distinct ordered send stamps required')
    rate = (y1-y0)/(x1-x0)
    return float(rate),float(y0-rate*x0)


def calibrate(data,count,reference=0):
    fits = {edge:fit_edge(samples) for edge,samples in data.items()}
    relative = {reference:1.}
    queue = deque([reference])
    while queue:
        i = queue.popleft()
        for (u,j),(rate,_) in fits.items():
            if u==i and j not in relative:
                if rate<=0:
                    raise ValueError('nonpositive fitted clock rate')
                relative[j] = relative[i]*rate
                queue.append(j)
    if len(relative)!=count:
        raise ValueError('not all clocks are reachable from the reference')
    residual = max(abs(relative[j]/relative[i]-rate)
                   for (i,j),(rate,_) in fits.items())
    if residual>1e-9:
        raise ValueError('constant-rate / constant-delay model is inconsistent')
    intercepts = {(i,j):intercept/relative[j] for (i,j),(_,intercept) in fits.items()}
    return relative,intercepts,float(residual)


def corrected_roundtrip(stamps,ratio):
    x1,y2,y3,x4 = stamps
    return float((x4-x1)-(y3-y2)/ratio)


def midpoint_offset(normalized_stamps):
    x1,y2,y3,x4 = normalized_stamps
    return float(((y2-x1)+(y3-x4))/2)


def slope_error_bound(ratio,span,send_error,receive_error):
    if span<=2*send_error:
        raise ValueError('timestamp span too short for this bound')
    return (2*receive_error+2*ratio*send_error)/(span-2*send_error)


def parameter_matrix(count,edges):
    """Rows: h_ij = delay_ij + beta_j - beta_i; beta_0 fixed."""
    edges = list(edges)
    matrix = np.zeros((len(edges),len(edges)+count-1))
    for row,(i,j) in enumerate(edges):
        matrix[row,row] = 1.
        if i:
            matrix[row,len(edges)+i-1] = -1.
        if j:
            matrix[row,len(edges)+j-1] = 1.
    return matrix


def fixture():
    delays = {}
    for k,(i,j) in enumerate(((0,1),(1,2),(2,3),(3,0),(0,2))):
        delays[i,j] = 1.+.2*k
        delays[j,i] = 2.+.15*k
    return AffineNetwork([1.25,.75,1.5,2.],[3.,-2.,4.,7.],delays)


def two_model_case():
    first = AffineNetwork([1,1],[0,4],{(0,1):1.,(1,0):5.})
    second = first.shifted([0,-2])
    stamps = first.exchange(0,1,10.,3.)
    return {'first_offsets':first.offsets.tolist(),
            'first_directional_delays':[first.delays[0,1],first.delays[1,0]],
            'second_offsets':second.offsets.tolist(),
            'second_directional_delays':[second.delays[0,1],second.delays[1,0]],
            'shared_four_timestamps':stamps.tolist(),
            'maximum_timestamp_difference':float(max(abs(stamps-second.exchange(0,1,10.,3.)))),
            'roundtrip':corrected_roundtrip(stamps,1.),
            'midpoint_estimate':midpoint_offset(stamps),
            'first_true_offset_difference':4.,'first_midpoint_bias':-2.}


def report():
    net = fixture()
    data = timestamp_data(net)
    relative,h,residual = calibrate(data,4)
    shifted = net.shifted([0,.2,-.3,.4])
    difference = max(float(np.max(abs(values-timestamp_data(shifted)[e])))
                     for e,values in data.items())
    matrix = parameter_matrix(4,net.delays)
    rtts = []
    for i,j in sorted(e for e in net.delays if e[0]<e[1]):
        ratio = relative[j]/relative[i]
        stamps = net.exchange(i,j,10.,7.)
        measured = corrected_roundtrip(stamps,ratio)/relative[i]
        rtts.append({'edge':[i,j],'roundtrip_reference_ticks':measured,
                     'intercept_sum':h[i,j]+h[j,i],
                     'half_roundtrip':measured/2})
    return {'round':289,
            'scope':'Conditional calibration of affine clocks on a known static network; no absolute time, one-way delay, light speed or spacetime derived.',
            'relative_rates':[relative[i] for i in range(4)],
            'rate_cycle_residual':residual,
            'roundtrips':rtts,
            'gauge_equivalent_timestamp_max_error':difference,
            'shifted_minimum_delay':min(shifted.delays.values()),
            'offset_delay_linear_system':{'observations':matrix.shape[0],
                                         'unknowns':matrix.shape[1],
                                         'rank':int(np.linalg.matrix_rank(matrix)),
                                         'nullity':matrix.shape[1]-int(np.linalg.matrix_rank(matrix))},
            'two_model_counterexample':two_model_case(),
            'nonstationary_counterexample':{
                'true_clock_rate_ratio':1.,'delay_law':'2 + 0.03 * send_time',
                'apparent_ratio_from_two_sends':1.03,
                'conclusion':'Delay drift can mimic clock-rate drift when stationarity is removed.'},
            'precision_example':{'true_ratio':1.2,'send_error_bound':.01,
                                 'receive_error_bound':.02,
                                 'slope_bounds_by_span':[
                                     {'span':s,'absolute_ratio_error_upper':slope_error_bound(1.2,s,.01,.02)}
                                     for s in (1.,10.,100.)]},
            'resource_ledger':{'directed_links':len(net.delays),
                               'calibration_one_way_probes':2*len(net.delays),
                               'local_timestamp_samples':4*len(net.delays),
                               'per_link_send_span_local_ticks':20.,
                               'extra_echo_transmissions':2*len(rtts),
                               'extra_echo_timestamp_samples':4*len(rtts),
                               'scope':'Counts probe transmissions and raw timestamp samples; encoding, delivery of records, channel occupancy and clock stability over the span remain explicit resources.'}}


class Audit(unittest.TestCase):
    def test_01_relative_rates_from_only_timestamps(self):
        net = fixture()
        relative,h,residual = calibrate(timestamp_data(net),4)
        self.assertLess(residual,1e-12)
        for i in range(4):
            self.assertAlmostEqual(relative[i],net.rates[i]/net.rates[0])
        for (i,j),value in h.items():
            expected = net.rates[0]*net.delays[i,j]+net.offsets[j]/relative[j]-net.offsets[i]/relative[i]
            self.assertAlmostEqual(value,expected)

    def test_02_wait_subtraction_and_reference_units(self):
        net = fixture()
        for i,j in net.delays:
            for hold in (0.,3.,100.):
                value = corrected_roundtrip(net.exchange(i,j,10.,hold),net.rates[j]/net.rates[i])
                self.assertAlmostEqual(value,net.rates[i]*(net.delays[i,j]+net.delays[j,i]))

    def test_03_all_local_transcripts_gauge_equivalent(self):
        net = fixture()
        other = net.shifted([0,.2,-.3,.4])
        for i,j in net.delays:
            for sent in (-10.,0.,123.):
                for hold in (0.,7.):
                    self.assertLess(np.max(abs(net.exchange(i,j,sent,hold)-
                                               other.exchange(i,j,sent,hold))),1e-12)

    def test_04_nonunique_offsets_and_midpoint_bias(self):
        row = two_model_case()
        self.assertEqual(row['maximum_timestamp_difference'],0.)
        self.assertEqual(row['shared_four_timestamps'],[10.,15.,18.,19.])
        self.assertEqual(row['roundtrip'],6.)
        self.assertEqual(row['midpoint_estimate'],2.)
        self.assertNotEqual(row['first_offsets'],row['second_offsets'])

    def test_05_positive_delay_nullspace(self):
        net = fixture()
        matrix = parameter_matrix(4,net.delays)
        self.assertEqual(matrix.shape[1]-np.linalg.matrix_rank(matrix),3)
        f = np.array([0.,.2,-.3,.4])
        direction = np.array([f[i]-f[j] for i,j in net.delays]+f[1:].tolist())
        self.assertLess(np.max(abs(matrix@direction)),1e-14)
        self.assertGreater(min(net.shifted(f).delays.values()),0)

    def test_06_rate_consistency_and_missing_links(self):
        data = timestamp_data(fixture())
        bad = {e:values.copy() for e,values in data.items()}
        bad[0,2][1,1] += 1.
        with self.assertRaises(ValueError):
            calibrate(bad,4)
        with self.assertRaises(ValueError):
            calibrate({(0,1):data[0,1]},4)

    def test_07_delay_drift_counterexample(self):
        measured = np.array([[x,x+2+.03*x] for x in (0.,10.)])
        ratio,intercept = fit_edge(measured)
        self.assertAlmostEqual(ratio,1.03)
        self.assertAlmostEqual(intercept,2.)
        self.assertNotAlmostEqual(ratio,1.)

    def test_08_bounded_timestamp_error(self):
        for span in (1.,10.,100.):
            bound = slope_error_bound(1.2,span,.01,.02)
            for e0,e1,f0,f1 in product((-1,1),repeat=4):
                measured = [[.01*e0,3+.02*f0],
                            [span+.01*e1,1.2*span+3+.02*f1]]
                estimate,_ = fit_edge(measured)
                self.assertLessEqual(abs(estimate-1.2),bound+1e-13)
        with self.assertRaises(ValueError):
            slope_error_bound(1.,.01,.01,.01)


if __name__ == '__main__':
    main(__name__,'affine_clock_audit',report)

