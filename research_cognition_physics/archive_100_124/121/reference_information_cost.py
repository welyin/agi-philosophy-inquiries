"""Round 121: sharp information cost in the fixed real-reference correction task.

This optimizes classical message channels and real reference corrections, not
all real quantum networks or their total physical resource cost.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import unittest
import numpy as np

from certified_intervals import Interval as I, SCALE, sin_interval
from joint_relation_records import interval_fields
from noisy_message_alignment import corrected_reference_by_channel, binary_symmetric_channel
from one_bit_real_network import records_for_reference, visibility_interval
from bell_network_statistics import endpoint_observables, bob_effects, network_table


def unit_log_interval(x, terms=70):
    """1<=x<=2, log x = 2 sum z^(2k+1)/(2k+1), z=(x-1)/(x+1)."""
    x=F(x)
    if not 1<=x<=2 or terms<1: raise ValueError("Use x in [1,2] and positive terms.")
    z=I.exact((x-1)/(x+1))
    power=z
    value=I.exact(0)
    square=z*z
    for k in range(terms):
        value+=2*power/(2*k+1)
        power*=square
    tail=2*power/((2*terms+1)*(1-square))
    return I(value.lo,value.hi+tail.hi)


@lru_cache(maxsize=1)
def ln2_interval():
    return unit_log_interval(2)


def log_interval(x):
    x=F(x)
    if x<=0: raise ValueError("Positive logarithm argument required.")
    exponent=0
    while x<1:
        x*=2
        exponent-=1
    while x>2:
        x/=2
        exponent+=1
    return unit_log_interval(x)+exponent*ln2_interval()


def entropy_interval(p):
    p=F(p)
    if not 0<=p<=1: raise ValueError("Use probability in [0,1].")
    if p in (0,1): return I.exact(0)
    if p==F(1,2): return I.exact(1)
    return (-p*log_interval(p)-(1-p)*log_interval(1-p))/ln2_interval()


def information_frontier(d):
    d=F(d)
    if not 0<=d<=1: raise ValueError("Use distinguishability in [0,1].")
    return 1-entropy_interval((1-d)/2)


def information_for_ratio_interval(ratio):
    if ratio.lo<0 or ratio.hi>SCALE:
        raise ValueError("The target exceeds the fixed alignment visibility.")
    lower=information_frontier(F(ratio.lo,SCALE))
    upper=information_frontier(F(ratio.hi,SCALE))
    return I(lower.lo,upper.hi)


def h2(p):
    if p<=0 or p>=1: return 0.
    return -p*math.log2(p)-(1-p)*math.log2(1-p)


def channel_statistics(channel):
    """Uniform binary input; finite output alphabet. Floating diagnostics."""
    w=np.asarray(channel,dtype=float)
    if w.ndim!=2 or len(w)!=2 or np.min(w)<0 or not np.allclose(w.sum(axis=1),1):
        raise ValueError("Two normalized channel rows required.")
    marginal=w.sum(axis=0)/2
    conditional_entropy=sum(mass*h2(w[0,j]/(2*mass)) for j,mass in enumerate(marginal) if mass>0)
    d=float(np.abs(w[0]-w[1]).sum()/2)
    return {"distinguishability":d,"mutual_information_bits":1-conditional_entropy,
            "sharp_lower_bound_diagnostic":1-h2((1-d)/2)}


def block_statistics(channel):
    """2^N equiprobable independent input flags, one shared classical message."""
    w=np.asarray(channel,dtype=float)
    n=int(math.log2(len(w)))
    if len(w)!=2**n or np.min(w)<0 or not np.allclose(w.sum(axis=1),1):
        raise ValueError("Use 2^N normalized rows.")
    marginal=w.mean(axis=0)
    conditional=0.
    for j,mass in enumerate(marginal):
        if mass:
            posterior=w[:,j]/(len(w)*mass)
            conditional-=mass*sum(p*math.log2(p) for p in posterior if p>0)
    distances=[]
    for bit in range(n):
        rows0=[i for i in range(2**n) if not (i>>bit)&1]
        rows1=[i for i in range(2**n) if (i>>bit)&1]
        distances.append(channel_statistics([w[rows0].mean(axis=0),w[rows1].mean(axis=0)])["distinguishability"])
    return n-conditional,distances


@lru_cache(maxsize=1)
def target_intervals():
    alpha=4*sin_interval(I.rational(1,4))
    nu=I.rational(1,5)
    for _ in range(5):
        nu=nu*(1+alpha)/(1+alpha*nu**2)
    gamma=visibility_interval(5)
    target=nu**2
    return nu,gamma,target,target/gamma


@lru_cache(maxsize=1)
def cost_certificate():
    nu,gamma,target,ratio=target_intervals()
    info=information_for_ratio_interval(ratio)
    error=(1-ratio)/2
    rows=[]
    for d in (F(0),F(1,2),F(9,10),F(99,100),F(1)):
        rows.append({"required_message_distinguishability":str(d),
                     "minimum_information_bits":interval_fields(information_frontier(d))})
    return {"private_complex_bias":interval_fields(nu),"real_alignment_visibility":interval_fields(gamma),
            "target_reference_correlation":interval_fields(target),
            "required_message_distinguishability":interval_fields(ratio),
            "max_symmetric_channel_error":interval_fields(error),
            "minimum_message_information_bits":interval_fields(info),
            "independent_1000_trial_information_bits_lower":str(F(1000*info.lo,SCALE)),
            "independent_1000_trial_fixed_length_bits_necessary":math.ceil(F(1000*info.lo,SCALE)),
            "frontier_examples":rows,"all_network_statistics_compared":288,
            "total_cost_real_greater_than_complex_proved":False,
            "scope":"fixed sources, alignment, target effects; arbitrary classical channels and real reference CPTP corrections"}


class ReferenceInformationCostTests(unittest.TestCase):
    def test_log_certificate_contains_independent_high_precision_values(self):
        with localcontext() as ctx:
            ctx.prec=85
            for x in (F(1),F(2),F(1,10**8),F(17,19),F(37,5)):
                result=log_interval(x)
                actual=(Decimal(x.numerator)/Decimal(x.denominator)).ln()
                self.assertLessEqual(Decimal(result.lo)/Decimal(SCALE),actual)
                self.assertGreaterEqual(Decimal(result.hi)/Decimal(SCALE),actual)

    def test_entropy_endpoints_symmetry_and_independent_decimal_values(self):
        self.assertEqual(entropy_interval(F(1,2)).lo,SCALE)
        self.assertEqual(entropy_interval(0).hi,0)
        with localcontext() as ctx:
            ctx.prec=85
            for p in (F(1,4),F(1,10000),F(2,7)):
                actual_p=Decimal(p.numerator)/Decimal(p.denominator)
                actual=-(actual_p*actual_p.ln()+(1-actual_p)*(1-actual_p).ln())/Decimal(2).ln()
                result=entropy_interval(p)
                self.assertLessEqual(Decimal(result.lo)/Decimal(SCALE),actual)
                self.assertGreaterEqual(Decimal(result.hi)/Decimal(SCALE),actual)
                symmetric=entropy_interval(1-p)
                self.assertLessEqual(result.lo,symmetric.hi)
                self.assertLessEqual(symmetric.lo,result.hi)

    def test_binary_symmetric_messages_reach_the_sharp_information_frontier(self):
        for error in (0,.001,.1,.25,.5):
            report=channel_statistics(binary_symmetric_channel(error))
            self.assertAlmostEqual(report["mutual_information_bits"],report["sharp_lower_bound_diagnostic"],places=14)

    def test_arbitrary_finite_messages_and_erasure_obey_the_bound(self):
        rng=np.random.default_rng(121)
        for alphabet in (2,3,9):
            for _ in range(20):
                w=rng.random((2,alphabet));w/=w.sum(axis=1,keepdims=True)
                report=channel_statistics(w)
                self.assertGreaterEqual(report["mutual_information_bits"]+1e-14,report["sharp_lower_bound_diagnostic"])
        erasure=channel_statistics([[.5,0,.5],[0,.5,.5]])
        self.assertAlmostEqual(erasure["mutual_information_bits"],.5)
        self.assertGreater(erasure["mutual_information_bits"],erasure["sharp_lower_bound_diagnostic"])

    def test_independent_flags_require_individual_information_even_with_joint_messages(self):
        for w in (np.eye(4),np.array([[1,0],[0,1],[0,1],[1,0]]),np.ones((4,1))):
            info,distances=block_statistics(w)
            self.assertGreaterEqual(info+1e-14,sum(1-h2((1-d)/2) for d in distances))
            self.assertGreaterEqual(info+1e-14,2*(1-h2((1-sum(distances)/2)/2)))
        info,distances=block_statistics(np.array([[1,0],[0,1],[0,1],[1,0]]))
        self.assertAlmostEqual(info,1)
        self.assertEqual(distances,[0,0])

    def test_matched_message_reproduces_all_288_finite_complex_statistics(self):
        nu,gamma,_,ratio=target_intervals()
        bias=sum(nu.floats())/2
        vis=sum(gamma.floats())/2
        d=sum(ratio.floats())/2
        corrected=corrected_reference_by_channel(binary_symmetric_channel((1-d)/2),alignment_count=5)
        real=records_for_reference(corrected,vis,vis,vis)
        alice,charlie=endpoint_observables(vis,bias)
        complex_records=network_table(alice,charlie,bob_effects(vis))
        self.assertEqual(len(real),288)
        for key in real:
            self.assertAlmostEqual(real[key],complex_records[key],places=13)

    def test_target_cost_is_nonzero_and_strictly_below_one_bit_of_information(self):
        report=cost_certificate()
        bound=report["minimum_message_information_bits"]
        self.assertGreater(F(bound["lower_exact"]),F(999,1000))
        self.assertLess(F(bound["upper_exact"]),1)
        error=report["max_symmetric_channel_error"]
        self.assertGreater(F(error["lower_exact"]),0)
        self.assertLess(F(error["upper_exact"]),F(3,100000))
        self.assertFalse(report["total_cost_real_greater_than_complex_proved"])

    def test_infeasible_targets_and_invalid_probabilities_are_rejected(self):
        with self.assertRaises(ValueError): information_for_ratio_interval(I.rational(101,100))
        with self.assertRaises(ValueError): entropy_interval(F(-1,2))
        with self.assertRaises(ValueError): log_interval(0)
        with self.assertRaises(ValueError): channel_statistics([[.2,.2],[.5,.5]])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceInformationCostTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":121,"certificate":cost_certificate(),
            "quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("reference_information_cost_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
