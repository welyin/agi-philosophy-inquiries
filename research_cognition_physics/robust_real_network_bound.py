"""Round 72: sharpen the real bound and quantify departures from real separability.

Let d be the trace distance of the outer premeasurement marginal to the real
separable set. The proof gives T <= 6 sqrt(2) - max(0,1-2 sqrt(2)d)^2 /
(34+35 sqrt(2)). The analytic proof, not finite sampling, establishes the bound.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bell_network_statistics import BELL_SIGNS, SETTINGS
from certified_intervals import Interval as I, SCALE
from finite_network_separation import finite_separation_certificate
from independent_source_alignment import keep_systems, pair_reference
from one_bit_real_network import real_feedback_score
from real_network_analytic_bound import branch_audit, pair_directions, random_density, random_reflection
from sourcewise_real_network_optimum import outer_axes


# Every listed pair denotes (a+b*sqrt(2))/2.
SHARP_COEFFICIENTS = ((4,1),(2,5),(2,2),(6,2),(2,1),(2,1))


def robust_bound_interval(distance_upper=Fraction(0)):
    distance_upper=Fraction(distance_upper)
    if not 0<=distance_upper<=1: raise ValueError("Use a trace-distance upper bound in [0,1].")
    root=I.exact(2).sqrt()
    deficit=1-2*root*distance_upper
    positive=I(max(0,deficit.lo),max(0,deficit.hi))
    return 6*root-positive**2/(34+35*root)


def required_distance_interval(score_lower):
    score_lower=Fraction(score_lower)
    root=I.exact(2).sqrt()
    score_loss=6*root-score_lower
    if score_loss.lo<0: raise ValueError("Score lower bound must be below the quantum maximum.")
    value=(1-((34+35*root)*score_loss).sqrt())/(2*root)
    return I(max(0,value.lo),max(0,value.hi))


def sharp_certificate():
    root=I.exact(2).sqrt()
    upper=robust_bound_interval()
    old=6*root-1/(108+95*root)
    score_lower=Fraction(finite_separation_certificate()["network_score_lower_exact"])
    required=required_distance_interval(score_lower)
    return {"real_upper_exact":"6*sqrt(2)-1/(34+35*sqrt(2))",
            "real_upper_interval":upper.floats(),"real_upper_rational":str(Fraction(upper.hi,SCALE)),
            "strict_improvement_over_round_68":upper.hi<old.lo,
            "robust_formula":"6*sqrt(2)-max(0,1-2*sqrt(2)*d)^2/(34+35*sqrt(2))",
            "d_definition":"Trace distance to real separable states of the outer premeasurement marginal; local dimensions arbitrary finite",
            "round_69_score_requires_distance_at_least_interval":required.floats(),
            "round_69_required_distance_lower_exact":str(Fraction(required.lo,SCALE)),
            "distance_tolerance_rows":[{"distance_upper_exact":str(d),"score_upper_interval":robust_bound_interval(d).floats()}
                                       for d in (Fraction(0),Fraction(1,100),Fraction(1,20),Fraction(1,10),Fraction(1,5),Fraction(1,2))]}


class RobustRealNetworkBoundTests(unittest.TestCase):
    def test_half_integer_coefficients_give_the_exact_improved_constant(self):
        rational=sum(Fraction(a*a+2*b*b,4) for a,b in SHARP_COEFFICIENTS)
        radical=sum(Fraction(2*a*b,4) for a,b in SHARP_COEFFICIENTS)
        self.assertEqual((rational,radical),(35,17))
        self.assertEqual((2*radical,rational),(34,35))

    def test_paired_charlie_product_is_half_a_commutator_with_norm_at_most_one(self):
        rng=np.random.default_rng(72)
        for dimension in (2,3,4,7):
            charlie=[random_reflection(rng,dimension) for _ in range(6)]
            pairs=pair_directions(charlie)
            d1,d2=pairs[0]
            np.testing.assert_allclose(d1@d2,(charlie[1]@charlie[0]-charlie[0]@charlie[1])/2,atol=1e-15)
            self.assertLessEqual(np.linalg.norm(d1@d2,2),1+2e-15)
            jc=d1@d2@pairs[1][1]
            self.assertLessEqual(np.linalg.norm(jc,2),math.sqrt(2)+2e-15)

    def test_sharper_branch_difference_and_weighted_error_bounds(self):
        rng=np.random.default_rng(172)
        for da,dc in ((2,3),(4,4),(5,6)):
            alice=[random_reflection(rng,da) for _ in range(3)]
            charlie=[random_reflection(rng,dc) for _ in range(6)]
            for signs in BELL_SIGNS:
                psi=rng.normal(size=da*dc)/math.sqrt(da*dc)
                row=branch_audit(alice,charlie,signs,psi)
                self.assertLessEqual(row["difference"],math.sqrt(row["mass"]/2)*row["reversal"]+1e-12)
                weighted=sum((a+b*math.sqrt(2))/2*e for (a,b),e in zip(SHARP_COEFFICIENTS,row["residuals"]))
                self.assertLessEqual(row["transfer_bound"]+row["reversal_bound"]/math.sqrt(2),weighted+1e-12)

    def test_real_separable_mixtures_annihilate_the_skew_witness(self):
        rng=np.random.default_rng(272)
        weights=rng.dirichlet(np.ones(4))
        state=sum(w*np.kron(random_density(rng,3),random_density(rng,4)) for w in weights)
        a=[random_reflection(rng,3) for _ in range(3)]
        ja=a[0]@a[1]@a[2]
        ka=(ja-ja.T)/2
        pairs=pair_directions([random_reflection(rng,4) for _ in range(6)])
        jc=pairs[0][0]@pairs[0][1]@pairs[1][1]
        self.assertAlmostEqual(np.trace(state@np.kron(ka,jc)),0.,places=15)

    def test_trace_distance_perturbation_controls_the_symmetric_witness_part(self):
        rng=np.random.default_rng(372)
        a=[random_reflection(rng,4) for _ in range(3)]
        ja=a[0]@a[1]@a[2]
        ka=(ja-ja.T)/2
        pairs=pair_directions([random_reflection(rng,4) for _ in range(6)])
        jc=pairs[0][0]@pairs[0][1]@pairs[1][1]
        witness=np.kron(ka,(jc-jc.T)/2)
        np.testing.assert_allclose(witness,witness.T,atol=0)
        self.assertLessEqual(np.linalg.norm(witness,2),math.sqrt(2)+2e-15)
        sep=np.kron(random_density(rng,4),random_density(rng,4))
        state=.9*sep+.1*random_density(rng,16)
        distance=np.abs(np.linalg.eigvalsh(state-sep)).sum()/2
        expectation=np.trace(state@witness)
        self.assertLessEqual(abs(expectation),2*math.sqrt(2)*distance+1e-15)
        self.assertAlmostEqual(expectation,np.trace(state@np.kron(ka,jc)),places=15)

    def test_feedback_reference_distance_is_exact_and_exposes_the_changed_marginal(self):
        axes=outer_axes()
        ja=axes[0]@axes[1]@axes[2]
        for corr in (0.,.2,.7,1.):
            state=keep_systems(np.kron(np.eye(4)/4,pair_reference(corr)),(0,2,1,3),4)
            distance=np.abs(np.linalg.eigvalsh(state-np.eye(16)/16)).sum()/2
            self.assertAlmostEqual(distance,corr/2,places=15)
            self.assertAlmostEqual(np.trace(state@np.kron(ja,ja)).real,-corr,places=15)
            # The YY witness proves the matching lower bound against every real
            # separable state; the identity matrix supplies this upper bound.
            upper=robust_bound_interval(Fraction(str(corr))/2).floats()[1]
            self.assertLessEqual(real_feedback_score(corr),upper+3e-15)

    def test_zero_distance_bound_is_stronger_and_remains_nontrivial_for_small_errors(self):
        report=sharp_certificate()
        self.assertTrue(report["strict_improvement_over_round_68"])
        self.assertLess(Fraction(report["real_upper_rational"]),Fraction(84734,10000))
        rows=report["distance_tolerance_rows"]
        self.assertEqual(sorted(r["score_upper_interval"][0] for r in rows),[r["score_upper_interval"][0] for r in rows])
        self.assertLess(robust_bound_interval(Fraction(1,20)).floats()[1],8.48)

    def test_inverse_certificate_for_round_69_requires_a_substantial_real_separable_distance(self):
        report=sharp_certificate()
        self.assertGreater(Fraction(report["round_69_required_distance_lower_exact"]),Fraction(3,10))
        self.assertLess(Fraction(report["round_69_required_distance_lower_exact"]),Fraction(31,100))
        self.assertEqual(required_distance_interval(Fraction(6)).lo,0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RobustRealNetworkBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":72,**sharp_certificate(),"bound_is_claimed_tight":False,
            "proof":"research_note_72.md","source_sdp_reproduced":False,
            "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("robust_real_network_bound_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
