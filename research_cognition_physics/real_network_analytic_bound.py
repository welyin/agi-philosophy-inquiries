"""Round 68: a conservative dimension-independent analytic real-network bound.

The proof is in research_note_68.md. Tests audit its identities and constants;
random matrices are not the proof. Assumptions: real matrices, ordinary tensor
composition, two independent sources, local measurements, all four Bob outputs.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bell_network_statistics import AXES, BELL_SIGNS, PAIRS, endpoint_observables
from certified_intervals import Interval as I, SCALE


# Coefficients a+b*sqrt(2), in residual order 01:0,1; 02:0,2; 12:1,2.
COEFFICIENTS = ((3,1),(2,4),(2,2),(4,2),(2,1),(2,1))


def bound_certificate():
    root = I.exact(2).sqrt()
    gap = 1/(108+95*root)
    upper = 6*root-gap
    return {"complex_ideal_exact":"6*sqrt(2)",
            "real_upper_exact":"6*sqrt(2)-1/(108+95*sqrt(2))",
            "gap_exact":"1/(108+95*sqrt(2))",
            "gap_interval":gap.floats(),"real_upper_interval":upper.floats(),
            "real_upper_rational":str(Fraction(upper.hi,SCALE)),
            "gap_lower_rational":str(Fraction(gap.lo,SCALE)),
            "squared_coefficient_sum":"95+54*sqrt(2)","interval_bits":100}


def random_reflection(rng, dimension):
    orthogonal, _ = np.linalg.qr(rng.normal(size=(dimension,dimension)))
    signs = np.ones(dimension)
    signs[:dimension//2] = -1
    return (orthogonal*signs)@orthogonal.T


def pair_directions(charlie):
    return tuple(((charlie[2*p]+charlie[2*p+1])/math.sqrt(2),
                  (charlie[2*p]-charlie[2*p+1])/math.sqrt(2)) for p in range(3))


def branch_audit(alice, charlie, signs, psi):
    """psi may be unnormalized; Alice and Charlie act on the two tensor factors."""
    da, dc = alice[0].shape[0], charlie[0].shape[0]
    aa = [np.kron(a,np.eye(dc)) for a in alice]
    dd = [(np.kron(np.eye(da),di),np.kron(np.eye(da),dj)) for di,dj in pair_directions(charlie)]
    residuals = []
    anticomm = []
    score = 0.
    for p,(i,j) in enumerate(PAIRS):
        di,dj = dd[p]
        residuals.extend((np.linalg.norm((signs[i]*aa[i]-di)@psi), np.linalg.norm((signs[j]*aa[j]-dj)@psi)))
        anticomm.append(np.linalg.norm((aa[i]@aa[j]+aa[j]@aa[i])@psi))
        score += math.sqrt(2)*np.vdot(psi,(signs[i]*aa[i]@di+signs[j]*aa[j]@dj)@psi).real
    ja = aa[0]@aa[1]@aa[2]
    jc = dd[0][0]@dd[0][1]@dd[1][1]
    # Adjoint equals transpose for the real case. This also audits the complex
    # negative control with a skew-Hermitian K, whose local trace need not vanish.
    ka = (ja-ja.conj().T)/2
    reversal = np.linalg.norm((ja+ja.conj().T)@psi)
    transfer = np.linalg.norm((ja@jc+np.eye(da*dc))@psi)
    e = residuals
    transfer_bound = e[0]+math.sqrt(2)*e[1]+2*e[3]
    reversal_bound = anticomm[2]+math.sqrt(2)*anticomm[1]+anticomm[0]+2*e[1]
    difference = abs(np.vdot(psi,(ja-ka)@jc@psi))
    mass = np.vdot(psi,psi).real
    weighted_error = sum((a+b*math.sqrt(2))*v for (a,b),v in zip(COEFFICIENTS,e))
    return {"mass":float(mass),"score":float(score),"residuals":e,"residual_square_sum":float(np.dot(e,e)),
            "anticommutators":anticomm,"transfer":float(transfer),"transfer_bound":float(transfer_bound),
            "reversal":float(reversal),"reversal_bound":float(reversal_bound),
            "difference":float(difference),"difference_bound":math.sqrt(2*mass)*reversal,
            "weighted_error":float(weighted_error),"skew_expectation":np.vdot(psi,ka@jc@psi)}


def real_naimark_dilation(contraction):
    values, vectors = np.linalg.eigh(contraction)
    root = (vectors*np.sqrt(np.maximum(0.,1-values**2)))@vectors.T
    return np.block([[contraction,root],[root,-contraction]])


def random_density(rng, dimension):
    raw = rng.normal(size=(dimension,dimension))
    density = raw@raw.T
    return density/np.trace(density)


def random_network_assemblage(rng, da=3, dc=4):
    left, right = random_density(rng,da*2), random_density(rng,2*dc)
    tensor = np.kron(left,right).reshape(da,4,dc,da,4,dc)
    orthogonal, _ = np.linalg.qr(rng.normal(size=(4,4)))
    effects = [np.outer(orthogonal[:,b],orthogonal[:,b]) for b in range(4)]
    assemblage = [np.einsum("kj,ajclkd->acld",e,tensor).reshape(da*dc,da*dc) for e in effects]
    marginal_a = np.einsum("abcb->ac",left.reshape(da,2,da,2))
    marginal_c = np.einsum("abac->bc",right.reshape(2,dc,2,dc))
    return assemblage,marginal_a,marginal_c


class RealNetworkAnalyticBoundTests(unittest.TestCase):
    def test_exact_coefficient_squares_give_the_stated_constant(self):
        rational = sum(a*a+2*b*b for a,b in COEFFICIENTS)
        radical = sum(2*a*b for a,b in COEFFICIENTS)
        self.assertEqual((rational,radical),(95,54))
        self.assertEqual((2*radical,rational),(108,95))

    def test_uniform_gap_has_a_strict_positive_interval_certificate(self):
        report = bound_certificate()
        self.assertGreater(Fraction(report["gap_lower_rational"]),Fraction(4,1000))
        self.assertLess(Fraction(report["real_upper_rational"]),Fraction(84812,10000))
        self.assertGreater(Fraction(report["real_upper_rational"]),Fraction(8481,1000))

    def test_pair_directions_have_exact_anticommutator_and_squared_sum_identities(self):
        rng = np.random.default_rng(681)
        for dimension in (2,3,6):
            charlie = [random_reflection(rng,dimension) for _ in range(6)]
            for di,dj in pair_directions(charlie):
                np.testing.assert_allclose(di@dj+dj@di,0,atol=3e-15)
                np.testing.assert_allclose(di@di+dj@dj,2*np.eye(dimension),atol=3e-15)
                self.assertLessEqual(np.linalg.norm(di,2),math.sqrt(2)+2e-15)

    def test_noncommuting_word_reversal_identity(self):
        rng = np.random.default_rng(682)
        for _ in range(6):
            a,b,c = [rng.normal(size=(5,5)) for _ in range(3)]
            rhs = a@(b@c+c@b)-(a@c+c@a)@b+c@(a@b+b@a)
            np.testing.assert_allclose(a@b@c+c@b@a,rhs,atol=3e-14)

    def test_sum_of_squares_and_each_transfer_bound_on_real_branches(self):
        rng = np.random.default_rng(683)
        for da,dc in ((2,2),(3,4),(4,6)):
            alice = [random_reflection(rng,da) for _ in range(3)]
            charlie = [random_reflection(rng,dc) for _ in range(6)]
            for signs in BELL_SIGNS:
                psi = rng.normal(size=da*dc)/math.sqrt(da*dc)
                row = branch_audit(alice,charlie,signs,psi)
                self.assertAlmostEqual(row["residual_square_sum"],math.sqrt(2)*(6*math.sqrt(2)*row["mass"]-row["score"]),places=12)
                for p,anti in enumerate(row["anticommutators"]):
                    self.assertLessEqual(anti,(1+math.sqrt(2))*sum(row["residuals"][2*p:2*p+2])+1e-12)
                for name in ("transfer","reversal","difference"):
                    self.assertLessEqual(row[name],row[name+"_bound"]+1e-12)
                self.assertLessEqual(row["transfer_bound"]+math.sqrt(2)*row["reversal_bound"],row["weighted_error"]+1e-12)

    def test_real_product_marginal_forces_the_skew_expectation_to_vanish(self):
        rng = np.random.default_rng(684)
        assemblage,ra,rc = random_network_assemblage(rng)
        np.testing.assert_allclose(sum(assemblage),np.kron(ra,rc),atol=5e-17)
        alice = [random_reflection(rng,3) for _ in range(3)]
        charlie = [random_reflection(rng,4) for _ in range(6)]
        ja = alice[0]@alice[1]@alice[2]
        ka = (ja-ja.T)/2
        dd = pair_directions(charlie)
        jc = dd[0][0]@dd[0][1]@dd[1][1]
        np.testing.assert_array_equal(ka.T,-ka)
        self.assertAlmostEqual(np.trace(ra@ka),0.,places=16)
        self.assertAlmostEqual(sum(np.trace(s@np.kron(ka,jc)) for s in assemblage),0.,places=16)

    def test_full_real_network_obeys_weighted_sos_and_cauchy_steps(self):
        rng = np.random.default_rng(685)
        assemblage,_,_ = random_network_assemblage(rng,4,4)
        alice = [random_reflection(rng,4) for _ in range(3)]
        charlie = [random_reflection(rng,4) for _ in range(6)]
        rows = []
        # Spectral components replace purification; every inequality is homogeneous
        # and Cauchy also applies with (b, eigenvector) as the branch label.
        for state,signs in zip(assemblage,BELL_SIGNS):
            values,vectors = np.linalg.eigh(state)
            for value,vector in zip(values,vectors.T):
                self.assertGreaterEqual(value,-1e-16)
                rows.append(branch_audit(alice,charlie,signs,math.sqrt(max(value,0))*vector))
        score = sum(r["score"] for r in rows)
        residual = sum(r["residual_square_sum"] for r in rows)
        weighted = sum(math.sqrt(r["mass"])*r["weighted_error"] for r in rows)
        self.assertAlmostEqual(sum(r["mass"] for r in rows),1.,places=14)
        self.assertAlmostEqual(sum(r["skew_expectation"] for r in rows),0.,places=14)
        self.assertGreaterEqual(weighted,1.)
        self.assertLessEqual(weighted,math.sqrt((95+54*math.sqrt(2))*residual)+1e-12)
        self.assertAlmostEqual(residual,math.sqrt(2)*(6*math.sqrt(2)-score),places=12)
        self.assertLess(score,float(Fraction(bound_certificate()["real_upper_rational"])))

    def test_general_real_binary_povms_have_real_reflection_dilations(self):
        rng = np.random.default_rng(686)
        for dimension in (2,3,5):
            contraction = random_density(rng,dimension)*1.5-.5*np.eye(dimension)
            dilation = real_naimark_dilation(contraction)
            np.testing.assert_allclose(dilation,dilation.T,atol=2e-16)
            np.testing.assert_allclose(dilation@dilation,np.eye(2*dimension),atol=2e-15)
            np.testing.assert_array_equal(dilation[:dimension,:dimension],contraction)

    def test_complex_ideal_is_a_negative_control_for_the_real_trace_step(self):
        alice,charlie = endpoint_observables()
        phi = np.array([1,0,0,1])/math.sqrt(2)
        row = branch_audit(alice,charlie,BELL_SIGNS[0],phi)
        self.assertLess(row["residual_square_sum"],1e-28)
        self.assertAlmostEqual(row["score"],6*math.sqrt(2),places=13)
        self.assertAlmostEqual(row["skew_expectation"],-1.,places=14)
        ja = AXES[0]@AXES[1]@AXES[2]
        np.testing.assert_array_equal(ja,1j*np.eye(2))
        self.assertNotEqual(np.trace(np.eye(2)/2@ja),0.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RealNetworkAnalyticBoundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round":68,**bound_certificate(),"scope":"Arbitrary finite local dimensions, real matrix states/effects, ordinary tensor product, two independent sources, no interparty measurement messages, all Bob outcomes",
              "bound_is_tight":False,"proof":"research_note_68.md: SOS, product transfer, real antisymmetric trace, Cauchy-Schwarz",
              "known_network_separation_source":"https://arxiv.org/pdf/2101.10873",
              "source_sdp_reproduced":False,"quantum_theory_derived_from_cognition":False,
              "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("real_network_analytic_bound_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
