"""Round 155: sharp conditional recovery costs on independent encoded sources.

Each matrix entry is a worst external trace distance for one true label and
one chosen decoder. Averaging these separately sharp costs is an upper bound,
not generally a sharp minimax value for the complete recovery channel.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import GROUPS, centered_decoder, purification
from certified_intervals import Interval as I
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics, negative_witness
from common_orientation_structure import encode_state, orientation
from flagged_subject_joining import joining_kraus, random_state
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary
from quantum_interface_audit import PAULI_X as X, PAULI_Z as Z


NAMES=("A_X","A_Y","A_Z","B_X","B_Y","B_Z")
# A cost is this linear combination of five radicals, divided by eight.
COEFFICIENTS={"0":(0,0,0,0,0),"1":(8,0,0,0,0),"s":(0,4,0,0,0),
              "d":(1,0,0,1,0),"g":(1,0,1,0,0),"b":(1,0,0,0,1)}


def radical_basis():
    return (I.exact(1),I.exact(3).sqrt(),I.exact(33).sqrt(),
            (I.exact(25)-16*I.exact(2).sqrt()).sqrt(),
            (I.exact(25)+16*I.exact(2).sqrt()).sqrt())


def symbols(coarse=False):
    if coarse:
        return (("d","g","b","g"),("g","d","g","b"),
                ("1","s","0","s"),("s","1","s","0"))
    matrix=[["0" if r==g else "s" for g in range(6)] for r in range(6)]
    for r,g in ((0,1),(1,2),(3,4),(4,5)): matrix[r][g]=matrix[g][r]="1"
    return tuple(tuple(row) for row in matrix)


def cost_interval(symbol):
    if symbol not in COEFFICIENTS: raise ValueError("Unknown cost symbol.")
    return sum((n*b for n,b in zip(COEFFICIENTS[symbol],radical_basis())),I.exact(0))/8


def weighted_coefficients(coarse=False):
    """Integer coefficients for prior(true)*cost(true,report), divided by 48."""
    return tuple(tuple(tuple(n*(2 if coarse and true<2 else 1) for n in COEFFICIENTS[s])
                       for s in row) for true,row in enumerate(symbols(coarse)))


def conditional_kraus(true,report,coarse=False):
    labels=4 if coarse else 6
    if any(type(v) is not int or v not in range(labels) for v in (true,report)):
        raise ValueError("True and reported labels must belong to the selected alphabet.")
    fine=GROUPS[true] if coarse else (true,)
    decoder=centered_decoder(report) if coarse else conditional_unitary(report+1).T
    return tuple(decoder @ conditional_unitary(r+1)/np.sqrt(len(fine)) for r in fine)


def attaining_states(true,report,coarse=False):
    conditional_kraus(true,report,coarse)
    states=[np.eye(2)/2,np.eye(2)/2]; kind=symbols(coarse)[true][report]
    if coarse and kind in ("1","b"):
        states[true%2]=(np.eye(2)-(X+Z)/np.sqrt(2))/2
    elif not coarse and kind=="1":
        axis=Z if {true%3,report%3}=={0,1} else X
        states[true//3]=(np.eye(2)-axis)/2
    return tuple(states)


def error_on_states(true,report,states,coarse=False):
    return pure_extension_metrics(conditional_kraus(true,report,coarse),
                                  purification(independent_encoding(states)))["trace_error"]


class HistoryRecoveryLossMatrixTests(unittest.TestCase):
    def test_negative_sector_has_the_declared_encoded_form_and_reference_symmetry(self):
        rng=np.random.default_rng(155); km=joining_kraus((2,2))[(0,1)]
        for _ in range(6):
            a,b=random_state(rng,2),random_state(rng,2)
            block=km @ independent_encoding((a,b)) @ km.T
            np.testing.assert_allclose(block,encode_state(np.kron(a,b.conj()))/2,atol=2e-16)
            np.testing.assert_allclose(block @ orientation(4),orientation(4) @ block,atol=2e-16)

    def test_fine_matrix_has_sharp_attaining_witnesses_for_every_entry(self):
        for true in range(6):
            for report in range(6):
                cost=sum(cost_interval(symbols()[true][report]).floats())/2
                self.assertAlmostEqual(error_on_states(true,report,attaining_states(true,report)),cost)

    def test_fine_conditional_errors_follow_the_exact_pure_state_overlap(self):
        rng=np.random.default_rng(255)
        for _ in range(5):
            states=(random_state(rng,2),random_state(rng,2)); omega=independent_encoding(states)
            for true in range(6):
                for report in range(6):
                    k=conditional_kraus(true,report)[0]; overlap=np.trace(omega @ k)
                    expected=np.sqrt(max(0.,1-abs(overlap)**2))
                    actual=error_on_states(true,report,states)
                    self.assertAlmostEqual(actual,expected,places=7)
                    self.assertLessEqual(actual,sum(cost_interval(symbols()[true][report]).floats())/2+2e-14)

    def test_coarse_matrix_has_sharp_attaining_witnesses_for_every_entry(self):
        for true in range(4):
            for report in range(4):
                cost=sum(cost_interval(symbols(True)[true][report]).floats())/2
                self.assertAlmostEqual(error_on_states(true,report,attaining_states(true,report,True),True),cost)

    def test_pair_to_single_maximum_respects_concavity_and_swap_symmetry(self):
        rng=np.random.default_rng(355)
        for _ in range(16):
            a=random_state(rng,2); x=np.trace(a @ X).real; z=np.trace(a @ Z).real
            swapped=(x+z)/2; averaged=(np.eye(2)+swapped*(X+Z))/2
            first=error_on_states(0,2,(a,np.eye(2)/2),True)
            symmetric=error_on_states(0,2,(averaged,np.eye(2)/2),True)
            self.assertLessEqual(first,symmetric+1e-14)
            overlap=(1+swapped)/2
            self.assertAlmostEqual(symmetric,(1+np.sqrt(49-64*overlap**2))/8)

    def test_all_coarse_costs_bound_complex_external_purifications_of_unknown_sources(self):
        rng=np.random.default_rng(455)
        for _ in range(6):
            states=(random_state(rng,2),random_state(rng,2))
            c=purification(independent_encoding(states))
            raw=rng.normal(size=(19,16))+1j*rng.normal(size=(19,16)); isometry=np.linalg.qr(raw)[0]
            for true in range(4):
                for report in range(4):
                    actual=pure_extension_metrics(conditional_kraus(true,report,True),c @ isometry.T)["trace_error"]
                    self.assertLessEqual(actual,sum(cost_interval(symbols(True)[true][report]).floats())/2+2e-14)

    def test_source_promise_is_essential_for_the_reduced_cross_subject_cost(self):
        outside=pure_extension_metrics(conditional_kraus(0,4),negative_witness())["trace_error"]
        self.assertAlmostEqual(outside,1.)
        self.assertGreater(outside,sum(cost_interval("s").floats())/2)
        self.assertNotEqual(symbols(True)[0][2],symbols(True)[2][0])

    def test_integer_weighted_costs_and_input_guards(self):
        for coarse in (False,True):
            for true,row in enumerate(weighted_coefficients(coarse)):
                prior=F(1,3) if coarse and true<2 else F(1,6)
                for report,coeffs in enumerate(row):
                    value=sum((n*b for n,b in zip(coeffs,radical_basis())),I.exact(0))/48
                    expected=prior*cost_interval(symbols(coarse)[true][report])
                    self.assertLessEqual(value.lo,expected.hi); self.assertGreaterEqual(value.hi,expected.lo)
        with self.assertRaises(ValueError): conditional_kraus(6,0)
        with self.assertRaises(ValueError): conditional_kraus(0,4,True)
        with self.assertRaises(ValueError): cost_interval("unknown")


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(HistoryRecoveryLossMatrixTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":155,"matrix_orientation":"rows=true label, columns=chosen decoder",
        "six_label_names":NAMES,"six_label_cost_symbols":symbols(),"four_label_cost_symbols":symbols(True),
        "constants":{s:interval_fields(cost_interval(s)) for s in COEFFICIENTS},
        "constant_formulas":{"s":"sqrt(3)/2","d":"(1+sqrt(25-16sqrt(2)))/8",
                             "g":"(1+sqrt(33))/8","b":"(1+sqrt(25+16sqrt(2)))/8"},
        "every_matrix_entry_sharp_on_independent_encoded_sources_and_external_extensions":True,
        "sum_of_separately_worst_costs_claimed_sharp":False,
        "all_decoder_matrices_and_source_promise_unchanged":True,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("history_recovery_loss_matrix_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
