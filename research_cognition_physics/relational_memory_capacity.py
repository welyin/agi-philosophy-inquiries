"""Round 96: repeated evidence versus independently variable relations."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from coherent_memory_compression import compressor, vector


def exact_rank(rows):
    work = [list(map(Fraction,row)) for row in rows]
    rank = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((r for r in range(rank,len(work)) if work[r][column]),None)
        if pivot is None: continue
        work[rank],work[pivot] = work[pivot],work[rank]
        scale=work[rank][column]
        work[rank]=[x/scale for x in work[rank]]
        for r in range(rank+1,len(work)):
            scale=work[r][column]
            work[r]=[x-scale*y for x,y in zip(work[r],work[rank])]
        rank+=1
        if rank==len(work): break
    return rank


def grouped_gram(groups):
    """Each group records ONE binary label; labels vary freely between groups."""
    groups=tuple(tuple(Fraction(c) for c in group) for group in groups)
    if not groups or any(not g for g in groups) or any(not 0<=c<=1 for g in groups for c in g):
        raise ValueError("Use nonempty groups of exact overlaps in [0,1].")
    cs=tuple(math.prod(g) for g in groups)
    labels=tuple(product((0,1),repeat=len(groups)))
    gram=[[math.prod(c if s!=t else Fraction(1) for c,s,t in zip(cs,left,right))
           for right in labels] for left in labels]
    return cs,labels,gram


def hypothesis_columns(groups):
    _,labels,_=grouped_gram(groups)
    columns=[]
    for label in labels:
        state=np.ones(1)
        for sign,group in zip(label,groups):
            for c in group:
                state=np.kron(state,vector(float(c)) if sign else [1.,0])
        columns.append(state)
    return np.column_stack(columns)


def capacity_certificate(groups):
    cs,_,gram=grouped_gram(groups)
    active=sum(c<1 for c in cs)
    rank=exact_rank(gram)
    expected=2**active
    if rank!=expected: raise AssertionError("Exact Gram rank disagrees with tensor proof.")
    determinant=math.prod((1-c*c)**(2**(len(cs)-1)) for c in cs)
    return {"group_overlaps_exact":[str(c) for c in cs],
        "number_of_weak_memories":sum(map(len,groups)),
        "number_of_independently_variable_binary_labels":len(groups),
        "active_nontrivial_labels":active,"gram_rank_exact":rank,
        "gram_determinant_exact":str(determinant),
        "minimum_coherent_output_dimension":rank,
        "minimum_coherent_qubits":active,
        "scope":"Reversible compression of the full candidate span with fixed, uncorrelated discarded systems",
        "minimal_memory_for_one_fixed_classical_decision_claimed":False}


def off_promise_scratch_probabilities(a,b):
    u=compressor(a,b)
    rows=[]
    for left,right in product((0,1),repeat=2):
        state=np.kron(vector(a) if left else [1.,0],vector(b) if right else [1.,0])
        output=(u@state).reshape(2,2)
        rows.append({"labels":[left,right],
                     "scratch_one_probability":float(np.abs(output[:,1])@np.abs(output[:,1]))})
    return rows


class RelationalMemoryCapacityTests(unittest.TestCase):
    def test_exact_gram_matches_actual_original_weak_memory_vectors(self):
        groups=((Fraction(4,5),)*2,(Fraction(3,5),),(Fraction(12,13),))
        _,_,gram=grouped_gram(groups)
        columns=hypothesis_columns(groups)
        np.testing.assert_allclose(columns.T@columns,np.array(gram,dtype=float),atol=5e-16)

    def test_many_observations_of_one_label_keep_two_dimensional_span(self):
        for n in (1,3,20,100):
            row=capacity_certificate(((Fraction(4,5),)*n,))
            self.assertEqual(row["gram_rank_exact"],2)
            self.assertEqual(row["minimum_coherent_qubits"],1)

    def test_independent_nontrivial_labels_require_exponentially_larger_dimension(self):
        for q in range(1,6):
            row=capacity_certificate(tuple((Fraction(4,5),) for _ in range(q)))
            self.assertEqual(row["gram_rank_exact"],2**q)
            self.assertEqual(row["minimum_coherent_qubits"],q)
            self.assertGreater(Fraction(row["gram_determinant_exact"]),0)

    def test_groupwise_compressors_attain_the_dimension_bound_and_keep_coherences(self):
        groups=((Fraction(4,5),)*2,(Fraction(3,5),)*2)
        before=hypothesis_columns(groups)
        u=np.kron(compressor(.8,.8),compressor(.6,.6))
        after=u@before
        expected=[]
        for left,right in product((0,1),repeat=2):
            a=vector(.8**2) if left else np.array([1.,0])
            b=vector(.6**2) if right else np.array([1.,0])
            expected.append(np.kron(np.kron(np.kron(a,[1.,0]),b),[1.,0]))
        expected=np.column_stack(expected)
        np.testing.assert_allclose(after,expected,atol=6e-16)
        coefficients=np.array([1,1j,-.3,.2j])
        np.testing.assert_allclose(u@(before@coefficients),expected@coefficients,atol=8e-16)

    def test_using_one_relation_compressor_for_independent_relations_leaves_garbage(self):
        rows=off_promise_scratch_probabilities(.8,.8)
        self.assertLess(rows[0]["scratch_one_probability"],1e-28)
        self.assertLess(rows[3]["scratch_one_probability"],1e-28)
        self.assertGreater(rows[1]["scratch_one_probability"],.2)
        self.assertGreater(rows[2]["scratch_one_probability"],.2)
        self.assertEqual(capacity_certificate(((Fraction(4,5),),(Fraction(4,5),)))["gram_rank_exact"],4)

    def test_completely_uninformative_groups_do_not_count_as_new_dimensions(self):
        row=capacity_certificate(((Fraction(1),)*3,(Fraction(0),),(Fraction(4,5),)*2))
        self.assertEqual(row["gram_rank_exact"],4)
        self.assertEqual(row["gram_determinant_exact"],"0")
        self.assertEqual(row["minimum_coherent_qubits"],2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RelationalMemoryCapacityTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":96,"dimension_theorem":"d_min=rank(Gram)=2^q_active for the declared product family",
        "one_relation_repeated_three_times":capacity_certificate(((Fraction(4,5),)*3,)),
        "three_independent_relations":capacity_certificate(((Fraction(4,5),),)*3),
        "two_relations_each_repeated_twice":capacity_certificate(((Fraction(4,5),)*2,)*2),
        "single_relation_compressor_off_promise":off_promise_scratch_probabilities(.8,.8),
        "whole_disappears_during_ideal_compression":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("relational_memory_capacity_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
