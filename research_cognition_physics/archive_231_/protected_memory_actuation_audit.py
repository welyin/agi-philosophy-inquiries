"""Round 460: minimum raw locality for a protected three-spin memory actuator.

The higher-body generators here are additional model inputs, not derived
from round 429. Completeness uses invariant Pauli tensors and exact integer
commutator Gram matrices, not selected coupling scans. NumPy only.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as old

TARGET = Path(__file__).with_name('protected_memory_actuation_audit_results.json')
OBS = {}
PAULI = [np.eye(2)] + old.PAULI
EDGES = list(itertools.combinations(range(4), 2))
TRIPLES = list(itertools.combinations(range(4), 3))


def kron(*items):
    out = np.ones((1, 1), complex)
    for a in items:
        out = np.kron(out, a)
    return out


def exact_rank(a):
    rows = [[Fraction(int(x)) for x in row] for row in a]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((k for k in range(rank, len(rows)) if rows[k][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [v/scale for v in rows[rank]]
        for k in range(len(rows)):
            if k != rank:
                scale = rows[k][column]
                rows[k] = [x-scale*y for x, y in zip(rows[k], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def operators():
    swaps = {e: old.swap(4, *e) for e in EDGES}
    chis = {(a,b,c): 1j*(swaps[a,b]@swaps[b,c]-swaps[b,c]@swaps[a,b])
            for a,b,c in TRIPLES}
    i = np.eye(16)
    p3 = 3*i-swaps[0,1]-swaps[0,2]-swaps[1,2]
    k = swaps[0,3]+swaps[1,3]+swaps[2,3]
    c = chis[0,1,3]-chis[0,2,3]+chis[1,2,3]
    h3 = k+c/(2*math.sqrt(3))
    record = (i-swaps[0,1])/2
    h4 = record@(i+swaps[2,3])
    return swaps, chis, p3, k, c, h3, h4


def gram(items):
    out = np.array([[np.vdot(a,b) for b in items] for a in items])
    assert np.max(np.abs(out.imag)) == 0
    assert np.max(np.abs(out.real-np.rint(out.real))) == 0
    return out.real.astype(np.int64)


def word_expansion(a):
    out = {}
    for word in itertools.product(range(4), repeat=4):
        value = np.trace(kron(*(PAULI[k] for k in word))@a)/16
        if abs(value) > 1e-11:
            assert abs(value.imag) < 1e-11
            out[''.join('IXYZ'[k] for k in word)] = float(value.real)
    return out


def input_embedding():
    singlet = np.array([0,1,-1,0], complex)/math.sqrt(2)
    return np.kron(old.encoding(), singlet[:,None])


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=4e-12):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_complete_invariant_local_basis(self):
        swaps, chis, _, _, _, _, _ = operators()
        basis = [np.eye(16)]+list(swaps.values())+list(chis.values())
        self.assertEqual(exact_rank(gram(basis)), 11)
        for mu in range(1,4):
            total = sum(kron(*(PAULI[mu] if a == b else PAULI[0]
                               for a in range(4))) for b in range(4))
            for h in basis:
                self.close(h@total, total@h)
        for h in basis:
            self.assertLessEqual(max(sum(x!='I' for x in w)
                                     for w in word_expansion(h)), 3)
        # Pauli tensor invariants: weight 0 scalar, no weight 1 invariant,
        # weight 2 delta_ab, weight 3 epsilon_abc; one per support.
        OBS['complete_basis'] = dict(dimension=11, exact_gram_rank=11,
            weight_dimensions=[1,0,6,4],
            completeness_source='SU(2) invariant tensors of 0,1,2,3 Pauli vectors',
            locality_is_raw_Pauli_support=True, auxiliary_carriers_allowed=False)

    def test_02_exact_complete_preservation_classification(self):
        swaps, chis, p3, _, _, _, _ = operators()
        cross = [swaps[a,3] for a in range(3)]
        chiral = [chis[t] for t in [(0,1,3),(0,2,3),(1,2,3)]]
        even_comm = [h@p3-p3@h for h in cross]
        odd_comm = [h@p3-p3@h for h in chiral]
        g2, g3 = gram(even_comm), gram(odd_comm)
        self.assertTrue(np.array_equal(g2, 24*np.array([[2,-1,-1],[-1,2,-1],[-1,-1,2]])))
        self.assertTrue(np.array_equal(g3, 72*np.array([[2,1,-1],[1,2,1],[-1,1,2]])))
        self.assertEqual(exact_rank(g2),2)
        self.assertEqual(exact_rank(g3),2)
        self.assertTrue(np.array_equal(g2@np.ones(3,dtype=int),np.zeros(3,dtype=int)))
        self.assertTrue(np.array_equal(g3@np.array([1,-1,1]),np.zeros(3,dtype=int)))
        # Hermitian coefficients are real. Real and imaginary commutator
        # entries cannot cancel; the real quadratic cross term is zero.
        for a in even_comm:
            for b in odd_comm:
                self.assertEqual(float(np.vdot(a,b).real),0)
        for key in [(0,1),(0,2),(1,2)]:
            self.close(swaps[key]@p3,p3@swaps[key])
        self.close(chis[0,1,2]@p3,p3@chis[0,1,2])
        OBS['preservation_classification'] = dict(pair_integer_gram=g2.tolist(),
            chiral_integer_gram=g3.tolist(), projector_scale=3,
            cross_pair_kernel=[1,1,1], cross_chiral_kernel=[1,-1,1],
            exact_ranks=[2,2], arbitrary_real_coefficients_covered=True,
            raw_two_body_memory_control_excluded_under_fixed_code=True)

    def test_03_actual_chiral_actuator_intertwining(self):
        swaps, chis, p3, k, c, h3, _ = operators()
        v = np.kron(old.encoding(),np.eye(2))  # G,L,B column order
        expected_c = math.sqrt(3)*sum(kron(q,old.PAULI[1],q) for q in old.PAULI)
        self.close(c@v,v@expected_c)
        self.close(k@v,v@(np.eye(8)+old.swap(3,0,2)))
        y = kron(np.eye(2),old.PAULI[1],np.eye(2))
        sgb = old.swap(3,0,2)
        expected_h = np.eye(8)-y/2+(np.eye(8)+y)@sgb
        self.close(h3@v,v@expected_h)
        self.close(h3@p3,p3@h3)
        y_raw = chis[0,1,2]/math.sqrt(3)
        self.close(h3@y_raw,y_raw@h3)
        rng = np.random.default_rng(46003)
        rho = old.density(8*3,rng)
        for t in [0.13,math.pi/4,1.1]:
            u, ueff = old.evolve(h3,t),old.evolve(expected_h,t)
            self.close(u@v,v@ueff)
            vr = np.kron(v,np.eye(3))
            ur = np.kron(u,np.eye(3))
            expected = vr@np.kron(ueff,np.eye(3))@rho@np.kron(ueff.conj().T,np.eye(3))@vr.conj().T
            self.close(ur@vr@rho@vr.conj().T@ur.conj().T,expected)
        OBS['chiral_actuator'] = dict(raw_qubits=4, code_dimension=8,
            generator='K + (chi013 - chi023 + chi123)/(2 sqrt(3))',
            effective_generator='I - Y_L/2 + (I+Y_L) S_GB',
            exact_all_time_code_and_Y_population_preservation=True,
            unknown_G_L_B_reference_scope=True, numerical_reference_dimension=3,
            raw_three_body_chirality_is_additional_input=True)

    def test_04_reference_complete_actual_response(self):
        _,_,_,_,_,h3,h4 = operators()
        v = input_embedding()  # source GL -> raw A,B,C with supplied BC singlet
        effect = (np.eye(32)-old.swap(5,3,4))/2
        records = [(h3,old.PAULI[1],2,math.pi/4),
                   (h4,old.PAULI[2],1,math.pi/2)]
        rows=[]
        for h, axis, rate, center in records:
            for t in [0.19,center,center+0.07]:
                w=old.evolve(np.kron(h,np.eye(2)),t)@v
                fl=np.eye(2)-0.75*math.sin(rate*t)**2*(np.eye(2)+axis)/2
                self.close(w.conj().T@effect@w,np.kron(np.eye(2),fl))
                self.close(w.conj().T@w,np.eye(4))
            rows.append(dict(generator='H3' if rate==2 else 'H4',
                center='pi/4' if rate==2 else 'pi/2',
                inactive_probability=1,active_probability=0.25,contrast=0.75))
        # Effect identity covers every GLR, including unknown correlated G.
        rng=np.random.default_rng(46004)
        rho=old.density(12,rng)
        w=old.evolve(np.kron(h3,np.eye(2)),math.pi/4)@v
        wr=np.kron(w,np.eye(3))
        out=wr@rho@wr.conj().T
        e=np.kron(effect,np.eye(3))
        actual=old.partial(e@out@e,(32,3),(1,))
        fl=np.eye(2)-.75*(np.eye(2)+old.PAULI[1])/2
        desired=old.partial(np.kron(np.kron(np.eye(2),fl),np.eye(3))@rho,
                            (2,2,3),(2,))
        self.close(actual,desired)
        OBS['actual_response'] = dict(rows=rows,
            full_unknown_GL_reference_effect_identity=True,
            effect='singlet of actual B and internal reference C',
            independent_pure_BC_resource_is_input=True,
            single_B_output_claimed_gauge_independent=False)

    def test_05_time_reversal_and_even_four_body_alternative(self):
        swaps,_,p3,k,c,h3,h4 = operators()
        record=(np.eye(16)-swaps[0,1])/2
        self.close(h4,record@k)
        self.close(h4@p3,p3@h4)
        self.close(h4@record,record@h4)
        v=np.kron(old.encoding(),np.eye(2))
        target=kron(np.eye(2),np.diag([1,0]),np.eye(2))@(np.eye(8)+old.swap(3,0,2))
        self.close(h4@v,v@target)
        theta_unitary=kron(*(1j*old.PAULI[1] for _ in range(4)))
        self.close(theta_unitary@c.conj()@theta_unitary.conj().T,-c)
        self.close(theta_unitary@h4.conj()@theta_unitary.conj().T,h4)
        self.assertGreater(np.linalg.norm(theta_unitary@h3.conj()@theta_unitary.conj().T-h3),1)
        w3,w4=word_expansion(h3),word_expansion(h4)
        self.assertEqual(max(sum(x!='I' for x in w) for w in w3),3)
        self.assertEqual(max(sum(x!='I' for x in w) for w in w4),4)
        self.assertAlmostEqual(w4['ZZZZ'],-1/8)
        OBS['locality_and_time_reversal'] = dict(
            minimum_raw_locality_without_T_restriction=3,
            minimum_raw_locality_with_T_even_restriction=4,
            H4='(I-S01)(I+S23)/2', H4_ZZZZ_coefficient='-1/8',
            H3_chiral_component_T_odd=True,H4_T_even=True,
            T_even_is_additional_comparison_not_cognitive_axiom=True,
            auxiliary_or_effective_or_time_dependent_implementations_not_excluded=True)

    def test_06_window_and_quantum_memory_disturbance(self):
        _,_,_,_,_,h3,h4=operators()
        gap_bound=Fraction(3,4)*(1-Fraction(1,4)**2)
        self.assertEqual(gap_bound,Fraction(45,64))
        self.assertEqual((1+gap_bound)/2,Fraction(109,128))
        v=input_embedding()
        # L=|0_Z> is a coherent superposition of the H3 Y records.
        rho=np.kron(np.eye(2)/2,np.diag([1,0]))
        out=old.evolve(np.kron(h3,np.eye(2)),math.pi/4)@v
        code=np.kron(old.encoding().conj().T,np.eye(4))@out
        evolved=code@rho@code.conj().T
        reduced=old.partial(evolved,(2,2,2,2),(1,))
        self.close(np.sort(np.linalg.eigvalsh(reduced)),np.array([.25,.75]))
        self.assertAlmostEqual(float(np.trace(reduced@reduced).real),5/8)
        for center,width,rate in [(math.pi/4,1/8,2),(math.pi/2,1/4,1)]:
            for delta in [-width,0,width]:
                self.assertGreaterEqual(.75*math.sin(rate*(center+delta))**2,float(gap_bound)-1e-14)
        OBS['finite_window_and_memory_scope'] = dict(
            H3_window_center='pi/4',H3_window_half_width='1/8',
            H4_window_center='pi/2',H4_window_half_width='1/4',
            analytic_gap_lower='45/64',equal_prior_success_lower='109/128',
            proof='sin(x)^2 <= x^2, so cos(rate*delta)^2 >= 15/16',
            inactive_branch_record_population_preserved=True,
            all_quantum_memory_marginals_preserved=False,
            coherent_record_output_eigenvalues=[.25,.75],
            coherent_record_output_purity='5/8',
            full_unknown_information_preserved_in_global_isometry=True)


def run():
    OBS.clear()
    output=io.StringIO()
    result=unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=460,baseline_round=459,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(complete_three_local_invariant_classification=True,
            protected_private_memory_actuation_exact=True,
            reference_complete_actual_response=True,
            higher_body_input_derived_from_429=False,
            arbitrary_quantum_memory_undisturbed=False,
            full_recursive_cognitive_functional_closure_proved=False,
            spatial_dimension_derived=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
