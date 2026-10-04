"""Round 450: compositional agreement versus fixed-state pair contracts.
No new cognitive axiom or physical source is inferred from a weaker contract.
"""
import argparse
from functools import lru_cache
import io
import itertools as it
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import agreement_exchange_audit as old
import internal_partner_response_audit as local_record
import joint_marker_dynamics_audit as complete_record

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'composite_agreement_audit_results.json'
OBS = {}


def short(x):
    return float(f'{float(x):.12g}')


def hermitian_basis(d):
    out = []
    for i in range(d):
        a = np.zeros((d, d), dtype=complex)
        a[i, i] = 1
        out.append(a)
    for i, j in it.combinations(range(d), 2):
        a = np.zeros((d, d), dtype=complex)
        a[i, j] = a[j, i] = 1
        out.append(a)
        b = np.zeros_like(a)
        b[i, j], b[j, i] = -1j, 1j
        out.append(b)
    return out


def phi(hodd, x):
    d = x.shape[0]
    return old.marginals(hodd@np.kron(np.eye(d), x), d)[0]


def derivative_difference(h, x):
    joint = np.kron(x, x)
    a, b = old.marginals(-1j*(h@joint-joint@h), x.shape[0])
    return a-b


@lru_cache(None)
def derivative_constraints(d):
    local, domain = hermitian_basis(d), hermitian_basis(d*d)
    probes = local+[a+b for i, a in enumerate(local) for b in local[i+1:]]
    matrix = np.vstack([np.column_stack([derivative_difference(h, x).reshape(-1)
                                        for h in domain]) for x in probes])
    realified = np.vstack((matrix.real, matrix.imag)).astype(np.int64)
    return local, domain, probes, realified


def swap_sites(dimensions, first, second):
    assert dimensions[first] == dimensions[second]
    words = tuple(it.product(*(range(d) for d in dimensions)))
    lookup = {w: i for i, w in enumerate(words)}
    out = np.zeros((len(words), len(words)), dtype=int)
    for col, word in enumerate(words):
        v = list(word)
        v[first], v[second] = v[second], v[first]
        out[lookup[tuple(v)], col] = 1
    return out


def product(operators):
    out = np.array([[1]], dtype=complex)
    for op in operators:
        out = np.kron(out, op)
    return out


class Audit(unittest.TestCase):
    def test_01_exact_weak_contract_kernel(self):
        reports = []
        for d in (2, 3):
            _, domain, probes, constraint = derivative_constraints(d)
            rank, pivots = old.modular_rank(constraint.astype(complex))
            expected = (d**4-d*d)//2
            self.assertEqual(rank, expected)
            s = old.swap(d)
            # All symmetric Hamiltonians are in the derivative kernel exactly.
            for h in domain:
                even = h+s@h@s
                for x in probes:
                    np.testing.assert_array_equal(derivative_difference(even, x), np.zeros((d, d)))
            # Constraints for [H,S]=0 independently have the same codimension.
            c = np.column_stack([(h@s-s@h).reshape(-1) for h in domain])
            c = np.vstack((c.real, c.imag)).astype(np.int64)
            swap_rank, _ = old.modular_rank(c.astype(complex))
            combined_rank, _ = old.modular_rank(np.vstack((constraint, c)).astype(complex))
            self.assertEqual((swap_rank, combined_rank), (expected, expected))
            reports.append(dict(d=d, Hermitian_unknowns=d**4,
                independent_derivative_rank_mod1009=rank,
                swap_commutator_rank_mod1009=swap_rank, combined_rank_mod1009=combined_rank,
                weak_contract_real_dimension=(d**4+d*d)//2,
                strong_429_real_dimension=2, nonzero_pivots=pivots))
        OBS['weak_contract_classification'] = dict(exact_finite_certificates=reports,
            all_dimensions_proof_analytic=True,
            condition='all identical independent inputs retain equal marginals iff [H,SWAP]=0',
            time_zero_derivative_equality_already_sufficient=True,
            closure_and_time_independence_required=True)

    def test_02_antisymmetric_map_derivative_identity(self):
        reports = []
        for d in (2, 3):
            basis = hermitian_basis(d)
            rng = np.random.default_rng(450+d)
            raw = rng.integers(-3, 4, size=(d*d, d*d))+1j*rng.integers(-3, 4, size=(d*d, d*d))
            h = raw+raw.conj().T
            s = old.swap(d)
            a = (h-s@h@s)/2
            for x in basis:
                np.testing.assert_array_equal(phi(a, x), phi(a, x).conj().T)
                for y in basis:
                    self.assertEqual(np.trace(y@phi(a, x)), -np.trace(x@phi(a, y)))
            for x in basis+[u+v for i,u in enumerate(basis) for v in basis[i+1:]]:
                p = phi(a, x)
                np.testing.assert_array_equal(derivative_difference(h, x), -2j*(p@x-x@p))
            reports.append(dict(d=d, integer_derivative_and_skew_pairings_exact=True))
        OBS['necessity_identity'] = reports

    def test_03_weak_agreement_allows_joint_update(self):
        q = np.diag([2, 1, 1, 2])
        s = old.swap(2)
        np.testing.assert_array_equal(q@s, s@q)
        plus = np.full((2, 2), 0.5)
        u = old.unitary(-q, np.pi/2)
        a, b = old.marginals(u@np.kron(plus, plus)@u.conj().T, 2)
        np.testing.assert_allclose(a, np.eye(2)/2, atol=1e-14)
        np.testing.assert_allclose(b, a, atol=1e-14)
        self.assertAlmostEqual(old.distance(a, plus), 0.5)
        worst = 0.0
        for d in (2, 3):
            rng = np.random.default_rng(45030+d)
            raw = rng.normal(size=(d*d,d*d))+1j*rng.normal(size=(d*d,d*d))
            h = raw+raw.conj().T
            sw = old.swap(d)
            h = h+sw@h@sw
            u = old.unitary(h, 0.31)
            for _ in range(3):
                rho = old.state(d, rng)
                first, second = old.marginals(u@np.kron(rho, rho)@u.conj().T, d)
                worst = max(worst, float(np.linalg.norm(first-second)))
        self.assertLess(worst, 1e-12)
        OBS['weak_example'] = dict(symmetric_diagonal_pair_reward_admitted=True,
            each_original_marginal_trace_distance=0.5,
            marginal_equality_error=short(worst),
            unique_partial_SWAP_rule_selected=False,
            weak_condition_adopted_as_cognitive_axiom=False)

    def test_04_composite_internal_correlations_and_microscopic_exchange(self):
        dims = (2,2,2,2)  # A1,A2,B1,B2
        s1, s2 = swap_sites(dims, 0, 2), swap_sites(dims, 1, 3)
        block = old.swap(4)
        np.testing.assert_array_equal(s1@s2, block)
        h = s1+s2
        np.testing.assert_array_equal(h@block, block@h)
        bell = np.array([1,0,0,1])/np.sqrt(2)
        rho = np.outer(bell, bell)
        # Exact dyadic input and the exact t=pi/4 primitive product unitary.
        psi = np.zeros(16)
        psi[[0,3,12,15]] = 0.5
        chi = s1@psi
        np.testing.assert_array_equal(s1@psi, s2@psi)
        np.testing.assert_array_equal(block@psi, psi)
        cross, _ = old.marginals(np.outer(psi, chi), 4)
        rho_exact = np.zeros((4,4))
        rho_exact[np.ix_([0,3],[0,3])] = 0.5
        np.testing.assert_array_equal(cross, rho_exact/2)
        u_quarter = (np.eye(16)-1j*s1)@(np.eye(16)-1j*s2)/2
        out = u_quarter@psi
        joint = np.outer(out, out.conj())
        a, b = old.marginals(joint, 4)
        np.testing.assert_array_equal(a, np.eye(4)/4)
        np.testing.assert_array_equal(b, a)
        self.assertEqual(old.distance(a, rho_exact), 0.75)
        worst = 0.0
        for t in (0.173, np.pi/8, np.pi/2):
            u = old.unitary(h, t)
            out = u@psi
            a, b = old.marginals(np.outer(out, out.conj()), 4)
            target = np.cos(2*t)**2*rho_exact+np.sin(2*t)**2*np.eye(4)/4
            worst = max(worst, float(np.linalg.norm(a-target)), float(np.linalg.norm(b-target)))
            for state in (a, b):
                for one in old.marginals(state, 2):
                    np.testing.assert_allclose(one, np.eye(2)/2, atol=2e-14)
        # Unknown internally correlated mixed complete states, no product-over-components promise.
        rng = np.random.default_rng(45004)
        rho_unknown = old.state(4, rng)
        u = old.unitary(0.7*s1-0.2*s2, 0.6)
        a, b = old.marginals(u@np.kron(rho_unknown,rho_unknown)@u.conj().T, 4)
        np.testing.assert_allclose(a, b, atol=2e-14)
        OBS['composition'] = dict(subjects=2, micro_qubits_per_subject=2,
            primitive_terms_each_429_partial_SWAP=True,
            complete_subject_initial_state_independent_Bell_copies=True,
            reduced_formula='cos^2(2*t)*Bell + sin^2(2*t)*I4/4',
            exact_t_pi_over_4_old_state_trace_distance='3/4',
            complete_subjects_remain_equal=True,
            all_single_micro_qubit_marginals_unchanged=True,
            initial_cross_subject_correlations_used=False,
            full_state_not_recovered_from_single_component_marginals=True,
            numerical_formula_error=short(worst),
            general_component_exchange_weak_contract_closed=True,
            arbitrary_initial_shared_reference_assumed_independent_copies=False)

    def test_05_microscopic_support_boundary_for_exact_strong_contract(self):
        reports = []
        for local_dims in ((2,2), (2,3), (2,2,2)):
            dimensions = local_dims+local_dims
            k, d = len(local_dims), int(np.prod(local_dims))
            swaps = [swap_sites(dimensions, i, k+i) for i in range(k)]
            full = np.eye(d*d, dtype=int)
            for sw in swaps:
                full = full@sw
            np.testing.assert_array_equal(full, old.swap(d))
            traceless = [np.diag([1,-1]+[0]*(n-2)) for n in local_dims]
            witness = product(traceless+traceless)
            exact_trace = np.trace(witness@full)
            self.assertEqual(exact_trace, 2**k)
            # Every elementary pair term lacks at least one traceless witness factor.
            for sw in swaps:
                self.assertEqual(np.trace(witness@sw), 0)
            self.assertEqual(np.trace(witness), 0)
            reports.append(dict(local_dimensions=list(local_dims), total_micro_factors=2*k,
                trace_witness_full_SWAP=int(exact_trace.real),
                trace_witness_identity_and_each_primitive_SWAP=0))
        OBS['strong_support_boundary'] = dict(exact_tensor_witnesses=reports,
            conclusion='r<2*k raw microscopic terms plus exact strong full-composite contract imply scalar H',
            all_time_static_raw_generator_scope=True,
            endpoint_SWAP_or_encoded_auxiliary_simulation_excluded=False,
            finite_precision_effective_interactions_excluded=False)

    def test_06_raw_record_potential_and_comparison_dictionary(self):
        local_words = tuple(it.product((0,1,-1),repeat=2))
        lookup = {w:i for i,w in enumerate(local_words)}
        dim = len(local_words)
        global_words = tuple(x+y for x in local_words for y in local_words)
        record = np.diag([complete_record.record_energy(w) for w in global_words])
        active = np.diag([local_record.raw_energy(w) for w in global_words])
        s = old.swap(dim)
        r1 = lookup[(1,-1)]
        r0 = lookup[(0,-1)]
        initial = np.eye(dim*dim)[:, r1*dim+r0]
        self.assertEqual(float(np.vdot(initial, record@initial).real), 0)
        self.assertEqual(float(np.vdot(s@initial, record@(s@initial)).real), 2)
        self.assertGreater(np.linalg.norm(record@s-s@record), 1)
        plus = np.zeros(dim)
        plus[r0] = plus[r1] = 1/np.sqrt(2)
        rho = np.outer(plus,plus)
        u = old.unitary(record, np.pi/4)
        a,b = old.marginals(u@np.kron(rho,rho)@u.conj().T, dim)
        self.assertAlmostEqual(old.distance(a,b), 0.5)
        # Change the fixed identification: relabel 0<->1 inside both registers on one side.
        f = np.zeros((dim,dim),dtype=int)
        for col, w in enumerate(local_words):
            renamed = tuple(1-v if v>=0 else -1 for v in w)
            f[lookup[renamed],col] = 1
        twisted_swap = np.kron(f,f)@s
        np.testing.assert_array_equal(record@twisted_swap, twisted_swap@record)
        np.testing.assert_array_equal(active@twisted_swap, twisted_swap@active)
        local_swap = np.zeros((dim,dim),dtype=int)
        for col, w in enumerate(local_words):
            local_swap[lookup[w[::-1]],col] = 1
        local = np.kron(local_swap,np.eye(dim))+np.kron(np.eye(dim),local_swap)
        full = 3*record+active+local+s
        np.testing.assert_array_equal(full@twisted_swap, twisted_swap@full)
        np.testing.assert_allclose(a, f@b@f.T, atol=2e-14)
        OBS['record_interface'] = dict(raw_local_dimension=9, raw_pair_dimension=81,
            existing_fixed_numeric_dictionary_weak_contract_satisfied=False,
            h_R_energy_before_packet_SWAP=0, h_R_energy_after_packet_SWAP=2,
            same_input_numeric_marginal_difference_at_pi_over_4='1/2',
            changed_identification_F_swaps_0_and_1=True,
            same_raw_model_commutes_with_twisted_SWAP=True,
            changed_comparison_marginals_equal=True,
            global_multi_subject_identification_consistency_proved=False,
            passive_relabeling_confused_with_physical_SWAP=False,
            comparison_dictionary_origin_still_input=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=450, baseline_round=449, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(weak_agreement_Hamiltonian_iff_classification_proved=True,
            microscopic_exchange_weak_contract_closed_under_composition=True,
            strong_contract_failure_for_independent_composite_Bell_inputs_proved=True,
            exact_all_time_raw_support_boundary_explicit=True,
            raw_record_comparison_dictionary_obstruction_and_two_subject_alternative_verified=True,
            internal_and_cross_subject_correlations_distinguished=True,
            fixed_identification_source_not_erased=True,
            weak_agreement_adopted_as_new_cognitive_axiom=False,
            strong_429_theorem_invalidated=False,
            new_record_potentials_derived_from_429=False,
            all_encodings_or_effective_low_body_implementations_excluded=False,
            multi_subject_reference_composition_selected=False,
            full_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))
