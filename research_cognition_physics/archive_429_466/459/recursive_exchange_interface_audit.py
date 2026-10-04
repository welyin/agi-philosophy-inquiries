"""Round 459: exact recursive exchange interfaces and inherited internal readout.

The scientific baseline is 457. No result from the independent round 458 is
used. Collective-spin subsystem theory is reused; the new task is the lifting
of an actual continuous reader, with all private multiplicities and references
retained. This is a conditional interface construction, not a preparation,
complete cognitive architecture, spatial dimension, or GR derivation.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import exchange_relation_audit as old
import exchange_reader_audit as reader

TARGET = Path(__file__).with_name("recursive_exchange_interface_audit_results.json")
OBS = {}


def short(value):
    return float(f"{float(value):.12g}")


def kron_all(items):
    out = np.ones((1, 1), complex)
    for item in items:
        out = np.kron(out, item)
    return out


@lru_cache(None)
def permutation(n, a, b):
    indices = np.arange(2**n, dtype=np.int64)
    bit_a, bit_b = 1 << (n-1-a), 1 << (n-1-b)
    different = ((indices & bit_a) != 0) != ((indices & bit_b) != 0)
    return np.where(different, indices ^ (bit_a | bit_b), indices)


def action(n, edges, state, shift=0):
    result = -shift*state
    for a, b, weight in edges:
        result = result + weight*state[permutation(n, a, b)]
    return result


def collective(n, mu):
    return sum(kron_all([old.PAULI[mu] if a == b else np.eye(2)
                         for a in range(n)]) for b in range(n))/2


@lru_cache(None)
def half_encoding(n):
    """Full total-j=1/2 sector; column order G, private multiplicity.

    The n=3 basis is precisely round 430, not a new relabelling of its G/L.
    For other odd n the highest-weight kernel constructs every multiplicity.
    """
    assert n >= 1 and n % 2 == 1
    if n == 1:
        return np.eye(2, dtype=complex)
    if n == 3:
        return old.encoding()
    k = (n-1)//2
    current = [i for i in range(2**n) if i.bit_count() == k]
    higher = [i for i in range(2**n) if i.bit_count() == k-1]
    where = {value: row for row, value in enumerate(higher)}
    raising = np.zeros((len(higher), len(current)))
    for column, value in enumerate(current):
        for bit in range(n):
            if value & (1 << bit):
                raising[where[value ^ (1 << bit)], column] = 1
    _, singular, vh = np.linalg.svd(raising, full_matrices=True)
    rank = int(np.count_nonzero(singular > 1e-11))
    high = np.zeros((2**n, len(current)-rank), complex)
    high[current, :] = vh[rank:, :].T
    lowering = collective(n, 0)-1j*collective(n, 1)
    low = lowering@high
    return np.column_stack((high, low))


def grouped_embedding(sizes):
    """Columns all exposed G first, then all private M factors (including 1D)."""
    encodings = [half_encoding(n) for n in sizes]
    multiplicities = [v.shape[1]//2 for v in encodings]
    raw = kron_all(encodings)
    dimensions = list(itertools.chain.from_iterable((2, d) for d in multiplicities))
    order = list(range(0, 2*len(sizes), 2))+list(range(1, 2*len(sizes), 2))
    indices = np.arange(math.prod(dimensions)).reshape(dimensions).transpose(order).ravel()
    return raw[:, indices], multiplicities


def lifted_edges(sizes, effective_edges):
    starts = np.cumsum([0]+list(sizes))
    edges = [(i, j, weight) for a, b, weight in effective_edges
             for i in range(starts[a], starts[a+1])
             for j in range(starts[b], starts[b+1])]
    shift = sum(weight*(sizes[a]*sizes[b]-1)/2 for a, b, weight in effective_edges)
    return edges, shift


def raw_evolve(n, edges, state, time, shift, steps=32, degree=24):
    """Direct raw-space action, with a priori exponential-series norm error.

    The ignored shift contributes only an overall phase. No compression of
    the Hamiltonian is used in this numerical evolution.
    """
    dt = time/steps
    current = state.astype(complex).copy()
    for _ in range(steps):
        term = current.copy()
        result = current.copy()
        for k in range(1, degree+1):
            term = (-1j*dt/k)*action(n, edges, term, shift)
            result += term
        current = result
    bound = sum(abs(w) for _, _, w in edges)+abs(shift)
    x = Fraction(str(abs(time)))*Fraction(str(bound))/steps
    # e^x < 3^ceil(x), using e < 3. Error telescopes about exact unitaries.
    single = Fraction(3**math.ceil(x))*x**(degree+1)/math.factorial(degree+1)
    total = steps*single*(1+single)**(steps-1)
    return current, total


def total_singlet_projection(n):
    assert n % 2 == 0
    jsquared = sum(collective(n, mu)@collective(n, mu) for mu in range(3))
    projector = np.eye(2**n, dtype=complex)
    for j in range(1, n//2+1):
        projector = projector@(np.eye(2**n)-jsquared/(j*(j+1)))
    return projector


class Audit(unittest.TestCase):
    def close(self, left, right, tolerance=7e-11):
        self.assertLess(float(np.linalg.norm(left-right)), tolerance)

    def test_01_full_odd_sector_and_unknown_reference(self):
        rows = []
        for n in (1, 3, 5, 7):
            v = half_encoding(n)
            multiplicity = math.comb(n, (n-1)//2)-math.comb(n, (n-3)//2) if n > 1 else 1
            self.assertEqual(v.shape, (2**n, 2*multiplicity))
            self.close(v.conj().T@v, np.eye(2*multiplicity))
            errors = [np.linalg.norm(collective(n, mu)@v-v@np.kron(old.PAULI[mu]/2, np.eye(multiplicity))) for mu in range(3)]
            self.assertLess(max(errors), 1e-10)
            rows.append(dict(raw_qubits=n, private_dimension=multiplicity,
                             code_dimension=2*multiplicity))
        for n, m in ((1, 1), (3, 3), (3, 5), (5, 1), (7, 1)):
            w, ds = grouped_embedding((n, m))
            edges, constant = lifted_edges((n, m), [(0, 1, 1)])
            encoded = np.kron(old.swap(2, 0, 1)+constant*np.eye(4), np.eye(math.prod(ds)))
            self.close(action(n+m, edges, w), w@encoded)
        w, ds = grouped_embedding((3, 3))
        edges, constant = lifted_edges((3, 3), [(0, 1, 1)])
        hraw = action(6, edges, np.eye(64, dtype=complex))
        rng = np.random.default_rng(45901)
        psi = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
        psi /= np.linalg.norm(psi)
        unitary = old.evolve(hraw, .37)
        expected = np.exp(-1j*constant*.37)*np.kron(old.evolve(old.swap(2, 0, 1), .37), np.eye(4))@psi
        self.close(unitary@w@psi, w@expected)
        private_before = psi.reshape(4, 4, 3).transpose(1, 2, 0).reshape(12, 4)
        private_after = expected.reshape(4, 4, 3).transpose(1, 2, 0).reshape(12, 4)
        self.close(private_before@private_before.conj().T, private_after@private_after.conj().T)
        OBS['odd_sector_intertwining'] = dict(full_sectors=rows,
            independently_checked_pairs=[[1, 1], [3, 3], [3, 5], [5, 1], [7, 1]],
            unknown_reference_dimension=3, arbitrary_G_private_reference_correlations_allowed=True,
            formula='K_AB = ((n_A*n_B-1)/2) I + S_GA_GB, tensored with I_private',
            previous_three_by_three_formula_is_reused_not_new=True)

    def test_02_network_full_raw_intertwining_and_effects(self):
        sizes = (3, 3, 3, 1, 1)
        w, ds = grouped_embedding(sizes)
        h, effect, *_ = reader.setup()
        edges, constant = lifted_edges(sizes, [(i, i+1, 1) for i in range(4)])
        self.assertEqual(len(edges), 22)
        self.assertEqual(constant, 9)
        encoded = np.kron(h+constant*np.eye(32), np.eye(8))
        self.close(action(11, edges, w), w@encoded, 3e-10)
        e2 = (np.eye(4)-old.swap(2, 0, 1))/2
        actual_effect_w = np.einsum('ab,rbc->rac', e2, w.reshape(512, 4, 256)).reshape(2048, 256)
        self.close(actual_effect_w, w@np.kron(effect, np.eye(8)), 3e-10)
        # If the reader blocks are themselves three-qubit subjects, total
        # singlet projection of their six actual qubits is a legal effect.
        v, _ = grouped_embedding((3, 3))
        macro_effect = total_singlet_projection(6)
        self.close(macro_effect@macro_effect, macro_effect)
        self.close(macro_effect@v, v@np.kron(e2, np.eye(4)))
        self.assertAlmostEqual(float(np.trace(macro_effect).real), 5)
        OBS['network_and_actual_effect'] = dict(raw_block_sizes=list(sizes),
            microscopic_constant_edges=22, discarded_global_phase_generator=9,
            original_effective_dimension=32, private_dimension=8,
            full_raw_dimension=2048, full_encoded_columns_checked=256,
            compression_without_leakage_check_used=False,
            two_composite_reader_effect='total-j=0 projection on six actual qubits',
            six_qubit_effect_rank=5, both_exposed_qubits_and_actual_effect_lift=True)

    def test_03_actual_continuous_readout_with_private_memory(self):
        w, _ = grouped_embedding((3, 3, 3, 1, 1))
        edges, constant = lifted_edges((3, 3, 3, 1, 1), [(i, i+1, 1) for i in range(4)])
        h, effect, *_ = reader.setup()
        singlet = np.array([0, 1, -1, 0], complex)/math.sqrt(2)
        v3 = old.encoding()
        pin = np.kron(v3, singlet[:, None])
        small_u = old.evolve(h, 1.5)
        pulled = pin.conj().T@small_u.conj().T@effect@small_u@pin
        self.close(pulled, np.kron(np.eye(2), pulled[:2, :2]))
        rng = np.random.default_rng(45903)
        private_reference = rng.normal(size=(8, 2))+1j*rng.normal(size=(8, 2))
        private_reference /= np.linalg.norm(private_reference)
        probabilities = []
        for sign in (1, -1):
            logical = np.array([1, sign*1j])/math.sqrt(2)
            source = v3@np.kron([1, 0], logical)
            effective = np.kron(source, singlet)
            initial = w@np.kron(effective[:, None], private_reference)
            actual, error = raw_evolve(11, edges, initial, 1.5, constant)
            expected = w@np.kron((small_u@effective)[:, None], private_reference)
            self.close(actual, expected, 2e-11)
            reshaped = actual.reshape(512, 4, 2)
            reader_amplitude = np.einsum('b,abc->ac', singlet.conj(), reshaped)
            probability = float(np.linalg.norm(reader_amplitude)**2)
            ideal = float(np.vdot(small_u@effective, effect@small_u@effective).real)
            self.assertAlmostEqual(probability, ideal, places=11)
            probabilities.append(short(probability))
        inherited = json.loads(Path(reader.TARGET).read_text(encoding='utf-8'))['observations']['exact_transfer_certificate']
        inherited_low = Fraction(inherited['contrast_lower_rational'])
        self.assertGreater(inherited_low, Fraction(2, 3))
        self.assertEqual(Fraction(2, 3)-2*Fraction(1, 12), Fraction(1, 2))
        self.assertGreater(probabilities[1]-probabilities[0], 2/3)
        self.assertLess(error, Fraction(1, 10**17))
        OBS['inherited_readout'] = dict(model_time='3/2', actual_raw_probabilities=probabilities,
            raw_contrast=short(probabilities[1]-probabilities[0]),
            all_old_private_qubits_retained=3, numerical_reference_dimension=2,
            generalized_unknown_reference_covered_by_operator_intertwining=True,
            gauge_independent_source_readout=True,
            source_readout_certificate_round=431,
            source_results_sha256=hashlib.sha256(Path(reader.TARGET).read_bytes()).hexdigest(),
            inherited_contrast_strict_lower='2/3', inherited_window=['17/12', '19/12'],
            inherited_window_contrast_strict_lower='1/2',
            raw_evolution_steps=32, raw_evolution_series_degree=24,
            raw_series_vector_error_upper=short(float(error)),
            round431_numeric_signal_is_reused_not_new=True,
            initial_source_and_singlet_resources_are_inputs=True,
            complete_readout_instrument_or_singlet_supply_constructed=False)

    def test_04_two_level_recursive_code_and_private_role_boundary(self):
        children, _ = grouped_embedding((3, 3, 3))
        v = children@np.kron(old.encoding(), np.eye(8))
        self.assertEqual(v.shape, (512, 32))
        self.close(v.conj().T@v, np.eye(32))
        for mu in range(3):
            self.close(collective(9, mu)@v, v@np.kron(old.PAULI[mu]/2, np.eye(16)), 1e-10)
        # Fixed internal operations at each level decompose into one local
        # operator on each private tree-node L. They do not couple those Ls.
        upper_edges, constant = lifted_edges((3, 3, 3), [(0, 1, 1), (1, 2, 2), (0, 2, 3)])
        actual_edges = upper_edges+[(0, 1, 1), (3, 4, 2), (6, 7, 3)]
        parent_h = -math.sqrt(3)/2*old.PAULI[0]+1.5*old.PAULI[2]
        children_h = -sum((i+1)*kron_all([old.PAULI[2] if i == j else np.eye(2) for j in range(3)]) for i in range(3))
        private_h = np.kron(parent_h, np.eye(8))+np.kron(np.eye(2), children_h)
        expected = constant*np.eye(32)+np.kron(np.eye(2), private_h)
        self.close(action(9, actual_edges, v), v@expected, 3e-10)
        # Cross-interface to one new qubit acts on the root G alone; all
        # sixteen private dimensions remain, not just selected private states.
        joint = np.kron(v, np.eye(2))
        cross = [(a, 9, 1) for a in range(9)]
        effective_cross = 4.5*np.eye(64)+sum(np.kron(np.kron(pauli, np.eye(16)), pauli)/2 for pauli in old.PAULI)
        self.close(action(10, cross, joint), joint@effective_cross, 3e-10)
        OBS['recursive_two_level_interface'] = dict(raw_qubits=9,
            root_exposed_dimension=2, root_private_dimension=16,
            full_code_dimension=32, new_parent_private_qubits=1,
            inherited_child_private_qubits=3,
            arbitrary_private_states_and_references_preserved_by_external_rule=True,
            private_tree_node_dynamics_factorize_under_this_hierarchy=True,
            private_memory_drives_exposed_exchange_under_this_hierarchy=False,
            current_L_communication_requires_access_to_constituent_G_interfaces=True,
            all_ancestor_codes_claimed_preserved_during_readout=False)

    def test_05_unknown_total_spin_is_not_the_same_qubit_port(self):
        edges = [(a, 3, 1) for a in range(3)]
        h = action(4, edges, np.eye(16, dtype=int))
        identity = np.eye(16, dtype=int)
        self.close((h+identity)@h@(h-2*identity)@(h-3*identity), np.zeros((16, 16)), 1e-15)
        values, counts = np.unique(np.round(np.linalg.eigvalsh(h), 10), return_counts=True)
        self.assertEqual(values.tolist(), [-1, 0, 2, 3])
        self.assertEqual(counts.tolist(), [3, 2, 6, 5])
        u = old.evolve(h, math.pi/2)
        low = np.kron(old.encoding()[:, 0], [0, 1])
        high = np.eye(16)[:, 1]  # |000> in j=3/2, partner |1>
        target_one = np.kron(np.eye(8), np.diag([0, 1]))
        low_probability = float(np.vdot(u@low, target_one@u@low).real)
        high_probability = float(np.vdot(u@high, target_one@u@high).real)
        self.assertAlmostEqual(low_probability, 0)
        self.assertAlmostEqual(high_probability, 1)
        # On the entire j=3/2 block plus partner, U is just i I.
        pa = np.eye(8)-old.operators()[3]
        high_projector = np.kron(pa, np.eye(2))
        self.close(u@high_projector, 1j*high_projector)
        OBS['unprepared_sector_boundary'] = dict(raw_cross_generator_spectrum=[-1, 0, 2, 3],
            degeneracies=[3, 2, 6, 5], time='pi/2',
            partner_one_probability_in_prepared_half_sector=0,
            partner_one_probability_in_highest_weight_three_halves_sector=1,
            whole_three_halves_sector_evolution='i I',
            all_odd_block_states_claimed_to_define_qubit_port=False,
            collective_J_squared_sector_is_required_input=True)

    def test_06_recursive_capacity_and_explicit_local_budget(self):
        rows = []
        private_qubits = 0
        for depth in range(6):
            n = 3**depth
            self.assertEqual(private_qubits, (n-1)//2)
            code_dimension = 2**(private_qubits+1)
            rank_fraction = Fraction(code_dimension, 2**n)
            # For the five-block chain and a per-raw-site absolute coefficient
            # budget kappa=1, equal microscopic weights obey 2*n*J <= 1.
            largest_uniform_strength = Fraction(1, 2*n)
            read_time = Fraction(3, 2)/largest_uniform_strength
            self.assertEqual(read_time, 3*n)
            rows.append(dict(depth=depth, raw_qubits_per_block=n,
                private_qubits_per_block=private_qubits,
                total_code_dimension=str(code_dimension),
                hierarchical_code_fraction_of_raw_space=str(rank_fraction),
                raw_pairs_per_isolated_contact=n*n,
                five_block_chain_raw_pairs=4*n*n,
                max_uniform_coefficient_for_kappa_one=str(largest_uniform_strength),
                model_read_time_for_kappa_one=str(read_time)))
            private_qubits = 1+3*private_qubits
        OBS['scale_and_resource_contract'] = dict(rows=rows,
            recursion='q_(l+1)=1+3*q_l; n_l=3^l; q_l=(n_l-1)/2',
            general_raw_site_coefficient_load='sum_b n_b*abs(J_ab)',
            isolated_pair_exchange_time_lower='abs(theta)*max(n_A,n_B)/kappa',
            equal_scale_reader_time_at_max_allowed_strength='3*n/kappa',
            budget_is_explicit_additional_model_input=True,
            coefficient_budget_identified_with_thermodynamic_work=False,
            finite_depth_only=True, autonomous_code_preparation_claimed=False,
            physical_fractal_dimension_inferred_from_branching=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=459, baseline_round=457, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(exact_recursive_exchange_port_intertwining=True,
            actual_existing_continuous_reader_lifted=True,
            private_multiplicities_and_unknown_references_retained=True,
            collective_spin_and_noiseless_subsystem_theory_are_prior_tools=True,
            new_signal_value_or_new_exchange_encoding_discovery_claimed=False,
            G_and_old_L_roles_identified=False,
            full_recursive_cognitive_functional_closure_proved=False,
            arbitrary_unprepared_raw_states_accepted=False,
            internal_preparation_and_contact_pattern_derived=False,
            independent_from_round_458=True,
            physical_space_dimension_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
