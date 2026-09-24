"""Round 433: operator-valued participation with explicit new primitives.

Inspired by dynamic graph approaches, this is a finite audit model, not a
reproduction of their Hamiltonians or a derivation from exchange alone.
"""
import argparse
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('quantum_participation_audit_results.json')
OBS = {}
EDGES = list(itertools.combinations(range(3), 2))


def one_factor(q, index, count=3):
    out = np.array([[1.]], complex)
    for j in range(count):
        out = np.kron(out, q if j == index else np.eye(2))
    return out


def model(omega=1., detuning=1., coupling=1.):
    swaps = [core.swap(3, *edge) for edge in EDGES]
    occupancies = [one_factor(np.diag([0., 1.]), j) for j in range(3)]
    flips = [one_factor(core.PAULI[0], j) for j in range(3)]
    h = sum(np.kron(np.eye(8), omega*x+detuning*n)+coupling*np.kron(s, n)
            for s, n, x in zip(swaps, occupancies, flips))
    blank = np.zeros((8, 8), complex)
    blank[0, 0] = 1
    return h, swaps, occupancies, blank


def factor_permutation(destinations):
    count = len(destinations)
    out = np.zeros((2**count, 2**count), complex)
    for col in range(2**count):
        old = [(col >> (count-1-j)) & 1 for j in range(count)]
        new = [0]*count
        for j, dest in enumerate(destinations):
            new[dest] = old[j]
        row = sum(value << (count-1-j) for j, value in enumerate(new))
        out[row, col] = 1
    return out


def joint_relabel(pi):
    edge_destinations = [EDGES.index(tuple(sorted((pi[a], pi[b])))) for a, b in EDGES]
    return np.kron(factor_permutation(pi), factor_permutation(edge_destinations))


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=6e-12):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_identical_rule_and_joint_relabeling(self):
        h, swaps, occ, blank = model()
        errors = []
        for pi in itertools.permutations(range(3)):
            transform = joint_relabel(pi)
            errors.append(float(np.linalg.norm(transform@h@transform.conj().T-h)))
            self.close(transform@h@transform.conj().T, h)
        one = core.evolve(.2*core.PAULI[0]+.3*core.PAULI[1]-.1*core.PAULI[2], 1.)
        data_rotation = np.kron(np.kron(one, one), one)
        rotation = np.kron(data_rotation, np.eye(8))
        self.close(rotation@h@rotation.conj().T, h)
        OBS['same_rule'] = dict(data_qubits=3, relation_qubits=3,
            all_pair_parameters=dict(omega=1, detuning=1, coupling=1),
            joint_relabelings_checked=6, maximum_relabeling_error=max(errors),
            initial_relation_state='|000><000|', prescribed_on_edge_pattern=False)

    def test_02_exact_short_time_operator_response(self):
        h, swaps, occ, blank = model()
        # Integer matrices; projections onto relation vacuum select every eighth row.
        hi = h.real.astype(np.int64)
        certificates = []
        for s, n in zip(swaps, occ):
            ad = np.kron(np.eye(8), n).real.astype(np.int64)
            projected = []
            for k in range(5):
                projected.append((1j**k)*ad[::8, ::8])
                ad = hi@ad-ad@hi
            for k in (0, 1, 3):
                self.assertTrue(np.array_equal(projected[k], np.zeros((8, 8))))
            self.assertTrue(np.array_equal(projected[2], 2*np.eye(8)))
            self.assertTrue(np.array_equal(projected[4], -12*np.eye(8)-4*s))
            certificates.append(dict(order_two='2 I', order_four='-12 I - 4 SWAP'))
        OBS['short_time_response'] = dict(exact_integer_certificates=certificates,
            probability_expansion='t^2 - (1/2 + <SWAP_e>/6) t^4 + O(t^5)',
            asymmetry_for_logical_z='p_12 - p_13 = t^4/4 + O(t^5)',
            coupling_table_replaced_by_internal_state_response=True)

    def test_03_participation_and_back_action(self):
        h, swaps, occ, blank = model()
        a, b, c, p, x, y, z = core.operators()
        source = (p+z)/4
        initial = np.kron(source, blank)
        rows = []
        for t in (.2, .5, .8):
            u = core.evolve(h, t)
            out = u@initial@u.conj().T
            probabilities = [float(np.trace(np.kron(np.eye(8), n)@out).real) for n in occ]
            self.assertGreater(probabilities[0], probabilities[1])
            self.assertAlmostEqual(probabilities[1], probabilities[2])
            disturbance = core.distance(core.partial(out, [8, 8], (0,)), source)
            self.assertGreater(disturbance, 0)
            rows.append(dict(time=t, occupation_probabilities=probabilities,
                source_trace_distance_change=disturbance))
        OBS['dynamic_participation'] = dict(source='rho(logical z=1)', rows=rows,
            unequal_activation_from_initial_data_relations=True,
            no_measurement_or_target_dependent_switch_used=True,
            actual_geometry_or_sparse_limit_derived=False)

    def test_04_affine_evolution_and_unknown_reference(self):
        h, swaps, occ, blank = model()
        _, _, _, p, x, y, z = core.operators()
        first, second = (p+x)/4, (p+y)/4
        u = core.evolve(h, .5)
        mix = u@np.kron((first+second)/2, blank)@u.conj().T
        average = (u@np.kron(first, blank)@u.conj().T+u@np.kron(second, blank)@u.conj().T)/2
        affinity_error = float(np.linalg.norm(mix-average))
        self.close(mix, average)
        rng = np.random.default_rng(43304)
        rho_dr = core.density(16, rng)
        initial = np.kron(rho_dr, blank).reshape(8, 2, 8, 8, 2, 8).transpose(0, 2, 1, 3, 5, 4).reshape(128, 128)
        ur = np.kron(u, np.eye(2))
        after = ur@initial@ur.conj().T
        inverse_error = float(np.linalg.norm(ur.conj().T@after@ur-initial))
        self.close(ur.conj().T@after@ur, initial)
        self.close(core.partial(after, [8, 8, 2], (2,)), core.partial(initial, [8, 8, 2], (2,)))
        OBS['complete_quantum_process'] = dict(affinity_error=affinity_error,
            arbitrary_reference_dimension=2, inverse_error=inverse_error,
            state_dependent_hamiltonian_inserted=False)

    def test_05_new_primitives_are_not_bare_exchange(self):
        h, swaps, occ, blank = model()
        collective_z = sum(one_factor(core.PAULI[2], j, 6) for j in range(6))
        defect = float(np.linalg.norm(h@collective_z-collective_z@h))
        self.assertGreater(defect, 1)
        # Every pair swap on these six raw qubits commutes with collective Z.
        for i, j in itertools.combinations(range(6), 2):
            s = core.swap(6, i, j)
            self.close(s@collective_z, collective_z@s)
        # A genuine weight-three Pauli term Z_edge12 X_data1 X_data2 is present.
        paulis = [core.PAULI[0], core.PAULI[0], np.eye(2), core.PAULI[2], np.eye(2), np.eye(2)]
        word = np.array([[1.]], complex)
        for q in paulis:
            word = np.kron(word, q)
        coefficient = float(np.trace(word@h).real/64)
        self.assertAlmostEqual(coefficient, -.25)
        OBS['primitive_boundary'] = dict(raw_six_qubit_collective_z_commutator_norm=defect,
            weight_three_pauli_coefficient=coefficient,
            implementable_by_bare_six_qubit_exchange_hamiltonian=False,
            encoded_exchange_implementation_ruled_out=False,
            added_relation_flip_and_controlled_exchange_are_model_inputs=True)

    def test_06_no_relation_update_without_new_generator(self):
        h, swaps, occ, blank = model(omega=0.)
        _, _, _, p, x, y, z = core.operators()
        source = (p+z)/4
        initial = np.kron(source, blank)
        self.close(h@initial, np.zeros((64, 64)))
        u = core.evolve(h, .8)
        self.close(u@initial@u.conj().T, initial)
        for n in occ:
            self.close(h@np.kron(np.eye(8), n), np.kron(np.eye(8), n)@h)
        OBS['update_source'] = dict(omega_zero_conserves_every_edge_occupation=True,
            vacuum_relations_with_no_flip_give_no_data_evolution=True,
            fresh_dynamic_edges_do_not_follow_from_controlled_swap_alone=True)


def run():
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=433, baseline_round=432, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(dynamic_graph_literature_is_bridge_not_derivation=True,
            operator_participation_generates_state_dependent_activation=True,
            fixed_global_linear_evolution_and_reference_preserved=True,
            new_relation_resources_and_generators_explicit=True,
            derived_from_429_exchange_alone=False,
            full_cognitive_implementation_completed=False,
            three_dimensional_space_unconditionally_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
