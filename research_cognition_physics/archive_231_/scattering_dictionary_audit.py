"""Round 405: a fixed scattering dictionary and its physical access obstruction.

Reuses only the frozen round-404 sparse evolution functions. No old report or
old experiment is executed. Universal claims have proofs in research_note_405.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np
import quiescent_program_sector_audit as q

TARGET = Path(__file__).with_name("scattering_dictionary_audit_results.json")


def free(state, time):
    return {tuple(sorted((p+time*q.SPEEDS[k], k) for p, k in c)): v.copy()
            for c, v in state.items()}


def radius(state):
    return max((abs(p) for c in state for p, k in c), default=0)


def wave_at(state, time):
    return free(q.evolve(state, time), -time)


def inverse_wave_at(state, time):
    return q.evolve(free(state, time), -time)


def wave(state):
    return wave_at(state, 4*(2*radius(state)+1))


def inverse_wave(state):
    return inverse_wave_at(state, 2*radius(state)+1)


def vector_distance(a, b):
    return float(np.sqrt(sum(np.linalg.norm(a.get(c, 0)-b.get(c, 0))**2
                             for c in a.keys() | b.keys())))


def pure_trace_distance(a, b):
    # Stable phase-aligned norm formula for two normalized pure states.
    z = sum(np.vdot(a[c], b[c]) for c in a.keys() & b.keys())
    phase = np.conj(z)/abs(z) if abs(z) else 1.
    d = vector_distance(a, {c: phase*v for c, v in b.items()})
    return float(d*np.sqrt(max(0., 1.-d*d/4)))


def depleted(config):
    positions = [p for p, k in config]
    speeds = [q.SPEEDS[k] for p, k in config]
    return len(set(positions)) == len(positions) and speeds == sorted(speeds)


def bad_weight(state):
    return float(sum(np.vdot(v, v).real for c, v in state.items()
                     if not depleted(c)))


def coherent_case(n, rng):
    keys = {((0, 0), (0, 1), (0, 5))}
    for _ in range(5):
        keys.add(tuple((p, k) for p in range(-n, n+1) for k in range(6)
                       if rng.random() < .14))
    state = {c: rng.normal(size=3)+1j*rng.normal(size=3) for c in sorted(keys)}
    scale = np.sqrt(q.norm2(state))
    return {c: v/scale for c, v in state.items()}


def gate_config(length, data, control):
    c = [(3*length, 1), (4*length, 0)]
    if data:
        c.append((0, 3))
    if control:
        c.append((-length, 4))
    return tuple(sorted(c))


def gate_matrix(length):
    basis = [gate_config(length, a, b) for a in (0, 1) for b in (0, 1)]
    # Four orthogonal reference components audit all columns at once.
    state = {c: np.eye(4, dtype=complex)[j] for j, c in enumerate(basis)}
    out = wave_at(state, length+2)
    matrix = np.array([out.get(c, np.zeros(4, complex)) for c in basis])
    leakage = sum(np.vdot(v, v).real for c, v in out.items() if c not in basis)
    return matrix, float(leakage), vector_distance(
        out, wave_at(state, length+7))


def data_program_state(length, tail_probability):
    state = {}
    controllers = tuple(sorted([(-length, 4), (3*length, 1), (4*length, 0)]))
    for a in (0, 1):
        data = ((0, 3),) if a else ()
        v = np.eye(2, dtype=complex)[a]/np.sqrt(2)
        state[data] = np.sqrt(1-tail_probability)*v
        state[tuple(sorted(data+controllers))] = np.sqrt(tail_probability)*v
    return {c: v for c, v in state.items() if np.linalg.norm(v) > 0}


def moving_data_reference(state, time):
    # Passive readout of the v=+1 qubit at its ballistic position.
    blocks = {}
    target = (time, 3)
    for c, v in state.items():
        bit = int(target in c)
        outside = tuple(item for item in c if item != target)
        blocks.setdefault(outside, np.zeros(4, complex))[2*bit:2*bit+2] += v
    return sum(np.outer(v, v.conj()) for v in blocks.values())


@lru_cache(None)
def report():
    rng = np.random.default_rng(405)
    finite = []
    for n in (1, 2):
        state = coherent_case(n, rng)
        tau, kappa = 4*(2*n+1), 2*n+1
        w = wave_at(state, tau)
        inv = inverse_wave_at(state, kappa)
        finite.append(dict(
            input_radius=n, depletion_time=tau, inverse_time=kappa,
            interacting_bad_weight=bad_weight(q.evolve(state, tau)),
            free_bad_weight=bad_weight(free(state, kappa)),
            forward_stabilization_error=vector_distance(w, wave_at(state, tau+9)),
            inverse_stabilization_error=vector_distance(inv, inverse_wave_at(state, kappa+9)),
            inverse_after_wave_error=vector_distance(inverse_wave(w), state),
            wave_after_inverse_error=vector_distance(wave(inv), state),
            intertwining_error=vector_distance(wave(q.step(state)), free(w, 1)),
            reference_error=float(np.linalg.norm(q.reference_state(w)-q.reference_state(state))),
            norm_error=abs(q.norm2(w)-1),
            input_configurations=len(state), output_configurations=len(w)))
    x = np.array([[0., 1.], [1., 0.]], complex)
    ideal = np.diag([1., 1., 1., np.exp(1j*np.pi/4)])
    gates = []
    for length in (2, 5, 11, 23):
        matrix, leakage, stabilization = gate_matrix(length)
        o = matrix.conj().T @ np.kron(x, np.eye(2)) @ matrix
        remote = np.kron(np.eye(2), x)
        comm = float(np.linalg.norm(o@remote-remote@o, 2))
        gates.append(dict(
            separation=length, latest_program_position=4*length,
            complete_gate_error=float(np.linalg.norm(matrix-ideal, 2)),
            leakage_weight=leakage, stabilization_error=stabilization,
            remote_commutator_norm=comm,
            predicted_commutator_norm=float(2*np.sin(np.pi/8)),
            local_approximation_lower_bound=comm/2))
    future = []
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    bell_rho = np.outer(bell, bell)
    for length in (9, 17):
        for tail in (.0001, .04, .36):
            initial = data_program_state(length, tail)
            start = q.evolve(initial, 4)  # P_0 program truncation, tau_0=4.
            for k in (0, length-4, length-3, 2*length):
                actual = q.evolve(start, k)
                comparison = free(start, k)
                distance = pure_trace_distance(actual, comparison)
                phase = np.pi/4 if 4+k > length else 0.
                predicted_squared = 2*tail*(1-tail/2)*np.sin(phase/2)**2
                moving = moving_data_reference(actual, 4+k)
                moving_free = moving_data_reference(comparison, 4+k)
                moving_distance = float(np.abs(np.linalg.eigvalsh(moving-moving_free)).sum()/2)
                future.append(dict(
                    delay=length, tail_probability=tail, initial_time=4, future_steps=k,
                    whole_reference_trace_distance=distance,
                    exact_squared_distance=predicted_squared,
                    analytic_complete_channel_bound=min(1., 2*np.sqrt(tail)),
                    moving_receiver_trace_distance=moving_distance,
                    reference_error=float(np.linalg.norm(q.reference_state(actual)-np.eye(2)/2)),
                    norm_error=abs(q.norm2(actual)-1)))
    # Boundary illustration, included in the same future/readout test:
    # fixed cell is empty while the moving passive qubit keeps all Bell data.
    ballistic = data_program_state(9, 0.)
    moved = q.evolve(ballistic, 50)
    moving_bell_error = float(np.linalg.norm(moving_data_reference(moved, 50)-bell_rho))
    return dict(
        round=405,
        scope="Specified six-layer QCA in the normal quiescent Hilbert sector: a fixed exact unitary scattering dictionary exists, but fails operator-norm quasilocality; late free-future approximation covers arbitrary passive readouts, including moving ones. No physical position or spatial dimension is generated.",
        finite_support_scattering_checks=finite,
        remote_controlled_phase_checks=gates,
        complete_future_tail_checks=future,
        ballistic_fixed_cell_nonvacuum=q.nonvacuum_probability(moved),
        ballistic_moving_bell_error=moving_bell_error,
        exact_fixed_global_unitary_dictionary_proved=True,
        dictionary_operator_norm_quasilocal=False,
        normal_program_future_bound_uniform_over_future_waits=True,
        passive_moving_readouts_covered=True,
        newly_inserted_active_interventions_covered=False,
        global_information_erased=False,
        convergence_uniform_over_all_programs=False,
        actual_subject_position_generated=False,
        spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False)


class Audit(unittest.TestCase):
    def test_depletion_on_coherent_inputs_with_reference(self):
        for row in report()["finite_support_scattering_checks"]:
            self.assertLess(row["interacting_bad_weight"], 1e-12)
            self.assertLess(row["free_bad_weight"], 1e-12)

    def test_forward_and_inverse_limits_stabilize(self):
        for row in report()["finite_support_scattering_checks"]:
            self.assertLess(row["forward_stabilization_error"], 1e-11)
            self.assertLess(row["inverse_stabilization_error"], 1e-11)

    def test_both_unitary_inverses_and_unknown_reference(self):
        for row in report()["finite_support_scattering_checks"]:
            for key in ("inverse_after_wave_error", "wave_after_inverse_error",
                        "reference_error", "norm_error"):
                self.assertLess(row[key], 1e-11)

    def test_fixed_dictionary_intertwines_the_full_update(self):
        for row in report()["finite_support_scattering_checks"]:
            self.assertLess(row["intertwining_error"], 1e-11)

    def test_remote_collision_exactly_implements_controlled_phase(self):
        for row in report()["remote_controlled_phase_checks"]:
            for key in ("complete_gate_error", "leakage_weight", "stabilization_error"):
                self.assertLess(row[key], 1e-11)

    def test_nondecaying_operator_norm_access_obstruction(self):
        for row in report()["remote_controlled_phase_checks"]:
            self.assertAlmostEqual(row["remote_commutator_norm"],
                                   row["predicted_commutator_norm"], places=11)
            self.assertAlmostEqual(row["local_approximation_lower_bound"],
                                   float(np.sin(np.pi/8)), places=11)

    def test_whole_future_bound_and_moving_passive_readout(self):
        for row in report()["complete_future_tail_checks"]:
            distance = row["whole_reference_trace_distance"]
            self.assertAlmostEqual(distance**2, row["exact_squared_distance"], places=11)
            self.assertLessEqual(distance, row["analytic_complete_channel_bound"]+1e-11)
            self.assertLessEqual(row["moving_receiver_trace_distance"], distance+1e-11)
            self.assertLess(row["reference_error"], 1e-11)
            self.assertLess(row["norm_error"], 1e-11)
        self.assertLess(report()["ballistic_fixed_cell_nonvacuum"], 1e-12)
        self.assertLess(report()["ballistic_moving_bell_error"], 1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result["checks"] = dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors))
    result["runtime"] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

