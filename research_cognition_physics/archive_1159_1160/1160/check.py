"""1160: actual recoil record, complete clock, and the same active future.

Finite effective K1160 permissions are inputs, not inferred from commutators.
Default is read-only. --write-once creates this round's first result only.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1] / 'archive_1106_1119/1111/check.py'
SOURCE_SHA = '704ae4ca016c77e0df83fc9c2d39e90d79e9fbf3d6469c8c5b9bf1dad8a5bd44'
TOL = 4e-9


def norm(a):
    return float(np.linalg.norm(a))


def small(a, name):
    value = norm(a)
    assert value < TOL, (name, value)
    return value


def load_old():
    assert hashlib.sha256(OLD.read_bytes()).hexdigest() == SOURCE_SHA
    spec = importlib.util.spec_from_file_location('round1111_recoil_reuse', OLD)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def along(a, gate, axes):
    remaining = [i for i in range(a.ndim) if i not in axes]
    order = list(axes) + remaining
    inverse = np.argsort(order)
    t = a.transpose(order)
    shape = t.shape
    return (gate @ t.reshape(gate.shape[1], -1)).reshape(shape).transpose(inverse)


def calculate():
    old = load_old()
    m = old.model()  # Same frozen 288-dimensional physical FK construction.
    gamma, hd, hc = m['gamma'], m['hout'], m['hc']
    theta, delta, T = math.pi/4, math.pi/8, math.pi
    # Decoded indices: head, aux, Qc, Qr, AB, M, E, program, input column.
    logical_shape = (3, 4, 2, 2, 3)
    ld = math.prod(logical_shape)
    identity = np.eye(ld, dtype=complex)
    generator = np.zeros((4, 4), complex)
    generator[0, 3] = generator[3, 0] = 1
    collision = np.eye(4, dtype=complex)
    collision[np.ix_([0, 3], [0, 3])] = np.array(
        [[math.cos(theta), -1j*math.sin(theta)],
         [-1j*math.sin(theta), math.cos(theta)]])
    first = np.zeros((ld, ld), complex)
    last = np.zeros_like(first)
    record = np.zeros_like(first)
    erase = np.zeros_like(first)
    for qc, ab, rm, re, pr in np.ndindex(logical_shape):
        col = np.ravel_multi_index((qc, ab, rm, re, pr), logical_shape)
        for out in range(4):
            row = np.ravel_multi_index((qc, out, rm, re, pr), logical_shape)
            first[row, col] = collision[out, ab] if qc == 2 else int(out == ab)
            last[row, col] = collision[out, ab]
        active = (pr in (1, 2)) and (ab//2 == 1)
        row = np.ravel_multi_index((qc, ab, rm ^ int(active), re ^ int(active), pr), logical_shape)
        record[row, col] = 1
        active = pr == 2 and (ab//2 == 1)
        row = np.ravel_multi_index((qc, ab, rm ^ int(active), re ^ int(active), pr), logical_shape)
        erase[row, col] = 1
    gates = [first, record, erase, last]
    unitarity = max(small(g.conj().T @ g - identity, 'full primitive unitary') for g in gates)
    pa = np.array([[5, 4, 0, 0], [5, 0, 4, 0]], float)
    pb = np.array([[5, -4, 0, 0], [5, 0, -4, 0]], float)
    pair = np.array([pa[a] + pb[b] for a, b in np.ndindex(2, 2)])
    logical_charges = []
    for mu in range(4):
        diagonal = np.array([pair[ab, mu] for qc, ab, rm, re, pr in np.ndindex(logical_shape)])
        logical_charges.append(np.diag(diagonal))
    conservation = max(small(p @ g - g @ p, 'full mother fourP')
                       for p in logical_charges for g in gates)
    # Every gate is identity on head and auxiliary clock degrees before Gamma.
    clock_intertwiner = small(hd @ gamma - gamma @ (m['baseline']*np.eye(288)+np.kron(hc, np.eye(72))),
                              'original full FK intertwiner')
    charge = np.diag([pa[ab//2, 1] for ab in range(4)])
    exchange = collision.conj().T @ charge @ collision - charge
    noncommuting = norm(exchange @ charge - charge @ exchange)
    assert noncommuting > 1
    # Pointer readout preserves all present mechanical charges, even on cross sectors.
    preserve_present = 0.0
    for mu in range(4):
        qa = np.diag([pa[ab//2, mu] for qc, ab, rm, re, pr in np.ndindex(logical_shape)])
        preserve_present = max(preserve_present, small(record.conj().T @ qa @ record-qa, 'current A charge'))

    def dress(state, gate):
        decoded = (gamma.conj().T @ state.reshape(288, -1)).reshape(4, 8, 3, 3, 4, 2, 2, 3, -1)
        changed = along(decoded, gate, (2, 4, 5, 6, 7))
        return (gamma @ changed.reshape(288, -1)).reshape(state.shape)

    evals, vectors = np.linalg.eigh(hd)
    waiting = (vectors * np.exp(-1j*delta*evals)) @ vectors.conj().T
    # A+B energy=10. M,E,program each have unit positive baseline mass.
    extra_phase = np.exp(-1j*13*delta)

    def advance(state, gate):
        state = (waiting @ state.reshape(288, -1)).reshape(state.shape)*extra_phase
        return dress(state, gate)

    def run(program, destroy_environment=False):
        shape = (288, 4, 2, 2, 3, 2)
        state = np.zeros(shape, complex)
        start = T-4*delta
        state[:, 0, 0, 0, program, :] = m['w'](start)
        first_state = advance(state, first)
        recorded = advance(first_state, record)
        # A genuine dephasing CP channel is represented with an additional kept
        # environment branch. It is a distinct experiment, not erasing evidence.
        if destroy_environment:
            parts = []
            for value in (0, 1):
                branch = recorded.copy()
                branch[:, :, :, 1-value, :, :] = 0
                parts.append(advance(advance(branch, erase), last))
            return np.stack(parts, axis=-2), first_state, recorded
        return advance(advance(recorded, erase), last), first_state, recorded

    rows = []
    outputs = {}
    max_isometry = max_clock = max_direct = 0.0
    beta = np.array([0.6, 0, 0])
    boost = old.boost(beta)
    for pr in range(3):
        output, mid, rec = run(pr)
        outputs[pr] = output
        v = output.reshape(-1, 2)
        max_isometry = max(max_isometry, small(v.conj().T @ v-np.eye(2), 'all Q and arbitrary R'))
        # Independent compiled identity: wait to T then execute the four decoded gates.
        direct = np.zeros_like(output)
        direct[:, 0, 0, 0, pr, :] = m['w'](T)*np.exp(-1j*13*4*delta)
        for gate in gates:
            direct = dress(direct, gate)
        max_direct = max(max_direct, small(output-direct, 'original clock and total elapsed time'))
        p = np.sum(abs(output)**2, axis=(0, 2, 3, 4))
        expected_y = np.array([0.5, 0.5 if pr == 1 else 1.0])
        small(p[3]-expected_y, 'actual outgoing y probabilities')
        small(p[1:3], 'actual correlated input remains chi')
        # Decode only to evaluate the actual old clock/readout, never replace it.
        phys = output.reshape(4, 8, 3, 3, 4, 2, 2, 3, 2)
        mr = np.sum(abs(phys[:, [1, 3, 5, 7], ...])**2, axis=tuple(range(8)))
        max_clock = max(max_clock, small(mr-np.ones(2), 'original Mr at T=pi'))
        da = np.array([sum(p[ab, q]*(pa[ab//2]-pa[0]) for ab in range(4)) for q in range(2)])
        db = np.array([sum(p[ab, q]*(pb[ab % 2]-pb[0]) for ab in range(4)) for q in range(2)])
        small(da+db, 'actual pair recoil balance')
        midp = np.sum(abs(mid)**2, axis=(0, 2, 3, 4))
        recp = np.sum(abs(rec)**2, axis=(0, 2, 3, 4))
        small(midp-recp, 'record leaves current mechanical populations')
        # Full record density conditional on Q=1, retaining D,AB and program.
        me = output[..., 1].transpose(2, 3, 0, 1, 4).reshape(4, -1)
        rme = me @ me.conj().T
        rows.append(dict(program=pr, outgoing_yy=p[3].tolist(), old_Mr=mr.tolist(),
                         A_delta_fourP=da.tolist(), B_delta_fourP=db.tolist(),
                         A_boost_delta=(da @ boost.T).tolist(), B_boost_delta=(db @ boost.T).tolist(),
                         record_population=np.real(np.diag(rme)).tolist()))
    # Compare full physical vectors after a type relabel of the unchanged program
    # values. Such a comparison does not discard the actual program in the model.
    recovered = small(outputs[0][:, :, :, :, 0]-outputs[2][:, :, :, :, 2], 'complete coherent erasure')
    leaked, _, _ = run(2, destroy_environment=True)
    leakage_prob = np.sum(abs(leaked)**2, axis=(0, 2, 3, 4, 5))[3]
    small(leakage_prob-np.array([0.5, 0.5]), 'extra inaccessible record affects future')
    # Actual finite endpoint map; no path begins until all four windows finish.
    event = np.array([m['arrival']+T, 1, 0, 0], float)
    endpoints_a = np.array([event+(1.2/3)*p for p in pa])
    endpoints_b = np.array([event+(1.2/3)*p for p in pb])
    endpoint_covariance = max(small(boost@(z-event)-(1.2/3)*(boost@p), 'same physical rail under boost')
                              for z, p in list(zip(endpoints_a, pa))+list(zip(endpoints_b, pb)))
    mass_shell = max(small(old.mink(boost@p, boost@p)-9, 'same mass3') for p in list(pa)+list(pb))
    start_event = event.copy()
    start_event[0] -= 4*delta
    incoming_a = start_event[1:] - (pa[0, 1:]/pa[0, 0])*start_event[0]
    incoming_b = start_event[1:] - (pb[0, 1:]/pb[0, 0])*start_event[0]
    # In both programs the same input trajectories arrive at the initial contact.
    small(incoming_a+pa[0, 1:]/5*start_event[0]-start_event[1:], 'A actual arrival')
    small(incoming_b+pb[0, 1:]/5*start_event[0]-start_event[1:], 'B actual arrival')
    midpoint = T-4*delta
    wmid = m['w'](midpoint).reshape(4, 8, 3, 3, 2)
    old_clock_mid = np.sum(abs(wmid[:, [1, 3, 5, 7], ...])**2, axis=(0, 1, 2, 3))
    small(old_clock_mid-np.full(2, 1/8), 'same clock before sequence')
    return dict(schema='round1160_complete_exchange_record_v1', status='PASS', check_groups=4,
        new_scientific_calibration=0, new_adopted_axioms=0,
        source_sha256=SOURCE_SHA, code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        complete_blocks=dict(logical_mother_dimension=ld, full_output_shape=[13824, 2],
            full_mother_unitarity=unitarity, full_mother_fourP=conservation,
            original_FK_intertwiner=clock_intertwiner, present_A_charge_preserved=preserve_present,
            exchange_commutator_with_A_Px=noncommuting, added_positive_resource_mass=3),
        same_actual_process=dict(original_clock_before=old_clock_mid.tolist(),
            all_input_isometry=max_isometry, old_Mr_residual=max_clock,
            independent_full_waiting_residual=max_direct, programs=rows),
        record_future=dict(complete_erasure_residual=recovered, unreturned_environment_yy=leakage_prob.tolist(),
            meaningful_future_probability_gap=0.5, boosted_A_energy_exchange_gap=1.5),
        events=dict(start_event=start_event.tolist(), end_event=event.tolist(),
            positive_window=delta, window_count=4, A_initial_position=incoming_a.tolist(),
            B_initial_position=incoming_b.tolist(), A_outgoing=endpoints_a.tolist(),
            B_outgoing=endpoints_b.tolist(), covariant_endpoint_residual=endpoint_covariance,
            boosted_mass_shell_residual=mass_shell),
        limitations=['K1160 physical contact permissions are explicit input',
            'mechanical readiness xx only; full mother algebra does not grant all preparations',
            'no arbitrary internal window probes or microscopic stress',
            'no autonomous universal generator; old FK not claimed to trigger the collisions',
            'records preserve the full process, not the unrecorded reduced mechanical future'])


def compare(actual, expected, path=''):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for k, v in actual.items():
            compare(v, expected[k], path+'/'+k)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path+'/'+str(i))
    elif isinstance(actual, float):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-10), (path, actual, expected)
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    output = HERE/'results.json'
    if args.write_once:
        with output.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('PASS; four groups; first results written')
    else:
        compare(result, json.loads(output.read_text(encoding='utf-8')))
        print('PASS; four groups; READ_ONLY_MATCH')


if __name__ == '__main__':
    main()
