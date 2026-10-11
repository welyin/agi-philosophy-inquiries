"""Admission only: complete free-history phases of one typed port candidate.

No normal localized collision, loading duration or autonomous controller is
certified. Default verifies saved results without changing files.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def run():
    length = 5*math.pi/18  # fixed geometry; nominal clock phase is -pi/6
    duration = 17*length/15
    uc, ua = np.array([17/8, 15/8]), np.array([1., 0.])
    q = np.array([[2., 0.], [25/8, 15/8], [17/4, 15/4]])
    masses = np.array([2., 2.5, 2.])
    velocities = q[:, 1]/q[:, 0]
    starts = np.array([length, 8*length/25, 0.])
    ends = starts+velocities*duration
    displacements = np.column_stack([np.full(3, duration), length-starts])
    action = q[:, 0]*displacements[:, 0]-q[:, 1]*displacements[:, 1]
    assert np.max(abs(ends-length)) < 1e-14
    assert np.max(abs(action-masses*masses*duration/q[:, 0])) < 1e-14
    tau_c, tau_a = 8*length/15, duration
    phi = tau_c-tau_a
    assert abs(phi+math.pi/6) < 1e-14
    # Known preparation is declared before fresh Q; it does not use its state.
    reference_in = np.exp(1j*action)/math.sqrt(3)
    reference_out = np.exp(-1j*action)*reference_in
    clock_in = np.array([0., 1., 1., 0.])/math.sqrt(2)
    proper_free = np.diag([np.exp(-1j*((2+c)*tau_c+(2+a)*tau_a))
                          for c in range(2) for a in range(2)])
    clock_out = proper_free @ clock_in
    bt = np.zeros((12, 12), complex)
    for n in range(2):
        bt[4+n, 6+n], bt[6+n, 4+n] = 1j, -1j
    effect = (np.eye(12)+bt)/2
    probability = lambda c, r: float(np.vdot(np.kron(c, r), effect @ np.kron(c, r)).real)
    calibrated = probability(clock_out, reference_out)
    uncalibrated = probability(clock_out, np.exp(-1j*action)/math.sqrt(3))
    predicted = .5+sum(np.sin(phi-np.diff(action)))/6
    assert abs(calibrated-1/3) < 1e-14 and abs(uncalibrated-predicted) < 1e-14
    # A different experiment: entire fixed momentum directory waits globally.
    # Its invariant effect cannot measure that uniform external parameter.
    pca = np.array([(2+c)*uc+(2+a)*ua for c in range(2) for a in range(2)])
    total_energy = (pca[:, 0, None]+q[None, :, 0]).reshape(-1)
    incoming = np.kron(clock_in, np.ones(3)/math.sqrt(3))
    global_out = np.exp(-1j*total_energy*duration)*incoming
    global_prob = float(np.vdot(global_out, effect @ global_out).real)
    comm = np.linalg.norm(effect @ np.diag(total_energy)-np.diag(total_energy) @ effect)
    assert comm < 1e-14 and abs(global_prob-.5) < 1e-14
    # The calibrated local port preparation pulled back to one momentum basis.
    # It is different from the uniform reference in the control above.
    reference_global_in = reference_in*np.exp(-1j*q[:, 1]*starts)
    joint_global_in = np.kron(clock_in, reference_global_in)
    joint_global_out = np.exp(-1j*total_energy*duration)*joint_global_in
    calibrated_initial_global = float(np.vdot(joint_global_in, effect @ joint_global_in).real)
    calibrated_final_global = float(np.vdot(joint_global_out, effect @ joint_global_out).real)
    assert abs(calibrated_initial_global-calibrated) < 1e-14
    assert abs(calibrated_final_global-calibrated) < 1e-14
    phase_identity = action-q[:, 1]*starts-(q[:, 0]*duration-q[:, 1]*length)
    assert np.max(abs(phase_identity)) < 1e-14
    # The same free segments and all momenta transform together.
    beta = -.3
    gamma = 1/math.sqrt(1-beta*beta)
    lam = gamma*np.array([[1., beta], [beta, 1.]])
    qp, dxp = q @ lam.T, displacements @ lam.T
    transformed_action = qp[:, 0]*dxp[:, 0]-qp[:, 1]*dxp[:, 1]
    assert np.max(abs(action-transformed_action)) < 1e-14
    return dict(status='PASS_TYPED_FREE_HISTORY_ONLY', numerical_cases=1,
                length=length, duration=duration,
                reference_start_positions=starts.tolist(), velocities=velocities.tolist(),
                arrival_position_residual=float(np.max(abs(ends-length))),
                complete_reference_actions=action.tolist(),
                proper_times=[tau_c, tau_a], clock_phase=phi,
                calibrated_plus=calibrated, uncalibrated_plus=uncalibrated,
                general_phase_formula_plus=float(predicted),
                global_wait_plus=global_prob, global_wait_commutator=float(comm),
                calibrated_initial_global_plus=calibrated_initial_global,
                calibrated_final_global_plus=calibrated_final_global,
                preparation_phase_identity_residual=float(np.max(abs(phase_identity))),
                standalone_dynamic_clock_certificate=False,
                transported_action_residual=float(np.max(abs(action-transformed_action))),
                actual_localized_loading_and_collision=False,
                complete_internal_timing_and_control=False,
                new_completed_rounds=0, scientific_calibration_increment=0,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def equal(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a: equal(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): equal(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, abs_tol=1e-12, rel_tol=1e-10), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    destination = HERE/'clock_port_phase_results.json'
    if args.write_results:
        with destination.open('x', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write('\n')
    equal(data, json.loads(destination.read_text(encoding='utf-8')))
    print(json.dumps(dict(status=data['status'], saved_results_match=True,
                          new_completed_rounds=0)))
