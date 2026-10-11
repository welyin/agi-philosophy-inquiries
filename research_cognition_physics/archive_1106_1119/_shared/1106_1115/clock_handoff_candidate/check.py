"""Candidate only: finite-window handoff of the unchanged round-925 material.

No formal research round or new scientific calibration is claimed.
Default: recompute and compare saved results, writing nothing.
--write-results: create the initial results once; never overwrite it.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "archive_923_934/925/record_resource_exchange.py"


def op(a):
    return float(np.linalg.norm(a, 2))


def td(a, b):
    d = a - b
    return float(np.sum(np.abs(np.linalg.eigvalsh((d + d.conj().T) / 2))) / 2)


def run():
    spec = importlib.util.spec_from_file_location("handoff_source_925", SOURCE)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    z = old.model()
    h60, idx = z["h"], z["idx"]
    w = np.zeros((60, 6), complex)
    for q in range(3):
        w[idx(0, q, 3), 2*q] = 1
        w[idx(1, q, 3+q), 2*q+1] = 1
    h = w.conj().T @ h60 @ w
    x = np.array([[0, 1], [1, 0]], complex)
    es = np.array([5.2, 5.9, 6.7])
    explicit_h = np.kron(np.diag(es), np.eye(2)) + np.kron(np.eye(3), .4*x)
    spectrum = np.linalg.eigvalsh(h60)
    grid_error = float(np.max(np.abs(10*spectrum - np.rint(10*spectrum))))
    assert op(h60 @ w - w @ h) < 1e-13
    assert op(h - explicit_h) < 1e-13 and grid_error < 1e-10
    g1 = dict(full_single_material_dimension=60, closed_single_material_dimension=6,
              closed_joint_dimension=72, full_joint_dimension=7200,
              exact_subspace_intertwining_error=op(h60 @ w - w @ h),
              inherited_h_formula_error=op(h - explicit_h),
              original_spectrum_tenth_grid_error=grid_error)

    d = 6
    ident = np.eye(d)
    swap = np.zeros((d*d, d*d), complex)
    for a in range(d):
        for b in range(d):
            swap[b*d+a, a*d+b] = 1
    h0 = np.kron(np.kron(h, ident) + np.kron(ident, h), np.eye(2))
    hx = np.kron(np.kron(h, ident), np.eye(2))
    hy = np.kron(np.kron(ident, h), np.eye(2))
    vin = np.zeros((d, 3), complex)
    for q in range(3):
        vin[2*q, q] = 1
    t925 = np.pi / (2*z["g"])
    vout = old.unitary(h, t925) @ vin
    blank = np.zeros(d, complex)
    blank[0], blank[1] = 1/np.sqrt(2), -1/np.sqrt(2)
    c0, c1 = np.array([1, 0], complex), np.array([0, 1], complex)
    initial = np.column_stack([np.kron(np.kron(vout[:, q], blank), c0)
                               for q in range(3)])
    kappa = 1/40
    ts = np.pi/(2*kappa)
    k = kappa*np.kron(swap, x)
    total_h = h0+k
    uh = old.unitary(h, ts)
    target = -1j*np.column_stack([np.kron(np.kron(uh@blank, uh@vout[:, q]), c1)
                                  for q in range(3)])
    actual = old.unitary(total_h, ts) @ initial
    dx = actual.conj().T@hx@actual-initial.conj().T@hx@initial
    dy = actual.conj().T@hy@actual-initial.conj().T@hy@initial
    expected_dx = np.diag(4.8-es)
    assert op(total_h@h0-h0@total_h) < 1e-12
    assert op(actual-target) < 1e-10
    assert op(dx-expected_dx) < 1e-10 and op(dx+dy) < 1e-10
    assert op(initial.conj().T@k@k@initial-kappa*kappa*np.eye(3)) < 1e-12
    g2 = dict(kappa=kappa, handoff_proper_time=ts,
              energy_commutator_error=op(total_h@h0-h0@total_h),
              all_input_handoff_isometry_error=op(actual-target),
              source_energy_transfer_diagonal=np.diag(dx).real.tolist(),
              receiver_energy_transfer_diagonal=np.diag(dy).real.tolist(),
              transfer_operator_sum_error=op(dx+dy),
              initial_contact_compression=op(initial.conj().T@k@initial),
              final_contact_compression=op(actual.conj().T@k@actual),
              contact_variance=kappa*kappa,
              contact_variance_identity_error=op(initial.conj().T@k@k@initial-kappa*kappa*np.eye(3)),
              minimum_closed_block_energy=float(np.min(np.linalg.eigvalsh(total_h))))

    windows = []
    for j in (1, 10):
        kappa = 1/(40*j)
        ts = 20*np.pi*j
        k = kappa*np.kron(swap, x)
        total_h = h0+k
        uh = old.unitary(h, ts)
        target = -1j*np.column_stack([np.kron(np.kron(uh@blank, uh@vout[:, q]), c1)
                                      for q in range(3)])
        for s in (.1, t925):
            a = old.unitary(total_h, ts+s)@initial
            b = old.unitary(h0, s)@target
            angle = kappa*s
            cross_error = op(b.conj().T@a-np.cos(angle)*np.eye(3))
            isometry_distance = op(a-b)
            assert cross_error < 1e-9
            assert abs(isometry_distance-2*abs(np.sin(angle/2))) < 1e-9

            # Complete Choi input: source code maximally entangled with a passive qutrit.
            amp = (a/np.sqrt(3)).reshape(6, 6, 2, 3)
            actual_yr = np.einsum("xycr,xzcq->yrzq", amp, amp.conj()).reshape(18, 18)
            receiver_v = old.unitary(h, ts+s)@vout
            ideal_vector = (receiver_v/np.sqrt(3)).reshape(-1)
            ideal_yr = np.outer(ideal_vector, ideal_vector.conj())
            blank_s = old.unitary(h, ts+s)@blank
            replace_yr = np.kron(np.outer(blank_s, blank_s.conj()), np.eye(3)/3)
            mixture = np.cos(angle)**2*ideal_yr+np.sin(angle)**2*replace_yr
            mixture_error = op(actual_yr-mixture)
            assert mixture_error < 1e-9

            # q=1 is orthogonal to the q=0 blank, giving the receiver diamond lower witness.
            amp1 = a[:, 1].reshape(6, 6, 2)
            receiver1 = np.einsum("xyc,xzc->yz", amp1, amp1.conj())
            ideal1 = np.outer(receiver_v[:, 1], receiver_v[:, 1].conj())
            witness = td(receiver1, ideal1)
            assert abs(witness-np.sin(angle)**2) < 1e-9
            flag1 = float(np.sum(abs(amp[:, :, 1, :])**2))
            assert abs(1-flag1-np.sin(angle)**2) < 1e-9
            windows.append(dict(period_multiplier=j, kappa=kappa, handoff_proper_time=ts,
                continuation_proper_time=s,
                whole_60d_h_period_error=op(old.unitary(h60, ts)-np.eye(60)),
                all_input_cross_gram_error=cross_error,
                complete_output_isometry_distance=isometry_distance,
                exact_complete_half_diamond=float(abs(np.sin(angle))),
                exact_receiver_half_diamond=float(np.sin(angle)**2),
                receiver_q1_trace_distance_witness=witness,
                receiver_choi_mixture_identity_error=mixture_error,
                actual_flag_failure_probability=float(1-flag1)))

    return dict(candidate="fixed_925_finite_window_handoff", formal_round=None,
        scientific_calibration_increment=0, numerical_groups=3,
        status="PASS_FINITE_WINDOW_ONLY",
        source={"path":str(SOURCE.relative_to(ROOT)).replace("\\", "/"),
                "sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
        group1_material_subspace=g1, group2_handoff_and_energy=g2,
        group3_finite_windows=windows,
        scope=dict(arbitrary_passive_reference_protected=True,
                   full_contact_hamiltonian_retained=True,
                   no_postselection=True,
                   autonomous_continuation_is_original_925_h=True,
                   actual_fresh_input_loading_implemented=False,
                   original_CO_meter_and_flag_readout_implemented=False,
                   permanent_isolated_takeover_implemented=False,
                   common_SR_model_completed=False,
                   formal_1106_completed=False))


def equal_saved(actual, saved, path="root"):
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and actual.keys() == saved.keys(), path
        for key in actual:
            equal_saved(actual[key], saved[key], path+"."+key)
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(actual) == len(saved), path
        for i, (a, b) in enumerate(zip(actual, saved)):
            equal_saved(a, b, path+"["+str(i)+"]")
    elif isinstance(actual, float):
        assert math.isfinite(actual) and math.isclose(actual, saved, rel_tol=1e-9, abs_tol=1e-10), path
    else:
        assert actual == saved, path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    data = run()
    target = HERE / "results.json"
    if args.write_results:
        assert not target.exists(), "Refusing to overwrite saved candidate results"
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    equal_saved(data, json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps({"candidate":data["candidate"], "status":data["status"],
                      "groups":data["numerical_groups"], "saved_results_match":True,
                      "formal_round":None, "scientific_calibration_increment":0}, ensure_ascii=False))
