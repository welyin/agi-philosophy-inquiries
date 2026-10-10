"""Finite calibration for the explicitly assumed 1103 causal implementation.

No empirical data, proof search, universal preparation audit, or new physics.
Default compares saved results; --write refuses to overwrite existing results.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = [
    "research_cognition_physics/archive_1086_/1090/proof.md",
    "research_cognition_physics/archive_1086_/1098/proof.md",
    "research_cognition_physics/archive_1086_/1102/axioms.md",
    "research_cognition_physics/archive_1086_/1102/results.json",
    "research_cognition_physics/archive_1086_/_admission/after_1102_atomic_macro_scope.md",
]
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
ETA = np.diag([1., -1., -1., -1.])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual(a):
    return float(np.max(np.abs(a)))


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b))) / 2)


def psqrt(a):
    v, w = np.linalg.eigh(a)
    assert float(v.min()) > -1e-12
    return (w * np.sqrt(np.maximum(v, 0))) @ w.conj().T


def rails(n, d=2):
    """Tensor digits: 0 = vacuum, 1..d = occupied memory."""
    j = np.zeros(((d+1)**n, n*d), complex)
    for site in range(n):
        for m in range(d):
            index = (m+1)*(d+1)**(n-1-site)
            j[index, site*d+m] = 1
    return j


def original_instrument(points):
    out = [np.zeros((2*len(points), 2*len(points)), complex) for _ in range(2)]
    for k, x in enumerate(points):
        f = x / np.sqrt(1+np.dot(x, x))
        e = .5*I2 + .25*(f[0]*X+f[1]*Y+f[2]*Z)
        for a, block in zip(out, (psqrt(e), psqrt(I2-e))):
            a[2*k:2*k+2, 2*k:2*k+2] = block
    return out


def instrument_with_reference(ks, rho, rd):
    """Explicit classical record plus full conditional quantum state."""
    dim = ks[0].shape[0]*rd
    out = np.zeros((len(ks)*dim, len(ks)*dim), complex)
    for i, k in enumerate(ks):
        kr = np.kron(k, np.eye(rd))
        out[i*dim:(i+1)*dim, i*dim:(i+1)*dim] = kr @ rho @ kr.conj().T
    return out


def clock_state(t, omega=1):
    v = np.exp(-.5j*omega*t*np.array([1., -1.])) / np.sqrt(2)
    return np.outer(v, v.conj())


def boost(v):
    beta2 = float(v @ v)
    g = 1/np.sqrt(1-beta2)
    a = np.eye(4)
    a[0, 0] = g
    a[0, 1:] = -g*v
    a[1:, 0] = -g*v
    if beta2:
        a[1:, 1:] += (g-1)*np.outer(v, v)/beta2
    return a


def compute():
    groups = {}
    points = np.array([[-1., 0, 0], [1., 0, 0], [0., 2., 0]])
    ks = original_instrument(points)
    j = rails(3)
    p = j @ j.conj().T
    physical = [j @ k @ j.conj().T for k in ks]+[np.eye(27)-p]
    vphys = np.vstack(physical)
    target = np.vstack([j @ k for k in ks]+[np.zeros_like(j)])
    eps = residual(vphys @ j-target)
    complete = residual(vphys.conj().T @ vphys-np.eye(27))
    embed = np.eye(6, dtype=complex)[:, :4]
    smallks = original_instrument(points[:2])
    natural = max(residual(ks[i] @ embed-embed @ smallks[i]) for i in range(2))
    rng = np.random.default_rng(1103)
    a = rng.normal(size=(12, 12))+1j*rng.normal(size=(12, 12))
    rho = a @ a.conj().T
    rho /= np.trace(rho)
    jr = np.kron(j, I2)
    actual = instrument_with_reference(physical, jr @ rho @ jr.conj().T, 2)
    logical = instrument_with_reference(ks+[np.zeros_like(ks[0])], rho, 2)
    recorded_j = np.kron(np.eye(3), jr)
    reference_error = residual(actual-recorded_j @ logical @ recorded_j.conj().T)
    assert max(eps, complete, natural, reference_error) < 1e-12
    groups["coherent_rails_full_instrument_and_common_window"] = {
        "logical_dim": 6, "complete_physical_dim": 27,
        "reference_dim_calibrated": 2,
        "analytic_reference_scope": "all finite R by isometry identity",
        "isometry_residual": residual(j.conj().T @ j-np.eye(6)),
        "complete_instrument_residual": complete,
        "compiled_isometry_residual": eps,
        "same_hardware_catalog_embedding_residual": natural,
        "explicit_record_and_reference_residual": reference_error,
    }

    psi = np.zeros(6, complex)
    psi[0] = psi[3] = 1/np.sqrt(2)
    pure = np.outer(psi, psi.conj())
    dephased = np.zeros_like(pure)
    for i in range(3):
        e = np.zeros((6, 6))
        e[2*i:2*i+2, 2*i:2*i+2] = I2.real
        dephased += e @ pure @ e
    w = vphys @ j
    loss = trace_distance(w @ pure @ w.conj().T, w @ dephased @ w.conj().T)
    assert abs(loss-.5) < 1e-12
    groups["measuring_location_is_not_coherent_routing"] = {
        "full_stinespring_output_trace_distance": loss,
        "location_was_not_measured_in_positive_construction": True,
    }

    # CO4 writes a function of orthogonal labels, retaining both inputs.
    coords = [-1., 1.]
    midpoints = [-1., 0., 1.]
    m = np.zeros((12, 4))
    for a in range(2):
        for b in range(2):
            k = midpoints.index((coords[a]+coords[b])/2)
            m[(2*a+b)*3+k, 2*a+b] = 1
    assert residual(m.T @ m-np.eye(4)) == 0
    groups["midpoint_retains_original_joint_labels"] = {
        "isometry_residual": residual(m.T @ m-np.eye(4)),
        "output_midpoints": midpoints,
        "copies_unknown_payload": False,
    }

    # A complete finite event schedule, including a separately fresh Q input.
    inputs = np.array([[0., -1., 0, 0], [0., 1., 0, 0], [.25, 0, 0, 0]])
    gather = np.array([2., 0, 0, 0])
    gate_end = np.array([3., 0, 0, 0])
    outputs = np.array([[5., -1., 0, 0], [5., 1., 0, 0]])
    differences = np.vstack([gather-inputs, gate_end-gather, outputs-gate_end])
    interval = np.einsum("ni,ij,nj->n", differences, ETA, differences)
    assert np.all(differences[:, 0] > 0) and np.all(interval > 0)
    lam = boost(np.array([.6, .0, .0]))
    transformed = differences @ lam.T
    transformed_interval = np.einsum("ni,ij,nj->n", transformed, ETA, transformed)
    invariant_error = residual(transformed_interval-interval)
    assert invariant_error < 1e-12 and np.all(transformed[:, 0] > 0)
    proper = np.sqrt(interval)
    clock_error = max(residual(clock_state(t)-clock_state(s))
                      for t,s in zip(proper, np.sqrt(transformed_interval)))
    outprime = outputs @ lam.T
    transformed_frame = lam
    intrinsic = -outprime @ ETA @ transformed_frame[:, 1:]
    cuttime = outprime @ ETA @ transformed_frame[:, 0]
    assert residual(intrinsic-outputs[:, 1:]) < 1e-12
    assert residual(cuttime-5) < 1e-12
    groups["causal_schedule_cut_dictionary_and_live_clocks"] = {
        "fresh_Q_input_time": .25,
        "first_joint_gate_start": 2.,
        "gate_duration": 1.,
        "output_time": 5.,
        "minimum_timelike_interval_squared": float(interval.min()),
        "proper_time_invariance_residual": invariant_error,
        "proper_clock_state_residual": clock_error,
        "boosted_coordinate_output_times": outprime[:, 0].tolist(),
        "transported_cut_time_residual": residual(cuttime-5),
        "original_intrinsic_position_residual": residual(intrinsic-outputs[:, 1:]),
        "whole_clock_law_assumed_not_derived": True,
    }

    old, extended = clock_state(5), clock_state(9)
    td = trace_distance(old, extended)
    assert abs(td-abs(np.sin(2))) < 1e-12
    pad_u = np.diag(np.exp(-.5j*4*np.array([1., -1.])))
    padded = pad_u @ old @ pad_u.conj().T
    assert residual(padded-extended) < 1e-12
    # Wrongly realigning the boosted outputs adds unequal physical proper times.
    target_t = outprime[:, 0].max()
    extra_proper = (target_t-outprime[:, 0])/lam[0, 0]
    wrong_slice_difference = max(trace_distance(clock_state(5),
                                                clock_state(5+dt))
                                 for dt in extra_proper)
    assert wrong_slice_difference > .5
    groups["deadline_and_reslicing_are_observable"] = {
        "old_deadline": 5., "larger_envelope_deadline": 9.,
        "clock_trace_distance_when_deadlines_differ": td,
        "same_window_with_real_padding_residual": residual(padded-extended),
        "extra_proper_times_if_outputs_resliced": extra_proper.tolist(),
        "clock_difference_if_reslicing_ignored": wrong_slice_difference,
    }
    return {
        "round": 1103, "status": "PASS_FINITE_CONTEXT_CALIBRATION",
        "groups": groups,
        "source_sha256": {p: sha(ROOT/p) for p in SOURCES},
        "new_universal_axioms": 0, "new_scientific_experiments": 0,
        "scientific_count_increment": 0, "scientific_count_total": 3860,
        "joint_model_for_entire_SR_package_certified": False,
        "Lorentz_structure_derived_from_cognition": False,
        "uniform_deadline_for_unbounded_catalogs_claimed": False,
        "empirical_data": False,
    }


def compare(expected, actual):
    if isinstance(expected, dict):
        assert expected.keys() == actual.keys()
        for k in expected:
            compare(expected[k], actual[k])
    elif isinstance(expected, list):
        assert len(expected) == len(actual)
        for a, b in zip(expected, actual):
            compare(a, b)
    elif isinstance(expected, float):
        assert np.isclose(expected, actual, atol=1e-12, rtol=1e-12)
    else:
        assert expected == actual, (expected, actual)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = compute()
    path = HERE/"results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(json.loads(path.read_text(encoding="utf-8")), result)
    print(json.dumps({"status": result["status"], "groups": len(result["groups"]),
                      "saved_results_match": not args.write}))
