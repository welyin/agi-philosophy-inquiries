"""1131: two normal source preparations and the original bare B effect.

Default is read-only. --write-results writes only this round's results.json.
Finite CAR represents the exact restriction of the old sea to four modes;
it is not a continuum field simulation or a constructed gradient QFT.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOL = 2e-10
FROZEN = {
    "archive_1120_/1127/proof.md": "d1c8827cc12563b996c951da6039e041e08aa711de18e84a56bb240660660c27",
    "archive_1120_/1129/proof.md": "c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa",
    "archive_1120_/1129/results.json": "69ed793bd5847d31003eb05b18205ecf72af34c26f9ad2919ebe548e7a25a02c",
    "archive_1120_/1130/proof.md": "9e4c2c64f062db6ef61777f532753a2a9cf3b0ab13df4280dcc3c6c35519a6ee",
    "archive_1120_/1130/results.json": "e8cc423302b5fda96d509b5a3f04cbe572874cb22955877d68be5602e211ba62",
    "archive_554_584/research_note_564.md": "3823a4afa317b562a2794d0eaf8163f0798c40d739a6908dafaa684693151955",
    "archive_585_628/research_note_606.md": "4cfbece6aa871e78d8b2d5e8ae615fb8674ca2c3b1ced94c1f15b32cd29ad0f1",
    "archive_1009_1043/research_note_1028.md": "1a28a50f0836b490d376231bfd483a58c473a8cd5b6b87f3174cb2f20b867503",
}


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e)) @ v.conj().T


def sources():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes == FROZEN
    field = json.loads((ROOT/"archive_1120_/1129/results.json").read_text(encoding="utf-8"))
    clock = json.loads((ROOT/"archive_1120_/1130/results.json").read_text(encoding="utf-8"))
    assert field["round"] == 1129 and field["status"] == "PASS"
    assert clock["round"] == 1130 and clock["status"].startswith("PASS")
    assert clock["source_identity"]["frozen_old_clock_dimension"] == 288
    assert clock["exact_budgets"]["fractions"]["time"] == "5/4"
    assert clock["exact_budgets"]["fractions"]["material_error"] == "11/8750"
    return dict(source_sha256=hashes, old_devices_reexecuted=False,
                same_autonomous_H=True, same_amplitude_parameter=True,
                old_complete_clock_dimension=288, same_original_bare_B_effect=True,
                scope="Frozen normal CAR/source and complete old clock imported; no old author file is modified.")


def budgets():
    eps, amp, lam, time = F(1, 10**5), F(10**34), F(11, 7), F(5, 4)
    q, c = 8/eps**3, 8*(8/eps**3)*524
    material = 80*lam*eps
    src = lam*c*512*3*time/(eps*amp)
    e1, e2 = material+src, material/2+src
    lower = F(1, 4)-e2
    assert 5+time+eps < 7 and 2+time+eps < 4
    assert lower > F(249370416742399, 10**15)
    assert lower > F(1, 100)+F(1, 100)
    items = dict(epsilon=eps, fixed_amplitude=amp, time=time, q_upper=q, C_upper=c,
                 material_error_s1=material, material_error_s2=material/2,
                 source_joint_error=src, effect_error_s1=e1, effect_error_s2=e2,
                 probability_s1_lower=1-e1, probability_s1_upper=F(1),
                 probability_s2_lower=F(3, 4)-e2, probability_s2_upper=F(3, 4)+e2,
                 two_actual_preparations_gap_lower=F(1, 4)-e1-e2,
                 mismatch_lower_before_kappa=lower,
                 illustrative_match_budget=F(1, 100), illustrative_gradient_budget=F(1, 100),
                 source_energy_upper_s1=1073*amp**2,
                 source_energy_upper_s2=4292*amp**2,
                 old_contact_work_deviation_upper=27*lam*512*c*q*time/amp**2)
    return dict(scope="Exact rational bounds using pi < 22/7; candidate kappa_grad remains unassigned.",
                fractions={k: str(v) for k, v in items.items()},
                decimal_views={k: float(v) for k, v in items.items()},
                rate_menu=[1, 2], source_energy_ratio=4,
                new_material_permission="M1131: two independent normal coherent preparations",
                candidate_gradient_error_proved=False, kappa_grad=None,
                same_W_commutes_with_B_effect=True,
                extra_clock_error_needed_for_B_effect=False,
                s2_old_Mr_reset_to_sigma0_claimed=False)


def annihilator(mode, modes=4):
    out = np.zeros((2**modes, 2**modes), complex)
    for bits in range(2**modes):
        if (bits >> mode) & 1:
            parity = (bits & ((1 << mode)-1)).bit_count()
            out[bits ^ (1 << mode), bits] = (-1)**parity
    return out


def normal_car_reference():
    # Modes A1,A2,B1,B2. Every occupation sector is retained.
    ca = [annihilator(i) for i in range(4)]
    ident = np.eye(16)
    car = max(norm(ca[i] @ ca[j].conj().T+ca[j].conj().T @ ca[i]
                   -(ident if i == j else 0*ident)) for i in range(4) for j in range(4))
    n = [x.conj().T @ x for x in ca]
    sea = ident/16
    empty_a = (ident-n[0]) @ (ident-n[1])
    success = float(np.trace(empty_a @ sea).real)
    cond = empty_a @ sea @ empty_a/success
    assert abs(success-.25) < TOL
    assert max(abs(np.trace(n[i] @ cond).real-.5) for i in (2, 3)) < TOL
    # Exact restricted input channel and a purification of its B sea remainder.
    maps = []
    joint = np.zeros((64, 2), complex)
    for b in range(4):
        ready = np.eye(16)[:, b << 2]
        v = np.column_stack([ca[0].conj().T @ ready, ca[1].conj().T @ ready])
        maps.append(v)
        for i in range(2):
            joint[:, i] += np.kron(v[:, i], np.eye(4)[:, b])/2
    assert norm(joint.conj().T @ joint-np.eye(2)) < TOL
    hop = sum(ca[i].conj().T @ ca[i+2]+ca[i+2].conj().T @ ca[i] for i in range(2))
    parity = np.diag([(-1)**bits.bit_count() for bits in range(16)])
    assert norm(hop @ parity-parity @ hop) < TOL
    assert norm(n[2] @ parity-parity @ n[2]) < TOL
    rows, max_units, max_reference, max_heisenberg = [], 0., 0., 0.
    for alpha in [0., np.pi/7, np.pi/4, np.pi/2]:
        unitary = exp_h(hop, alpha)
        pullback = unitary.conj().T @ n[2] @ unitary
        compressed = sum(v.conj().T @ pullback @ v for v in maps)/4
        expected = .5*np.cos(alpha)**2*np.eye(2)+np.sin(alpha)**2*np.diag([1., 0.])
        # Annihilators are antilinear in one-particle modes.
        cv = np.cos(alpha)*ca[2]-1j*np.sin(alpha)*ca[0]
        heis = norm(pullback-cv.conj().T @ cv)
        max_heisenberg = max(max_heisenberg, heis)
        assert heis < TOL and norm(compressed-expected) < TOL
        for i in range(2):
            for j in range(2):
                eij = np.zeros((2, 2), complex)
                eij[i, j] = 1
                rho = sum(v @ eij @ v.conj().T for v in maps)/4
                observed = np.trace(n[2] @ unitary @ rho @ unitary.conj().T)
                residual = float(abs(observed-np.trace(expected @ eij)))
                max_units = max(max_units, residual)
        # Bell input on Q/R with a retained pure B-environment, 16*4*2 dimensions.
        bell = np.zeros((2, 2), complex)
        bell[0, 0] = bell[1, 1] = 1/np.sqrt(2)
        state = (np.kron(unitary, np.eye(4)) @ joint @ bell).reshape(16, 4, 2)
        effect_state = np.einsum('ab,ber->aer', n[2], state)
        ref_record = np.einsum('aer,aes->rs', effect_state, effect_state.conj())
        expected_ref = expected.T/2
        reference_error = norm(ref_record-expected_ref)
        max_reference = max(max_reference, reference_error)
        assert reference_error < TOL
        rows.append(dict(angle=alpha, rail1_probability=float(compressed[0, 0].real),
                         rail2_probability=float(compressed[1, 1].real),
                         complete_effect_residual=norm(compressed-expected),
                         reference_record_residual=reference_error))
    original = .25*np.eye(2)+.5*np.diag([1., 0.])
    grad = np.diag([1., 0.])
    assert abs(norm(grad-original)-.25) < TOL
    assert car < TOL and max_units < TOL and max_reference < TOL
    return dict(scope="Exact four-mode restriction with B sea density I/4 and retained purifier; all occupation sectors, matrix units and Bell reference checked. Not a continuum cutoff simulation.",
                CAR_modes=4, full_Fock_dimension=16, retained_sea_purifier_dimension=4,
                passive_reference_dimension=2, CAR_residual=car,
                original_A_empty_projection_probability=success,
                initial_B_occupation=.5, all_matrix_units_residual=max_units,
                arbitrary_angle_Heisenberg_residual=max_heisenberg,
                Bell_reference_residual=max_reference, rows=rows,
                s2_original_probability=.75, gradient_endpoint_probability=1.,
                bare_effect_gap=.25, incorrect_empty_B_prediction=.5,
                complete_effect_norm_gap=norm(grad-original))


def integrate(fun, lo, hi, order=160):
    z, w = np.polynomial.legendre.leggauss(order)
    x = lo+(hi-lo)*(z+1)/2
    return float((hi-lo)/2*(w @ np.array([fun(t) for t in x])))


def gradient_geometry():
    def bump(x):
        return float(np.exp(-1/(x*(1-x)))) if 0 < x < 1 else 0.
    normalization = integrate(bump, 0., 1., 256)
    def gate(theta):
        return 2*bump(2*theta-.5)/normalization
    def primitive(theta):
        if theta <= .25:
            return 0.
        if theta >= .75:
            return 1.
        return integrate(gate, .25, theta)
    area = integrate(gate, .25, .75)
    assert abs(area-1) < TOL
    rows = []
    for rate in (1., 2.):
        base = integrate(lambda t: gate(rate*t), .25/rate, .75/rate)
        pure = integrate(lambda t: rate*gate(rate*t), .25/rate, .75/rate)
        normalized = integrate(lambda t: rate/np.sqrt(rate**2)*gate(rate*t),
                               .25/rate, .75/rate)
        endpoint = primitive(rate*1.25)-primitive(0)
        angles = np.pi/2*np.array([base, pure, normalized, endpoint])
        probabilities = 1-.5*np.cos(angles)**2
        assert abs(base-1/rate) < TOL and abs(pure-endpoint) < TOL
        assert abs(normalized-base) < TOL
        rows.append(dict(source_rate=rate, original_angle=float(angles[0]),
                         pure_gradient_angle=float(angles[1]), normalized_gradient_angle=float(angles[2]),
                         primitive_endpoint_angle=float(angles[3]),
                         original_B_probability=float(probabilities[0]),
                         pure_gradient_B_probability=float(probabilities[1]),
                         normalized_gradient_B_probability=float(probabilities[2])))
    # Nonuniform timelike test jet theta=t(1+c*x). This is a classical identity
    # diagnostic, not an added actual M1131 preparation or a third no-go witness.
    t, x, c, d, step = .4, .2, .3, .4, 2e-6
    def data(tt, xx):
        theta = tt*(1+c*xx)
        v = np.array([1+c*xx, c*tt])
        X = v[0]**2-v[1]**2
        return theta, v, X
    def pure_form(tt, xx):
        theta, v, _ = data(tt, xx)
        return gate(theta)*v
    def normalized_form(tt, xx):
        theta, v, X = data(tt, xx)
        return gate(theta)*v/np.sqrt(X)
    def eta_form(tt, xx):
        return (1+d*xx)*pure_form(tt, xx)
    def completed_eta_form(tt, xx):
        theta, _, _ = data(tt, xx)
        return eta_form(tt, xx)+np.array([0., d*primitive(theta)])
    def curl(form):
        dt = (form(t+step, x)[1]-form(t-step, x)[1])/(2*step)
        dx = (form(t, x+step)[0]-form(t, x-step)[0])/(2*step)
        return float(dt-dx)
    theta, v, X = data(t, x)
    dX = np.array([-2*c*c*t, 2*c*(1+c*x)])
    norm_exact = -gate(theta)*(dX[0]*v[1]-dX[1]*v[0])/(2*X**1.5)
    eta_exact = -d*gate(theta)*v[0]
    pure_curl, norm_curl = curl(pure_form), curl(normalized_form)
    eta_curl, completed_curl = curl(eta_form), curl(completed_eta_form)
    assert abs(pure_curl) < 2e-7 and abs(completed_curl) < 2e-7
    assert abs(norm_curl-norm_exact) < 2e-7 and abs(norm_exact) > .1
    assert abs(eta_curl-eta_exact) < 2e-7 and abs(eta_exact) > .1
    return dict(scope="Actual bump areas plus independent finite-difference checks of classical gradient/curl identities. No normal gradient QFT is constructed.",
                actual_gate_integral=area, rate_rows=rows,
                nonuniform_test_jet=dict(theta=theta, timelike_X=X,
                    pure_gradient_curl=pure_curl,
                    normalized_gradient_curl=norm_curl, normalized_gradient_curl_analytic=float(norm_exact),
                    eta_times_gradient_curl=eta_curl, eta_times_gradient_curl_analytic=float(eta_exact),
                    completed_eta_gradient_curl=completed_curl),
                normalized_gradient_escapes_two_rate_witness=True,
                nonuniform_jet_added_to_actual_preparation_menu=False,
                pure_gradient_field_redefinition_used_as_mature_tool=True)


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual, expected, atol=1e-10, rtol=1e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = dict(round=1131, schema=1, status="PASS_CONDITIONAL_PURE_GRADIENT_MATCHING_OBSTRUCTION",
                  new_adopted_axioms=0, scientific_calibration=0, stage="active",
                  numerical_check_groups=4, source_identity=sources(), exact_budgets=budgets(),
                  normal_CAR_and_reference=normal_car_reference(), gradient_geometry=gradient_geometry())
    target = HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print("PASS: four groups; only 1131/results.json registered explicitly.")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        print("PASS: four groups; eight frozen source hashes and recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
