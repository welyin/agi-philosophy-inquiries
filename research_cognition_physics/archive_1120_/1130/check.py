"""1130: original clock / continuous-source calibration and common identity.

Default is read-only. --write-results registers this round's results.json only.
No old 288-dimensional device is executed. The Gaussian diagnostic is an
independent continuous canonical coordinate, not the actual scalar field.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = (
    "archive_1106_1119/1106/proof.md", "archive_1106_1119/1106/results.json",
    "archive_1106_1119/1119/proof.md", "archive_1106_1119/1119/results.json",
    "archive_1120_/1127/proof.md", "archive_1120_/1129/proof.md",
    "archive_1120_/1129/results.json",
)
ETA = np.diag([1., -1., -1., -1.])
TOL = 2e-10


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def clock():
    hc = 1.5 * np.eye(4)
    for j in range(3):
        hc[j, j+1] = hc[j+1, j] = .5 * np.sqrt((j+1)*(3-j))
    return hc, (15+2*np.sqrt(19))*np.eye(4)+hc


def sources():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    old = json.loads((ROOT/SOURCES[1]).read_text(encoding="utf-8"))
    affine = json.loads((ROOT/SOURCES[3]).read_text(encoding="utf-8"))
    field = json.loads((ROOT/SOURCES[6]).read_text(encoding="utf-8"))
    assert old["round"] == 1106 and old["status"].startswith("PASS")
    assert affine["round"] == 1119 and affine["status"].startswith("PASS")
    assert field["round"] == 1129 and field["status"] == "PASS"
    mat = old["group3_independently_defined_local_material"]
    assert mat["local_dimension"] == 288 and mat["data_dimension"] == 72
    assert mat["full_72_data_endpoint_residual"] < TOL
    assert mat["launch_spectrum"] == [0., 1., 2., 3.]
    assert mat["launch_probabilities"] == [.125, .375, .375, .125]
    row = old["group5_clock_record_and_real_receiver"]["clock_record_rows"][1]
    assert abs(row["tau"]-np.pi/2) < TOL
    assert abs(row["M_probability"]-.125) < TOL
    assert row["all_input_effect_residual"] < TOL
    assert F(field["exact_budget"]["fractions"]["source_error"]) == F(34584, 42724609375)
    return dict(source_sha256=hashes, frozen_old_clock_dimension=288,
                old_gamma_full_data_residual=mat["full_72_data_endpoint_residual"],
                old_actual_record_probability=row["M_probability"],
                old_actual_record_all_input_residual=row["all_input_effect_residual"],
                old_devices_reexecuted=False,
                scope="Frozen Gamma and original Mr identity imported; the four-dimensional block is a diagnostic only.")


def budgets():
    eps, amplitude, lam = F(1, 10**5), F(10**34), F(11, 7)
    time, duration, gamma, bound = F(5, 4), F(1, 4), F(5, 4), F(27)
    q = 8/eps**3
    c = 8*q*524
    material = 80*lam*eps
    source = lam*c*512*3*time/(eps*amplitude)
    back = lam*512*c*q*time**2/(2*amplitude**2)
    kappa = 3/(eps*amplitude)+back
    calibration = bound*kappa
    joint = material+source+calibration
    pi_back = lam*512*c*q*time/amplitude**2
    work_error = bound*pi_back
    work_quadratic = q*bound**2/(2*amplitude**2)
    derivative_kernel_product = 102/eps**4
    impulse = bound*lam*512*c*derivative_kernel_product*time**2/(2*amplitude**2)
    tail = lam**10/F(factorial(10))/(1-lam/11)
    causal = 4*tail+2*source
    d = (gamma-1)*time
    record_gap = ((1+d-d**3/6)/2)**3-F(1, 8)-2*calibration
    assert 5+time+eps < 7 and 2+time+eps < 4
    assert duration > 0 and 2 > time+20*eps
    assert calibration < F(86, 10**29)
    assert joint < F(1, 100) and causal+2*calibration < F(1, 1000)
    assert record_gap > F(154350301450, 10**12)
    items = dict(epsilon=eps, amplitude=amplitude, time=time, calibration_duration=duration,
                 old_clock_norm_upper=bound, q_upper=q, C_upper=c,
                 material_error=material, source_error=source,
                 normalized_field_backreaction=back, pointer_rms_bound=kappa,
                 calibration_joint_error=calibration, full_handoff_error=joint,
                 momentum_pointer_backreaction=pi_back, calibration_work_deviation=work_error,
                 calibration_work_quadratic_upper=work_quadratic,
                 derivative_kernel_product_upper=derivative_kernel_product,
                 calibration_spatial_impulse_component_upper=impulse,
                 finite_causal_menu=causal, joint_clock_causal_menu=causal+2*calibration,
                 wrong_old_clock_residual_proper_time=d,
                 actual_record_gap_lower=record_gap)
    return dict(scope="Exact rational upper and lower bounds; pi < 22/7. No quadrature certifies these budgets.",
                fractions={k: str(v) for k, v in items.items()},
                decimal_views={k: float(v) for k, v in items.items()},
                finite_positive_window=True, autonomous_microcontact_not_claimed=True)


def original_clock():
    hc, hd = clock()
    e, v = np.linalg.eigh(hd)
    ready = np.eye(4)[:, 0]
    sigma0, t, gamma = np.pi/2, 1.25, 1.25
    state0 = exp_h(hd, sigma0) @ ready
    after_wait = exp_h(hd, t) @ state0
    after_cal = exp_h(hd, -t) @ after_wait
    wrong = exp_h(hd, -t) @ exp_h(hd, gamma*t) @ state0
    correct_p = float(abs(after_cal[3])**2)
    wrong_p = float(abs(wrong[3])**2)
    closed_wrong = float(((1+np.sin(5/16))/2)**3)
    full_block_residual = norm(exp_h(hd, -t) @ exp_h(hd, t)-np.eye(4))
    spectrum_residual = norm(np.diag(e)-np.diag(15+2*np.sqrt(19)+np.arange(4)))
    weights = abs(v.conj().T @ ready)**2
    assert full_block_residual < TOL and spectrum_residual < TOL
    assert max(abs(weights-np.array([1, 3, 3, 1])/8)) < TOL
    assert abs(correct_p-.125) < TOL and abs(wrong_p-closed_wrong) < TOL
    assert np.pi/4 < sigma0 < sigma0+(gamma-1)*t < np.pi
    assert sigma0+gamma*t < np.pi
    expectation = float(np.vdot(state0, hd @ state0).real)
    second = float(np.vdot(state0, hd @ hd @ state0).real)
    return dict(scope="Original Mr compressed by frozen Gamma on the declared old preparation; not a new head detector.",
                old_clock_mass_spectrum=e.tolist(), spectral_residual=spectrum_residual,
                ready_spectral_weights=weights.tolist(), sigma0=sigma0,
                source_proper_time=t, positive_contact_duration=.25,
                initial_record_probability=float(abs(state0[3])**2),
                correct_record_probability=correct_p, wrong_record_probability=wrong_p,
                actual_record_nominal_gap=wrong_p-correct_p,
                full_clock_spectral_block_echo_residual=full_block_residual,
                prepared_state_echo_residual=float(np.linalg.norm(after_cal-state0)),
                complete_clock_energy_mean=expectation, complete_clock_energy_second_moment=second,
                no_intermediate_record_measurement=True,
                source_time_held_fixed_in_wrong_dictionary=True)


def boost_and_recipe():
    speed, gamma, t = .6, 1.25, 1.25
    lam = np.eye(4)
    lam[0, 0] = lam[3, 3] = gamma
    lam[0, 3] = lam[3, 0] = gamma*speed
    n = np.array([1., 0., 0., 0.])
    moved_n = lam @ n
    z = t*n
    moved_z = lam @ z
    lorentz_residual = norm(lam.T @ ETA @ lam-ETA)
    tangent_gram = -lam[:, 1:].T @ ETA @ lam[:, 1:]
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    alpha = [np.kron(sx, p) for p in (sx, sy, sz)]
    beta = np.kron(sz, np.eye(2))
    rapidity = np.arctanh(speed)
    spin = np.cosh(rapidity/2)*np.eye(4)+np.sinh(rapidity/2)*alpha[2]
    flux = gamma*(np.eye(4)-speed*alpha[2])
    flux_residual = norm(spin.conj().T @ flux @ spin-np.eye(4))
    scalar_residual = norm(spin.conj().T @ beta @ spin-beta)
    _, ev = np.linalg.eigh(alpha[2])
    psi = ev[:, -1]
    moved = spin @ psi
    current = np.array([np.vdot(psi, psi).real]+[np.vdot(psi, a @ psi).real for a in alpha])
    moved_current = np.array([np.vdot(moved, moved).real]+[np.vdot(moved, a @ moved).real for a in alpha])
    current_residual = float(np.linalg.norm(moved_current-lam @ current))
    normal_density = float(moved_n @ ETA @ moved_current)
    wrong_lab_density = float(moved_current[0])
    intrinsic_points = [np.array([0., .2, .1, .4]), np.array([0., -.3, .2, .1])]
    kernel_geometry_error = max(abs(x @ ETA @ x-(lam @ x) @ ETA @ (lam @ x)) for x in intrinsic_points)
    normalized_radius, wrong_radius = .4, gamma*.4
    kernel = lambda r: float(np.exp(-1/(1-r*r)))
    hc, hd = clock()
    e, p = gamma*hd, gamma*speed*hd
    dx, dt = moved_z[3], moved_z[0]
    port_identity = norm(exp_h(p, dx) @ exp_h(hd, t)-exp_h(e, dt))
    old = json.loads((ROOT/SOURCES[1]).read_text(encoding="utf-8"))
    events = old["group1_material_and_two_port_events"]["events"][0]["events_o_a_b_f"]
    anchor = np.array([events[1][0]+np.pi/2, events[1][1], 0., 0.])
    old_history = [np.array([x[0], x[1], 0., 0.])-anchor for x in events[:3]]
    moved_history = [lam @ x for x in old_history]
    assert max(lorentz_residual, flux_residual, scalar_residual, current_residual,
               kernel_geometry_error, port_identity, norm(tangent_gram-np.eye(3))) < TOL
    assert abs(moved_n @ ETA @ moved_z-t) < TOL
    assert abs(normal_density-1) < TOL and abs(wrong_lab_density-2) < TOL
    assert old_history[-1][0] < 0
    return dict(scope="Finite boost checks the independently stated flux, kernel and original-port recipe; analytic proof covers Fock and arbitrary inputs.",
                lorentz_residual=lorentz_residual, spinor_flux_residual=flux_residual,
                spinor_scalar_residual=scalar_residual, current_vector_residual=current_residual,
                boosted_normal=moved_n.tolist(), same_completion_event=moved_z.tolist(),
                proper_elapsed=float(moved_n @ ETA @ moved_z), lab_elapsed=float(moved_z[0]),
                correct_normal_current_density=normal_density, incorrectly_reused_lab_time_density=wrong_lab_density,
                spin_not_ordinary_euclidean_unitary=norm(spin.conj().T @ spin-np.eye(4)),
                induced_kernel_metric_residual=float(kernel_geometry_error),
                correct_intrinsic_bump=kernel(normalized_radius), wrong_lab_resliced_bump=kernel(wrong_radius),
                original_port_phase_identity_residual=port_identity,
                original_clock_history_relative_to_field_start=[x.tolist() for x in old_history],
                same_boosted_old_history=[x.tolist() for x in moved_history],
                complete_fock_phase_and_source_data_must_transport=True,
                actual_predeployment_is_declared_material_permission=True,
                interacting_total_four_momentum_not_constructed=True)


def gaussian_contact():
    # A genuine L2(R) Gaussian at the contact: q = A*T + x, Var(x)=1/2,
    # momentum mean A. No finite-dimensional CCR representation is used.
    hc, hd = clock()
    energy, modes = np.linalg.eigh(hd)
    input_clock = exp_h(hd, np.pi/2) @ np.eye(4)[:, 0]
    amplitudes = modes.conj().T @ input_clock
    weights = abs(amplitudes)**2
    nodes, quadrature = np.polynomial.hermite.hermgauss(96)
    quadrature = quadrature/np.sqrt(np.pi)
    mean_energy = float(weights @ energy)
    second_energy = float(weights @ (energy**2))
    rows = []
    for amplitude in (40., 80., 160.):
        phases = np.exp(1j*np.outer(nodes/amplitude, energy))
        kernel_numeric = np.einsum("r,rn,rm->nm", quadrature, phases, phases.conj())
        kernel_exact = np.exp(-(energy[:, None]-energy[None, :])**2/(4*amplitude**2))
        gram_error = float(np.max(abs(kernel_numeric-kernel_exact)))
        prepared_error_sq_numeric = float(quadrature @ (abs(phases-1)**2 @ weights))
        prepared_error_sq_exact = float(weights @ (-2*np.expm1(-energy**2/(4*amplitude**2))))
        full_isometry_error = float(np.sqrt(max(-2*np.expm1(-energy**2/(4*amplitude**2)))))
        bound = float(max(energy)/(np.sqrt(2)*amplitude))
        reduced = modes @ (np.outer(amplitudes, amplitudes.conj())*kernel_numeric) @ modes.conj().T
        mr_probability = float(reduced[3, 3].real)
        # Maximal-reference output in the energy basis checks all matrix units.
        choi = np.zeros((16, 16), complex)
        choi_exact = np.zeros_like(choi)
        for i in range(4):
            for j in range(4):
                choi[5*i, 5*j] = kernel_numeric[i, j]/4
                choi_exact[5*i, 5*j] = kernel_exact[i, j]/4
        reference_error = norm(choi-choi_exact)
        work = mean_energy+second_energy/(2*amplitude**2)
        m0 = 15+2*np.sqrt(19)
        bare_energy = energy-m0
        bare_work = float(weights @ bare_energy+weights @ (bare_energy**2)/(2*amplitude**2))
        assert gram_error < TOL and reference_error < TOL
        assert abs(prepared_error_sq_numeric-prepared_error_sq_exact) < TOL
        assert full_isometry_error <= bound+TOL
        assert abs(np.trace(choi)-1) < TOL and min(np.linalg.eigvalsh(choi)) > -TOL
        rows.append(dict(amplitude=amplitude, complete_prepared_isometry_error=float(np.sqrt(prepared_error_sq_exact)),
                         full_clock_isometry_error=full_isometry_error, rms_lipschitz_bound=bound,
                         analytic_vs_hermite_gram_error=gram_error, all_matrix_units_reference_error=reference_error,
                         old_record_probability=mr_probability, complete_contact_work=work,
                         incorrectly_dropped_mass_work=bare_work,
                         source_momentum_shift_mean=mean_energy/amplitude))
    assert rows[-1]["complete_prepared_isometry_error"] < rows[0]["complete_prepared_isometry_error"]
    assert all(r["complete_contact_work"] > 25 for r in rows)
    return dict(scope="Continuous canonical Gaussian confirms Weyl sign, full-clock Gram and work mechanism only; not actual three-dimensional scalar-field certification.",
                coordinate_variance=.5, mean_complete_clock_energy=mean_energy,
                mean_complete_clock_energy_squared=second_energy, rows=rows,
                finite_CCR_cutoff_used=False, old_full_clock_replaced=False,
                tiny_state_error_does_not_make_contact_work_free=True)


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
    result = dict(round=1130, schema=1, status="PASS_CONDITIONAL_OLD_CLOCK_CONTINUUM_CALIBRATION",
                  new_adopted_axioms=0, scientific_calibration=0, stage="active", numerical_check_groups=5,
                  source_identity=sources(), exact_budgets=budgets(), original_clock=original_clock(),
                  hyperplane_recipe=boost_and_recipe(), continuous_gaussian_contact=gaussian_contact())
    target = HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print("PASS: five groups; only 1130/results.json registered explicitly.")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        print("PASS: five groups; frozen source hashes and recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
