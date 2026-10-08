"""1019: positive scalar CAR, separated smearings, and a finite signal witness.

Analytic support/kernel bounds and spectral-tail bounds are separate from
non-certified Gauss--Legendre quadrature calibration. The four-mode CAR
algebra is an exact finite subalgebra, not an invariant energy truncation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "spin_statistics_locality_results.json"
TAU = .25
A_RADIUS = 1.
PI2 = math.pi**2


def ball_transform(z):
    z = np.asarray(z, float)
    result = np.empty_like(z)
    small = np.abs(z) < .05
    s = z[small]
    result[small] = sum((-1)**k*3*(2*k+2)/math.factorial(2*k+3)*s**(2*k)
                        for k in range(6))
    t = z[~small]
    result[~small] = 3*(np.sin(t)-t*np.cos(t))/t**3
    return result


def quadrature(cutoff, panel_width, order, separation):
    count = int(math.ceil(cutoff/panel_width))
    edges = np.linspace(0., cutoff, count+1)
    z, w = leggauss(order)
    half = (edges[1:]-edges[:-1])/2
    mid = (edges[1:]+edges[:-1])/2
    p = (mid[:, None]+half[:, None]*z).ravel()
    weights = (half[:, None]*w).ravel()
    common = ball_transform(p)**4*np.sinc(p*TAU/(2*math.pi))**4/(4*PI2)
    radial_cross = np.sinc(p*separation/math.pi)
    moments = [float(weights @ (p**(k+1)*common)) for k in range(3)]
    cross = [float(weights @ (p**(k+1)*common*radial_cross)) for k in range(3)]
    tails = [324/(PI2*(6-k)*cutoff**(6-k)) for k in range(3)]
    return dict(cutoff=cutoff, panel_width=panel_width, gauss_order=order,
                points=len(p), moments= moments, cross_moments=cross,
                overlap=cross[0]/moments[0], analytic_absolute_tail_bounds=tails,
                quadrature_error_certified=False)


def analytical_bounds(separation):
    n_low = float(F(9, 10)**4*F(383, 384)**4)/(8*PI2)
    n_up = 1/PI2
    s_low = 1/(4*PI2*(separation+4)**2)
    s_up = 1/(4*PI2*((separation-4)**2-(2*TAU)**2))
    r_low, r_up = s_low/n_up, min(s_up/n_low, 1.)
    assert 0 < r_low < r_up < 1
    d_low = r_low*math.sqrt(1-r_up*r_up)/math.sqrt(2)
    theta = d_low/4
    signal = d_low*d_low/8
    cutoff = 2000
    tail = 54/(PI2*cutoff**6)
    epsilon = tail/n_low
    probability_error = 2*(1+abs(theta))*math.sqrt(epsilon)
    assert probability_error < signal/2
    return dict(n_lower=n_low, n_upper=n_up, s_lower=s_low, s_upper=s_up,
                r_lower=r_low, r_upper=r_up, d_lower=d_low,
                fixed_theta=theta, exact_continuum_probability_lower_bound=signal,
                inequalities_derived_without_numerical_quadrature=True,
                numerical_values_of_analytic_bounds_only=True,
                cutoff_control=dict(cutoff=cutoff, normalization_tail_bound=tail,
                    normalized_lost_norm_squared_upper_bound=epsilon,
                    projection_operator_norm_error_upper_bound=math.sqrt(epsilon),
                    same_state_delta_probability_error_upper_bound=probability_error,
                    error_over_signal_bound=probability_error/signal,
                    error_less_than_half_signal=True,
                    direct_quadrature_to_cutoff_performed=False,
                    cutoff_operations_strictly_local=False,
                    all_states_assumed_bandlimited=False))


def annihilators(n):
    result = []
    for mode in range(n):
        op = np.zeros((2**n, 2**n), complex)
        for state in range(2**n):
            if state & (1 << mode):
                target = state ^ (1 << mode)
                sign = (-1)**((state & ((1 << mode)-1)).bit_count())
                op[target, state] = sign
        result.append(op)
    return result


def commutator(a, b):
    return a @ b-b @ a


def norm(a):
    return float(np.linalg.norm(a, 2))


def expectation(state, operator):
    value = np.vdot(state, operator @ state)
    assert abs(value.imag) < 2e-11
    return float(value.real)


def probability_witness(P, Q, omega, theta):
    C = 1j*commutator(P, Q)
    d = norm(C @ omega)
    assert abs(expectation(omega, C)) < 2e-14 and d > 0
    candidates = [(omega+sign*C @ omega/d)/math.sqrt(2) for sign in (1, -1)]
    derivatives = [expectation(state, C) for state in candidates]
    index = int(np.argmax(np.abs(derivatives)))
    state = candidates[index]
    U = np.eye(P.shape[0])+(np.exp(-1j*theta)-1)*P
    after = U @ state
    initial_p, final_p = expectation(state, Q), expectation(after, Q)
    difference = final_p-initial_p
    derivative = derivatives[index]
    error = abs(difference-theta*derivative)
    assert norm(U.conj().T @ U-np.eye(P.shape[0])) < 3e-14
    assert abs(np.vdot(state, state)-1) < 3e-14
    assert abs(abs(derivative)-d) < 3e-14
    assert 0 <= initial_p <= 1 and 0 <= final_p <= 1
    assert error <= 2*theta*theta+3e-15
    assert abs(difference) >= d*theta/2-3e-15
    return dict(theta=theta, d=d, candidate_C_expectations=derivatives,
                chosen_sign=1 if index == 0 else -1,
                initial_probability=initial_p, final_probability=final_p,
                probability_difference=difference, absolute_difference=abs(difference),
                linear_response=theta*derivative, measured_Taylor_remainder=error,
                Taylor_remainder_upper_bound=2*theta*theta,
                finite_signal_lower_bound_using_computed_d=d*theta/2,
                postselection_used=False), state, after


def one_particle_energy_matrices(moments, cross):
    n = moments[0]
    r = cross[0]/n
    d = math.sqrt(1-r*r)
    transform = np.array([[1., -r/d], [0., 1/d]])
    matrices = []
    for k in range(3):
        gram = np.array([[moments[k], cross[k]], [cross[k], moments[k]]])/n
        matrices.append(transform.T @ gram @ transform)
    assert norm(matrices[0]-np.eye(2)) < 3e-14
    leak = matrices[2]-matrices[1] @ matrices[1]
    assert np.min(np.linalg.eigvalsh(matrices[1])) > 0
    assert np.min(np.linalg.eigvalsh(leak)) > 0
    return matrices[1], matrices[2], leak


def second_quantized(matrix, operators):
    result = np.zeros_like(operators[0])
    for i in range(len(operators)):
        for j in range(len(operators)):
            result += matrix[i, j]*operators[i].conj().T @ operators[j]
    return result


def energy_check(operators, h1, h2, leak, states):
    block1 = np.zeros((4, 4)); block2 = np.zeros((4, 4)); block_leak = np.zeros((4, 4))
    for start in (0, 2):
        block1[start:start+2, start:start+2] = h1
        block2[start:start+2, start:start+2] = h2
        block_leak[start:start+2, start:start+2] = leak
    H_compressed = second_quantized(block1, operators)
    # On Fock(P), compression of H^2 includes excursions out of the packet
    # subspace: P H^2 P = (dGamma(P omega P))^2 + dGamma(P omega (1-P) omega P).
    full_H2_compressed = H_compressed @ H_compressed+second_quantized(block_leak, operators)
    rows = []
    for name, state in states:
        mean = expectation(state, H_compressed)
        second = expectation(state, full_H2_compressed)
        assert mean >= -2e-14 and second >= mean*mean-2e-12
        rows.append(dict(state=name, mean_energy=mean, energy_second_moment=second,
                         energy_variance=second-mean*mean))
    return dict(one_particle_H_compression=h1.tolist(),
                one_particle_H_squared_compression=h2.tolist(),
                leakage_matrix=leak.tolist(),
                leakage_eigenvalues=np.linalg.eigvalsh(leak).tolist(),
                states=rows, positive_total_H="dGamma(|p|) on the full particle-antiparticle Fock space",
                two_wavepacket_space_energy_invariant=False,
                local_Klein_Gordon_energy_density_assumed_positive=False,
                finite_moments_are_analytic_tail_consequences=True,
                exact_energy_values_are_quadrature_calibrations=True)


def finite_car_calibration(overlap, quadrature_row, bounds):
    ops = annihilators(4)
    eye = np.eye(16)
    max_car = 0.
    for i, ai in enumerate(ops):
        for j, aj in enumerate(ops):
            max_car = max(max_car, norm(ai @ aj+aj @ ai),
                          norm(ai @ aj.conj().T+aj.conj().T @ ai-(eye if i == j else 0)))
    assert max_car == 0
    u = np.array([1., 0.])
    v = np.array([overlap, math.sqrt(1-overlap*overlap)])
    fields = [sum(vector[i]*(ops[i]+ops[i+2].conj().T) for i in range(2))
              for vector in (u, v)]
    P, Q = [field.conj().T @ field/2 for field in fields]
    charge = sum((1 if i < 2 else -1)*op.conj().T @ op for i, op in enumerate(ops))
    parity = np.diag([(-1)**i.bit_count() for i in range(16)])
    projector_error = max(norm(P @ P-P), norm(Q @ Q-Q))
    charge_error = max(norm(commutator(P, charge)), norm(commutator(Q, charge)))
    parity_error = max(norm(commutator(P, parity)), norm(commutator(Q, parity)))
    assert max(projector_error, charge_error, parity_error) < 3e-14
    omega = np.zeros(16, complex); omega[0] = 1
    bracket = commutator(P, Q)
    vacuum_norm_squared = norm(bracket @ omega)**2
    operator_norm = norm(bracket)
    assert abs(vacuum_norm_squared-overlap**2*(1-overlap**2)/2) < 3e-14
    assert abs(operator_norm-overlap*math.sqrt(1-overlap**2)) < 3e-14
    d = math.sqrt(vacuum_norm_squared)
    adapted, state_a, after_a = probability_witness(P, Q, omega, d/4)
    fixed, state_f, after_f = probability_witness(P, Q, omega, bounds["fixed_theta"])
    assert fixed["absolute_difference"] >= bounds["exact_continuum_probability_lower_bound"]
    states = [("psi_adapted", state_a), ("U_adapted_psi", after_a),
              ("psi_fixed", state_f), ("U_fixed_psi", after_f)]
    for name, state in [("before", state_f), ("after", after_f)]:
        for label, effect in [("yes", Q), ("no", eye-Q)]:
            branch = effect @ state
            probability = float(np.vdot(branch, branch).real)
            assert probability > 0
            states.append((f"Q_{label}_{name}", branch/math.sqrt(probability)))
    h1, h2, leak = one_particle_energy_matrices(quadrature_row["moments"],
                                              quadrature_row["cross_moments"])
    return dict(matrix_dimension=16, number_of_CAR_modes=4, overlap=overlap,
                maximum_CAR_error=max_car, projector_error=projector_error,
                charge_zero_commutator_error=charge_error,
                even_parity_commutator_error=parity_error,
                vacuum_commutator_norm_squared=vacuum_norm_squared,
                expected_vacuum_commutator_norm_squared=overlap**2*(1-overlap**2)/2,
                commutator_operator_norm=operator_norm,
                expected_commutator_operator_norm=overlap*math.sqrt(1-overlap**2),
                computed_d=d, adapted_theta_witness=adapted,
                fixed_analytic_theta_witness=fixed,
                finite_energy=energy_check(ops, h1, h2, leak, states))


def cutoff_projection_algebra():
    ops = annihilators(2)
    rows = []
    for epsilon in (1e-6, 1e-3, .25, .8):
        a = ops[0]
        b = math.sqrt(1-epsilon)*ops[0]+math.sqrt(epsilon)*ops[1]
        P, Q = a.conj().T @ a, b.conj().T @ b
        observed = norm(P-Q)
        assert abs(observed-math.sqrt(epsilon)) < 2e-14
        rows.append(dict(lost_norm_squared=epsilon, projection_difference_norm=observed,
                         expected=math.sqrt(epsilon)))
    return dict(samples=rows, exact_two_CAR_mode_identity="||c†c-d†d||=sqrt(1-|{c,d†}|^2)",
                samples_are_identity_calibrations_not_analytic_proof=True)


def compare(fresh, saved, path="root"):
    if isinstance(fresh, float):
        assert isinstance(saved, (int, float)) and math.isclose(fresh, saved, rel_tol=3e-10, abs_tol=2e-13), path
    elif isinstance(fresh, dict):
        assert isinstance(saved, dict) and fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert isinstance(saved, list) and len(fresh) == len(saved), path
        for i, (a, b) in enumerate(zip(fresh, saved)):
            compare(a, b, path+f"[{i}]")
    else:
        assert fresh == saved, (path, fresh, saved)


def run():
    rows = []
    for separation in (6, 8):
        assert separation > 4*A_RADIUS+2*TAU
        qs = [quadrature(cutoff, width, order, separation)
              for cutoff, width, order in ((64., 1., 16), (128., 1., 32), (256., .5, 32))]
        bounds = analytical_bounds(separation)
        final = qs[-1]
        r = final["overlap"]
        assert bounds["n_lower"] < final["moments"][0] < bounds["n_upper"]
        assert bounds["s_lower"] < final["cross_moments"][0] < bounds["s_upper"]
        assert bounds["r_lower"] < r < bounds["r_upper"]
        differences = [dict(from_cutoff=qs[i]["cutoff"], to_cutoff=qs[i+1]["cutoff"],
                            maximum_moment_difference=max(abs(a-b) for a, b in
                                zip(qs[i]["moments"]+qs[i]["cross_moments"],
                                    qs[i+1]["moments"]+qs[i+1]["cross_moments"])),
                            overlap_difference=abs(qs[i]["overlap"]-qs[i+1]["overlap"]))
                       for i in range(2)]
        assert differences[-1]["maximum_moment_difference"] < 1e-9
        rows.append(dict(separation=separation, strict_spacelike_support_gap=separation-4-2*TAU,
                         quadrature=qs, quadrature_refinement_differences=differences,
                         quadrature_refinement_is_not_rigorous_error_certificate=True,
                         analytical_bounds=bounds,
                         CAR_calibration=finite_car_calibration(r, final, bounds)))
    sources = [BASE/"archive_1009_/research_note_1018.md",
               BASE/"archive_1009_/1018/gravity_ideal_selection_results.json",
               BASE/"archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
               BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return dict(round=1019, new_calibration_groups=1, cumulative_test_groups=3797,
                new_cognitive_axioms=0, all_scientific_calibrations_passed=True,
                scope="positive-Fock CAR quantization of a free complex massless scalar; 3+1 Minkowski; bounded even globally charge-neutral double-smeared observables",
                microcausality_generated_from_FUCP=False,
                spin_statistics_complete_theorem_numerically_proved=False,
                QED_local_gauge_invariance_claimed=False,
                physical_device_implementation_claimed=False,
                no_signaling_failure_conditional_on_declared_local_operation_menu=True,
                failure_only_at_infinite_energy_claimed=False,
                norm_cutoff_equals_physical_minimum_scale=False,
                smearing=dict(temporal="w_tau(t)=max(1-|t|/tau,0)/tau",
                    spatial="convolution of two normalized radius-a ball indicators",
                    a=A_RADIUS, tau=TAU, spatial_support_radius=2*A_RADIUS,
                    temporal_support_halfwidth=TAU, total_integral=1,
                    on_shell_Fourier="U(ap)^2 sinc(p*tau/2)^2",
                    normalized_ball_Fourier="U(z)=3(sin(z)-z*cos(z))/z^3",
                    sinc_convention="sinc(z)=sin(z)/z",
                    distribution_test_smoothing_limit_requires_operator_norm_continuity=True),
                separated_support_examples=rows,
                cutoff_projection_algebra=cutoff_projection_algebra(),
                strict_tail_formula="|tail_k(P)| <= 324/[pi^2*(6-k)*P^(6-k)], k=0,1,2; P>=1",
                tail_bounds_rigorous_analytic=True, quadrature_errors_certified=False,
                whole_state_energy_moments_include_subspace_leakage=True,
                preparation_of_entangled_witness_is_explicit_extra_menu_assumption=True,
                analytic_obligations=[
                    "positive CAR Fock and scalar covariance do not themselves supply local observables",
                    "strictly separated nonnegative normalized smearings have 0<s<n by spacelike kernel bounds",
                    "even globally neutral P,Q have a nonzero commutator and a finite-probability witness",
                    "the same ideal operation menu violates spacelike noninterference on a finite-energy state",
                    "analytic cutoff operator error is below half a strictly positive continuum lower bound",
                    "bounded CAR continuity extends the nonsmooth compact smearing from smooth local approximants"],
                historical_source_sha256={str(p.relative_to(BASE)).replace("\\", "/"):
                    hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1019, passed=result["all_scientific_calibrations_passed"],
        examples=[dict(L=row["separation"], r=row["CAR_calibration"]["overlap"],
            signal=row["CAR_calibration"]["fixed_analytic_theta_witness"]["absolute_difference"],
            certified_continuum_lower_bound=row["analytical_bounds"]["exact_continuum_probability_lower_bound"],
            cutoff_error=row["analytical_bounds"]["cutoff_control"]["same_state_delta_probability_error_upper_bound"])
            for row in result["separated_support_examples"]]), indent=2))
