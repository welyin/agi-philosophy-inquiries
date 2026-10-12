"""1146 fixed packet checks; finite diagnostics do not simulate the continuous field.

Default is read-only. --write explicitly saves packet_results.json.
run() returns a JSON-serializable dict with exactly three check groups.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
import math
from pathlib import Path

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
P1 = np.diag([0.0, 1.0]).astype(complex)
I4 = np.eye(4, dtype=complex)
I8 = np.eye(8, dtype=complex)
BETA = np.kron(Z, I2)
ALPHA = [np.kron(X, pauli) for pauli in (X, Y, Z)]


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    assert norm(h - h.conj().T) < 1e-10
    e, w = np.linalg.eigh(h)
    return (w * np.exp(-1j * t * e)) @ w.conj().T


def fq(x):
    return {"numerator": x.numerator, "denominator": x.denominator, "float": float(x)}


def bump(s):
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    mask = np.abs(s) < 1
    out[mask] = np.exp(-1 / (1 - s[mask] ** 2))
    return out


def source_b(s):
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    mask = (s > 0) & (s < 1)
    out[mask] = np.exp(-1 / (s[mask] * (1 - s[mask])))
    return out


def group_1():
    # Analytic constants use rational inequalities, not quadrature observations.
    e_lower = Q(8, 3)
    qsecond = 2 * 2**2 / e_lower**2 + 8 * 3**3 / e_lower**3 + 4 * 4**4 / e_lower**4
    assert qsecond == Q(2097, 64) < 33
    assert 33 * Q(3, 2) < 50  # sqrt(2) < 3/2
    assert 3 * 50 * 4 == 600
    # e < 11/4 implies exp(8/3) < 16, and hence ||q||^2 > 1/16.
    assert Q(11, 4) ** 8 < 16**3

    lam = Q(11, 7)  # pi/2 < 11/7
    d = Q(1, 1000)
    m = 10**10
    epsilon = Q(1, 10**16)
    tc = Q(1, 4)
    k1, k2 = Q(20), Q(1000)
    B = d + 16 * lam + Q(3, 2) * (d*d + 2*d*lam + 16*lam*lam)
    packet = (B + (1+d+lam)*k1 + k2/2) / m
    kernel = 2*lam*epsilon*(Q(3*m, 4)+k1)
    early = (d+Q(3, 2)*tc*d*d+(1+tc*d)*k1+tc*k2/2)/m
    assert packet < Q(6358625, 10**14)
    assert kernel < Q(2357142864, 10**15)
    assert packet+kernel < Q(2420729109, 10**15)
    assert packet+kernel < Q(3, 10**6)
    assert early < Q(14501, 10**12)
    neg_amplitude = (k1+d)/(2*m)
    assert 1 - 2*neg_amplitude**2 > Q(99, 100)

    # Independent envelope diagnostics; no numerical value is used in certificates.
    points, weights = np.polynomial.legendre.leggauss(256)
    yy = 1-points*points
    qv = np.exp(-1/yy)
    qp = -2*points*qv/(yy*yy)
    qpp = (-2/yy**2-8*points**2/yy**3+4*points**2/yy**4)*qv
    nq = float(np.dot(weights, qv*qv))
    nq1 = float(np.dot(weights, qp*qp))
    nq2 = float(np.dot(weights, qpp*qpp))
    cross = float(np.dot(weights, qv*qpp))
    k1_num = math.sqrt(3*nq1/nq)
    k2_num = math.sqrt(3*nq2/nq+6*(cross/nq)**2)
    assert abs(cross+nq1) < 1e-9
    assert k1_num < 20 and k2_num < 600
    xx = (points+1)/2
    ww = weights/2
    ib = float(np.dot(ww, source_b(xx)))
    # F(1/2) is its exact single maximum.
    fmax = 2*math.exp(-4)/ib
    f_sq_integral = 2*float(np.dot(ww, source_b(xx)**2))/(ib*ib)
    assert fmax < 16 and f_sq_integral < 16 and 2*fmax < 32
    return {
        "pass": True,
        "kind": "exact_rational_certificate_plus_independent_envelope_diagnostic",
        "fixed_inputs": {
            "m": m, "p0_over_m": "3/4", "delta": "1/1000",
            "T": 1, "epsilon": "1/10^16", "lambda": "pi/2",
            "kappa1_upper": 20, "kappa2_used_upper": 1000,
            "kappa2_analytic_stronger_upper": 600,
        },
        "qsecond_sup_upper": fq(qsecond),
        "packet_upper": fq(packet),
        "kernel_upper": fq(kernel),
        "packet_plus_kernel_upper": fq(packet+kernel),
        "early_packet_upper": fq(early),
        "initial_negative_frequency_amplitude_upper": fq(neg_amplitude),
        "two_mode_empty_probability_lower": fq(1-2*neg_amplitude**2),
        "diagnostic_only": {
            "quadrature_points": 256, "q_norm_squared": nq,
            "kappa1": k1_num, "kappa2": k2_num,
            "integration_by_parts_residual": abs(cross+nq1),
            "F_maximum": fmax, "F_squared_integral": f_sq_integral,
            "F_total_variation": 2*fmax,
        },
    }


def band_data(momentum, mass):
    h = mass*BETA
    for ai, pi in zip(ALPHA, momentum):
        h = h + pi*ai
    e = math.sqrt(float(np.dot(momentum, momentum))+mass*mass)
    plus = (I4+h/e)/2
    minus = (I4-h/e)/2
    return h, e, plus, minus


def rotating_test(momentum, source_type):
    # A genuinely noncommuting, smoothly rotating source with exactly soluble
    # 8x8 propagator. Moderate diagnostic mass avoids catastrophic rest-energy
    # subtraction. Fixed production m=10^10 is certified only analytically.
    mass, delta, lam, omega, terminal = 20.0, 0.13, 0.71, 1.1, 0.8
    h, e, pp, pm = band_data(momentum, mass)
    rr = mass/e
    qq = h/e
    ph, mh = np.kron(pp, I2), np.kron(pm, I2)
    h0 = np.kron(h, I2)
    rotgen = np.kron(I4, omega*Z/2)
    if source_type == "scalar":
        wconst = np.kron(BETA, delta*P1+lam*X)
        dconst = h0+np.kron(rr*qq, delta*P1+lam*X)
    else:
        wconst = np.kron(BETA, delta*P1)+np.kron(I4, lam*X)
        dconst = h0+np.kron(rr*qq, delta*P1)+np.kron(I4, lam*X)
    hconst = h0+wconst

    def propagator(t, diagonal=False):
        rot = np.kron(I4, exp_h(omega*Z/2, t))
        return rot@exp_h((dconst if diagonal else hconst)-rotgen, t)

    t = 0.31
    c, s = math.cos(omega*t), math.sin(omega*t)
    vt = lam*(c*X+s*Y)
    vdot = lam*omega*(-s*X+c*Y)
    if source_type == "scalar":
        ws = delta*P1+vt
        vsdot = vdot
        w = np.kron(BETA, ws)
        d1 = np.kron(rr*qq, ws)
    else:
        ws = delta*P1
        vsdot = np.zeros((2, 2), complex)
        w = np.kron(BETA, ws)+np.kron(I4, vt)
        d1 = np.kron(rr*qq, ws)+np.kron(I4, vt)
    oo = w-d1
    ss = (mh@oo@ph-ph@oo@mh)/(2*e)
    odot = np.kron(BETA-rr*qq, vsdot)
    sdot = (mh@odot@ph-ph@odot@mh)/(2*e)
    comm_res = norm(h0@ss-ss@h0+oo)
    assert comm_res < 2e-14
    assert norm(ss+ss.conj().T) < 1e-15
    assert norm(ph@oo@ph)+norm(mh@oo@mh) < 1e-14

    # Independently differentiate the exact endpoint expression in (11).
    def exact_ss(tt):
        vv = lam*(math.cos(omega*tt)*X+math.sin(omega*tt)*Y)
        wws = delta*P1+vv if source_type == "scalar" else delta*P1
        ooo = np.kron(BETA-rr*qq, wws)
        return (mh@ooo@ph-ph@ooo@mh)/(2*e)

    def endpoint(tt):
        backwards = propagator(terminal)@propagator(tt).conj().T
        return backwards@exact_ss(tt)@propagator(tt, True)

    dt = 1e-6
    derivative = (endpoint(t+dt)-endpoint(t-dt))/(2*dt)
    rhs = (propagator(terminal)@propagator(t).conj().T @
           (sdot-1j*oo+1j*(d1@ss-ss@d1)+1j*oo@ss) @
           propagator(t, True))
    derivative_residual = norm(derivative-rhs)
    assert derivative_residual < 2e-8
    actual_gap = norm(propagator(terminal)-propagator(terminal, True))
    mu = delta+lam
    # A common safe bound (including an unnecessary derivative term in vector case).
    gap_upper = (2*mu+lam*omega*terminal+3*mu*mu*terminal)/(2*e)
    assert actual_gap <= gap_upper+1e-12
    noncomm = norm((delta*P1+lam*X)@(delta*P1+lam*Y) -
                   (delta*P1+lam*Y)@(delta*P1+lam*X))
    assert noncomm > 0.5
    return {
        "momentum": list(momentum), "source_type": source_type,
        "full_eight_dimensional_gap": actual_gap,
        "analytic_bound_for_this_diagnostic": gap_upper,
        "commutator_identity_residual": comm_res,
        "integration_by_parts_derivative_residual": derivative_residual,
        "noncommuting_internal_source_norm": noncomm,
    }


def group_2():
    mass = 20.0
    p0 = np.array([0., 0., 15.])
    _, e0, pp0, _ = band_data(p0, mass)
    # Positive center spinor, retained as a constant spinor away from center.
    chi = np.array([1., 0.], dtype=complex)
    u0 = np.concatenate([
        math.sqrt((e0+mass)/(2*e0))*chi,
        (sum((p0[i]*s for i, s in enumerate((X, Y, Z))), np.zeros((2,2), complex))
         @chi)/math.sqrt(2*e0*(e0+mass))
    ])
    assert abs(np.vdot(u0, u0)-1) < 1e-14
    assert norm((pp0-np.eye(4))@u0[:, None]) < 1e-14
    samples = [
        p0,
        p0+np.array([0.7, -0.4, 0.3]),
        p0+np.array([-5.0, 2.0, -9.0]),
    ]
    rows = []
    for momentum in samples:
        h, e, pp, pm = band_data(momentum, mass)
        dq = float(np.linalg.norm(momentum-p0))
        projection_error = norm(pp-pp0)
        negative_amplitude = float(np.linalg.norm(pm@u0))
        remainder = abs(e-e0-(3/5)*(momentum[2]-p0[2]))
        assert projection_error <= dq/(2*mass)+1e-14
        assert negative_amplitude <= dq/(2*mass)+1e-14
        assert remainder <= dq*dq/(2*mass)+1e-13
        assert norm(pp@BETA@pp-(mass/e)*pp) < 1e-14
        assert norm(pm@BETA@pm+(mass/e)*pm) < 1e-14
        row = {
            "momentum": momentum.tolist(), "projector_error": projection_error,
            "projector_lipschitz_bound": dq/(2*mass),
            "retained_negative_band_amplitude": negative_amplitude,
            "linear_dispersion_remainder": remainder,
            "hessian_bound": dq*dq/(2*mass),
        }
        rows.append(row)
    assert rows[-1]["retained_negative_band_amplitude"] > 0.05
    rotations = [
        rotating_test(samples[1], "scalar"),
        rotating_test(samples[1], "vector"),
    ]
    return {
        "pass": True,
        "kind": "full_two_band_8x8_noncommuting_exact_time_evolution_diagnostic",
        "diagnostic_mass": mass,
        "production_mass_not_floating_point_diagonalized": True,
        "fixed_momentum_geometry_samples": rows,
        "exact_rotating_source_samples": rotations,
        "continuous_field_simulation": False,
    }


def pulse_internal(delta, source_type, steps):
    pts, weights = np.polynomial.legendre.leggauss(256)
    ib = float(np.dot(weights/2, source_b((pts+1)/2)))
    dt = 1/steps
    times = (np.arange(steps)+0.5)*dt
    fs = 2*source_b(2*times-0.5)/ib
    out = I2.copy()
    source_coeff = 4/5 if source_type == "scalar" else 1.0
    for ft in fs:
        hh = (4/5)*delta*P1+source_coeff*(math.pi/2)*ft*X
        out = exp_h(hh, dt)@out
    return out


def group_3():
    delta, tc, rr = Q(1, 1000), Q(1, 4), Q(4, 5)
    x, y = delta*tc, rr*delta*tc
    clock_lower = (x-y-x**3/6)/2
    assert clock_lower > Q(24999998697, 10**15)
    # sqrt(5) < 223607/100000 is a strict rational upper bound.
    root_upper = Q(223607, 100000)
    assert root_upper**2 > 5
    source_gap_lower = (3-root_upper)/8-2*delta
    assert source_gap_lower > Q(934912, 10**7)
    plus = np.array([1,1], complex)/math.sqrt(2)
    minus_y_effect = (I2-Y)/2
    clock_vector = exp_h(float(rr*delta)*P1, float(tc))@plus
    clock_probability = float(np.vdot(clock_vector, minus_y_effect@clock_vector).real)
    clock_formula = (1+math.sin(float(y)))/2
    assert abs(clock_probability-clock_formula) < 3e-16
    wrong_clock_probability = (1+math.sin(float(x)))/2
    assert wrong_clock_probability-clock_probability > float(clock_lower)-1e-15

    initial = np.array([1.,0.], complex)
    zero_s = exp_h(float(rr)*(math.pi/2)*X, 1)@initial
    zero_v = exp_h((math.pi/2)*X, 1)@initial
    ps0, pv0 = abs(zero_s[1])**2, abs(zero_v[1])**2
    exact_zero_gap = (3-math.sqrt(5))/8
    assert abs((pv0-ps0)-exact_zero_gap) < 2e-15

    # The actual smooth pulse is checked at two fixed time resolutions only.
    # These observations do not establish the analytic Duhamel certificate.
    data = {}
    for kind in ("scalar", "vector"):
        u512 = pulse_internal(float(delta), kind, 512)
        u1024 = pulse_internal(float(delta), kind, 1024)
        convergence = norm(u1024-u512)
        assert convergence < 2e-8
        value = float(abs((u1024@initial)[1])**2)
        baseline = ps0 if kind == "scalar" else pv0
        assert abs(value-baseline) < float(delta)
        data[kind] = {
            "smooth_pulse_probability_diagnostic": value,
            "delta_zero_exact_probability": float(baseline),
            "512_to_1024_operator_difference": convergence,
            "delta_perturbation_probability_difference": abs(value-baseline),
        }
    observed_gap = (data["vector"]["smooth_pulse_probability_diagnostic"] -
                    data["scalar"]["smooth_pulse_probability_diagnostic"])
    assert observed_gap > float(source_gap_lower)
    return {
        "pass": True,
        "kind": "same_actual_flavor_effects_early_clock_and_terminal_source_readouts",
        "clock_ideal_gap_lower": fq(clock_lower),
        "source_ideal_gap_lower": fq(source_gap_lower),
        "clock_correct_probability": clock_probability,
        "clock_wrong_dictionary_probability": wrong_clock_probability,
        "clock_effect_sign_checked_by_matrix": True,
        "mass_split_not_deleted": True,
        "same_clock_coefficient_both_sources": "4/5",
        "scalar_source_coefficient": "4/5",
        "vector_source_coefficient": "1",
        "actual_F_internal_diagnostics_only": data,
        "diagnostic_source_gap": observed_gap,
        "full_CAR_and_quantum_source_errors_must_be_subtracted_elsewhere": True,
    }


def run():
    results = [group_1(), group_2(), group_3()]
    return {
        "round": 1146,
        "check_groups": 3,
        "all_passed": all(x["pass"] for x in results),
        "scope": "packet theorem; full CAR and quantum-source bridge is in main proof",
        "groups": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    text = json.dumps(data, ensure_ascii=False, indent=2)+"\n"
    destination = Path(__file__).with_name("packet_results.json")
    if args.write:
        if destination.exists():
            old = json.loads(destination.read_text(encoding="utf-8"))
            if old != data:
                raise RuntimeError("Existing results differ; refuse silent replacement")
        else:
            destination.write_text(text, encoding="utf-8")
    print(text, end="")
