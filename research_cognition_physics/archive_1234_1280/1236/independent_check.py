"""Independent round-1236 checks; NumPy and exact fractions, no author imports.

Default replay is read-only. --write creates this reviewer's bound result.
The finite instrument contract is a hypothesis, not a device certification.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def adm_worldline():
    h = np.array([[1.4, .2, .1], [.2, 1.1, .15], [.1, .15, .9]])
    hi = np.linalg.inv(h)
    sqrt_h = np.sqrt(np.linalg.det(h))
    N, mass, omega, action = 1.2, 1.0, 1.0, .2
    max_identity, max_gradient = 0.0, 0.0
    minimum_source_loss = float("inf")
    rows = []
    for p in (np.array([.17, -.23, .31]), np.zeros(3)):
        for shift in (np.array([.12, -.08, .04]), np.zeros(3)):
            E = np.sqrt(mass**2 + p @ hi @ p)
            velocity = N * (hi @ p) / E - shift
            tau_dot = np.sqrt(N*N - (velocity + shift) @ h @ (velocity + shift))
            phase_dot = omega * tau_dot
            g = np.empty((4, 4))
            g[0, 0] = -N*N + shift @ h @ shift
            g[0, 1:] = g[1:, 0] = h @ shift
            g[1:, 1:] = h
            u = np.r_[1.0, velocity] / tau_dot
            u_cov = g @ u
            normal_cov = np.r_[-N, np.zeros(3)]
            # At fixed coordinate t: dτ/dt=1/u0, sqrt(-g)=N sqrt(h).
            stress = mass / (N * sqrt_h * u[0]) * np.outer(u, u)
            normal_energy = normal_cov @ stress @ normal_cov
            momentum_density = -(normal_cov @ stress @ g)[1:]
            ham = N * E - shift @ p
            lag = action * phase_dot - mass * tau_dot
            legendre = p @ velocity + action * phase_dot - lag
            discrepancies = [abs(u @ g @ u + 1), abs(tau_dot - N*mass/E),
                             abs(legendre - ham), abs(normal_energy*sqrt_h - E),
                             np.max(np.abs(mass*u_cov[1:] - p)),
                             np.max(np.abs(momentum_density*sqrt_h - p))]
            max_identity = max(max_identity, *discrepancies)

            def H(momentum, internal_action=action, metric=h):
                M = .8 + omega * internal_action
                return N*np.sqrt(M*M + momentum @ np.linalg.inv(metric) @ momentum) - shift @ momentum

            step = 1e-5
            gradp = np.array([(H(p + step*unit) - H(p - step*unit))/(2*step)
                              for unit in np.eye(3)])
            gradI = (H(p, action+step)-H(p, action-step))/(2*step)
            max_gradient = max(max_gradient, np.max(np.abs(gradp-velocity)), abs(gradI-phase_dot))
            # A symmetric metric variation, with all off-diagonal terms retained.
            dh = np.array([[.1, -.2, .07], [-.2, .05, .03], [.07, .03, -.08]])
            metric_grad = (H(p, metric=h+step*dh)-H(p, metric=h-step*dh))/(2*step)
            expected_metric_grad = -N * (hi @ p) @ dh @ (hi @ p) / (2*E)
            max_gradient = max(max_gradient, abs(metric_grad-expected_metric_grad))
            # Delete only clock energy in the source; preserve actual phase/trajectory.
            wrong_mass = mass - omega*action
            wrong_normal_energy = (wrong_mass / (N*sqrt_h*u[0])
                                   * (normal_cov @ u)**2)
            minimum_source_loss = min(minimum_source_loss, normal_energy-wrong_normal_energy)
            rows.append({"p_norm": float(np.linalg.norm(p)), "shift_norm": float(np.linalg.norm(shift)),
                         "normal_energy_density": float(normal_energy), "coordinate_H": float(ham),
                         "phase_rate": float(phase_dot), "proper_time_rate": float(tau_dot)})
    assert max_identity < 4e-15 and max_gradient < 2e-10
    assert minimum_source_loss > 0
    # Positive normal energy does not imply positivity for arbitrary coordinate shift.
    p = np.array([.17, -.23, .31])
    E = np.sqrt(mass**2+p @ hi @ p)
    large_shift = 10*p
    coordinate_H = N*E-large_shift @ p
    assert coordinate_H < 0 and E > 0
    # Physical (not doubled) source normalization in the CT canonical convention.
    kappa, Lambda_phys, R3 = .73, .17, .16
    K = .1*h+np.diag([.03,-.02,.04])
    shape = hi @ K
    kin = np.trace(shape @ shape)-np.trace(shape)**2
    B = R3-kin  # contracted Gauss identity, no matter field equation used
    a, b, Lambda_H = 2*kappa,1/(2*kappa),Lambda_phys/(2*kappa)
    rho_phys = E/sqrt_h
    C = kin/a-b*R3+2*Lambda_H+rho_phys
    geometric_C = -B/(2*kappa)+Lambda_phys/kappa+rho_phys
    assert abs(C-geometric_C) < 3e-16
    assert abs((C+rho_phys)-geometric_C) > .1
    return {"identity_max_error": float(max_identity), "derivative_max_error": float(max_gradient),
            "minimum_deleted_source_density": float(minimum_source_loss), "moving_and_shift_rows": rows,
            "large_shift_coordinate_H": float(coordinate_H),
            "physical_source_constraint_normalization_error": float(abs(C-geometric_C)),
            "physical_source_is_not_doubled": True,
            "gravitational_pi_independence": "analytic: matter H has only h,N,beta,p,I"}


def dust_return(t, H=2/3, distance=.03):
    z = np.cbrt(1+1.5*H*t)
    return t + 2*distance*z*z + 2*H*distance**2*z + 2/3*H*H*distance**3


def determinant(emissions, returns):
    return ((emissions[1]-emissions[0])*(returns[2]-returns[0])
            -(emissions[2]-emissions[0])*(returns[1]-returns[0]))


def flat_return(t, b, velocity, delay=0.0):
    # Solve (s-t)^2=|b+v*s|^2, choosing the future reflection event.
    vv = velocity @ velocity
    q = t+b @ velocity
    s = (q+np.sqrt(q*q-(1-vv)*(t*t-b @ b)))/(1-vv)
    return s+delay+np.linalg.norm(b+velocity*(s+delay))


def flat_convexity_and_finite_comparison():
    times = np.array([0., .331, .728])
    rows = []
    for b, velocity in ((np.array([.4, .2, .3]), np.array([.12, -.07, .05])),
                        (np.array([.4, 0., 0.]), np.array([.13, 0., 0.]))):
        for delay in (0.0, .07):
            ret = [flat_return(t, b, velocity, delay) for t in times]
            D = determinant(times, ret)
            assert D >= -3e-15
            # Independent differentiable norm calculation, including nonzero delay.
            s = .6
            def radial(k):
                x = b+velocity*k
                r = np.linalg.norm(x)
                rp = (x @ velocity)/r
                rpp = ((velocity @ velocity)*r*r-(x @ velocity)**2)/(r**3)
                return rp, max(0.0, rpp)
            rp, rpp = radial(s)
            rp2, rpp2 = radial(s+delay)
            second = (rpp2*(1-rp)+(1+rp2)*rpp)/(1-rp)**3
            assert second >= 0
            rows.append({"transverse": bool(np.linalg.norm(np.cross(b, velocity)) > 0),
                         "constant_coordinate_delay": delay, "D": float(D), "R_second": float(second)})

    ret = dust_return(times)
    dust_D = float(determinant(times, ret))
    assert abs(dust_D+.00044232) < 2e-16
    # Delete free-geodesic status: an accelerated, timelike flat mirror exactly mimics R.
    reflection_times = (ret+times)/2
    positions = (ret-times)/2
    mimicked = reflection_times+positions
    assert np.max(np.abs(reflection_times-positions-times)) < 3e-16
    assert np.max(np.abs(mimicked-ret)) < 3e-16
    H, distance = 2/3, .03
    y = 1+1.5*H*times
    rp = 1+2*H*distance*y**(-1/3)+H**2*distance**2*y**(-2/3)
    rpp = -H**2*distance*y**(-4/3)-H**3*distance**2*y**(-5/3)
    speeds = (rp-1)/(rp+1)
    accelerations = 4*rpp/(rp+1)**3
    assert np.max(np.abs(speeds)) < 1 and np.max(accelerations) < 0
    return {"flat_free_rows": rows, "dust_D": dust_D,
            "accelerated_flat_mirror_same_D": float(determinant(times, mimicked)),
            "accelerated_mirror_max_speed": float(np.max(np.abs(speeds))),
            "accelerated_mirror_max_acceleration": float(np.max(accelerations)),
            "delayed_reflection_scope": "constant proper delay on an inertial mirror"}


def exact_source_and_record_budget():
    e = [F(0), F(331,1000), F(728,1000)]
    r = [F(61208,1000000), F(404928,1000000), F(815848,1000000)]
    D = determinant(e, r)
    assert D == -F(5529,12500000)
    error_coefficient = 2*(e[1]-e[0]+e[2]-e[0]+r[1]-r[0]+r[2]-r[0])
    assert error_coefficient == F(431472,100000)
    epsilon = F(1,1000000)
    single_error = error_coefficient*epsilon+8*epsilon**2
    max_vertex_error = F(0)
    for signs in itertools.product((-1,1), repeat=6):
        noisy_e = [e[i]+signs[i]*epsilon for i in range(3)]
        noisy_r = [r[i]+signs[3+i]*epsilon for i in range(3)]
        delta = abs(determinant(noisy_e,noisy_r)-D)
        max_vertex_error = max(max_vertex_error, delta)
    assert max_vertex_error <= single_error
    # Two models' noisy record boxes can overlap only if ideals differ by <=2 eps.
    flat_epsilon = F(1,100000)
    flat_overlap_error = error_coefficient*(2*flat_epsilon)+8*(2*flat_epsilon)**2
    flat_exclusion_margin = abs(D)-flat_overlap_error
    assert flat_exclusion_margin > 0
    assert float(max(r)+flat_epsilon) < 2*np.pi
    # Multi-affinity in each of the six independent record errors makes vertices sufficient.
    a,b,d = e[1],e[2],F(3,100)
    H_lower, gap_lower = F(1118,1875), F(211,3000)
    assert F(8944,10000)**2 < F(4,5) < F(8945,10000)**2
    derivative_lower = 2*H_lower*d*F(5,6)**7+3*H_lower**2*d**2*F(5,6)**8
    source_lower = a*b*(b-a)/2*gap_lower*derivative_lower
    assert source_lower == F(36049397777742167,1049760000000000000000)
    source_margin = source_lower-2*single_error
    assert source_margin > F(25711,1000000000)
    full = [dust_return(float(t)) for t in e]
    dropped = [dust_return(float(t), H=2*np.sqrt(.8)/3) for t in e]
    measured_difference = determinant(list(map(float,e)), dropped)-determinant(list(map(float,e)), full)
    assert measured_difference > float(source_lower)
    # The deleted source uses its own constraint-compatible H, not the full H.
    H_full, H_dropped = 2/3,2*np.sqrt(.8)/3
    assert abs(H_full**2-(4/3)*1/3) < 1e-16
    assert abs(H_dropped**2-(4/3)*.8/3) < 1e-16
    return {"D_exact": str(D), "D": float(D),
            "all_record_error_coefficient": float(error_coefficient),
            "epsilon_per_record": float(epsilon), "single_branch_error": float(single_error),
            "exact_64_vertex_max_error": float(max_vertex_error),
            "flat_comparison_epsilon_per_record": float(flat_epsilon),
            "flat_two_box_overlap_error": float(flat_overlap_error),
            "flat_class_exclusion_margin": float(flat_exclusion_margin),
            "phase_winding_counter_required_in_this_menu": False,
            "source_gap_lower_exact": str(source_lower), "source_gap_lower": float(source_lower),
            "two_branch_certified_margin": float(source_margin),
            "source_gap_float_diagnostic": float(measured_difference),
            "source_parameters_fixed_except_H_from_constraint": True}


def run():
    return {"status": "PASS_CONDITIONAL_CLASSICAL_COMMON_IDENTITY", "diagnostic_groups": 3,
            "ADM_worldline_and_source": adm_worldline(),
            "finite_radar_and_deleted_free_endpoint": flat_convexity_and_finite_comparison(),
            "exact_source_and_six_record_budget": exact_source_and_record_budget(),
            "scope": {"actual_signal_record_device_certified": False,
                      "unknown_quantum_code_and_reference_certified": False,
                      "GR_derived_from_two_CT_ratios": False,
                      "instrument_error_numbers_are_conditional_contracts": True,
                      "fixed_energy_record_term_is_a_memory_implementation": False}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    output["author_binding"] = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                for name in ("proof.md","check.py","results.json","NEXT.md","admission.md")}
    destination = HERE/"independent_results.json"
    if args.write:
        if destination.exists():
            raise FileExistsError("independent_results.json already exists; default replay is read-only")
        destination.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        if json.loads(destination.read_text(encoding="utf-8")) != output:
            raise SystemExit("Independent replay differs from saved checks or author binding")
    print(json.dumps({"status": output["status"],"groups":3,
                      "mode":"write" if args.write else "read_only_match"}))


if __name__ == "__main__":
    main()
