"""Fixed-model scale and coarse CAR readout audit; stdout JSON only.

No author imports, temperature fit, larger boost cutoff, or change of Hamiltonian.
All decimal results are NumPy diagnostics, not interval certificates.
"""
import json
import math
import numpy as np

N_SITE = 65
A_SITE = 32
Z = 0.989835554078098
ELL = 16.0 / N_SITE
MASS = 2.0
EPSILON = 0.15
NU = 0.35
WIDTHS = (4, 8, 16)  # Declared powers of two; not selected to optimize error.


def fermi(t):
    t = np.asarray(t)
    result = np.empty_like(t)
    pos = t >= 0.0
    small = np.exp(-t[pos])
    result[pos] = small / (1.0 + small)
    small = np.exp(t[~pos])
    result[~pos] = 1.0 / (1.0 + small)
    return result


def original_matrices(mu):
    """The original transverse-zero 2-component block, not a new material."""
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
    sz = np.diag([1.0, -1.0]).astype(complex)
    k = 2.0 * np.pi * np.fft.fftfreq(N_SITE)
    energies = np.sqrt(mu * mu + np.sin(k) ** 2)
    symbols = 0.5 * (np.eye(2)[None, :, :] -
        (mu * sz[None, :, :] + np.sin(k)[:, None, None] * sx[None, :, :])
        / energies[:, None, None])
    covariance = np.empty((2*A_SITE, 2*A_SITE), complex)
    for i in range(A_SITE):
        for j in range(A_SITE):
            covariance[2*i:2*i+2, 2*j:2*j+2] = np.einsum(
                "k,kab->ab", np.exp(1j*k*(i-j)), symbols) / N_SITE
    b = np.zeros_like(covariance)
    for i in range(A_SITE):
        b[2*i:2*i+2, 2*i:2*i+2] = mu * (i + 0.5) * sz
        if i + 1 < A_SITE:
            b[2*i:2*i+2, 2*i+2:2*i+4] = -0.5j * (i+1) * sx
            b[2*i+2:2*i+4, 2*i:2*i+2] = 0.5j * (i+1) * sx
    ev, eu = np.linalg.eigh(b)
    thermal = (eu * fermi(2.0*np.pi*ev)) @ eu.conj().T
    return k, covariance, b, ev, eu, thermal


def fourier_weights(regional_modes):
    field = np.zeros((N_SITE, regional_modes.shape[1]), complex)
    field[:A_SITE] = regional_modes
    transform = np.fft.fft(field, axis=0) / np.sqrt(N_SITE)
    return np.sum(abs(transform) ** 2, axis=1)


def mode_stats(k, weights):
    return {
        "total_weight": float(weights.sum()),
        "central_cone_abs_k_below_pi_over_2": float(
            weights[abs(k) < np.pi/2].sum()),
        "zone_edge_cones_abs_k_above_pi_over_2": float(
            weights[abs(k) > np.pi/2].sum()),
        "central_abs_k_at_most_0p3": float(weights[abs(k) <= 0.3].sum()),
        "dimensionless_rms_k": float(np.sqrt(np.dot(weights, k*k))),
    }


def main():
    lapse = 1.0 + EPSILON*Z
    velocity = 1.0 + NU*Z
    mu = MASS*lapse*ELL/velocity
    strip = math.asinh(mu)
    xi = ELL/strip
    k, c, b, ev, eu, thermal = original_matrices(mu)
    selected = int(np.flatnonzero(ev > 0.0)[0])
    witness = eu[:, selected]
    witness_k = fourier_weights(witness.reshape(A_SITE, 2))

    rows = []
    for width in WIDTHS:
        profile = np.zeros(A_SITE)
        profile[:width] = math.sqrt(2.0/(width+1))*np.sin(
            np.pi*np.arange(1, width+1)/(width+1))
        isometry = np.kron(profile[:, None], np.eye(2))
        cw = isometry.conj().T @ c @ isometry
        bw = isometry.conj().T @ thermal @ isometry
        weights = fourier_weights(profile[:, None])
        p = float(cw[0, 0].real)
        q = float(bw[0, 0].real)
        # For these real envelopes both compressed matrices are diag(p,1-p).
        rho = np.array([p*(1-p), p*p, (1-p)*(1-p), p*(1-p)])
        sigma = np.array([q*(1-q), q*q, (1-q)*(1-q), q*(1-q)])
        trace_distance = float(np.abs(rho-sigma).sum()/2.0)
        trace_formula = abs(p-q)*(1.0+abs(1.0-p-q))
        sin2_exact = width/(width+1.0)*math.sin(math.pi/(width+1))**2
        max_energy = math.sqrt(mu*mu+1.0)
        signal_lower = sin2_exact/(2.0*max_energy*(max_energy+mu))
        # Static family derivative only, not time response or backreaction.
        mu_prime = MASS*ELL*(EPSILON-NU)/(velocity*velocity)
        p_prime = -0.5*mu_prime*np.sum(
            weights*np.sin(k)**2/(mu*mu+np.sin(k)**2)**1.5)
        rows.append({
            "support_sites": width,
            "support_coordinate_length": width*ELL,
            "support_over_xi_band": width*strip,
            "mean_coordinate_x": float(ELL*np.dot(profile*profile,
                                                     np.arange(A_SITE)+0.5)),
            "source_upper_occupation": p,
            "bare_boost_upper_occupation": q,
            "two_by_two_covariance_difference_opnorm": float(np.linalg.norm(cw-bw, 2)),
            "two_CAR_mode_reference_state_trace_distance": trace_distance,
            "source_occupation_variance": p*(1-p),
            "analytic_positive_occupation_lower_bound": signal_lower,
            "static_preparation_family_dp_dz_not_feedback": float(p_prime),
            "momentum": mode_stats(k, weights),
            "sanity": {
                "isometry_residual": float(np.linalg.norm(
                    isometry.conj().T@isometry-np.eye(2), 2)),
                "covariance_diagonal_form_residual": float(np.linalg.norm(
                    cw-np.diag([p, 1-p]), 2)),
                "boost_diagonal_form_residual": float(np.linalg.norm(
                    bw-np.diag([q, 1-q]), 2)),
                "trace_distance_formula_residual": abs(trace_distance-trace_formula),
                "sin_squared_moment_formula_residual": abs(
                    float(np.dot(weights, np.sin(k)**2))-sin2_exact),
                "positive_signal_bound_pass": bool(p >= signal_lower-1e-13),
            },
        })

    dispersion = []
    for mode in (1, 2, 3, 4):
        momentum = 2.0*np.pi*mode/16.0
        ka = momentum*ELL
        e_cont = math.sqrt((MASS*lapse)**2+(velocity*momentum)**2)
        e_latt = math.sqrt((MASS*lapse)**2+
                           (velocity*math.sin(ka)/ELL)**2)
        squared_bound = velocity**2*ELL**2*momentum**4/3.0
        dispersion.append({
            "central_momentum_index": mode, "coordinate_p": momentum, "ka": ka,
            "continuum_energy_comparator": e_cont, "lattice_energy": e_latt,
            "relative_energy_difference": (e_cont-e_latt)/e_cont,
            "kinetic_energy_squared_relative_difference": 1-(math.sin(ka)/ka)**2,
            "relative_group_velocity_difference": 1-(math.sin(2*ka)/(2*ka))*(e_cont/e_latt),
            "full_spinor_symbol_opnorm_difference": velocity*abs(momentum-math.sin(ka)/ELL),
            "symbol_error_upper_for_one_particle_Duhamel_per_time":
                velocity*ELL**2*abs(momentum)**3/6.0,
            "analytic_energy_squared_error_bound": squared_bound,
            "bound_pass": bool(0 <= e_cont**2-e_latt**2 <= squared_bound+1e-12),
        })
    edge_kinetic = math.sin(np.pi/N_SITE)/ELL
    edge_levels = [{
        "near_pi_coordinates": d,
        "lowest_energy_at_that_zone_corner": math.sqrt(
            (MASS*lapse)**2+d*(velocity*edge_kinetic)**2),
    } for d in range(4)]
    passed = all(row["sanity"]["positive_signal_bound_pass"] and
                 max(v for name, v in row["sanity"].items()
                     if name.endswith("residual")) < 1e-12 for row in rows)
    passed &= all(row["bound_pass"] for row in dispersion)
    output = {
        "status": "PASS_WITH_SCOPE" if passed else "CHECK_FAILED",
        "params": {"period": N_SITE, "fixed_region": A_SITE,
                   "z_reference": Z, "ell": ELL, "mass": MASS,
                   "lapse_N": lapse, "velocity_v": velocity, "mu": mu,
                   "boost_beta_fixed": 2.0*np.pi},
        "scales": {
            "complex_band_strip_width": strip,
            "auxiliary_longest_axis_correlation_scale_in_sites": 1.0/strip,
            "auxiliary_longest_axis_correlation_scale_coordinate": xi,
            "continuum_Compton_coordinate_comparator": ELL/mu,
            "box_over_xi_band": 16.0/xi,
            "proper_to_coordinate_space_factor": lapse/velocity,
            "proper_lattice_spacing": ELL*lapse/velocity,
            "proper_xi_band": xi*lapse/velocity,
            "balanced_minimum_error_for_a_over_L_and_L_over_xi": math.sqrt(strip),
        },
        "original_1196_witness_momentum": {
            "positive_B_eigenvalue": float(ev[selected]),
            "source_occupation": float(np.vdot(witness, c@witness).real),
            "boost_occupation": float(fermi(np.array([2*np.pi*ev[selected]]))[0]),
            **mode_stats(k, witness_k),
        },
        "fixed_macro_readout_families": rows,
        "central_cone_dispersion": dispersion,
        "finite_odd_box_zone_corner_levels": edge_levels,
        "scope": [
            "No larger boost matrix and no temperature adjustment.",
            "xi is the auxiliary infinite-lattice band scale, not a full-parent interacting correlation theorem.",
            "No parameter-separated a << L << xi window; this is not necessary for a massive wedge.",
            "Each macro entry restricts the fixed known references to two smeared CAR modes.",
            "It does not match the whole regional state, dynamics, instrument, entropy flux or geometric feedback.",
            "Transverse momentum is zero as in 1196; the mode is not a local three-dimensional ball.",
            "Static dp/dz compares preparations; it is not a retarded response.",
            "All entries are floating-point diagnostics, not new rational interval certificates.",
        ],
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
