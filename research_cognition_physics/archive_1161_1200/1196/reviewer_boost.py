"""Independent finite boost witness audit; stdout JSON only, no author imports.

This script checks the declared 65-site, transverse-zero, two-component
effective CAR block. It does not certify a Hartree state of the complete
finite-K parent or implement a physical occupation instrument.
"""
import json
import math
import numpy as np

Z_REFERENCE = 0.989835554078098
PERIOD = 65
REGION = 32
ETA = 0.4


def fermi(values):
    values = np.asarray(values, dtype=float)
    answer = np.empty_like(values)
    positive = values >= 0.0
    small = np.exp(-values[positive])
    answer[positive] = small / (1.0 + small)
    small = np.exp(values[~positive])
    answer[~positive] = 1.0 / (1.0 + small)
    return answer


def matrices():
    return (np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
            np.diag([1.0, -1.0]).astype(complex))


def periodic_covariance(mu, n, length):
    sx, sz = matrices()
    momenta = 2.0 * np.pi * np.arange(n) / n
    blocks = []
    for k in momenta:
        energy = math.sqrt(mu * mu + math.sin(k) ** 2)
        blocks.append(0.5 * (np.eye(2) -
                      (mu * sz + math.sin(k) * sx) / energy))
    result = np.empty((2 * length, 2 * length), dtype=complex)
    for i in range(length):
        for j in range(length):
            block = np.zeros((2, 2), dtype=complex)
            for k, symbol in zip(momenta, blocks):
                block += np.exp(1j * k * (i - j)) * symbol
            result[2*i:2*i+2, 2*j:2*j+2] = block / n
    return result


def periodic_hamiltonian(mu, n):
    sx, sz = matrices()
    h = np.zeros((2*n, 2*n), dtype=complex)
    for i in range(n):
        h[2*i:2*i+2, 2*i:2*i+2] = mu * sz
        j = (i + 1) % n
        h[2*i:2*i+2, 2*j:2*j+2] += -0.5j * sx
        h[2*j:2*j+2, 2*i:2*i+2] += 0.5j * sx
    return h


def boost(mu, length):
    sx, sz = matrices()
    b = np.zeros((2*length, 2*length), dtype=complex)
    for i in range(length):
        b[2*i:2*i+2, 2*i:2*i+2] = mu * (i + 0.5) * sz
        if i + 1 < length:
            t = -0.5j * (i + 1) * sx
            b[2*i:2*i+2, 2*i+2:2*i+4] = t
            b[2*i+2:2*i+4, 2*i:2*i+2] = t.conj().T
    return b


def alias_bound(mu, n, site_norms, eta):
    # Analytic strip of P_-(k): Re(mu^2 + sin(k+i eta)^2) > 0.
    gap_squared = mu * mu - math.sinh(eta) ** 2
    if gap_squared <= 0:
        raise ValueError("eta must lie strictly below asinh(mu)")
    symbol_bound = 0.5 * (
        1.0 + (mu + math.cosh(eta)) / math.sqrt(gap_squared))
    sites = np.arange(len(site_norms))
    distances = np.abs(sites[:, None] - sites[None, :])
    image_bound = symbol_bound * (
        np.exp(-eta * (n - distances)) +
        np.exp(-eta * (n + distances))) / (-math.expm1(-eta * n))
    return float(site_norms @ image_bound @ site_norms), symbol_bound


def main():
    lapse = 1.0 + 0.15 * Z_REFERENCE
    velocity = 1.0 + 0.35 * Z_REFERENCE
    ell = 16.0 / PERIOD
    mu = 2.0 * lapse * ell / velocity
    c = periodic_covariance(mu, PERIOD, REGION)
    h = periodic_hamiltonian(mu, PERIOD)
    h_values, h_vectors = np.linalg.eigh(h)
    negative_vectors = h_vectors[:, h_values < 0.0]
    c_direct = negative_vectors @ negative_vectors.conj().T
    fourier_residual = float(np.linalg.norm(
        c - c_direct[:2*REGION, :2*REGION], ord=2))

    b = boost(mu, REGION)
    values, vectors = np.linalg.eigh(b)
    selected = int(np.flatnonzero(values > 0.0)[0])
    psi = vectors[:, selected]
    psi = psi / np.linalg.norm(psi)
    rayleigh = float(np.vdot(psi, b @ psi).real)
    inside_residual = float(np.linalg.norm(b @ psi - rayleigh * psi))
    site_norms = np.linalg.norm(psi.reshape(REGION, 2), axis=1)
    outgoing_squared = (REGION / 2.0) ** 2 * float(site_norms[-1] ** 2)
    total_residual_squared = inside_residual ** 2 + outgoing_squared
    # Mean(B_M - rayleigh)=0. Taylor remainder bound for F(x)=f(2 pi x).
    extension_bound = (np.pi ** 2 / (3.0 * np.sqrt(3.0))) * total_residual_squared
    periodic_bound, symbol_bound = alias_bound(
        mu, PERIOD, site_norms, ETA)
    p_periodic = float(np.vdot(psi, c @ psi).real)
    p_boost = float(fermi(np.array([2.0*np.pi*rayleigh]))[0])
    raw_difference = p_periodic - p_boost

    extensions = []
    for length in (48, 64, 96):
        bm = boost(mu, length)
        ev, eu = np.linalg.eigh(bm)
        extended = np.zeros(2*length, dtype=complex)
        extended[:2*REGION] = psi
        weights = np.abs(eu.conj().T @ extended) ** 2
        p = float(np.dot(weights, fermi(2.0*np.pi*ev)))
        measured_residual_squared = float(np.linalg.norm(
            bm @ extended - rayleigh * extended) ** 2)
        extensions.append({
            "auxiliary_boost_length": length,
            "occupation_same_zero_extended_mode": p,
            "difference_from_finite_boost": abs(p - p_boost),
            "residual_squared": measured_residual_squared,
            "analytic_bound_with_fp_slack_pass":
                bool(abs(p - p_boost) <= extension_bound + 5e-12),
        })

    prefixes = []
    for prefix in (8, 12, 16, 24):
        tail = float(np.sum(site_norms[prefix:] ** 2))
        short = psi.copy()
        short[2*prefix:] = 0.0
        short /= np.linalg.norm(short)
        p_short = float(np.vdot(short, c @ short).real)
        prefixes.append({
            "prefix_sites": prefix,
            "discarded_mode_weight": tail,
            "one_model_occupation_change_bound": math.sqrt(tail),
            "normalized_prefix_occupation": p_short,
            "measured_probability_change": abs(p_short - p_periodic),
        })

    passed = (fourier_residual < 1e-12
              and all(x["analytic_bound_with_fp_slack_pass"] for x in extensions)
              and raw_difference > periodic_bound + extension_bound)
    output = {
        "status": "PASS_WITH_SCOPE" if passed else "CHECK_FAILED",
        "params": {
            "z_reference_decimal": Z_REFERENCE, "lapse_N": lapse,
            "velocity_v": velocity, "ell": ell, "mass": 2.0, "mu": mu,
            "period": PERIOD, "region_sites": REGION,
            "boost_beta_test": 2.0*np.pi, "analytic_strip_eta": ETA,
        },
        "original_covariance_crosscheck": {
            "negative_rank": int(np.sum(h_values < 0.0)),
            "fourier_vs_direct_negative_projector_opnorm": fourier_residual,
            "hamiltonian_min_abs_energy_dimensionless": float(np.min(np.abs(h_values))),
        },
        "selected_mode": {
            "positive_mode_index": selected, "rayleigh_lambda": rayleigh,
            "inside_eigen_residual": inside_residual,
            "last_site_weight": float(site_norms[-1]**2),
            "covariance_occupation": p_periodic,
            "geometric_boost_fermi_occupation": p_boost,
            "signed_difference": raw_difference,
        },
        "analytic_boundary_controls": {
            "all_finite_boost_extensions_bound": float(extension_bound),
            "periodic_vs_infinite_reference_bound": periodic_bound,
            "complex_symbol_norm_bound": symbol_bound,
            "combined_boundary_bound": float(periodic_bound + extension_bound),
            "difference_minus_boundary_bounds": float(
                raw_difference - periodic_bound - extension_bound),
        },
        "auxiliary_extension_checks": extensions,
        "local_read_support": prefixes,
        "scope": [
            "Independent floating-point deployment of explicit analytic inequalities, not interval certification.",
            "Original finite 65-site negative-energy covariance, no logarithm or clipping of covariance eigenvalues.",
            "The larger boost matrices are auxiliary boundary comparisons, not an enlarged physical preparation.",
            "Fixed Hartree/reference effective block only; no complete finite-K ground-state certification.",
            "Occupation-effect comparison only; no implemented measurement instrument, post-state or full source-insertion certificate.",
            "Mismatch of this fixed lattice candidate is not a no-go for Unruh, effective Lorentz symmetry, or general relativity.",
        ],
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
