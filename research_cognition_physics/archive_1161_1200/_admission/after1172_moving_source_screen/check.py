"""Unnumbered finite moving-source screen; default is read-only.

No field Fock cutoff: second-order vacuum contractions use the exact finite
oscillator kernel. This does not certify a continuum or an actual CO instrument.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = {
    "1129": ("archive_1120_1158/1129/proof.md",
             "c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa"),
    "1144": ("archive_1120_1158/1144/proof.md",
             "205dc68aec100de83b59f0f156a1afc8048dcdd87d45280dd2df242c89cde939"),
    "1146": ("archive_1120_1158/1146/proof.md",
             "0350230f9649aba387ff89aeeff5201300e1d999418ee3aa0639b79936ab8583"),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def models():
    n, mu, t = 9, 0.8, 2.0
    x = np.arange(n) - n // 2
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1.0, -1.0])
    masses = np.array([2.0, 2.4])
    m = np.diag(masses)
    p = np.zeros((n, n), complex)
    kg = (mu * mu + 2) * np.eye(n)
    for j in range(n):
        p[j, (j + 1) % n] = -0.5j
        p[j, (j - 1) % n] = 0.5j
        kg[j, (j + 1) % n] = -1
        kg[j, (j - 1) % n] = -1
    h = np.kron(p, np.kron(sx, np.eye(2))) + np.kron(
        np.eye(n), np.kron(sz, m))
    ev, evec = np.linalg.eigh(h)
    kv, ko = np.linalg.eigh(kg)
    omega = np.sqrt(kv)
    current_block = np.kron(sz, m)

    def packet(center, momentum):
        # Project each flavor onto its full positive band, then normalize it.
        # The subsequent dynamics and currents retain both bands.
        envelope = np.exp(-(x - center) ** 2 / (4 * 0.8 ** 2)
                          + 1j * momentum * x)
        columns = []
        for a in range(2):
            seed = np.zeros((n, 2, 2), complex)
            seed[:, 0, a] = envelope
            seed = seed.ravel()
            v = evec[:, ev > 0] @ (evec[:, ev > 0].conj().T @ seed)
            v /= np.linalg.norm(v)
            columns.append(v)
        return sum(columns) / math.sqrt(2), np.column_stack(columns)

    source0, js = packet(-2, 0)
    # Exact translation, NOT an independently cut Gaussian centered at +3.
    tr = np.roll(np.eye(n), 5, axis=0)
    tr = np.kron(tr, np.eye(4))
    source1 = tr @ source0
    probe, jp = packet(-1, 0.7)
    c = np.kron(np.eye(2), (np.eye(2) + sy) / 2)
    region = ((x >= 0) & (x <= 2)).astype(float)
    effects = {"clock": np.kron(np.eye(n), c),
               "endpoint": np.kron(np.diag(region), np.eye(4)),
               "joint": np.kron(np.diag(region), c)}
    return locals()


def evolve(v, times, ev, evec):
    return (np.exp(-1j * np.outer(times, ev))
            * (evec.conj().T @ v)[None, :]) @ evec.T


def coefficients(m, count):
    ts = np.linspace(0, m["t"], count)
    dt = ts[1] - ts[0]
    ev, evec = m["ev"], m["evec"]
    s0 = evolve(m["source0"], ts, ev, evec).reshape(count, 9, 4)
    s1 = evolve(m["source1"], ts, ev, evec).reshape(count, 9, 4)
    j0 = np.einsum("txi,ij,txj->tx", s0.conj(),
                   m["current_block"], s0).real
    j1 = np.einsum("txi,ij,txj->tx", s1.conj(),
                   m["current_block"], s1).real
    dj = j0 - j1
    modal = dj @ m["ko"]
    r = np.sin(ts[:, None] * m["omega"]) / m["omega"]
    length = 1 << (2 * count - 1).bit_length()
    convolution = np.fft.irfft(np.fft.rfft(r, n=length, axis=0)
                              * np.fft.rfft(modal, n=length, axis=0),
                              n=length, axis=0)[:count]
    response = (dt * (convolution - 0.5 * r * modal[0])) @ m["ko"].T
    pt = evolve(m["probe"], ts, ev, evec).reshape(count, 9, 4)
    result = {}
    for name, effect in m["effects"].items():
        final = pt[-1].ravel()
        back = evec @ (np.exp(1j * m["t"] * ev)
                       * (evec.conj().T @ (effect @ final)))
        ct = evolve(back, ts, ev, evec).reshape(count, 9, 4)
        comm = 2 * np.einsum("txi,ij,txj->tx", pt.conj(),
                            m["current_block"], ct).imag
        integrand = np.sum(comm * response, axis=1)
        coefficient = dt * (np.sum(integrand) -
                            0.5 * (integrand[0] + integrand[-1]))
        # Source-only terms vanish; probe self/noise terms are identical in
        # both source preparations. We do not remove them from the full U.
        result[name] = {
            "free_probability": float(np.vdot(final, effect @ final).real),
            "source_cross_coefficient": float(coefficient)}
    return result


def quadrature_bound(count):
    t, mass, mu = 2.0, 2.4, 0.8
    step = t / (count - 1)
    a = 2 * math.sqrt(mass * mass + 1)
    omega_max = math.sqrt(mu * mu + 4)
    cj = (2 * mass) ** 2
    b = cj * (omega_max + 2 * a + t * a * a)
    return step ** 2 * (t * t * b / 24 + t * (cj + t * b) / 12)


def compute():
    fingerprints = {}
    for key, (relative, expected) in SOURCES.items():
        actual = sha(ROOT / relative)
        if actual != expected:
            raise AssertionError(("source changed", key, actual))
        fingerprints[key] = {"path": relative, "sha256": actual}
    m = models()
    h, tr = m["h"], m["tr"]
    translated_energy = []
    for order in (1, 2, 3, 4):
        hp = np.linalg.matrix_power(h, order)
        translated_energy.append(float(abs(
            np.vdot(m["source0"], hp @ m["source0"])
            - np.vdot(m["source1"], hp @ m["source1"]))))
    base = {
        "translation_commutator_norm": float(np.linalg.norm(h @ tr - tr @ h, 2)),
        "source_energy_moment_differences_1_to_4": translated_energy,
        "source_norms": [float(np.linalg.norm(m[k]))
                         for k in ("source0", "source1")],
        "probe_norm": float(np.linalg.norm(m["probe"])),
        "source_code_isometry_error": float(np.linalg.norm(
            m["js"].conj().T @ m["js"] - np.eye(2), 2)),
        "probe_code_isometry_error": float(np.linalg.norm(
            m["jp"].conj().T @ m["jp"] - np.eye(2), 2)),
        "source_energy": float(np.vdot(m["source0"], h @ m["source0"]).real),
        "probe_energy": float(np.vdot(m["probe"], h @ m["probe"]).real),
    }
    coarse = coefficients(m, 2049)
    fine = coefficients(m, 4097)
    bound = quadrature_bound(4097)
    # The analytic truncation bound is rigorous real-arithmetic mathematics.
    # This reserve is explicitly an unvalidated floating point allowance.
    floating_reserve = 1e-10
    for name in fine:
        coefficient = fine[name]["source_cross_coefficient"]
        fine[name].update({
            "coarse_coefficient": coarse[name]["source_cross_coefficient"],
            "refinement_difference": abs(
                coefficient - coarse[name]["source_cross_coefficient"]),
            "analytic_trapezoid_error_bound": bound,
            "floating_allowance_not_interval_certified": floating_reserve,
            "screen_lower_coefficient": coefficient - bound - floating_reserve,
            "screen_upper_coefficient": coefficient + bound + floating_reserve})
    omega_inverse_diagonal = np.sum(m["ko"] ** 2 / m["omega"], axis=1)
    fmax = math.sqrt(float(max(omega_inverse_diagonal)) / 2)
    c = 2 * 2.4 * fmax
    coupling = 1e-6
    z = 2 * coupling * c * m["t"]
    eps = z ** 3 / (math.sqrt(6) * (1 - z / 2))
    state_to_polynomial_contrast_bound = 2 * (2 + eps) * eps
    # Vacuum parity cancels the odd terms of <D_2^* E D_2>. Its degree-four
    # term remains: each preparation contributes at most z^4/2.
    polynomial_quartic_contrast_bound = z ** 4
    probability_contrast_bound = (state_to_polynomial_contrast_bound
                                  + polynomial_quartic_contrast_bound)
    dyson = {
        "domain": "two first-quantized one-particle sectors and exact KG vacuum",
        "fmax": fmax, "current_number_bound_C": c,
        "coupling_example": coupling, "z": z,
        "prepared_second_order_norm_remainder_bound": eps,
        "state_to_Dyson_polynomial_contrast_bound":
            state_to_polynomial_contrast_bound,
        "Dyson_polynomial_quartic_contrast_bound":
            polynomial_quartic_contrast_bound,
        "two_preparation_probability_contrast_remainder_bound":
            probability_contrast_bound,
        "joint_screen_signal_lower": coupling ** 2 *
            fine["joint"]["screen_lower_coefficient"],
        "joint_screen_net_lower": coupling ** 2 *
            fine["joint"]["screen_lower_coefficient"] -
            probability_contrast_bound,
        "instrument_error_budget_certified": False,
        "CAR_normal_ordering_constant_included_in_this_C": False,
    }
    # Finite CAR normal-ordering data, no exponentially large Fock matrices.
    negative = m["evec"][:, m["ev"] < 0]
    sea = negative @ negative.conj().T
    kappa = np.array([np.trace(sea[4*j:4*j+4, 4*j:4*j+4]
                               @ m["current_block"]).real for j in range(9)])
    car = {
        "negative_band_rank_per_species": int(negative.shape[1]),
        "normal_order_energy_shift_per_species":
            float(-np.sum(m["ev"][m["ev"] < 0])),
        "local_scalar_normal_order_constants_per_species": kappa.tolist(),
        "constant_translation_variation": float(np.ptp(kappa)),
        "normal_Fock_completion_implemented": False,
        "comment": "Constants cancel in source contrast but enlarge a naive Dyson C.",
    }
    assert base["translation_commutator_norm"] < 1e-12
    assert max(translated_energy) < 1e-10
    assert all(fine[k]["screen_lower_coefficient"] > 0 for k in fine)
    assert dyson["joint_screen_net_lower"] > 0
    return {
        "status": "ADMISSION_SCREEN_ONLY_NOT_A_NUMBERED_ROUND",
        "source_fingerprints": fingerprints,
        "parameters": {"sites": 9, "spacing": 1, "periodic": True,
                       "mass_values": [2.0, 2.4], "KG_mass": 0.8,
                       "time": 2.0, "packet_width": 0.8,
                       "source_translation_sites": 5, "probe_p0": 0.7,
                       "endpoint_sites": [0, 1, 2],
                       "fine_grid": 4097, "coarse_grid": 2049},
        "same_energy_preparation": base,
        "second_order_effects": fine, "prepared_Dyson": dyson,
        "normal_CAR_responsibility": car,
        "not_certified": ["actual_CO_instrument", "positive_full_CAR_joint_model",
                          "continuum_limit", "strict_null_causality",
                          "common_dynamic_metric", "Einstein_equation"],
    }


def compare(new, old, path=""):
    if isinstance(new, dict):
        if set(new) != set(old):
            raise AssertionError(("keys", path))
        for key in new:
            compare(new[key], old[key], f"{path}/{key}")
    elif isinstance(new, list):
        if len(new) != len(old):
            raise AssertionError(("length", path))
        for i, (a, b) in enumerate(zip(new, old)):
            compare(a, b, f"{path}/{i}")
    elif isinstance(new, float):
        # The tiny net signal is recomputed from the much larger coefficient;
        # use its own relative scale instead of a blanket absolute tolerance.
        if not math.isclose(new, old, rel_tol=2e-9,
                            abs_tol=1e-24 if abs(new) < 1e-9 else 1e-11):
            raise AssertionError(("value", path, new, old))
    elif new != old:
        raise AssertionError(("value", path, new, old))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true",
                        help="Write results once; refuse existing file.")
    args = parser.parse_args()
    result = compute()
    result["check_sha256"] = sha(Path(__file__).resolve())
    target = HERE / "results.json"
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps({
        "status": "PASS_SCREEN",
        "write_once": args.write,
        "joint_coefficient": result["second_order_effects"]["joint"]
            ["source_cross_coefficient"],
        "coefficient_analytic_error": result["second_order_effects"]["joint"]
            ["analytic_trapezoid_error_bound"],
        "joint_net_example": result["prepared_Dyson"]["joint_screen_net_lower"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
