"""T2K rate projections and a constructive counterexample to true-energy inversion."""
import argparse
import itertools
import json
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PRIOR = runpy.run_path(str(HERE.parent / "1668/check.py"))
CONFIG = json.loads((HERE / "inputs.json").read_bytes())
DATA = CONFIG["t2k_2023"]
PROBE = CONFIG["competition"]


def coefficients(y):
    minus, zero, plus, pi = np.asarray(y, dtype=float)
    return np.array([(minus + zero + plus + pi) / 4, (zero - pi) / 2,
                     (plus - minus) / 2, (zero + pi - plus - minus) / 4])


def harmonics(delta):
    return np.array([1., np.cos(delta), np.sin(delta), np.cos(2 * delta)])


def pdg_mix(delta):
    s12, s13, s23 = np.sqrt(PROBE["mixing_reference_sin2_theta12_13_23"])
    c12, c13, c23 = np.sqrt(1 - np.array([s12, s13, s23])**2)
    p = np.exp(1j * delta)
    return np.array([[c12*c13, s12*c13, s13/p],
                     [-s12*c23-c12*s23*s13*p, c12*c23-s12*s23*s13*p, s23*c13],
                     [s12*s23-c12*c23*s13*p, -c12*s23-s12*c23*s13*p, c23*c13]])


def kernel(mix, d21, d31, polarity):
    """Pure muon-flavour appearance per incident selected unit, plus flat background."""
    energies = np.array(PROBE["true_energy_GeV"])
    geom = PRIOR["INPUT"]["propagation"]
    length = geom["baseline_km"] * geom["km_per_inverse_eV_conversion"]
    if polarity == "antinu":
        mix = mix.conj()
    vac = (mix * np.array([0., d21, d31])) @ mix.conj().T
    matter = geom["V_per_rhoYe_eV"] * PROBE["rhoYe"] * (1 if polarity == "nu" else -1)
    h = vac[None, :, :] / (2e9 * energies[:, None, None]) + np.diag([matter, 0., 0.])
    eig, modes = np.linalg.eigh(length * h)
    amp = np.sum(modes[:, 0, :] * modes[:, 1, :].conj() * np.exp(-1j * eig), axis=1)
    return np.r_[abs(amp)**2, 1.]


def compatible_kernels(polarity):
    """Extreme positive solutions of three harmonic moment equations on a fixed grid."""
    d21 = PROBE["dm21_eV2"]
    d31 = PROBE["dm32_eV2"] + d21
    phases = np.pi * np.array(DATA["phase_rad_in_pi"])
    samples = np.array([kernel(pdg_mix(x), d21, d31, polarity) for x in phases])
    mat = coefficients(samples)
    assert np.max(abs(mat[3])) < 1e-12
    published = np.array(DATA["templates"][polarity])
    target = coefficients(published)[:3]
    weights = []
    for active in itertools.combinations(range(samples.shape[1]), 3):
        sub = mat[:3, active]
        if abs(np.linalg.det(sub)) < 1e-12:
            continue
        v = np.linalg.solve(sub, target)
        if v.min() < -1e-9:
            continue
        full = np.zeros(samples.shape[1])
        full[list(active)] = v
        assert np.max(abs(samples @ full - published)) < .005001
        weights.append(full)
    assert weights
    return np.array(weights), samples[1]


def deviance(observed, expected):
    return float(2*np.sum(expected - observed + observed*np.log(observed/expected)))


def study():
    baseline = json.loads((HERE / CONFIG["parent"]).read_bytes())["branches"]
    charged, _ = PRIOR["OLD"]["charged_basis"](PRIOR["read"](PRIOR["INPUT"]["charged"]))
    phi = -np.pi/6
    branches = {}
    mixes = {}
    observed = np.array(list(DATA["one_ring_e_observed"].values()))
    for name, sign, R, d21 in PRIOR["INPUT"]["branches"]:
        phases = ([2*phi, 3*phi + np.angle(3 + 2*np.exp(4j*phi))] if sign == 1
                  else [10*phi + np.radians(11.5), 7*phi])
        mix = PRIOR["neutrino_mix"](charged, phases, sign, R)
        delta = np.pi * PRIOR["phase_delta"](mix)
        assert abs(delta - baseline[name]["witness"]["delta_over_pi"]*np.pi) < 1e-9
        rates = np.array([coefficients(y) @ harmonics(delta)
                          for y in DATA["templates"].values()])
        scaled = rates * observed.sum()/rates.sum()
        branches[name] = {"delta_rad":delta, "reference_slice_counts_nu_antinu":rates.tolist(),
                          "fixed_poisson_deviance":deviance(observed, rates),
                          "free_common_scale_deviance":deviance(observed, scaled),
                          "relative_nu_over_antinu_scale_for_exact_ratio":
                              (observed[0]/observed[1])/(rates[0]/rates[1])}
        mixes[name] = (mix, d21, R**2*d21)
    counter = {}
    comparative_scores = []
    for channel, polarity in enumerate(DATA["templates"]):
        w, rate_zero = compatible_kernels(polarity)
        band_counts = w[:, 0] * rate_zero[0]
        low, high = int(np.argmin(band_counts)), int(np.argmax(band_counts))
        assert band_counts[low] == 0 and band_counts[high] > 0
        rates_by_branch = {}
        channel_deviances = []
        for name, pars in mixes.items():
            predicted = w @ kernel(*pars, polarity)
            rates_by_branch[name] = [float(predicted.min()), float(predicted.max())]
            n = observed[channel]
            channel_deviances.append(2*(predicted-n+n*np.log(n/predicted)))
        # Each endpoint comparison uses identical weights for both hypotheses.
        differences = channel_deviances[1] - channel_deviances[0]
        comparative_scores.append([float(min(differences)), float(max(differences))])
        counter[polarity] = {
            "extreme_kernels_count":len(w),
            "band_count_at_reference_delta_zero":[band_counts[low], band_counts[high]],
            "example_weights_low_high":[w[low].tolist(), w[high].tolist()],
            "full_branch_witness_counts_over_competitors":rates_by_branch,
            "rate_preference_examples_weights_low_high":[
                w[int(np.argmin(differences))].tolist(), w[int(np.argmax(differences))].tolist()]}
    a,b = branches.values()
    attained = np.sum(comparative_scores, axis=0)
    assert attained[0] < 0 < attained[1]
    return {"reference_slice":branches,
            "attained_twice_log_L_third_vs_second_same_kernel_per_polarity":
                attained.tolist(),
            "twice_log_L_third_vs_second_reference_slice":b["fixed_poisson_deviance"]-a["fixed_poisson_deviance"],
            "twice_log_L_third_vs_second_common_scale_profiled":b["free_common_scale_deviance"]-a["free_common_scale_deviance"],
            "nonunique_response":counter}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = PRIOR["tidy"](study())
    target = HERE/"results.json"
    if args.write:
        target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    else:
        assert json.loads(target.read_bytes()) == result
    print(json.dumps(result, indent=2))
