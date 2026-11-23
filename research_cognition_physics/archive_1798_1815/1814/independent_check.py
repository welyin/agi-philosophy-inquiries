"""Sharp distance between the two normalized adjoint seesaw portals at a common readout."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
VEV = 246.22
CG2 = 3/5
EV_TO_GEV = 1e-9


def lower_distance(light_eV, singlet, triplet, gain):
    return math.sqrt(2*light_eV*EV_TO_GEV/VEV**2 * abs(singlet/CG2 - gain**2*triplet))


def rotation(t):
    c, s = np.cosh(t), np.sinh(t)
    return np.array([[c, -1j*s], [1j*s, c]])


def physical_test(masses, singlet, triplet, rng):
    # Rank-two PMNS embedding; all measured angles are allowed by the algebra.
    q, _ = np.linalg.qr(rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2)))
    roots = q @ np.diag(np.sqrt(np.array(masses)*EV_TO_GEV))
    alpha, beta = np.sqrt(2*np.array([singlet/CG2, triplet]))/VEV
    for gain in (.2, 1., 4.):
        scale = math.sqrt(abs(alpha**2-(gain*beta)**2))
        t = np.arctanh(min(alpha, gain*beta)/max(alpha, gain*beta))
        swap = np.eye(2) if alpha > gain*beta else np.array([[0., 1.], [1., 0.]])
        r = swap @ rotation(t)
        np.testing.assert_allclose(r@r.T, np.eye(2), atol=1e-13)
        x = roots@r
        np.testing.assert_allclose(x@x.T, roots@roots.T, rtol=1e-13, atol=1e-24)
        p, v = x[:, 0]*alpha, x[:, 1]*beta
        actual = np.linalg.norm(p-1j*gain*v)
        np.testing.assert_allclose(actual, scale*math.sqrt(masses[0]*EV_TO_GEV), rtol=1e-12)
        for _ in range(9):
            angle, phase = rng.uniform(-2, 2), rng.uniform(-math.pi, math.pi)
            r = rotation(angle) @ rotation(1j*rng.normal())
            cols = roots@r
            assert np.linalg.norm(alpha*cols[:, 0]-gain*np.exp(1j*phase)*beta*cols[:, 1]) >= actual*(1-1e-12)


def check_transport(rng):
    # Test the Duhamel bound with noncommuting positive matrices and a common projection.
    def exp_h(matrix):
        vals, vecs = np.linalg.eigh(matrix)
        return (vecs*np.exp(vals))@vecs.conj().T
    us = ut = np.eye(5)
    exposure = 0.
    for n in (5, 5, 3, 3):
        if us.shape[0] > n:
            us, ut = us[:n], ut[:n]
        f = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
        h = .003*(f@f.conj().T)
        f = rng.normal(size=(n, 1))+1j*rng.normal(size=(n, 1))
        common = .11*(f@f.conj().T)
        us = exp_h(1.5*h-common)@us
        ut = exp_h(-2.5*h-common)@ut
        exposure += np.linalg.norm(h, 2)
    assert np.linalg.norm(us-ut, 2) <= (4/1.5)*math.expm1(1.5*exposure)


def compute():
    cfg = json.loads((HERE.parent/"1370"/"inputs.json").read_bytes())
    # This is the explicit low-triplet gauge example, not a fit to the lepton sources.
    mx, triplet = 1.85e15, 1200.
    octet = mx/2
    singlet = (2*octet+3*triplet)/5
    length = math.log(singlet/triplet)
    cap, relative, overall = math.sqrt(4*math.pi), 4., 2.
    charged_cap = .02
    h = length*charged_cap**2/(16*math.pi**2)
    max_difference = (4/1.5)*math.expm1(1.5*h)
    rng = np.random.default_rng(1814)
    all_rows = {}
    for label in ("NO", "IO"):
        entry = cfg[label]
        solar = entry["solar"]
        atmospheric = entry["atm_abs"]
        if label == "NO":
            light = math.sqrt(solar[1])
            example = [math.sqrt(solar[0]), math.sqrt(atmospheric[0])]
            e_range = entry["s13"][1:]
            worst_mass_ratio = math.sqrt(solar[2]/atmospheric[1])
        else:
            light = math.sqrt(atmospheric[1]-solar[2])
            example = [math.sqrt(atmospheric[0]-solar[0]), math.sqrt(atmospheric[0])]
            e_range = [entry["s12"][1]*(1-entry["s13"][2]),
                       entry["s12"][2]*(1-entry["s13"][1])]
            worst_mass_ratio = math.sqrt(1-solar[1]/atmospheric[2])
        delta = lower_distance(light, singlet, triplet, 1.)
        wide_delta = lower_distance(light, singlet, triplet, relative)
        required_h = math.log1p(1.5*wide_delta/(4*overall*cap))/1.5
        max_singlet = CG2*(relative**2*triplet+
            VEV**2*(overall*cap*max_difference)**2/(2*light*EV_TO_GEV))
        angle_ratio = math.sqrt(CG2*triplet/singlet*worst_mass_ratio)
        max_e_change = angle_ratio/math.sqrt(1+angle_ratio**2)
        physical_test(example, singlet, triplet, rng)
        physical_test(example, .02*triplet, triplet, rng)
        all_rows[label] = {
            "lighter_nonzero_mass_floor_eV": light,
            "minimum_Clebsch_stripped_distance_unit_gain": delta,
            "lower_distance_over_relative_gain_at_most_four": wide_delta,
            "minimum_distance_solution_triplet_e_fraction_max_shift_from_heavier_PMNS_column": max_e_change,
            "minimum_distance_solution_e_fraction_range": [e_range[0]-max_e_change, e_range[1]+max_e_change],
            "different_map_norm_needed_under_scalar_and_parent_budgets": wide_delta/(overall*cap),
            "charged_Y_e_rms_necessary_if_only_differential_driver": math.sqrt(required_h*16*math.pi**2/length),
            "singlet_mass_ceiling_under_registered_charged_exposure_GeV": max_singlet}
    for _ in range(10):
        check_transport(rng)
    return {"mass_example_GeV":{"X":mx, "octet":octet, "triplet":triplet, "singlet":singlet},
            "Duhamel_test_budgets":{"parent_vector_norm":cap, "relative_scalar_factor_max":relative,
                "singlet_scalar_and_low_common_norm_gain_max":overall,
                "interval_log_length": length, "charged_Y_e_spectral_norm_cap":charged_cap,
                "dimensionless_exposure_upper":h, "two_map_norm_difference_upper":max_difference},
            "oscillation_cases":all_rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    path = HERE/"independent_results.json"
    raw = json.dumps(compute(), indent=2, allow_nan=False)+"\n"
    if parser.parse_args().check:
        assert json.loads(raw) == json.loads(path.read_bytes())
    else:
        path.write_text(raw, encoding="utf-8")
    print("1814 independent", hashlib.sha256(raw.encode()).hexdigest()[:12])
