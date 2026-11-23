"""Antineutrino-inclusive kaon widths. All masses and hadronic factors come from inputs.json."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent
INPUT = json.loads((BASE / "inputs.json").read_text(encoding="utf-8"))


def phase(mass, meson):
    return mass * (1 - (meson / mass) ** 2) ** 2


def shape(a, b, c, angle):
    # In RHS-neutrino row order (s,d). Neutron overall isospin sign arbitrary.
    t = math.tan(angle)
    return np.array([a, t * b]), np.array([c, -t * b])


def interval(width):
    # Exact vertex search: ratio=(c²+x)/(a²+x); x=tan²(theta)*b².
    boxes = []
    for key in ("us_d", "ud_s", "ds_u"):
        center, *errors = INPUT["W0_GeV2"][key]
        delta = width * sum(errors)
        boxes.append((max(0, abs(center) - delta), abs(center) + delta))
    theta, err = INPUT["CKM_matching"]["angle_rad"]
    boxes.append((theta - width * err, theta + width * err))
    values = []
    for a, b, c, angle in itertools.product(*boxes):
        p, n = shape(a, b, c, angle)
        values.append(float(n @ n / (p @ p)))
    return min(values), max(values)


def gram_envelope(p, n, lower):
    pn, nn = p / np.linalg.norm(p), n / np.linalg.norm(n)
    sin_angle = abs(float(pn[0]*nn[1]-pn[1]*nn[0]))
    x = (1 - lower) * sin_angle / (2 * math.sqrt(lower))
    ampl = math.sqrt(1 + x*x) + x
    ratio = float(n @ n / (p @ p))
    lo, hi = ratio / ampl**2, ratio * ampl**2
    for target, maximize in [(lo, False), (hi, True)]:
        # Optimize tr(G H) for H=nn^T-target pp^T. Active eigenvalues attain bounds.
        h = np.outer(n, n) - target*np.outer(p, p)
        w, v = np.linalg.eigh(h)
        positive = (w > 0) if maximize else (w < 0)
        g = (v * np.where(positive, 1, lower)) @ v.T
        assert abs((n @ g @ n) / (p @ g @ p) - target) < 1e-12
    return lo, hi


def run():
    mass = INPUT["mass_GeV"]
    q = phase(mass["n"], mass["K0"]) / phase(mass["p"], mass["K+"])
    theta = INPUT["CKM_matching"]["angle_rad"][0]
    p, n = shape(*(INPUT["W0_GeV2"][k][0] for k in ("us_d", "ud_s", "ds_u")), theta)
    central = q * float(n @ n / (p @ p))
    width_boxes = {str(i): [q*j for j in interval(i)] for i in (1, 2)}
    tau_p = INPUT["proton"]["partial_lifetime_lower_yr"]
    tau_n = INPUT["neutron"]["partial_lifetime_lower_yr"]
    # Start from the two general Wilson terms, not from their reduced norms.
    rng = np.random.default_rng(1483)
    for ratio in (0., 1., 150.):
        unitary, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))
        v1 = np.exp(1j*rng.uniform(-3., 3.))
        v4 = np.diag(np.exp(1j*rng.uniform(-3., 3., 3)))
        v = np.array([math.cos(theta), math.sin(theta)*np.exp(0.7j), 0j])
        def coefficient(alpha, beta):
            return v1*(v[alpha]*unitary[beta] + ratio*v4[beta, alpha]*(v@v4.conj().T@unitary))
        ca, cb = coefficient(0, 1), coefficient(1, 0)
        a, b, c = [INPUT["W0_GeV2"][k][0] for k in ("us_d", "ud_s", "ds_u")]
        physical = q*(np.linalg.norm(c*ca-b*cb)/np.linalg.norm(a*ca+b*cb))**2
        assert abs(physical-central) < 1e-12
    lo, hi = gram_envelope(p, n, INPUT["illustrative_common_Gram_floor"])
    eps = INPUT["extra_amplitude_fraction"]
    return {"ratio_Gamma_nK0_over_pKplus": central, "phase_space_ratio": q,
            "working_rectangles_ratio": width_boxes,
            "direct_90_lower_years": {"p": tau_p, "n": tau_n},
            "model_transfer_neutron_lower_years": {
                "central":tau_p/central,
                **{f"worst_{i}_width":tau_p/bounds[1] for i,bounds in width_boxes.items()}},
            "ratio_where_direct_neutron_limit_becomes_stronger":tau_p/tau_n,
            "common_Gram_sharp_central_ratio": [q*lo,q*hi],
            "independent_extra_amplitude_central_ratio":
                [central*((1-eps)/(1+eps))**2, central*((1+eps)/(1-eps))**2]}


def tidy(x):
    if isinstance(x, dict):
        return {k: tidy(i) for k,i in x.items()}
    if isinstance(x, list):
        return [tidy(i) for i in x]
    return float(f"{x:.11g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = tidy(run())
    path = BASE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
    elif args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
        print("1483: both kaon widths, common Gram and lifetime transfer verified")
    else:
        print(json.dumps(result,indent=2))
