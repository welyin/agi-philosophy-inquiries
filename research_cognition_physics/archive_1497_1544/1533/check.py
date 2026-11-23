"""Reconstruct neutral couplings, a soft-Z2 potential and independent rate windows."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

DIR = Path(__file__).resolve().parent
I = json.loads((DIR / "inputs.json").read_text(encoding="utf-8"))


def interval(row, n):
    return row[0] - n * row[1], row[0] + n * row[2]


def couplings(rate, angle, sign):
    beta2 = 1 - 4 * (I["tau_mass_GeV"] / I["higgs_mass_GeV"]) ** 2
    radians = math.radians(angle)
    cs, sn = math.cos(radians), math.sin(radians)
    length = sign * math.sqrt(rate / (cs ** 2 + sn ** 2 / beta2))
    u, w = length * cs, length * sn
    a = 1 / math.sqrt(1 + (1 - u) ** 2 + w ** 2)
    return {"Stau_a": u, "Ptau_a": w, "Stop_a": 2-u,
            "alpha_top_deg": math.degrees(math.atan2(w, 2-u)),
            "h": [a, a*(1-u), -a*w], "scalar_ttH_floor": (2-u)**2}


def soft_potential(h):
    a, x, y = h
    mh, v = I["higgs_mass_GeV"], I["vev_GeV"]
    dx, dy = I["witness_gap_squared_GeV2"]
    sx, sy = np.array([-x/a, 1, 0]), np.array([-y/a, 0, 1])
    mass = mh**2 * np.eye(3) + dx*np.outer(sx, sx) + dy*np.outer(sy, sy)
    q = 1 / math.sqrt(2)
    rotation = np.array([[q, q, 0], [-q, q, 0], [0, 0, 1]])
    weak = rotation.T @ mass @ rotation
    nu, charged_sq = 0.4*mh**2, mh**2 + dx
    l1, l2 = (2*weak[0, 0] - nu)/v**2, (2*weak[1, 1] - nu)/v**2
    l345 = (2*weak[0, 1]+nu)/v**2
    l5 = complex(nu-weak[2, 2], -2*mass[0, 2])/v**2
    l4 = 2*(nu-charged_sq)/v**2 - l5.real
    l3 = l345 - l4 - l5.real
    m12 = complex(nu/2, v**2*l5.imag/4)
    margins = [l1, l2, l3+math.sqrt(l1*l2),
               l3+l4-abs(l5)+math.sqrt(l1*l2)]
    reconstructed = np.array([
        [(l1*v**2+nu)/2, (l345*v**2-nu)/2, -l5.imag*v**2*q/2],
        [(l345*v**2-nu)/2, (l2*v**2+nu)/2, -l5.imag*v**2*q/2],
        [-l5.imag*v**2*q/2, -l5.imag*v**2*q/2, nu-v**2*l5.real]])
    masses, vectors = np.linalg.eigh(mass)
    assert min(margins)>0 and np.max(abs(weak-reconstructed)) < 1e-9
    assert np.linalg.norm(mass @ h - mh**2*np.array(h)) < 1e-8
    assert len(set(np.round(masses,5))) == 3 and charged_sq > 0
    row = vectors.T
    tau_even, top_even, odd = row[:,0]-row[:,1], row[:,0]+row[:,1], -row[:,2]
    assert max(abs(tau_even@odd), abs(top_even@odd), abs(row[:,0]@odd)) < 1e-12
    return {"neutral_masses_GeV": np.sqrt(masses).tolist(), "charged_mass_GeV": math.sqrt(charged_sq),
            "lambda1_to_4": [l1,l2,l3,l4], "lambda5_real_imag": [l5.real,l5.imag],
            "m12_squared_real_imag_GeV2": [m12.real,m12.imag], "BFB_margins": margins,
            "Im_m12sq_squared_lambda5star_over_v4": (m12**2*l5.conjugate()).imag/v**4}


def working_box(n):
    rmin = interval(I["cms_tau_over_z_sm"], n)[0] / max(I["rate_transfer"])
    center, width = I["cms_alpha_tau_deg"]
    angular_max = abs(center)+n*width+I["spin_offset_deg"]
    umin = couplings(rmin, angular_max, 1)["Stau_a"]
    data, sm = I["atlas_gamma_data_fb"], I["atlas_gamma_SM_fb"]
    top, vbf = interval(data["ttH"], n), interval(data["VBF"], n)
    top_sm, vbf_sm = interval(sm["ttH"], n), interval(sm["VBF"], n)
    atlas_range = [top[0]/vbf[1] * vbf_sm[0]/top_sm[1],
                   top[1]/vbf[0] * vbf_sm[1]/top_sm[0]]
    ideal = (2+umin)**2
    return {"min_abs_Stau_over_a": umin, "ideal_mirror_ttH_over_VBF_min": ideal,
            "transferred_mirror_min": I["min_hard_transfer"]*ideal,
            "ATLAS_range": atlas_range, "critical_hard_transfer": atlas_range[1]/ideal,
            "mirror_disjoint": bool(I["min_hard_transfer"]*ideal > atlas_range[1])}


def run():
    result = {}
    for sign in (1, -1):
        h = couplings(I["cms_tau_over_z_sm"][0], I["cms_alpha_tau_deg"][0], sign)
        h["soft_Z2_stationary_witness"] = soft_potential(h["h"])
        result["positive" if sign == 1 else "negative"] = h
    result["windows"] = {str(n):working_box(n) for n in (1,2)}
    assert all(x["mirror_disjoint"] for x in result["windows"].values())
    return result


def same(left, right):
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(same(left[k],right[k]) for k in left)
    return np.allclose(left, right, rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = cli.parse_args()
    path = DIR / "results.json"
    values = run()
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(values, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert same(values, json.loads(path.read_text(encoding="utf-8")))
    print("1533 OK", values["windows"], flush=True)
