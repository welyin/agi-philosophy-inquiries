"""Exact upper envelope from a muon ellipse and an untagged neutrino rate."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def mu_interval(experiment, sm, widths):
    center = experiment["branching_in_1e-9"]
    base = sm["branching_in_1e-9"]
    lo = center - widths * sum(experiment["lower_widths"])
    hi = center + widths * sum(experiment["upper_widths"])
    return [lo / (base + widths * sm["width"]), hi / (base - widths * sm["width"])]


def nu_interval(info, widths):
    delta = (info["total_width_in_1e-6"] + info["tau_width_in_1e-6"]) / info["reference_SD_in_1e-6"]
    return [(info["mu"] - widths*info["width"]) / (1+widths*delta),
            (info["mu"] + widths*info["width"]) / (1-widths*delta)]


def envelope(r, weight, k, sign):
    # z=x+iy, x^2+weight*y^2<=r; maximize |A-B*z|^2 over the ellipse.
    a, b = 1 + sign*k, sign*k
    radius = math.sqrt(r)
    vertex = -a*weight/(b*(1-weight)) if weight < 1 else -math.copysign(math.inf, a*b)
    x = max(-radius, min(radius, vertex))
    y2 = max(0., (r-x*x)/weight)
    amp = math.hypot(a-b*x, b*math.sqrt(y2))
    return amp, x


def evaluate():
    source = json.loads((HERE/"inputs.json").read_bytes())
    cms, lhcb = source["CMS"], source["LHCb"]
    sm, belle = source["sm_Bs"], source["BelleII"]
    lower_a = cms["accepted_light_to_heavy_weight_range"][0]
    ratio = source["sm_coefficient_ratio"]
    tau_h, tau_l = cms["heavy_ps"], cms["light_ps_example"]
    c0,c1 = cms["lifetime_BF_factor"]
    mid_weight = (tau_l/tau_h)*(c0+c1*tau_h)/(c0+c1*tau_l)
    assert lower_a <= mid_weight <= 1
    out = {}
    for experiment, data, weight in (("CMS",cms,lower_a),("LHCb_real_only",lhcb,1.0)):
        rows = {}
        for n in source["working_multipliers"]:
            mu = mu_interval(data,sm,n)
            nr = nu_interval(belle,n)
            opt = {}
            for name, sign in (("singlet",1),("triplet",-1)):
                # Norm is convex in k for each z: endpoints of the k interval suffice.
                amp, x, kval = max((*envelope(mu[1],weight,k,sign),k) for k in ratio["range"])
                # LHCb use of the disk is safe here: extrema are real and allowed.
                assert math.isclose(abs(x),math.sqrt(mu[1]))
                max_rate = (2+amp*amp)/3
                min_other = max(0., (math.sqrt(max(0.,3*nr[0]-2))-amp)/2)
                opt[name] = {"max_Rnu": max_rate, "saturating_z_real": x,
                             "saturating_k":kval, "lower_abs_other_Wilson_over_abs_CL":min_other}
            rows[str(n)] = {"accepted_mu_rate":mu,"Belle_Rnu":nr,"sources":opt}
        out[experiment]=rows
    # Finite 2-width singlet witness, real z, true SM central and SM Belle normalization.
    z, k = -1., ratio["center"]
    singlet_allowed = (2 + (1+k*(1-z))**2)/3
    assert 2 <= singlet_allowed <= 7.2
    assert out["CMS"]["1"]["sources"]["singlet"]["max_Rnu"] < out["CMS"]["1"]["Belle_Rnu"][0]
    assert out["CMS"]["2"]["sources"]["triplet"]["max_Rnu"] < out["CMS"]["2"]["Belle_Rnu"][0]
    # Pauli tensor contraction and change from weak to neutrino-mass basis.
    pauli = (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
    triplet = sum(np.kron(p,p) for p in pauli)
    assert np.allclose(np.diag(triplet), [1,-1,-1,1])
    assert triplet[1,2] == triplet[2,1] == 2
    v = np.array([1,1j,2])/math.sqrt(6)
    p = np.outer(v,v.conj())
    for alpha in (.7j, -2, 1.2+.8j):
        amp = np.eye(3)+alpha*p
        assert math.isclose(np.linalg.norm(amp,"fro")**2, 2+abs(1+alpha)**2)
    target = math.sqrt(3*belle["mu"]-2)
    # Set z=1 using s+t=0; w=target is then realized at any common k.
    return {"accepted_light_relative_weight_example":mid_weight,"comparisons":out,
            "two_width_real_singlet_witness":{"z":z,"Rnu":singlet_allowed},
            "two_coefficient_center_neutrino_amplitude_at_z_one":target,
            "opposite_Wilsons_abs_each_over_abs_CL_at_center":(target-1)/2}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key,value in a.items():
            compare(value,b[key])
    else:
        assert np.allclose(a,b,atol=1e-13,rtol=1e-11)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    value=evaluate()
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(value,f,indent=2)
            f.write("\n")
    else:
        compare(value,json.loads(path.read_bytes()))
    print(json.dumps(value))
