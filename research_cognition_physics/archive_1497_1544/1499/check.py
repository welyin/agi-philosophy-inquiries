"""Selector for 331 up alignment and charged-current to signed neutral-B holdout."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def interval(values, n):
    return [values[0]-n*sum(values[1:]), values[0]+n*sum(values[1:])]


def unitary(us, cb, rho, gamma):
    ub = rho * cb
    c13 = math.sqrt(1-ub**2)
    s12, s23 = us/c13, cb/c13
    c12, c23 = math.sqrt(1-s12**2), math.sqrt(1-s23**2)
    h = c12*s23*ub / (s12*c23)
    delta = gamma + math.asin(h*math.sin(gamma))
    z = np.exp(1j*delta)
    return np.array([
        [c12*c13, s12*c13, ub*z.conjugate()],
        [-s12*c23-c12*s23*ub*z, c12*c23-s12*s23*ub*z, s23*c13],
        [s12*s23-c12*c23*ub*z, -c12*s23-s12*c23*ub*z, c23*c13]
    ])


def selectors(charges):
    table = {}
    for row, scalar, conjugate, right in [
        ("Qplus", "eta", False, ["u_light_R","u_special_R"]),
        ("Qminus", "rho", True, ["u_light_R","u_special_R"]),
        ("Qplus", "rho", False, ["all_down_R"]),
        ("Qminus", "eta", True, ["all_down_R"]),
    ]:
        table[row + "_" + right[0]] = [
            name + "-" + fermion
            for name in [scalar+"D", scalar+"U"] for fermion in right
            if (charges[row] + (1 if conjugate else -1)*charges[name] - charges[fermion]) % 3 == 0
        ]
    return table


def samples(data, v):
    example = data["independent_mass_source_example_GeV"]
    u = example["v"]/math.sqrt(2)
    eta, rho = example["etaU"], example["rhoU"]
    eta_d = example["etaD"]
    rho_d = math.sqrt(u*u-eta*eta-rho*rho-eta_d*eta_d)
    # The singlet fixes both relative phases. Check all three discrete vacua.
    e1, e2, e3 = np.eye(3)
    for branch in range(3):
        phase = np.exp(2j * np.pi * branch / 3)
        singlet = 10 * phase
        eta_u_vec, rho_u_vec = eta * phase * e1, rho * phase * e2
        chi_vec = 1200 * (phase.conjugate()**2) * e3
        assert abs(singlet**2 - 10 * singlet.conjugate()) < 1e-12
        assert np.isclose(np.vdot(eta_d*e1, eta_u_vec), eta_d*eta*singlet/10)
        assert np.isclose(np.vdot(rho_d*e2, rho_u_vec), rho_d*rho*singlet/10)
        assert np.allclose(np.cross(eta_u_vec, rho_u_vec), eta*rho/1200*chi_vec.conjugate())
    ys = read(data["matched_masses"])["yukawa"]
    mu, md = (u*np.array([row[0] for row in ys[name]]) for name in ["up","down"])
    light_norm, special = max(mu[:2])/rho, mu[2]/eta
    assert eta * example["lower_special_yukawa"] > rho * example["upper_light_block_norm"]
    assert special > example["lower_special_yukawa"] and light_norm < example["upper_light_block_norm"]
    down = v @ np.diag(md)
    factors = [np.diag([1/eta_d,1/eta_d,0]), np.diag([0,0,1/rho_d])]
    couplings = [f @ down for f in factors]
    assert np.allclose(eta_d*couplings[0] + rho_d*couplings[1], down)
    p = np.diag([0,0,1])
    pd = v.conjugate().T @ p @ v
    assert np.allclose(v.conjugate().T @ couplings[1], pd @ np.diag(md)/rho_d)
    return {"rhoD_GeV":rho_d,"special_yukawa":special,"light_up_source_norm":light_norm,
            "moduli_K_Bd_Bs_D": [abs(pd[0,1]),abs(pd[0,2]),abs(pd[1,2]),0.0]}


def certificate():
    data = read("inputs.json")
    allowed = selectors(data["Z3_four_component_charges"])
    assert list(allowed.values()) == [["etaU-u_special_R"], ["rhoU-u_light_R"],
                                     ["rhoD-all_down_R"], ["etaD-all_down_R"]]
    tree, belle = read(data["tree"]), read(data["belle_and_gamma"])
    n = data["tree_moduli_widths"]
    kernel, form = interval(tree["K_product"], n), interval(tree["K_form"], n)
    us = [kernel[0]/form[1],kernel[1]/form[0]]
    cb = interval(belle["BelleII_Vcb_BGL"], n)
    rho = interval(tree["ub_over_cb"], n)
    ubmax = cb[1]*rho[1]
    c13min = math.sqrt(1-ubmax**2)
    s12max, s23max = us[1]/c13min, cb[1]/c13min
    cdlo = us[0]*math.sqrt(1-s23max**2)-s23max*ubmax
    cdhi = s12max+s23max*ubmax
    udlo = math.sqrt(1-us[1]**2-ubmax**2)
    rlo, rhi = rho[0]*udlo/cdhi, rho[1]/cdlo
    gmin,gmax = np.radians(belle["gamma_expanded_1_width"])
    # beta=arg(1 - Rb*exp(-i gamma)). Retain the charged-tree mirror at gamma+pi.
    sine_lo, sine_hi = math.sin(gmin), math.sin(gmax)
    cos_lo, cos_hi = math.cos(gmax), math.cos(gmin)
    beta = [math.atan(rlo*sine_lo/(1-rlo*cos_lo)),
            math.atan(rhi*sine_hi/(1-rhi*cos_hi))]
    mirror = [-math.atan(rhi*sine_hi/(1+rlo*cos_lo)),
              -math.atan(rlo*sine_lo/(1+rhi*cos_hi))]
    cp = data["signed_holdout"]
    sdata = interval(cp["S"], cp["widths"])
    cdata = interval(cp["C"], cp["widths"])
    mirror_gap = math.asin(sdata[0])-2*mirror[1]
    assert 2*mirror[0] > -math.pi/2 and sdata[0]>0
    point = data["frequency_witness"]
    v = unitary(tree["K_product"][0]/tree["K_form"][0], belle["BelleII_Vcb_BGL"][0],
                point["rho"], math.radians(point["gamma_deg"]))
    assert np.allclose(v @ v.conjugate().T, np.eye(3))
    actual_beta = np.angle(-v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate()))
    observed_ratio = point["Bd_MeV"]/point["Bs_MeV"]/point["xi"]**2*abs(v[2,0]/v[2,1])**2
    s_beta = math.sin(2*actual_beta)
    assert sdata[0] < s_beta < sdata[1] and point["observed_ratio"][0] < observed_ratio < point["observed_ratio"][1]
    return {"allowed_ordinary_yukawas":allowed,"Rb_interval":[rlo,rhi],
            "beta_direct_outer_degrees":np.degrees(beta).tolist(),
            "beta_mirror_outer_degrees":np.degrees(mirror).tolist(),
            "S_holdout_interval":sdata,"C_holdout_interval":cdata,
            "mirror_total_phase_minimum_degrees_lower_bound":math.degrees(mirror_gap),
            "one_common_point":{"sin_2beta":s_beta,"beta_degrees":math.degrees(actual_beta),
                                "delta_md_over_ms":observed_ratio, **samples(data, v)}}


def equal(a,b):
    if isinstance(a, dict):
        return a.keys()==b.keys() and all(equal(x,b[k]) for k,x in a.items())
    if isinstance(a, list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(a, str):
        return a==b
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    g=parser.add_mutually_exclusive_group()
    g.add_argument("--write",action="store_true")
    g.add_argument("--check",action="store_true")
    args=parser.parse_args()
    cert=certificate()
    result=ROOT/"results.json"
    if args.write:
        with result.open("x",encoding="utf-8") as dest:
            json.dump(cert,dest,indent=2)
            dest.write("\n")
    elif not equal(cert,read("results.json")):
        raise ValueError("Mismatch in protected flavor orientation certificate")
    print(json.dumps(cert,indent=2))
