"""SM gauge-singlet probe of the previously registered spacelike flavor responses."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
ETA = np.diag([1., -1., -1., -1.])


def read(name):
    return json.loads((HERE / name).read_bytes())


def dot(p, q):
    return p @ ETA @ q


def photon_polarizations(k):
    n = k[1:] / k[0]
    first = np.cross(n, np.array([0., 1., 0.]))
    first /= np.linalg.norm(first)
    return [np.r_[0., first], np.r_[0., np.cross(n, first)]]


def event(t, mass, cfg):
    e, eh, mh = (cfg[k] for k in ("incoming_energy", "higgs_energy", "higgs_mass"))
    ph = math.sqrt(eh**2 - mh**2)
    cos_theta = (eh - (mh**2 + t)/(2*e)) / ph
    assert abs(cos_theta) < 1
    h = np.array([eh, ph*math.sqrt(1-cos_theta**2), 0., ph*cos_theta])
    p = np.array([e, 0., 0., e])
    k = np.array([e, 0., 0., -e])
    recoil = p+k-h
    mr = math.sqrt(dot(recoil, recoil))
    w = (mr**2-mass**2)/(2*mr)
    direction = np.array([.48, .64, .60])
    beta = recoil[1:] / recoil[0]
    gamma = recoil[0] / mr
    longitudinal = np.dot(beta, w*direction)
    kp = np.r_[gamma*(w+longitudinal),
               w*direction + (gamma*w+(gamma**2/(gamma+1))*longitudinal)*beta]
    xp = recoil-kp
    np.testing.assert_allclose([dot(h,h),dot(kp,kp),dot(xp,xp),-dot(p-h,p-h)],
                               [mh**2,0.,mass**2,t],atol=8e-14)
    assert min(xp[0],kp[0],h[0]) > 0
    return p,k,h,xp,kp


def reduced_rate(form, t, mass, cfg, coupling):
    p,k,h,xp,kp = event(t,mass,cfg)
    tensor = dot(k,kp)*ETA - np.outer(kp,k)
    ward = max(np.max(abs(k@ETA@tensor)), np.max(abs(tensor@ETA@kp)))
    photon = sum(abs(e@ETA@tensor@ETA@ep)**2 for e in photon_polarizations(k)
                 for ep in photon_polarizations(kp)) / 2
    np.testing.assert_allclose(photon, dot(k,kp)**2)
    spin = 2*dot(p,xp)
    u = -dot(k-kp,k-kp)
    matrix_element_squared = (coupling*abs(form))**2 * spin * photon / (u+cfg['spectator_mass']**2)**2
    stripped = matrix_element_squared * (u+cfg['spectator_mass']**2)**2 / (spin*dot(k,kp)**2)
    return stripped, ward


def inverse_transfer(x, w, y, t):
    value = 0j
    for pole, weight in zip(x,w):
        phase = np.exp(.37j)
        b = np.array([[math.sqrt(pole), y*phase],[-y*phase.conjugate(),math.sqrt(pole)]])
        inv_flip = b @ np.linalg.inv(t*np.eye(2)+b.conj().T@b)
        value += weight*(1+pole)*inv_flip[0,1]
    return value


def certificate():
    cfg = read("inputs.json")
    kin = cfg["kinematics"]
    lower = read(cfg["bounds"])["single_lower_when_Q_at_least_largest_neutral_mass"]
    upper = read(cfg["bounds"])["positive_multipole_upper_at_the_same_scale"]
    positive = read(cfg["double"])
    single = read(cfg["single"])["witnesses"][0]
    r = math.sqrt(7)/5
    alpha = kin["low_fraction"]
    predictions = {}
    for name, q, x, weights in [
            ("single",1/math.sqrt(3),[single["x"]],[1.]),
            ("two_positive",.3,positive["masses_squared_over_g2"],positive["pole_fractions_at_reference"])]:
        low = alpha*(min(x)+min(r,q)**2)
        high = max(x)+r**2
        rates = []
        wards = []
        for t in [low,high]:
            row = []
            for y, mass, coupling in zip([q,r],kin["final_fermion_masses"],[.37,.81]):
                reduced, ward = reduced_rate(inverse_transfer(x,weights,y,t),t,mass,kin,coupling)
                row.append(reduced)
                wards.append(ward)
            rates.append(row)
        value = rates[1][0]*rates[0][1]/(rates[1][1]*rates[0][0])
        assert max(wards) < 5e-12
        if name == "single":
            assert value > lower/(1+alpha)**2
        else:
            assert value < upper*(1+alpha)**2
        predictions[name] = dict(Q2_low=round(low,10), Q2_high=round(high,10),
                                 four_rate_ratio=round(value,10), ward_pass=True)

    _,e2,e3 = np.eye(3)
    high_scale = read("../1397/inputs.json")["high_scale"]
    center,width = np.asarray(high_scale["down"]+high_scale["angles_theta12_theta13_theta23_delta"][:3]).T
    witness = read("../1678/independent_results.json")["joint_witness"]
    inputs = center+width*np.array(witness["quark_pulls"])
    basis = runpy.run_path(str(HERE.parent/"1678/independent_check.py"))
    _,target,_,_ = basis["get_matrices"](inputs,witness["quark_gamma_high"])
    z,e,f = target[1,2],target[2,1],target[2,2]
    recon = np.outer(z*e2+f*e3,e3)+np.outer(e*e3,e2)
    np.testing.assert_allclose(target,recon,atol=2e-15)
    assert np.linalg.matrix_rank(target)==2
    source = np.array([[0.,1.],[1.,3.],[1.,1.]])
    assert np.linalg.matrix_rank(np.diag([0,0,1])@source)==1
    for hypercharge,xright in [(-1,-1),(0,0)]:
        higgs = -0.5 if hypercharge==0 else 0.5
        assert .5+higgs+hypercharge == 0 and -hypercharge+xright == 0
    return dict(charged_minimum_rank=2,charged_minor_abs_in_frozen_texture=round(float(abs(z*e)),10),
                neutral_selected_rank=1,hypercharge_sums_zero=True,
                max_low_fraction=(lower/upper)**.25-1,
                low_fraction=alpha,universal_A_lower=lower/(1+alpha)**2,
                universal_B_upper=upper*(1+alpha)**2,on_shell_witnesses=predictions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args = parser.parse_args()
    values = certificate()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(values,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        assert read("results.json")==values
        print("PASS: port ranks, hypercharges, photon Ward identities, on-shell kinematics and four-rate gap")
