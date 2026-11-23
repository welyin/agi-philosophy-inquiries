"""Cross-run branching-ratio holdout and explicitly stable neutral mass matrix."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def extrema(values):
    numbers = list(values)
    return [min(numbers), max(numbers)]


def width(row, n):
    return [row[0] - n * row[1], row[0] + n * row[-1]]


def fit_box(cfg, n, tau_kind):
    mw, mz, angle = (width(cfg[k], n) for k in
                     ("CMS_MW_GeV", "LEP_MZ_GeV", "SLD_sin2_eff_e"))
    proxy = extrema(w*w / (z*z*(1-s))
                    for w,z,s in itertools.product(mw, mz, angle))
    clo, chi = cfg["rho_proxy_over_kinetic_rho"]
    rho = [max(.5, proxy[0]/chi), min(1., proxy[1]/clo)]
    assert rho[0] <= rho[1]
    eta = cfg["rate_transfer_over_tree_ratio_squared"]
    if tau_kind == "Run1":
        raw = width(cfg["run1_tau_over_Z_absolute"], n)
        sm = width(cfg["run1_tau_over_Z_SM"], n)
        obs_tau = [raw[0]/sm[1], raw[1]/sm[0]]
    else:
        obs_tau = width(cfg["run2_tau_over_Z_normalized"], n)
    obs_w = width(cfg["run2_W_over_Z_normalized"], n)
    tau_mag = [math.sqrt(obs_tau[0]/eta[1]), math.sqrt(obs_tau[1]/eta[0])]
    w_mag = [math.sqrt(obs_w[0]/eta[1]), math.sqrt(obs_w[1]/eta[0])]
    predictions = {}
    for sign in [1, -1]:
        signed_w = extrema((1+(2*r-1)*sign*t)/(2*r)
                           for r,t in itertools.product(rho,tau_mag))
        low, high = signed_w
        minimum_abs = 0 if low <= 0 <= high else min(abs(low), abs(high))
        predictions[str(sign)] = {
            "signed_W_amplitude": signed_w,
            "normalized_W_rate": [eta[0]*minimum_abs**2, eta[1]*max(abs(low),abs(high))**2]
        }
    residuals = {}
    for wsign, tsign in itertools.product([1,-1], repeat=2):
        rr = extrema(2*r*wsign*w - 1 - (2*r-1)*tsign*t
                     for r,w,t in itertools.product(rho,w_mag,tau_mag))
        residuals[f"W{wsign:+d}_tau{tsign:+d}"] = {
            "range":rr, "minimum_abs":max(0,rr[0],-rr[1])}
    assert residuals["W+1_tau+1"]["minimum_abs"] == 0
    assert all(row["minimum_abs"] > 0 for label,row in residuals.items()
               if label != "W+1_tau+1")
    # When lambda_tau=-a and positive lambda_W>=b, rho=(1+a)/(2*(a+b)).
    # This declines with b, and for b<1 it declines with a, so the largest rho
    # is obtained at the two lower endpoints. The wrong-branch small root
    # lambda_W<=-b instead requires a>=1+2b at rho<=1.
    a,b = tau_mag[0],w_mag[0]
    assert b < 1 and tau_mag[1] < 1+2*b
    wrong_rho_max = (1+a)/(2*(a+b))
    return {"rho_proxy":proxy, "rho_selected":rho, "tau_magnitude":tau_mag,
            "W_magnitude":w_mag, "W_observed":obs_w, "predictions":predictions,
            "residuals":residuals, "negative_tau_required_rho_at_most":wrong_rho_max,
            "negative_tau_required_C_at_least":proxy[0]/wrong_rho_max}


def couplings(vev, h):
    d,c,r = vev
    wvec, zvec = vev*np.array([1,2,4]), vev*np.array([1,4,0])
    vw, vz = wvec@vev, zvec@vev
    ref = math.sqrt(vw)
    return vw/vz, ref*(wvec@h)/vw, ref*(zvec@h)/vz, ref*h[0]/d


def witness(cfg):
    x = (2*math.sqrt(6)-3)/10
    c, r, d = math.sqrt(x), math.sqrt(x/2), math.sqrt(1-4*x)
    vev = np.array([d,c,r])
    h = np.array([-d, (1+d*d)/(4*c), (1+d*d)/(8*r)])
    assert abs(h@h-1) < 2e-14
    alt = np.cross(h, [1,0,0])
    alt /= np.linalg.norm(alt)
    rot = np.array([h,alt,np.cross(h,alt)])
    masses = [cfg["h_mass_GeV"],*cfg["witness_other_neutral_masses_GeV"]]
    target = rot.T @ np.diag(np.square(masses)) @ rot
    lock = np.array([[-2*c/d,1,0],[-2*r/d,0,1]])
    lock_hessian = cfg["witness_lock_mass_squared_GeV2"] * lock.T @ lock
    norm_hessian = target - lock_hessian
    assert np.linalg.eigvalsh(norm_hessian)[0] > 0
    out = couplings(vev,h)
    assert np.allclose(out,[1,1,1,-1])
    rho, kw, kz, kt = out
    residual = 2*rho*kw/kz - 1 - (2*rho-1)*kt/kz
    assert abs(residual - 8*r*(h[2]-r*h[0]/d)/(d*h[0]+4*c*h[1])) < 1e-12
    return {"vevs_in_vW":vev.tolist(),"direction":h.tolist(),
            "rho_kW_kZ_kTau":out,"minimum_norm_hessian_GeV2":np.linalg.eigvalsh(norm_hessian)[0],
            "neutral_masses_GeV":np.sqrt(np.linalg.eigvalsh(norm_hessian+lock_hessian)).tolist(),
            "identity_residual":residual}


def remaining_example():
    # Nonzero complex-triplet VEV, ratio chosen by the Run 2 tau central rate.
    rho, ktz = .995, math.sqrt(.8)
    d, c = math.sqrt(2*rho-1), math.sqrt((1-rho)/2)
    hratio = (1/(d*ktz)-d)/(4*c)
    h = np.array([1,hratio,0])/math.hypot(1,hratio)
    measured = couplings(np.array([d,c,0]),h)
    _,kw,kz,kt = measured
    assert abs(kt/kz-ktz) < 1e-12
    return {"rho":rho, "h":h.tolist(),"kW_kZ_kTau":measured[1:],
            "Rtau":(kt/kz)**2, "RW":(kw/kz)**2}


def rounded(v):
    if isinstance(v,dict):
        return {k:rounded(val) for k,val in v.items()}
    if isinstance(v,(tuple,list)):
        return [rounded(val) for val in v]
    return round(float(v),10)


def run():
    cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    outputs = {f"{t}_{n}_width":fit_box(cfg,n,t) for t in ("Run1","Run2") for n in (1,2)}
    outputs["signed_triplet_competitor"]=witness(cfg)
    outputs["positive_branch_example"]=remaining_example()
    return rounded(outputs)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    result=run()
    dst=HERE/"results.json"
    if args.write_new:
        with dst.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    elif args.check:
        assert json.loads(dst.read_text(encoding="utf-8"))==result
        print("1477: cross-run ratios, phases, matching thresholds and scalar potential checked")
    else:
        print(json.dumps(result,indent=2))
