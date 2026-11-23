"""LEP/SLD Z chiral couplings, a finite bottom mass matrix and a CMS holdout."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def project_ratio(mu, covariance, a, b, radius):
    """Extrema of (a.g)/(b.g) on a covariance ellipsoid with positive denominator."""
    am, bm = a @ mu, b @ mu
    aa, ab, bb = a @ covariance @ a, a @ covariance @ b, b @ covariance @ b
    assert bm > radius*math.sqrt(bb)
    quad = bm*bm-radius**2*bb
    center = am*bm-radius**2*ab
    discr = center**2-quad*(am*am-radius**2*aa)
    return [(center-sign*math.sqrt(max(0.,discr)))/quad for sign in (1,-1)]


def run():
    cfg = read(HERE/"inputs.json")
    baseline = read(HERE/cfg["bottom_and_threshold"])
    m = baseline["kinematic_GeV"]["b"]
    md, ms = cfg["example_doublet_hard_GeV"], cfg["additional_singlet_hard_GeV"]
    hdata = read(HERE/cfg["CMS_holdout"])["CMS_ratio_bb_ZZ_normalized"]
    mu, sig = np.array(cfg["central"]), np.array(cfg["error"])
    cov = np.outer(sig,sig)*np.array(cfg["correlation"])
    assert min(np.linalg.eigvalsh(cov))>0
    delta_r = np.array([0,0,-1/3,1])
    delta_l = np.array([-1,1,2/3,0])
    den = np.array([-1,0,1,0])  # k/2
    dr, dl, half_k = delta_r@mu, delta_l@mu, den@mu
    sigma_r = math.sqrt(delta_r@cov@delta_r)
    t = (m/md)**2
    numerator = delta_l+(1-2*t)*delta_r
    kappa = 1 - numerator@mu/half_k
    right_d, left_s = dr/half_k, (dl-2*t*dr)/half_k
    left_d, right_s = t*right_d, (m/ms)**2*left_s

    left = np.sqrt([1-left_s-left_d,left_s,left_d])
    right = np.sqrt([1-right_s-right_d,right_s,right_d])
    x_l = (m*right[1]-ms*left[1])/left[0]
    x_r = (m*left[2]-md*right[2])/right[0]
    m0 = (m*left[0]-x_l*right[1])/right[0]
    matrix = np.array([[m0,x_l,0],[0,ms,0],[x_r,0,md]])
    hard = np.diag([0,ms,md])
    assert max(np.linalg.norm(matrix@right-m*left),np.linalg.norm(matrix.T@left-m*right))<1e-12
    assert abs(left@(matrix-hard)@right/m-kappa)<1e-12
    masses = np.sort(np.linalg.svd(matrix,compute_uv=False))
    assert abs(masses[0]-m)<1e-12

    envelope = []
    for radius in (1,2):
        low, high = project_ratio(mu,cov,numerator,den,radius)
        kappas = [1-high,1-low]
        rates = [v*v for v in kappas]
        cms = [hdata[0]-radius*hdata[1],hdata[0]+radius*hdata[2]]
        smallest_den = half_k-radius*math.sqrt(den@cov@den)
        extra = (1+abs(1-2*t))*cfg["comparison_residual_budget_each"]/smallest_den
        envelope.append({"radius":radius, "min_signed_residual_R":dr-radius*sigma_r,
                         "kappa":kappas,"rate":rates,"CMS":cms,
                         "overlap": bool(max(cms[0],rates[0])<=min(cms[1],rates[1])),
                         "largest_allowed_rate_transfer":cms[1]/rates[0],
                         "kappa_outer_with_effective_residuals":[kappas[0]-extra,kappas[1]+extra]})
    return {"right_difference":dr, "right_difference_sigma":sigma_r, "zero_pull":dr/sigma_r,
            "left_difference":dl,"k":2*half_k,"right_doublet_probability":right_d,
            "left_singlet_probability":left_s,"kappa_b":kappa, "normalized_bb_over_ZZ":kappa**2,
            "triangular_matrix_GeV":matrix.tolist(), "physical_masses_GeV":masses.tolist(),
            "untouched_exotic_mass_GeV":md, "envelopes":envelope}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args = ap.parse_args()
    result = run()
    out = HERE/"results.json"
    if args.check:
        assert result==read(out)
        print("1633 PASS: Z-pole covariance, finite mass witness and CMS Higgs holdout")
    else:
        with out.open("x",encoding="utf-8") as file:
            file.write(json.dumps(result,indent=2)+"\n")
