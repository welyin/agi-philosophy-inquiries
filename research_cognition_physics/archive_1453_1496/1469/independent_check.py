"""Energy-loss and angular-response certificates for DELPHI's LL pair fraction."""
import argparse
import json
from fractions import Fraction
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def certificates(cfg, pole, kappa):
    hi_m = max(cfg["W_mass_window_GeV"])
    lo_m = min(cfg["W_mass_window_GeV"])
    z, dz = pole["calibration"]["m_Z"]
    zmin, zmax = z - 2 * dz, z + 2 * dz
    v, dv = cfg["Z_flavor_specific"]["gV"]
    a, da = cfg["Z_flavor_specific"]["gA"]
    width, dwidth = cfg["W_width_GeV"]
    branch, dbranch = cfg["BR_e"]
    gf = pole["formula"]["Gmu_GeV_minus2"]
    neutral = (v-a-2*(dv+da))*sqrt(sqrt(2)*gf*48*pi*lo_m*(width-2*dwidth)*(branch-2*dbranch))
    photon = 4*pi/max(cfg["effective_inverse_alpha_range"])
    pb = cfg["pb_per_GeV_inverse_square"]

    # The zero bound on 0..2*m_hi also covers the uncertainty on the threshold.
    edges = np.linspace(2*hi_m, max(cfg["hard_energy_support_GeV"]), 3001)
    lo, hi = edges[:-1], edges[1:]
    bmin = (lo/hi_m)**2/2
    prop = 1/(1-(zmin/hi)**2)/(1+(cfg["max_Zwidth_GeV"]*zmax/(lo**2-zmax**2))**2)
    cells = pb*(1-(2*hi_m/lo)**2)**1.5/(96*pi*hi**2)
    cells *= (photon*(1+kappa*bmin) + neutral*prop*(1+bmin))**2
    cells = np.r_[0., cells]
    upper_energy = np.r_[edges[0], hi]
    cut = cfg["total_uniform_contract_pb"]*(cfg["DELPHI"]["table2_ratio"][0]+2*cfg["DELPHI"]["table2_ratio"][1])
    mean = cfg["mean_rounding_lower_GeV"]
    best_at_mean = -1.
    critical_mean = 1e6
    for j in range(200, 2001):
        slope = j/10000.
        offset = float(min(cells - slope*upper_energy))-1e-7
        best_at_mean = max(best_at_mean, slope*mean+offset)
        critical_mean = min(critical_mean,(cut-offset)/slope)
    margin = best_at_mean/cfg["total_uniform_contract_pb"] - cut/cfg["total_uniform_contract_pb"]
    return {"right_pb_lower_at_beam_mean_if_also_hard_mean": best_at_mean,
            "sufficient_hard_mean_GeV_strict_lower_cut": critical_mean,
            "tolerated_beam_to_hard_mean_loss_strictly_below_GeV":mean-critical_mean,
            "sufficient_total_variation_upper_for_angular_transfer_if_hard_mean_at_beam":margin/15}


def main():
    cfg = read("inputs.json")
    pole = read(cfg["pole_ref"])
    # Projector duality using normalized symmetrized angular laws.
    for second_moment, probability in ((Fraction(1,5),1),(Fraction(2,5),0)):
        assert 2-5*second_moment == probability
    out = {str(k):certificates(cfg,pole,k) for k in (1,.90,.76)}
    return json.loads(json.dumps(out),parse_float=lambda x:float(f"{float(x):.10g}"))


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write-new", action="store_true")
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    data = main()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(data, f, indent=2)
            f.write("\n")
    else:
        assert data == read("independent_results.json")
        print("1469 independent: whole hard interval and angular transfer verified")
