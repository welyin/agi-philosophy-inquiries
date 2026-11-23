"""On-shell longitudinal W pairs and independent charge/weak-pole calibrations."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
POLE = json.loads((HERE/CFG["pole_ref"]).read_text(encoding="utf-8"))
GF = POLE["formula"]["Gmu_GeV_minus2"]
MZ = POLE["calibration"]["m_Z"][0]
MW = CFG["nominal_W_mass_GeV"]
GV, DGV = CFG["Z_flavor_specific"]["gV"]
GA, DGA = CFG["Z_flavor_specific"]["gA"]
PB = CFG["pb_per_GeV_inverse_square"]
BR, DBR = CFG["BR_e"]
WIDTH, DWIDTH = CFG["W_width_GeV"]


def amplitudes(energy, z, sign, kappa_photon=1.0):
    s = energy**2
    beta = math.sqrt(1-4*MW**2/s)
    sine = np.sqrt(1-z*z)
    left = MZ*math.sqrt(math.sqrt(2)*GF)*(GV+GA)
    right = MZ*math.sqrt(math.sqrt(2)*GF)*(GV-GA)
    charged = math.sqrt(48*math.pi*BR*WIDTH/MW)
    zww = sign*charged*MW/MZ
    charge2 = 4*math.pi/CFG["effective_inverse_alpha_nominal"]
    prop = s/(s-MZ**2)
    b = s/(2*MW**2)
    photon = -charge2*(1+b*kappa_photon)
    def neutral(coupling):
        return beta*sine*(photon + coupling*zww*prop*(1+b))
    t = MW**2-s/2*(1-beta*z)
    numerator = s/(2*MW**2)*sine*(s/2*(z-beta)-beta*MW**2)
    return neutral(right), neutral(left)+charged**2*numerator/(2*t)


def nominal(energy, sign):
    z, w = np.polynomial.legendre.leggauss(96)
    r,l = amplitudes(energy,z,sign)
    s = energy**2
    norm = PB*math.sqrt(1-4*MW**2/s)/(128*math.pi*s)
    return {"right_pb":float(norm*np.dot(w,r*r)), "left_pb":float(norm*np.dot(w,l*l))}


def robust_bound(kappa):
    # Compare an opposite sign of the entire YM Z tensor against the photon.
    # Positive factors have independent lower endpoints on each energy cell.
    energy = np.linspace(*CFG["hard_energy_support_GeV"],8193)
    lo,hi = energy[:-1],energy[1:]
    mlo,mhi = CFG["W_mass_window_GeV"]
    zlo,zhi = MZ-2*POLE["calibration"]["m_Z"][1],MZ+2*POLE["calibration"]["m_Z"][1]
    gw_times_mw = math.sqrt(48*math.pi*(WIDTH-2*DWIDTH)*(BR-2*DBR)*mlo)
    electron = GV-GA-2*(DGA+DGV)
    product_lo = gw_times_mw*math.sqrt(math.sqrt(2)*GF)*electron
    qlo = 4*math.pi/max(CFG["effective_inverse_alpha_range"])
    bw_real = 1/(1-zlo**2/hi**2)/(1+(zhi*CFG["max_Zwidth_GeV"]/(lo**2-zhi**2))**2)
    bmin = lo**2/(2*mhi**2)
    term = qlo*(1+kappa*bmin)+product_lo*bw_real*(1+bmin)
    sigma_lo = PB*(1-4*mhi**2/lo**2)**1.5*term**2/(96*math.pi*hi**2)
    # A global affine minorant of the lower rectangles remains valid at the
    # observed luminosity mean. A positive slope uses each upper cell edge.
    mean_min = CFG["mean_rounding_lower_GeV"]
    slopes = np.linspace(0,0.20,601)
    candidates = [float(np.min(sigma_lo-q*hi) + q*mean_min)-1e-8 for q in slopes]
    best = int(np.argmax(candidates))
    avg = candidates[best]
    arbitrary_weights = float(np.min(sigma_lo))-1e-8
    cap = CFG["total_uniform_contract_pb"]
    data_upper = sum([CFG["DELPHI"]["table2_ratio"][0],2*CFG["DELPHI"]["table2_ratio"][1]])
    def transfer(bound):
        return 1-data_upper/(bound/cap)
    return {"product_Z_lower":product_lo,"any_energy_mix_right_pb_lower":arbitrary_weights,
            "measured_mean_mix_right_pb_lower":avg,"minorant_slope_pb_per_GeV":float(slopes[best]),
            "two_width_ratio_upper":data_upper,"any_mix_ratio_lower":arbitrary_weights/cap,
            "mean_mix_ratio_lower":avg/cap,"mean_mix_fractional_loss_margin":transfer(avg),
            "max_total_for_two_width_test_pb":avg/data_upper,
            "incoherent_or_coherent_escape_amplitude_L2_min_sqrt_pb":
                max(0.,math.sqrt(avg)-math.sqrt(cap*data_upper))}


def compute():
    rows = CFG["LEP_table53_energy_total_pb_error"]
    lephi = max(r[1]+2*r[2] for r in rows)
    assert lephi < CFG["total_uniform_contract_pb"]
    combos = {str(sign):{str(e):nominal(e,sign) for e in (189,198.2,200,206)}
              for sign in (1,-1)}
    branches = {str(k):robust_bound(k) for k in (1.,.90,.76)}
    return {"nominal_at_each_energy":combos,"opposite_product_uniform_bounds":branches,
            "measured_points_max_two_width_pb":lephi,
            "nominal_weak_g_from_W":math.sqrt(48*math.pi*WIDTH*BR/MW),
            "nominal_Z_right_effective":MZ*math.sqrt(math.sqrt(2)*GF)*(GV-GA)}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for j in a:
            compare(a[j],b[j])
    else:
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-10),(a,b)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write",action="store_true")
    args=cli.parse_args()
    results=compute()
    dest=HERE/"results.json"
    if args.write:
        with dest.open("x",encoding="utf-8") as f:
            json.dump(results,f,ensure_ascii=False,indent=2)
            f.write("\n")
    else:
        compare(results,json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(results))
