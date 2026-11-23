"""Calibrate the absolute weak-coupling tradeoff without a CKM-unitarity fit."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def compute():
    inp = read("inputs.json")
    ew = read(inp["mass_and_GF"])
    spectrum = read(inp["charged_masses"])["yukawa"]
    low = read("../1381/inputs.json")
    shape = read("../1398/independent_results.json")
    mass_z, dz = ew["calibration"]["m_Z"]
    alpha_inv, gf = (ew["formula"][k] for k in ("alpha_inverse", "Gmu_GeV_minus2"))
    ratio_norm = sqrt(2) * gf * alpha_inv/pi
    r_high = inp["ew_working_factor"][1]
    rows = {}
    for n in inp["widths"]:
        def edge(q, sign):
            return q[0]+n*sign*sum(q[1:])
        ratios = {s: [(edge(m,1)/edge(spectrum[s][-1],-1))**2 for m in spectrum[s][:2]]
                  for s in ("up","down")}
        upper = 1 + sum(ratios["up"]) + ratios["down"][1]
        rlo = edge(inp["pole"]["converted_sum"], -1)
        rhi = edge(inp["pole"]["converted_sum"], 1)
        ell2_pole = upper/rlo
        ell2_low = (shape[str(n)]["max_inert_sigma"]/shape[str(n)]["min_light_sigma"])**2
        def four_norm(sign):
            return (edge(low["ud_beta"],sign)**2 + edge(low["charm_d"],sign)**2
                    + edge(low["charm_s"],sign)**2 +
                    (edge(low["K_product"],sign)/edge(low["K_form"],-sign))**2)
        masses = {}
        for label, pair in ew["measured"].items():
            wl, wh = edge(pair,-1), edge(pair,1)
            zlo, zhi = mass_z-n*dz,mass_z+n*dz
            def k(w,z):
                return ratio_norm*w*w*(1-(w/z)**2)
            kmin, kmax = k(wh,zlo), k(wl,zhi)
            mass_out = {}
            for tag,ell2 in (("pole_alone",ell2_pole),("with_low_shape",min(ell2_pole,ell2_low))):
                s2_upper = r_high*ell2/(ratio_norm*wl**2)
                rho_upper = (wh/zlo)**2/(1-s2_upper)
                gl2_min = 4*sqrt(2)*gf*wl**2/(r_high*ell2)
                mass_out[tag] = {"min_muon_factor_if_mass_angle": kmin/ell2,
                                 "angle_squared_sine_upper_if_bounded_factor": s2_upper,
                                 "rho_upper_if_bounded_factor": rho_upper,
                                 "lower_g_squared_over_16pi2": gl2_min/(16*pi*pi)}
            masses[label] = {"K_range": [kmin,kmax], "certificates": mass_out}
        rows[str(n)] = {"pole_ratio": [rlo,rhi], "inert_hadronic_strength_upper": upper,
                        "ell_upper_pole": sqrt(ell2_pole),
                        "ell_upper_low_shape": sqrt(ell2_low),
                        "apparent_bottom_column_squared_from_subtraction":
                            [rlo-four_norm(1),rhi-four_norm(-1)], "mass_calibrations": masses}
    return rows


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x,y in zip(a,b):
            compare(x,y)
    else:
        assert abs(a-b) < 1e-11*max(1,abs(b))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true")
    result = compute()
    if p.parse_args().check:
        compare(result, read("results.json"))
        print("1400 weak-ratio and absolute-calibration bounds reproduced")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
