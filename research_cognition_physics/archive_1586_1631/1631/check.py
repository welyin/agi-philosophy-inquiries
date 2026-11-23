"""Neutral triplet proper decay length with the gauge-threshold partner fixed."""
import argparse
import json
import math
from pathlib import Path

D = Path(__file__).resolve().parent
EV = 1e-9


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def rates_per_light_mass(mass, vev, masses):
    def vector(m):
        x = (m/mass)**2
        return (1-x)**2*(1+2*x) if x < 1 else 0.
    def scalar(m):
        x = (m/mass)**2
        return (1-x)**2 if x < 1 else 0.
    scale = mass**2/(16*math.pi*vev**2)
    return {"W":2*scale*vector(masses["W"]), "Z":scale*vector(masses["Z"]),
            "h":scale*scalar(masses["h"])}


def run():
    cfg=read(D/"inputs.json")
    Plo, Phi=read(D/cfg["gauge_result"])["whole_triplet_geometric_mean_GeV"]
    box=read(D/cfg["oscillation_input"])["box"]
    v=read(D/cfg["scalar_input"])["vev_GeV"]
    hbar_c=cfg["hbar_c_GeV_mm"]
    lo,hi=box["abs31_eV2"]
    rlo,rhi=box["solar_over_abs31"]
    light_bounds={
        "normal":(math.sqrt(lo*rlo),math.sqrt(hi),math.sqrt(hi)*(1-math.sqrt(rlo))),
        "inverted":(math.sqrt(lo),math.sqrt(hi*(1+rhi)),
                    math.sqrt(hi)*(math.sqrt(1+rhi)-1)) }
    cap=cfg["heavy_partner_Y_norm_ceiling"]
    rows=[]
    for mass in cfg["lighter_triplet_mass_GeV"]:
        assert mass<=Plo
        weights=rates_per_light_mass(mass,v,cfg["masses_GeV"])
        width_coeff=math.fsum(weights.values())
        for name,(min_a,max_b,max_diff) in light_bounds.items():
            weak_min=min_a*EV/mass
            weak_max=v*v*cap*cap/(2*Plo*Plo)+max_diff*EV/mass
            width_low=width_coeff*mass*weak_min
            width_high=width_coeff*mass*weak_max
            rows.append({"ordering":name, "light_triplet_GeV":mass,
                         "other_triplet_GeV":[Plo*Plo/mass,Phi*Phi/mass],
                         "phase_space_fraction_of_asymptotic":width_coeff/(mass**2/(4*math.pi*v*v)),
                         "min_a_eV":min_a, "max_b_eV":max_b, "max_pair_difference_eV":max_diff,
                         "weak_mixing_squared":[weak_min,weak_max],
                         "neutral_width_GeV":[width_low,width_high],
                         "ctau_mm":[hbar_c/width_high,hbar_c/width_low],
                         "real_angle_ctau_min_mm":hbar_c/(width_coeff*max_b*EV)})
    return {"bounds":rows, "isolated_gauge_term_in_mixing_upper":v*v*cap**2/(2*Plo**2)}


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=run()
    f=D/"results.json"
    if args.check:
        if result!=read(f):
            raise SystemExit("1631 values differ")
        print("1631 PASS: both orderings, two finite-mass widths and partner bounds")
    else:
        with f.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
