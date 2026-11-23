"""Solar two-level source bounds and coherent versus incoherent jump budgets."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def calculate():
    spec=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    reactor=json.loads((HERE/spec["reactor_inputs"]).read_text(encoding="utf-8"))
    solar=spec["sno"]
    day,juno=reactor["daya_bay"],reactor["juno_primary"]
    rows=[]
    for width in [0,1,2]:
        p=solar["P"]+width*(solar["stat"]+solar["syst"])
        z=(1-math.sqrt(1-day["sin_sq_2theta13"]-width*day["one_sigma"]))/2
        c=1-2*(juno["sin_sq_theta12"]-width*juno["one_sigma"])
        cmin=1-2*(juno["sin_sq_theta12"]+width*juno["one_sigma"])
        k=1+2*(z*z-p)/(1-z)**2
        assert c*c<k<cmin and cmin>0
        g_target=-k/c
        jump=(1-c*g_target-math.sqrt((1-c*c)*(1-g_target*g_target)))/2
        phi=math.acos(g_target)-math.acos(c)
        unitary=np.array([[math.cos(phi/2),-math.sin(phi/2)],
                          [math.sin(phi/2),math.cos(phi/2)]])
        initial=np.sqrt([(1+c)/2,(1-c)/2])
        final=unitary@initial
        actual_g=final[0]**2-final[1]**2
        assert abs(actual_g-g_target)<1e-12 and abs(unitary[1,0]**2-jump)<1e-13
        density_cap = c+k*math.sqrt((1-c*c)/(c*c-k*k))
        # Uniform-material necessary lower bound is decreasing in c on the full interval.
        derivative_ratio_min=k*cmin*(1-k*k)/(math.sqrt(1-cmin*cmin)*(c*c-k*k)**1.5)
        assert derivative_ratio_min>1
        assert abs((c-density_cap)/math.hypot(c-density_cap,math.sqrt(1-c*c))-g_target)<1e-13
        mass_fraction=(1-g_target)/2
        rows.append({"width":width,"mass_eigenstate_low_electron_overlap_fraction_lower":mass_fraction,
                     "adiabatic_fraction_aligned_min":(c-g_target)/(1+c),
                     "adiabatic_fraction_beyond_resonance_min":-g_target,
                     "minimum_mean_coherent_nonadiabatic_jump_if_sign_wrong":jump,
                     "minimum_mean_incoherent_nonadiabatic_jump_if_sign_wrong":(c-g_target)/(1+c),
                     "necessary_essential_sup_of_oriented_potential":density_cap})
    return {"round":1420,"rows":rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=parser.parse_args()
    data=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            f.write(data)
    elif json.loads(path.read_text(encoding="utf-8"))!=json.loads(data):
        raise SystemExit("independent results mismatch")
    print("1420 mass components, matter support and explicit coherent jump verified")
