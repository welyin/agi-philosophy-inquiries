"""Tree-level weak charges: analytic elimination and an independent pure-lepton read."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def get(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def run():
    data=get("inputs.json")
    results=get("results.json")
    past=get("../1385/inputs.json")
    # Actual Qweak target and Z=55,N=78; normalize each quark's response.
    charges=np.array([[2.,-1.],[2*55+78.,-(55.+2*78)]])
    inv=np.linalg.inv(charges)
    ang=np.linalg.solve(charges,[4,4*55])
    assert np.allclose(ang,[8/3,4/3])
    shift=.2397-.2381
    mol_width=.0010+.0008+.0006
    transport=.0003
    bounds={}
    for model, (value,stat,atom) in data["atomic_branches"].items():
        rows=[]
        for i,n in enumerate((1,2)):
            pp=past["proton"]
            pmax=pp["pv_electron_scattering_only"][0]-pp["SM"][0]+n*(pp["pv_electron_scattering_only"][1]+pp["SM"][1])
            cmin=value-data["baseline_Cs"][0]-n*(stat+atom+data["baseline_Cs"][1])
            x,y=inv@np.array([pmax,cmin])
            assert np.allclose([x,y],results[model][i]["attaining_powers"])
            assert math.isclose(math.sqrt(x*y),results[model][i]["abs_right_ud_bounds"][1],rel_tol=1e-11)
            zmax=shift+n*mol_width+transport
            u,d=np.array([x,y])+zmax*ang
            assert np.allclose(charges@[u,d]-zmax*np.array([4,220]),[pmax,cmin])
            rows.append(round(math.sqrt(u*d),10))
        bounds[model]=rows

    weak=get("../1460/results.json")["charged_windows"][1]
    old=get("../../archive_1046_1062/1061/public_inputs.json")
    X=float(old["x"]["mean"])+float(old["x"]["sigma"])+float(old["x"]["rounding_half_unit"])
    Y=float(old["y"]["mean"])+float(old["y"]["sigma"])+float(old["y"]["rounding_half_unit"])
    r=float(old["r"]["mean"])-math.hypot(*map(float,old["r"]["sigma_components"]))-float(old["r"]["rounding_half_unit"])
    # Independent bisection of 4S <= (X+t)^2+(Y+rt)^2.
    lo,hi=0.,2.
    for _ in range(80):
        mid=(lo+hi)/2
        if (X+mid)**2+(Y+r*mid)**2>=4*weak["light_ds_left_square_lower"]:
            hi=mid
        else:
            lo=mid
    rsmin=(r*hi-Y)/2
    assert rsmin > 0
    rare=get("../1461/results.json")["branches"]["broad_left_angle_contract"]["right_Z_sd_modulus_cap"]
    prod=weak["abs_RdRs_lower"]
    # Conservative upward padding for a stored rounded rare amplitude.
    f=rare+1e-14
    factor=2-(f/rsmin)**2
    phi=.0011+2*(.0045+.0003)
    xmax=phi/factor
    zmin=(factor*prod/f-phi)/4
    assert xmax*f < prod
    return {"atomic_max_abs_R_with_E158_and_transfer":bounds,
            "E158_shift_upper_one_two_width":[round(shift+i*mol_width+transport,8) for i in (1,2)],
            "one_width_Rus_min":round(rsmin,12),"two_width_proton_only_fixed_angle_RdRs_cap":round(xmax*f,15),
            "necessary_universal_angle_shift_proton_only":round(zmin,10)}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    choices=p.add_mutually_exclusive_group(required=True)
    choices.add_argument("--check",action="store_true")
    choices.add_argument("--write-new",action="store_true")
    opt=p.parse_args(); out=run()
    if opt.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2);f.write("\n")
    else:
        assert out==get("independent_results.json")
    print(json.dumps(out))
