"""Separate a Cabibbo holdout from other replacements of the measured-input mother."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def solve():
    joint=runpy.run_path(str(HERE/"joint.py"))
    cfg=read("inputs.json")
    setup=read("../1744/inputs.json")
    old=read("../1739/inputs.json")["nearby_measured_inputs_in_marginal_widths"][setup["parent_point"]]
    weak=runpy.run_path(str(HERE/"../1739/check.py"))["DIRECT"]
    mod={k:weak["V"+k][0]+old["weak"][k]*math.hypot(*weak["V"+k][1:]) for k in ("us","ub","cb")}
    shifts={k:old[k] for k in ("up","charm","top","g3")}
    baseline=joint["point"](HERE,mod["us"],mod["ub"],mod["cb"],.25,shifts)
    log={}
    values={}
    o=read("../1484/inputs.json")["DayaBay_sin2_2theta13"]
    for subset in itertools.product((0,1),repeat=3):
        strong,bottom,mother=subset
        tail=(cfg["inclusive_ub"][0],cfg["inclusive_cb"][0]) if bottom else (mod["ub"],mod["cb"])
        ratio=1/cfg["QCDSF_s_d"][0] if strong else baseline["down_over_s"]
        row=joint["shoot"](HERE,ratio,*tail,.25,None if mother else shifts)
        solar,reactor=row["angle_readings"][1],row["angle_readings"][0]
        values[subset]=np.array([row["low_moduli"]["Vus"],reactor,solar,row["sun_calibrated_dm31_eV2"]])
        log["".join(map(str,subset))]=dict(d_over_s=ratio,Vus=values[subset][0],
            reactor=reactor,reactor_pull=(reactor-o[0])/(o[1] if reactor<o[0] else o[2]),
            solar_angle=solar,dm31=row["sun_calibrated_dm31_eV2"])
    np.testing.assert_allclose(values[(0,0,0)][1:], [baseline["angle_readings"][0],baseline["angle_readings"][1],
        baseline["sun_calibrated_dm31_eV2"]],atol=1e-9)
    # Read the independently reported nominal scalar via a frozen observed endpoint, not as a fit input.
    exact=read("../1748/results.json")["outputs"][1]["angle_readings"][0]
    assert abs(values[0,0,0][1]-exact)<3e-6
    components=np.zeros((3,4))
    for order in itertools.permutations(range(3)):
        s=[0,0,0]
        for i in order:
            before=values[tuple(s)]
            s[i]=1
            components[i]+=(values[tuple(s)]-before)/6
    np.testing.assert_allclose(sum(components),values[1,1,1]-values[0,0,0],atol=5e-15)
    assert abs(values[1,1,1][1]-.100694)<1e-5
    data=read("../1369/inputs.json")["g3"]
    delta=4*math.pi*(data[0]**-2-(data[0]+shifts["g3"]*data[1])**-2)
    # At one loop g1=g2 and all principal thresholds stay fixed; only inverse g3 changes.
    return dict(old_weak_moduli=mod,old_mother_shifts=shifts,old_implied_d_over_s=baseline["down_over_s"],
                QCDSF_d_over_s=1/cfg["QCDSF_s_d"][0],effects_order=["d_over_s","ub_cb","up_charm_top_g3"],
                effects_quantities=["Vus","reactor","solar_angle","dm31_eV2"],cases=log,
                average_contrast_over_all_factor_orders=components.tolist(),
                spectrum_fixed_inverse_alpha3_minus_alpha2_at_old_crossing_after_g3_reset=delta)


def native(x):
    if isinstance(x,dict):
        return {k:native(v) for k,v in x.items()}
    if isinstance(x,list):
        return [native(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=native(solve())
    path=HERE/"factorial_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1755 eight refitted Cabibbo holdouts and frozen-threshold gauge residual verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
