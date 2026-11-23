"""Independent 1536 audit: scalar CP and finite-mass ggF/VBF at an actual ATLAS port."""
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent
ARCH=ROOT.parent
sys.dont_write_bytecode=True


def get(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def kernels(t):
    if t>=1:
        f=math.asin(1/math.sqrt(t))**2
    else:
        f=-(math.acosh(1/math.sqrt(t))-1j*math.pi/2)**2
    return 2*t*(1+(1-t)*f), 2*t*f


def calculate():
    config=get(ROOT/"inputs.json")
    parent=get(ARCH/"1533"/"inputs.json")
    model=get(ROOT/"results.json")
    soft_module=importlib.util.spec_from_file_location("soft",ARCH/"1533"/"check.py")
    soft=importlib.util.module_from_spec(soft_module)
    soft_module.loader.exec_module(soft)
    t=config["shared_odd_witness_ratio"]
    a=1/math.hypot(1,t)
    h=(a,0,-t*a)
    potential=soft.soft_potential(h)
    assert abs(potential["lambda1_to_4"][0]-potential["lambda1_to_4"][1])<1e-14
    beta_tau2=1-4*(parent["tau_mass_GeV"]/parent["higgs_mass_GeV"])**2
    assert math.isclose(model["full_witness"]["rates"]["tau"],1+t*t/beta_tau2)
    assert model["types"]["X"]["protected"] and not model["types"]["II"]["protected"]
    assert model["rate_boxes"]["2"]["full_vertex_angle_outer_deg"]<13
    historic=ARCH/"1373"/"check.py"
    spec=importlib.util.spec_from_file_location("mass_inputs",historic)
    prior=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    _,_,matter=prior.common_inputs(prior.read("inputs.json"))
    mh=parent["higgs_mass_GeV"]
    scalar=0j
    pseudo=0j
    for mass,ncolor,charge in matter:
        if ncolor==3:
            sc,ps=kernels((2*mass/mh)**2)
            scalar+=sc
            pseudo+=ps*(1 if charge>0 else -1)
            assert abs(sc-prior.loop_functions(mh,mass)[1])<1e-10
    pg=abs(pseudo/scalar)**2
    raw=get(ARCH/"1377"/"inputs.json")["atlas_four_leptons"]
    gg,vbf=raw["gg_plus_bb"],raw["vbf"]
    norm=raw["central_normalization"]
    fb=norm["inclusive_bbh_pb"]*norm["br_zz"]/norm["stage_gg_plus_bb_times_br_pb"]
    # The b-associated part contributes at least its scalar component at the CP-symmetric LO port.
    slope=(1-fb)*pg
    bounds={}
    for n in (1,2):
        upper=(gg["mu"]+n*math.hypot(*gg["width_plus"]))/(vbf["mu"]-n*math.hypot(*vbf["width_minus"]))
        bounds[str(n)]={}
        for transfer in (1,.9):
            max_t=math.sqrt(max(0,upper/transfer-1)/slope)
            bounds[str(n)][str(transfer)]={"abs_angle_cap_deg":math.degrees(math.atan(max_t)),
                                           "ratio_upper":upper}
    assert pg>2 and slope>2
    return {"witness_h":h,"charged_GeV":potential["charged_mass_GeV"],
            "equal_lambdas":potential["lambda1_to_4"][:2],
            "scalar_form_real_imag":[scalar.real,scalar.imag],
            "odd_form_real_imag":[pseudo.real,pseudo.imag],"pure_g_slope":pg,
            "bb_fraction_cap":fb,"inclusive_slope_floor":slope,
            "witness_gg_VBF_floor":1+slope*t*t,"ATLAS_working_bounds":bounds}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=calculate()
    path=ROOT/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        def equal(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
            return np.allclose(a,b,rtol=1e-10,atol=1e-13)
        assert equal(result,get(path))
        print("1536: compatible scalar potential, full-mass loops and published ATLAS windows verified")
