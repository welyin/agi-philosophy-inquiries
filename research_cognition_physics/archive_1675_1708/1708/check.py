"""Verify the same weak-scale seed against the quark and lepton sectors."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
def read(path): return json.loads((HERE/path).read_bytes())
cfg=read("inputs.json")
sim=runpy.run_path(str(HERE/cfg["neutrino"]))
full=sim["fl"]
initial=read(cfg["charged_source"])["direct_semileptonic_moduli"]
clock=runpy.run_path(str(HERE/cfg["independent_clocks"]))
measure=read(cfg["clock_measurements"])
qb=read(cfg["quark_comparisons"])


def quarks(gamma,angles):
    ckm=full["low_ckm"](gamma,angles)
    beta=np.angle(-ckm[1,0]*ckm[1,2].conjugate()/(ckm[2,0]*ckm[2,2].conjugate()))
    bs=np.angle(-ckm[2,1]*ckm[2,2].conjugate()/(ckm[1,1]*ckm[1,2].conjugate()))
    return {"beta_deg":float(np.degrees(beta)),"golden_S":float(np.sin(2*beta)),
            "phi_s_rad":float(-2*bs),
            "Vtd_over_Vts":float(abs(ckm[2,0]/ckm[2,1])),"gamma_deg":float(gamma)}


def trial(vub,dz,k,role,step=.20):
    value=initial["Vub"][0]+vub*np.linalg.norm(initial["Vub"][1:])
    # Only local in-memory inputs are updated. Reused high and low runners
    # both read this very dictionary when they propagate the three CKM angles.
    full["cfg"]["direct_semileptonic_moduli"]["Vub"][0]=float(value)
    angles=full["angles"]("direct")
    assert abs(np.sin(angles[1])-value)<1e-17
    joint=sim["solve_case"](k,dz,role,step)
    return {"Vub":value,"Vub_widths":vub,"yd_widths":dz,"k":k,"port":role,
            "quark":quarks(joint["charged_gamma_MZ_deg"],angles),
            "neutral_Higgs":joint["same_Higgs_calibration"],"third_modulus":joint["third_modulus"],
            "neutrino_roots":joint["roots"]}


def compute():
    cases={f"d{dz}_k{k}_{p}":trial(-1.,dz,k,p) for dz in cfg["down_Yukawa_displacements"]
           for k in cfg["phases_k"] for p in cfg["third_ports"]}
    golden=measure["one_width_golden_S_range"]
    qr=measure["one_width_clock_ratio_range"]
    alt=measure["alternative_FNAL_MILC_2plus1_with_omitted_charm_budget"]["one_width_ratio"]
    ga=qb["time_integrated_gamma_deg"]["68"]
    def within(x,r): return bool(r[0]<=x<=r[1])
    labels=[]
    for key,case in cases.items():
        q=case["quark"]
        tests={"time_integrated_gamma":within(q["gamma_deg"],ga),
               "Bd_golden":within(q["golden_S"],golden),
               "Bd_Bs_no_extra_residual":within(q["Vtd_over_Vts"],qr)}
        case["quark_tests"]=tests
        case["other_clock_conventions"]={
            "one_percent_independent_mixing_corrections":
                within(q["Vtd_over_Vts"],[qr[0]*math.sqrt(.99/1.01),qr[1]*math.sqrt(1.01/.99)]),
            "alternative_two_plus_one_lattice":within(q["Vtd_over_Vts"],alt)}
        for src,v in qb["phi_s_rad"].items():
            case["other_clock_conventions"][src+"_marginal_displacement"]=(q["phi_s_rad"]-v[0])/math.hypot(*v[1:])
        width=sim["fl"]["physics"]["bisect"](lambda t:clock["ratio_at"](t)[0]-q["Vtd_over_Vts"],-4.,4.)
        case["lattice_clock_witness"]={"common_signed_width_fraction":width,**clock["ratio_at"](width)[1]}
        for i,nu in enumerate(case["neutrino_roots"]):
            fc=nu["accelerator"]
            joint=all(tests.values()) and all(nu["windows_pass"]) and (fc["conf90_0"] or fc["conf90_1"])
            nu["all_registered_comparisons"]=bool(joint)
            if joint: labels.append([key,i])
    # This kinematic inversion uses only the two measured clocks and lattice ratio,
    # not the GJ model, charged shooting, lepton angles or CP comparison.
    th=full["angles"]("direct")
    gamma_gate=[sim["fl"]["physics"]["bisect"](
        lambda g:quarks(g,th)["Vtd_over_Vts"]-r,10.,160.) for r in qr]
    testscores=[]
    for tag in ["unit_third","shared_third"]:
        def hadronic(d):
            # Eight-argument shooting uses the same neutral normalization as each case.
            return sim["graph"](1,d,tag,.20)[1]["gamma_MZ_deg"]-gamma_gate[0]
        l,h=-1.,0.
        assert hadronic(l)>0 and hadronic(h)<0
        # Bisection locates the local crossing only; global monotonicity is not assumed.
        for _ in range(14):
            d=(l+h)/2
            if hadronic(d)>0: l=d
            else: h=d
        testscores.append({"port":tag,"critical_d_local_sign_bracket":[l,h]})
    assert labels==[["d-1_k1_unit_third",1],["d-1_k1_shared_third",1]]
    ratio=cases["d-1_k1_unit_third"]["Vub"]/initial["Vcb"][0]
    checks={key:(ratio-mode)/(math.hypot(*errs_high) if ratio>=mode else math.hypot(*errs_low))
            for key,(mode,errs_low,errs_high) in qb["ratio_ports"].items()}
    return {"predeclared_endpoints":cases,"joint_pass_labels":labels,
            "fixed_tree_ratio":ratio,"additional_tree_ratio_marginal_displacements":checks,
            "quark_only_gamma_gate_deg":gamma_gate,"quark_only_yd_crossings":testscores}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    out=sim["tidy"](compute())
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==out
        print("1708 verified")
    else:
        path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
