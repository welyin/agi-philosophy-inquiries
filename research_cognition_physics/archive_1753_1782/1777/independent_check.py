"""Sequential heavy-current production, prompt/long-lived coverage, and the same Yukawa limit."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE/"../1373/check.py"))
OLD_CFG = OLD["read"]("inputs.json")
OBS = json.loads((HERE/"../1377/results.json").read_bytes())
GF = json.loads((HERE/"../1404/inputs.json").read_bytes())["const_GeV"]["GF"]
HBARC_MM = 1.973269804e-13       # GeV mm
MW, MT, MB = 80.379, 172.5, 4.8   # tree-level width convention, GeV
OLD_FERMIONS = OLD["common_inputs"](OLD_CFG)[2]
UP_MASSES = [f[0] for f in OLD_FERMIONS[:3]]


def width(parent, child):
    x, w = (child/parent)**2, (MW/parent)**2
    assert child+MW < parent
    return GF*parent**3/(8*math.sqrt(2)*math.pi)*math.sqrt((1-x-w)**2-4*x*w)*(
        (1-x)**2+(1+x)*w-2*w**2)


def mass_test(v, known):
    budget = 4*math.pi*v**2/3
    sum_old = math.fsum(m*m for m in known)
    bounds = [math.sqrt(k*budget-sum_old) for k in (1,2)]
    return dict(v=v,ordinary_up_masses=known,up_fourth_bounds_pure_and_common_shift_GeV=bounds,
                differences_from_CMS_685=[top-685 for top in bounds],
                fractions_of_half_circle_at_685=[(685**2+sum_old)/(k*budget) for k in (1,2)],
                A09_approximate_CMS_margin_GeV=755-bounds[1],
                A09_threshold_pole_over_matching_ratio=755/bounds[1],
                A09_nonuniversal_partial_wave_norm_needed_at_755=max(0.,.25*(755**2+sum_old)/budget-.5))


def rare_log(x):
    return math.log(-math.expm1(-x))


def geometry(g, k):
    """Both prompt or both survive long radius, normalized by mean upper-quark length."""
    def derivative_log_ratio(x):
        a = (k*(1+g)-1)*x + rare_log(g*x)
        b = math.log(g) + (k*(1+g)-g)*x + rare_log(x)
        return float(np.logaddexp(a,b))-math.log(k*(1+g))

    lo, hi = 0., 2.
    for _ in range(100):
        x = (lo+hi)/2
        if derivative_log_ratio(x)>0:
            hi=x
        else:
            lo=x
    x=(lo+hi)/2
    prompt=math.expm1(-x)*math.expm1(-g*x)
    long=math.exp(-k*(1+g)*x)
    prob=prompt+long
    assert abs(derivative_log_ratio(x)) < 2e-12
    for z in (x*0.5,x*.99,x*1.01,x*2):
        assert math.expm1(-z)*math.expm1(-g*z)+math.exp(-k*(1+g)*z) >= prob
    return dict(x_min=x,mean_upper_length_in_inner_radii=1/x,
                both_prompt=prompt,both_long=long,ideal_all_lifetimes_floor=prob,
                produced_pairs_95pct=math.ceil(math.log(.05)/math.log1p(-prob)))


def detector_point(mass):
    gu, gd = width(mass,MB), width(mass,MT)
    g=gd/gu
    # In the rest frame Q=2 sqrt(2) m, each particle has beta*gamma = 1.
    rates=[]
    for log_mix in (-6,-12,-14,-16,-18):
        mixing=10.**log_mix
        x=mixing*gu/HBARC_MM
        rates.append(dict(one_minus_A=mixing,mean_U_mm=1/x,mean_D_mm=1/(g*x),
              both_decay_within_1mm=math.expm1(-x)*math.expm1(-g*x),
              W_UD_spectral_weight_relative_massless_lepton=3*(1-mixing)*(7/8)/math.sqrt(2)))
    envelope={str(k):geometry(g,k) for k in [10,100,1000,10000]}
    for item in envelope.values():
        item["one_minus_A_at_least_efficient_point"]=item["x_min"]*HBARC_MM/gu
    mh=125.
    base = complex(*OBS["SM_quark_loop_real_imag"])
    amp = OLD["loop_functions"](mh,mass)[1]
    hratio=abs(base+2*amp)**2/abs(base)**2
    return dict(mass_GeV=mass,off_shell_invariant_mass_GeV=2*math.sqrt(2)*mass,
                unmixed_widths_U_D_GeV=[gu,gd],width_ratio_D_over_U=g,
                one_mm_mean_criteria_one_minus_A=[HBARC_MM/gu,HBARC_MM/gd],
                rates=rates,ideal_prompt_or_long=envelope,
                Higgs_ggF_ratio_first_nonzero_order=hratio)


def calculate():
    cases={}
    for mass in (360.,420.,600.,685.):
        cases[str(int(mass))]=detector_point(mass)
    same=geometry(1.,10000)
    return dict(tree_partial_wave_examples={
                "archived_MZ_parameters_used_as_constant":mass_test(OLD_CFG["v_MZ_GeV"],UP_MASSES),
                "on_shell_top_no_uc_conservative_tree_comparison":mass_test((math.sqrt(2)*GF)**-.5,[MT])},
                stable_limit_angular_weight=3*(7/8)/math.sqrt(2), mass_points=cases,
                equal_width_example_10m_over_1mm=same,
                sequential_quark_Higgs_lower=OBS["min_gg_plus_bb_over_vbf_ratio"],
                atlas_two_marginal_error_upper=OBS["marginal_rectangles"]["2"]["upper_ratio"])


def clean(obj):
    if isinstance(obj,dict):
        return {k:clean(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):
        return [clean(x) for x in obj]
    return float(f"{obj:.10g}") if isinstance(obj,float) else obj


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=clean(calculate())
    payload=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
    else:
        path.write_text(payload,encoding="utf-8")
    print("1777 fourth-quark production and finite-lifetime coverage",hashlib.sha256(payload.encode()).hexdigest()[:12])
