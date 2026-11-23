"""Recompute the rare-K norm bound and transport it on the aligned singlet branch."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def annulus(blo, bhi, coefficient, backgrounds):
    """Amplitudes in per-species units where the tree contribution is the complex P."""
    average = sum(backgrounds)/3
    variance = sum(abs(a-average)**2 for a in backgrounds)/3
    if bhi < coefficient*variance:
        return None
    return -average, math.sqrt(max(0, blo/coefficient-variance)), math.sqrt(bhi/coefficient-variance)


def certificate():
    spec = read(HERE/"inputs.json")
    ch, last, history = [read(HERE/spec[k]) for k in ("charged","row_and_projection","row_boxes")]
    neutral = read(HERE/spec["prior_neutral_targets"])
    m = last["row_dsb"]
    s = sum(abs(z)**2 for z in m)
    d = 1-s
    vd, vs = m[0]/math.sqrt(s), m[1]/math.sqrt(s)
    hidden = last["branches"]["proton_null"]["ds_projector"]
    assert math.isclose(hidden, last["branches"]["proton_null"]["Y"]*abs(vd*vs))
    rate, e3, r, other = [spec[k] for k in ("B_rare","B_e3","r_conversion","other_amplitude")]
    blo, bhi = rate["central"]-2*rate["minus"], rate["central"]+2*rate["plus"]
    prior_cap = other["SM_rate_anchor"]+2*other["SM_width"]
    product, form = ch["K_product"], ch["K_form"]
    max_us = (product[0]+2*product[1])/(form[0]-2*form[1])
    min_e3 = e3["central"]-2*(e3["stat"]+e3["sys"])
    c0 = 1.5*e3["central"]*r["central"]/abs(m[1])**2
    cmin = 1.5*min_e3*r["sensitivity"][0]/max_us**2
    cap = (math.sqrt(bhi)+math.sqrt(prior_cap))/math.sqrt(cmin)
    btree_min = cmin*hidden**2
    other_required = (max(0, math.sqrt(btree_min)-math.sqrt(bhi)))**2

    def transported(deficit, prod_lower, d_weight_upper):
        ycap = min(1.,cap/prod_lower)
        assert ycap < deficit < 1
        xlower = (deficit-ycap)/(1-ycap)
        return {"Y_upper":ycap, "X_lower":xlower,
                "proton_shift_upper":-2*xlower+d_weight_upper*ycap,
                "Cs_shift_upper":-188*xlower+211*d_weight_upper*ycap}

    smallest_d = 1-history["bounds"]["1"]["official_mixed"]["row"][1]
    min_d_amp = ch["ud_beta"][0]-ch["ud_beta"][1]
    min_s_amp = (product[0]-product[1])/(form[0]+form[1])
    whole_box = transported(smallest_d,min_d_amp*min_s_amp,1.)
    whole_box["Cs_one_width_lower"] = neutral["data_minus_baseline"][1]-neutral["QCs_work_width"]
    zlower = (2/3)*smallest_d*min_d_amp*min_s_amp

    # Coherent alignment saturates the needed cancellation; a nonuniform
    # background leaves a positive orthogonal component even after the best choice of P.
    a0 = math.sqrt(rate["central"]/c0)
    tree_max_with_fixed_other = a0 + math.sqrt(prior_cap/c0)
    an = annulus(blo,bhi,c0,[1j*a0,-1j*a0,0j])
    assert an is not None and an[0]==0
    assert annulus(blo,bhi,c0,[2j*a0,-2j*a0,0j]) is None
    tuned = a0-hidden
    assert math.isclose(c0*abs(hidden+tuned)**2,rate["central"],rel_tol=1e-12)
    assert math.isclose(c0*(tree_max_with_fixed_other-math.sqrt(prior_cap/c0))**2,rate["central"])
    return {"C_central":c0,"C_conservative_lower":cmin,"two_width_rare_rate":[blo,bhi],
            "other_only_rate_budget":prior_cap, "projection_abs_upper":cap,
            "null_witness":{"P_abs":hidden,"tree_rate_central":c0*hidden**2,
                            "tree_rate_lower":btree_min,"minimum_compensating_rate":other_required,
                            "rate_ratio_to_other_budget":other_required/prior_cap,
                            "exact_central_cancellation_other_P_units":tuned},
            "central_aligned_after_rare":transported(d,abs(vd*vs),vd**2),
            "one_charged_box":{"D_lower":smallest_d,"proton_null_projection_lower":zlower,**whole_box},
            "annulus_test":{"heterogeneous_upper_radius":an[2],"large_orthogonal_background_allowed":False}}


def same(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=5e-12,abs_tol=1e-18) if isinstance(a,float) else a==b


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    out=certificate()
    path=HERE/"results.json"
    if parser.parse_args().write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert same(out,read(path))
        print("1387: rare-K annulus, conditional cancellation and neutral-charge bounds verified")
