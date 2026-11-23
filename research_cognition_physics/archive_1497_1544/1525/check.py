"""Charged Higgs tests with independent K and B semileptonic standards."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def box(row, width):
    error = width * sum(row[1:])
    return row[0] - error, row[0] + error


def extent(values):
    values = list(values)
    return [min(values), max(values)]


def ratio_enclosure(num, den):
    assert min(den) > 0 and min(num) >= 0
    return [num[0] / den[1], num[1] / den[0]]


def positive_ratio(v, w, u):
    return [v[0]*w[0]/u[1], v[1]*w[1]/u[0]]


def b_normalization(cfg, b, width):
    tau = b["charged_tau_control"]
    choices = [box(tau[k], width) for k in (
        "f_Bplus_GeV", "tau_Bplus_ps", "mass_Bplus_GeV")]
    mt = tau["mass_tau_GeV_window"] if width else [1.77686]
    choices += [mt, box(b["B0_lifetime_ps"], width),
                box(b["HPQCD_partial_massless_per_ps"], width)]
    def value(f, tp, mb, ml, t0, rat):
        width_GeV = cfg["GF_GeV_m2"]**2*f*f*mb*ml**2*(1-(ml/mb)**2)**2/(8*math.pi)
        return tp*width_GeV/cfg["hbar_GeV_ps"]/(t0*rat)
    k = extent(value(*trial) for trial in itertools.product(*choices))
    electron = np.sum([box(row, width) for row in b["branch_e_high_bins_1e4"]], axis=0)
    tau_counts = box(tau["branch_tau_1e4"], width)
    return k, ratio_enclosure(tau_counts, electron)


def kaon_roots(r_range, eta, ak, ap, xmax, masses):
    mk, mp = masses
    phi = (1-mp**2*xmax)/(1-mk**2*xmax)
    trange = [r_range[0]/(eta[1]*phi), r_range[1]/eta[0]]
    answer = {}
    for sign in (1, -1):
        vals = []
        for t, k, p in itertools.product(trange, ak, ap):
            assert k - sign*t*p > 0
            vals.append((1-sign*t)/(k-sign*t*p))
        hlow, hhigh = extent(vals)
        answer["normal" if sign == 1 else "mirror"] = [
            math.sqrt(max(0, hlow-xmax)), math.sqrt(max(0, hhigh))]
    return answer


def amplitude(mass, quark_ratio, x, y):
    return (1-mass**2*(y+x)/(1+quark_ratio))/(1-mass**2*x)


def calculate():
    cfg = read("inputs.json")
    b = read(cfg["refs"]["B_inputs"])
    quartic = read(cfg["refs"]["quartic"])
    g = read("../1523/inputs.json")["sm_dr_at_one_TeV"]
    dterm = (0.6*g["g1"][0]**2 + g["g2"][0]**2)/8
    mk, mp = cfg["masses_GeV"]["K"], cfg["masses_GeV"]["pi"]
    ratios = cfg["common_scheme_quark_ratio_boxes"]
    rd, rs, rb = [ratios[key] for key in ("u_over_d", "u_over_s", "u_over_b")]
    ak, ap = (extent(m*m/(1+r) for r in rr) for m,rr in [(mk,rs),(mp,rd)])
    xmax = cfg["minimum_charged_mass_GeV"]**-2
    ek, eb = cfg["assumed_transfers"].values()
    old = cfg["KLOE"]
    tau = b["charged_tau_control"]
    window = {}
    for width in (1, 2):
        published = box(old["published_R_with_unitarity"], width)
        separate = positive_ratio(*[box(old[k],width) for k in (
            "separate_Vus_over_Vud_mu2","nuclear_Vud","separate_Vus_kl3")])
        k, br = b_normalization(cfg, b, width)
        max_amp = math.sqrt(br[1]/(k[0]*eb[0]))
        min_amp = math.sqrt(br[0]/(k[1]*eb[1]))
        mb = box(tau["mass_Bplus_GeV"],width)
        qmax = math.sqrt((1+max_amp)*(1+rb[1])/mb[0]**2)
        # For each signed B root solve q²=(1+r)(1-F)/m²+x*((1+r)*F-1).
        # It is separately monotone in each variable at fixed values of the others.
        b_roots = {}
        for sign, label in [(1,"positive"),(-1,"negative")]:
            yb = extent((1+r)*(1-sign*f)/m**2+x*((1+r)*sign*f-1)
                        for f,r,m,x in itertools.product((min_amp,max_amp),rb,mb,(0,xmax)))
            b_roots[label] = [math.sqrt(max(0, bound)) for bound in yb]
        kres = {name:kaon_roots(r, ek, ak, ap, xmax, (mk,mp)) for name,r in
                [("KLOE_fit",published),("separate",separate)]}
        assert all(pair["mirror"][0] > qmax for pair in kres.values())
        historical = {}
        for name, estimate in cfg["BaBar_q_per_GeV"].items():
            if name == "source":
                continue
            q = box(estimate,width)
            qsquared = [v*v for v in q]
            # Factor the known common x from the kaon ratio, enclose h and phi separately.
            hrange = [qsquared[0], qsquared[1]+xmax]
            bare = extent((1-aa*h)/(1-bb*h) for aa,bb,h in
                          itertools.product(ak,ap,hrange))
            assert 1-max(ak)*max(hrange) > 0
            phi_max = (1-mp**2*xmax)/(1-mk**2*xmax)
            kpred = [bare[0],bare[1]*phi_max]
            bamps = extent(abs(amplitude(m, r, x, y)) for m,r,x,y in
                           itertools.product(mb,rb,(0,xmax),qsquared))
            assert amplitude(mb[0],rb[1],0,qsquared[0]) < 0
            historical[name] = dict(q=q, kaon_bare=kpred, B_abs_amp=bamps,
                min_K_transfer_free_CKM=separate[0]/kpred[1],
                largest_B_rate_transfer_consistent=br[1]/(k[0]*bamps[0]**2))
        window[str(width)] = dict(K_published_R=published, K_separate_R=separate,
            q_from_K=kres, B_ratio=br, K_B=k, B_amplitude_range=[min_amp,max_amp],
            q_from_B=b_roots, q_upper_B=qmax, historical=historical)
    k0, _ = b_normalization(cfg,b,0)
    wk = window["2"]
    assert wk["q_from_K"]["separate"]["normal"][1] < .4
    assert wk["historical"]["Dstar"]["largest_B_rate_transfer_consistent"] < eb[0]
    witnesses = []
    for tan, mass in [(1,1000),(20,1000),(20,20000)]:
        x,y = mass**-2,(tan/mass)**2
        fK = amplitude(mk,sum(rs)/2,x,y)/amplitude(mp,sum(rd)/2,x,y)
        fB = amplitude(tau["mass_Bplus_GeV"][0],sum(rb)/2,x,y)
        ratioB = k0[0]*fB**2
        tree = dterm*((1-tan*tan)/(1+tan*tan))**2
        stop_req = quartic["quartic"]-tree
        coef = quartic["coefficient"]
        stop_u = 6*(1-math.sqrt(1-stop_req/(3*coef))) if 0<=stop_req<=3*coef else None
        witnesses.append(dict(tan_beta=tan, charged_mass_GeV=mass, q=math.sqrt(y),
            R_K=fK, B_ratio=ratioB, B_tau_1e4=.38*ratioB,
            lambda_tree_D=tree, stop_needed=stop_req, stop_X2_over_MS2=stop_u))
        assert wk["K_separate_R"][0] <= fK <= wk["K_separate_R"][1]
        assert wk["B_ratio"][0] <= ratioB <= wk["B_ratio"][1]
        if stop_u is not None:
            assert abs(tree+coef*stop_u*(1-stop_u/12)-quartic["quartic"])<1e-12
    return dict(widths=window, witnesses=witnesses,
        K_finite_propagator_ratio_max=(1-mp**2*xmax)/(1-mk**2*xmax),
        D_term_unit=dterm, stop_positive_max=3*quartic["coefficient"])


def check(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            check(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for l,r in zip(a,b):
            check(l,r)
    elif a is None:
        assert b is None
    else:
        np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-14)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    flags = parser.parse_args()
    result = calculate()
    if flags.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        check(result,read("results.json"))
        print("1525 OK: both CKM schemes, signed K roots, B normalization, shared threshold")
