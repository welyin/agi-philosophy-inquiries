"""Joint leading gauge and flavour matching for one vectorlike 15 pair split by an adjoint."""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PI = math.pi


def read(relative):
    return json.loads((HERE/relative).read_bytes())


def zero(fun, low, high):
    assert fun(low)<0<fun(high)
    for _ in range(80):
        point = (low+high)/2
        if fun(point)<0:
            low = point
        else:
            high = point
    return (low+high)/2


def group_coefficients():
    # Columns T=(1,3,1), S=(6,1,-2/3), Z=(3,2,1/6); vectorlike pair of chiral superfields.
    reps = [(1,3,1,0,2), (6,1,-2/3,2.5,0), (3,2,1/6,.5,.5)]
    B = np.array([[2*3/5*c*w*y*y,2*c*iw,2*w*ic] for c,w,y,ic,iw in reps]).T
    differential = np.array([3.,-4.,1.])
    np.testing.assert_allclose(B.sum(axis=1),7.)
    np.testing.assert_allclose(-5*(B[0]-B[1])/54 + 20*(B[1]-B[2])/27,differential)
    # Check the 5bar component multiplicities with a normalized symmetric tensor: diagonal
    # term has squared vertex 2 for the same species, off-diagonal squared vertex 1.
    w = np.ones((5,5))+np.eye(5)
    tags = np.array([0,0,0,1,1])  # colour first, then weak
    weights = []
    for a in (0,4):
        weights.append([sum(w[a,b] for b in range(5) if tags[a]+tags[b]==tag) for tag in (2,0,1)])
    np.testing.assert_allclose(np.subtract(weights[1],weights[0]),differential)
    return B


def one_gauge_box(steps, ew, lattice, top, a3_base):
    d1,d2,da,dt = steps
    def inverse_shift(key,index):
        g,dg=ew[key]
        return 4*PI/(g+index*dg)**2-4*PI/g**2
    a1,a2 = inverse_shift("g1_GUT",d1),inverse_shift("g2",d2)
    a3 = 1/(lattice[0]+da*lattice[1])+math.log((top[0]+dt*top[1])/ew["mu_ref_GeV"])/(3*PI)-a3_base
    g12,g23 = 2*PI*(a1-a2),2*PI*(a2-a3)
    u0 = -(14*g23+5*g12)/54
    v0 = (10*g12+g23)/27
    eta = 10*u0+7*v0
    c0 = -3*u0-v0
    k=2*math.exp(eta/7)
    # For r<1, r**(10/7)(1-r) reaches its maximum at r=10/17.
    assert k>(10/17)**(10/7)*(7/17)
    r_lo = zero(lambda r:r**(10/7)*(1+r)-k,0.,5.)
    r_hi = zero(lambda r:r**(10/7)*(r-1)-k,1.,5.)
    ends=[(-eta-11*math.log(r))/7 for r in (r_hi,r_lo)]
    assert ends[0] < ends[1]
    ell=[(math.log(r)-u0)*9/14 for r in (r_lo,r_hi)]
    assert ell[0]<ell[1]
    return np.array([u0,v0,1+math.exp(u0)-2*math.exp(v0),eta,*ends,*ell,c0])


def gauge_sample(name, b, ew, conf, baseline, bottom, B):
    a3 = baseline["a3"]
    au, base = baseline["inverse_alpha_U"], baseline["U"]
    observed = conf["strong_lattice_5f"]["ALPHA_2026"]
    top = conf["top_MS_matching_GeV"]
    sample={}
    for n in (0,2):
        points=np.array([one_gauge_box(ind,ew,observed,top,a3)
                         for ind in itertools.product((-n,n),repeat=4)])
        lo,hi=points.min(axis=0),points.max(axis=0)
        magnitude=max(abs(lo[4]),abs(hi[5]))
        fixed_margin_upper=1+math.exp(hi[0])-2*math.exp(lo[1])
        sample[str(n)] = {
            "fixed_scale_min_max_complex_mass_sum_margin_among_corners": [lo[2],hi[2]],
            "fixed_scale_entire_box_margin_upper":fixed_margin_upper,
            "gauge_invariant_eta_interval":[lo[3],hi[3]],
            "floating_scale_leading_flavour_log_range":[lo[4],hi[5]],
            "max_relative_leg_log_if_Ytop_squared_le_4pi":magnitude/(4*PI),
            "max_frozen_generator_gain_minus_one":math.expm1(magnitude/(4*PI)),
            "Ytop_at_least_for_b_only_goals":[math.sqrt(16*PI**2*math.log1p(v)/magnitude) for v in bottom]}
        if n==0:
            # At the first allowable crossing all masses can be real and aligned. Fix one
            # common mass for an elementary constructive gauge example (this is a free scale).
            h=lo[6]
            r,q=math.exp(lo[0]+14*h/9),math.exp(lo[1]-20*h/9)
            assert abs(2*q-1-r)<2e-13
            uu=base*math.exp(h)
            masses=1e13*np.array([r,1.,q])
            ell=np.log(uu/masses)
            a_new=au-(b[0]*h+B[0]@ell)/(2*PI)
            residual = (a_new-au)+(b*h+B@ell)/(2*PI)
            target3=1/observed[0]+math.log(top[0]/ew["mu_ref_GeV"])/(3*PI)-a3
            np.testing.assert_allclose(residual,[0,0,target3],atol=7e-13)
            assert a_new>1 and min(masses)>1e7
            sample[str(n)].update({
                "fixed_parent_required_MT_over_MS_MZ_over_MS":[math.exp(lo[0]),math.exp(lo[1])],
                "unification_scale_relative_allowed":[math.exp(lo[6]),math.exp(lo[7])],
                "lowest_scale_positive_mass_witness_GeV_T_S_Z":masses.tolist(),
                "witness_U_GeV_and_inverse_alpha_U":[uu,a_new]})
    return {"baseline":name,"source_U_GeV":base,"source_inverse_alpha_U":au,
            "degenerate_15_floor_if_alphaU_le_one_and_frozen_crossing_GeV":base*math.exp(-2*PI*(au-1)/7),
            "experimental_working_boxes":sample}


def run():
    B = group_coefficients()
    main=read("results.json")
    histories=[v["stop"] for v in main["two_pure_II_histories"].values()]
    bottom=[p["fixed_direction_delta_b_lower"] for p in histories]
    old=read("../1698/results.json")["shared_witness"]
    new=read("../1769/results.json")["at_225_floor_relax_high_pair"]
    ew=read("../1369/inputs.json")
    conf=read("../1769/inputs.json")
    samples={}
    for tag,n,point in [
        ("1698_two_pairs",2,dict(U=old["U_GeV"],inverse_alpha_U=old["inverse_alpha_U"],a3=4*PI/ew["g3"][0]**2)),
        ("1769_four_pairs",4,dict(U=new["crossing_GeV"],inverse_alpha_U=new["inverse_alpha_U"],a3=new["a3_SM"]))]:
        samples[tag]=gauge_sample(tag,np.array([33/5,1.,-3.])+n,ew,conf,point,bottom,B)
    flav=[]
    for p in histories:
        probabilities=np.array(p["full_PMNS_moduli"])**2
        mass_squared=np.array(p["neutral_singular_ratios"])**2
        fractions=probabilities@mass_squared/mass_squared.sum()
        assert abs(sum(fractions)-1)<8e-10
        flav.append(dict(normalized_leptonic_triplet_single_tag_fractions=fractions.tolist(),
                         tau_minus_e_diag_over_top_F_squared=(probabilities[2]-probabilities[0])@mass_squared))
    return dict(step_beta_rows_123_columns_T_S_Z=B.tolist(),threshold_cases=samples,
                isolated_bottom_boost_targets=bottom,pure_typeII_parent_tagged_decay_weights=flav)


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):
        return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    answer=clean(run())
    raw=json.dumps(answer,indent=2,allow_nan=False)+"\n"
    target=HERE/"independent_results.json"
    if args.check:
        assert json.loads(target.read_bytes())==answer
    else:
        target.write_text(raw,encoding="utf8")
    print("1800 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])
