"""All-scalar-partner alternative for the 5-colour plus 4-weak gauge thresholds."""
import argparse
import json
import itertools
import math
from pathlib import Path
import numpy as np


HERE = Path(__file__).resolve().parent
SM = np.array([41/10, -19/6, -7.])
MSSM = np.array([33/5, 1., -3.])
COLOR = np.array([2/5, 0., 1.])
WEAK = 1-COLOR
ADJ = np.array([0., 2., 3.])


def read(file):
    return json.loads((HERE/file).read_bytes())


def totals(data, old, bare, newleptons, oldleptons):
    mz, soft = data[3], old["soft_GeV"]
    measured = 4*np.pi/np.array(data[:3])**2
    a = measured - ADJ/(12*np.pi) - (SM-MSSM)*np.log(soft/mz)/(2*np.pi)
    projection = np.array([5, -3, -2])
    length = 2*np.pi*(projection@a)/(projection@MSSM)
    upper = mz*np.exp(length)
    target = 2*np.pi*(a[1]-a[2])-(MSSM[1]-MSSM[2])*length
    def effective_log(m):
        return (2*np.log(m)+np.log(np.hypot(m,soft)))/3
    old_colors = np.array(old["colored_fermions_GeV"])
    delta = np.sum(effective_log(old_colors))-np.sum(effective_log(oldleptons))
    value = (target-delta+np.log(upper)+np.sum(effective_log(newleptons))-2*effective_log(bare))/3
    low, high = np.log(bare), np.log(upper)
    assert effective_log(np.exp(low))<value<effective_log(np.exp(high))
    for _ in range(80):
        mid = (low+high)/2
        if effective_log(np.exp(mid)) < value:
            low = mid
        else:
            high = mid
    third = math.exp((low+high)/2)
    weak_all = np.concatenate((newleptons,oldleptons))
    inverse = a[1]-(MSSM[1]*length+sum(np.log(upper)-effective_log(weak_all)))/(2*np.pi)
    events = [(soft,MSSM-SM,-ADJ/(12*np.pi))]
    for mats, coeff in [(np.r_[old_colors, [bare]*2, [third]*3], COLOR),(weak_all, WEAK)]:
        for m in mats:
            events.append((m,2*coeff/3,0))
            events.append((np.hypot(m,soft),coeff/3,0))
    ai, b, last = measured.copy(), SM.copy(), mz
    min_inverse = np.min(ai)
    for threshold, added, scheme in sorted(events,key=lambda row:row[0]):
        assert mz <= threshold <= upper
        ai -= b*np.log(threshold/last)/(2*np.pi)
        min_inverse = min(min_inverse, np.min(ai))
        ai += scheme
        min_inverse = min(min_inverse, np.min(ai))
        b += added
        last = threshold
    ai -= b*np.log(upper/last)/(2*np.pi)
    np.testing.assert_allclose(ai, np.full(3,inverse),atol=1e-10)
    return dict(R=third/bare, three_upper_color_poles_GeV=third, unified_scale_GeV=upper,
                inverse_alpha_U=inverse, minimum_inverse_alpha_along_RGE=min(min_inverse,ai.min()))


def full_soft(kappa_ratio, yb, K, bare, soft):
    # r=soft squared/heavy squared; Q=soft. All eight ordinary and added squark soft terms equal.
    # EFT ordinary down contribution is subtracted: both scenarios use the same light Yukawas.
    def f(r):
        return (1+r)*math.log1p(r)-r-r*math.log(r) if r else 0.
    R = kappa_ratio
    g2 = yb**2/(1-1/R**2)
    a = (soft/bare)**2
    return 3*bare**2*g2/(8*math.pi**2)*((3-K)*f(a)+K*f(a/R**2))


def supertrace_check(K, bare, soft):
    spectrum=read("../1713/results.json")["tagged_yukawa_values_desc"]
    yb,ys,yd=spectrum
    s=np.array([-yd/yb,ys/yb,-1.])
    def f(vals, scale):
        vals = vals[vals>0]
        return np.sum(vals**2*(np.log(vals/scale)-1.5))
    for R in [1.1,1.4,2.8]:
        u=math.sqrt(R*R-1)
        a=(soft/bare)**2
        g2=yb**2/(1-1/R**2)
        def pot(z):
            poles=[]
            for k in s:
                w=np.sqrt((1+np.array([1,-1])*k)/2)
                used=w>0
                v=w[used]
                mat=np.zeros((len(v)+1,len(v)+1))
                mat[0,1:]=np.sqrt(z)*v
                mat[1:,0]=-u*v
                mat[1:,1:]=np.diag(np.array([1,-1])[used])
                poles.extend(np.linalg.svd(mat,compute_uv=False)**2)
            eig=np.array(poles)
            light=s*s*z*u**2/R**2
            return f(eig+a,a)-f(eig,a)-f(light+a,a)+f(light,a)
        delta=2e-4
        fd=(4*pot(delta)-pot(2*delta)-3*pot(0))/(2*delta)*3*bare**2*g2/(16*np.pi**2)
        assert abs(fd/full_soft(R,yb,K,bare,soft)-1)<2e-6


def compute():
    gauge=read("../1369/inputs.json")
    old=read("../1698/results.json")["shared_witness"]
    mz=gauge["mu_ref_GeV"]
    data=np.array([gauge[k] for k in ("g1_GUT","g2","g3")])
    cen=[*data[:,0],mz]
    extremes=[[*row,mz] for row in itertools.product(*(data[:,0,None]+2*data[:,1,None]*np.array([-1,1])))]
    bare=read("../1713/inputs.json")["mediator_matching_GeV"]
    new_weak=np.array(read("../1714/results.json")["same_1712_point"]["heavy_masses_GeV"])
    low_weak=np.full(4,1.e7)
    old_weak=np.array(old["weak_fermions_GeV"])
    physical=totals(cen,old,bare,new_weak,old_weak)
    ranges={}
    for name, pair in {"old_weak_fixed":old_weak,"old_weak_any_above_MZ":np.full(2,mz)}.items():
        rows=[totals(c,old,bare,low_weak,pair) for c in extremes]
        ranges[name]={"two_width_minimum_necessary_R":min(x["R"] for x in rows),
                      "two_width_minimum_upper_color_pole_GeV":min(x["three_upper_color_poles_GeV"] for x in rows),
                      "two_width_minimum_inverse_alpha_U":min(x["inverse_alpha_U"] for x in rows)}
    R=physical["R"]
    rad=np.sqrt(R*R-1)
    delta=1/(rad*(R+rad))
    yb,ys,yd=read("../1713/results.json")["tagged_yukawa_values_desc"]
    K=1+(ys/yb)**2+(yd/yb)**2
    first=read("../1719/results.json")
    budget=(first["target_central_and_two_width_abs_upper_GeV2"][1]+
            first["one_parameter_neutral_rescaling"]["maximum_positive_neutral_c"])
    pref=3*bare**2*yb**2/(8*np.pi**2)
    cap=np.exp((budget/pref-1.5)/(2*K))
    assert ranges["old_weak_any_above_MZ"]["two_width_minimum_necessary_R"]>cap
    soft=old["soft_GeV"]
    supertrace_check(K,bare,soft)
    active=full_soft(R,yb,K,bare,soft)
    unshifted=read("../1698/independent_results.json")["frozen_30TeV_universal_soft"]["soft_scale"]
    retune=unshifted["required_other_and_tree_coefficient_GeV2"]-active
    return {"actual_R1714_weak_spectrum_central_all_superpartners":physical,
            "central_weak_lower_bound_all_new_equal_10PeV":totals(cen,old,bare,low_weak,old_weak),
            "full_marginal_boxes_all_new_weak_above_10PeV":ranges,
            "central_g_over_yb_minus_one":delta,
            "M_times_kappa_GeV":bare*rad,
            "hypothetical_R_cap_if_colored_pure_fermion_curvature_imposed":cap,
            "same_colored_spectrum_yukawa_soft_subsector_at_30TeV":{
                "added_threshold_c_GeV2_after_subtracting_ordinary_down_squarks":active,
                "original_required_other_local_GeV2":unshifted["required_other_and_tree_coefficient_GeV2"],
                "required_other_local_GeV2":retune,
                "pure_colored_fermion_at_R1719_Q0_GeV2":-pref*(1.5+2*K*np.log(R))/(1-1/R**2)}}


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    return float(f"{x:.11g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read("independent_results.json")
        print("1720 independent complete multiplet running and color supertrace verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

