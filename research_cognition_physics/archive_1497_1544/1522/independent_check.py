"""A physical stop-pair readout on the one-loop gauge-threshold flat direction."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def run_sm(scale, steps):
    inp = read("../1369/inputs.json")
    quartic = read("../1396/inputs.json")
    z = np.array([inp[key][0] for key in ("g1_GUT","g2","g3","yt","yb","ytau")]
                 + [quartic["lambda_paper"][0]/quartic["lambda_paper_factor"]])
    h = math.log(scale/inp["mu_ref_GeV"])/steps
    def deriv(y):
        g1,g2,g3,t,b,e,lam=y
        gy2=3*g1*g1/5
        tr=3*t*t+3*b*b+e*e
        return np.array([4.1*g1**3,-19*g2**3/6,-7*g3**3,
                         t*(tr+1.5*(t*t-b*b)-17*gy2/12-9*g2*g2/4-8*g3*g3),
                         b*(tr+1.5*(b*b-t*t)-5*gy2/12-9*g2*g2/4-8*g3*g3),
                         e*(tr+1.5*e*e-15*gy2/4-9*g2*g2/4),
                         24*lam*lam+(4*tr-9*g2*g2-3*gy2)*lam
                         -2*(3*t**4+3*b**4+e**4)+3/8*(2*g2**4+(g2*g2+gy2)**2)])/LOOP
    for _ in range(steps):
        k1=deriv(z)
        k2=deriv(z+h*k1/2)
        k3=deriv(z+h*k2/2)
        k4=deriv(z+h*k3)
        z+=h*(k1+2*k2+2*k3+k4)/6
    return z


def scalar_loop(tau):
    assert 0<tau<0.1
    # Positive Feynman-parameter moments; no subtractive cancellation.
    term=1/3
    answer=term
    for k in range(100):
        term*=4*tau*(k+2)**2/((2*k+6)*(2*k+5))
        answer+=term
        if term<1e-20:
            break
    return answer


def amplitudes(mass, top, higgs, x2):
    shift=top*mass*math.sqrt(x2)
    common=mass*mass+top*top
    tau_top=(higgs/(2*top))**2
    fermion=2*(tau_top+(tau_top-1)*math.asin(math.sqrt(tau_top))**2)/tau_top**2
    ratios=[]
    masses=[]
    for sign in (-1,1):
        eigen=common+sign*shift
        assert eigen>higgs*higgs/4
        tau=higgs*higgs/(4*eigen)
        ratios.append((top*top+sign*shift/2)/eigen*scalar_loop(tau)/fermion)
        masses.append(math.sqrt(eigen))
    return sum(ratios),masses


def invert_gg(mass, top, higgs, target):
    lo,hi=0.,8.
    assert amplitudes(mass,top,higgs,lo)[0]>target>amplitudes(mass,top,higgs,hi)[0]
    for _ in range(64):
        mid=(lo+hi)/2
        if amplitudes(mass,top,higgs,mid)[0]>target:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def mixing(x2):
    return x2*(1-x2/12)


def exact_curvature_factor(mass, top, x2):
    """CW stop-only radial curvature after removing the corresponding tadpole."""
    r=top/mass
    terms=[]
    for s in (-math.sqrt(x2),math.sqrt(x2)):
        u=1+r*r+r*s
        terms.append(2*math.log(u)*(2*r+s)**2-2*u*(math.log(u)-1)*s/r)
    return sum(terms)/(16*r*r)


def certificate():
    ms=read("../1518/inputs.json")["susy_threshold_GeV"]
    assert ms == 1000 and read("../1521/inputs.json")["tan_beta"]==1
    mz=read("../1369/inputs.json")["mu_ref_GeV"]
    g1,g2,g3,yt,yb,ytau,lam=run_sm(ms,256)
    np.testing.assert_allclose(run_sm(ms,64),run_sm(ms,256),rtol=2e-10)
    exp=read("../1394/inputs.json")
    vev=(math.sqrt(2)*exp["formula"]["Gmu_GeV_minus2"])**-0.5
    mh=exp["calibration"]["m_H"][0]
    mt=yt*vev/math.sqrt(2)
    pref=3*yt**4/LOOP
    blind=invert_gg(ms,mt,mh,0.)
    xhi=invert_gg(ms,mt,mh,-.01)
    xlo=invert_gg(ms,mt,mh,.01)
    assert 2*(ms**2+mt**2-mh**2/4)/ms**2 < blind < 2*(ms**2+mt**2)/ms**2
    assert 0<xlo<blind<xhi<6
    p=(5,-12,7)
    weights=np.array([[1/30,1/2,1/3],[4/15,0,1/6],[1/5,0,0]])
    np.testing.assert_allclose(weights.sum(axis=0),[.5]*3)
    # Scalar-10 common masses cancel from differences. Permuting the same
    # three masses among generations leaves all gauge running identical.
    comparison={}
    for masses in [(ms/2,ms,2*ms),(2*ms,ms,ms/2)]:
        slopes=[]
        for mu in [ms/3,.75*ms,1.5*ms,3*ms]:
            slopes.append((sum(m<mu for m in masses)*weights.sum(axis=0)).tolist())
        logs=sum(math.log(m/ms) for m in masses)
        assert abs(logs)<1e-14
        comparison["_".join(str(m) for m in masses)]=dict(active_scalar_10_beta=slopes,
                top_leading_shift_vs_common=2*pref*math.log(masses[2]/ms))
    first,second=comparison.values()
    assert first["active_scalar_10_beta"]==second["active_scalar_10_beta"]
    assert abs(sum(w*b for w,b in zip(p,weights.sum(axis=0))))<1e-14
    max_amp,physical=amplitudes(ms,mt,mh,6.)
    mb=yb*vev/math.sqrt(2)
    squark_log=math.log(math.prod(physical)/ms**2)+math.log1p((mb/ms)**2)/2
    return {
        "scale_GeV":ms, "v_from_Gmu_GeV":vev,"matched_top_yukawa":yt,
        "leading_top_mass_for_loop_GeV":mt, "external_Higgs_mass_GeV":mh,
        "target_lambda":lam, "max_dimension4_mixing_threshold":3*pref,
        "remaining_at_max":lam-3*pref, "max_mixing_stop_masses_GeV":physical,
        "color_inverse_alpha_shift_at_high_scale_from_broken_third_Q_U_masses":squark_log/(12*math.pi),
        "max_mixing_finite_gg_amplitude_over_top":max_amp,
        "extra_quartic_from_finite_vev_curvature_at_max":pref*(exact_curvature_factor(ms,mt,6)-3),
        "zero_gg_x_squared":blind, "quartic_at_zero_gg":pref*mixing(blind),
        "fraction_of_max_quartic_at_zero_gg":mixing(blind)/3,
        "one_percent_amplitude_window_x_squared":[xlo,xhi],
        "max_quartic_at_one_percent":pref*mixing(xhi),
        "small_momentum_blind_x_squared":2*(1+(mt/ms)**2),
        "threshold_run_test_mu_GeV":[ms/3,.75*ms,1.5*ms,3*ms],
        "gauge_equivalent_family_permutations":comparison,
        "quartic_difference_of_examples":first["top_leading_shift_vs_common"]-second["top_leading_shift_vs_common"],
        "mu_ref_GeV":mz}


def compare(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(compare(i,j) for i,j in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write",action="store_true")
    args=parser.parse_args()
    numbers=certificate()
    if args.write:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as stream:
            json.dump(numbers,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    elif not compare(numbers,read("independent_results.json")):
        raise SystemExit("Independent stage-1522 certificate differs.")
    print("Verified threshold, finite stop loop, and piecewise gauge degeneracy.")
