"""Positive two-particle cut, pole derivative and finite recoil response, in GeV."""
import argparse
import json
from itertools import product
from math import pi, sqrt
from pathlib import Path

import numpy as np


V, MH = 248.401, 125.11
MZ, ECM = 91.1876, 2000.
QLO,QHI = 1100.,1500.
NODES, WEIGHTS = np.polynomial.legendre.leggauss(80)


def integrate(fun, lo, hi):
    return float((hi-lo)/2 * (WEIGHTS @ fun((hi-lo)/2*NODES+(hi+lo)/2)))


def charged_trace(g1,g2,q,e):
    return (6*g2**4+3*(g2**2+.6*g1**2)**2)/16-3*np.sum(q**4)-np.sum(e**4)


def known_trace_values():
    q=np.array([7.04e-6,.00356,.967,1.54e-5,.000306,.01630])
    dq=np.array([.15e-6,.00006,.004,.02e-5,.000004,.00009])
    ell=np.array([2.77713e-6,.000585042,.0099378])
    dell=np.array([.00036e-6,.000000075,.0000014])
    return [charged_trace(.461228+k*.000026,.65096+k*.00004,q-k*dq,ell-k*dell)
            for k in (-1,0,1)]


def d_target(v,mh,c):
    return (8*pi**2*(mh/v)**2-c)*v**4


def single_delta(m, v, mh):
    # Dispersion: s'=4m²/(1-u²), cut Im(Pi)=g²u/(32 pi).
    g=2*m*m/v
    value = g*g/(32*pi**2)*integrate(
        lambda u: 8*m*m*u**2/(4*m*m-mh*mh*(1-u**2))**2, 0,1)
    parameter_check = g*g/(32*pi**2)*integrate(
        lambda u: u*(1-u)/(m*m-mh*mh*u*(1-u)), 0,1)
    assert np.isclose(value,parameter_check,rtol=2e-13,atol=1e-15)
    return value


def zh_kinematic(s):
    sr=ECM**2
    z=MZ**2/sr
    lam=(1-z-s/sr)**2-4*z*s/sr
    return np.sqrt(lam)*(lam+12*z)


def recoil(v,mh,weights_m4):
    # (rho ds) gives the spectrum divided at each s by the known unit
    # Higgsstrahlung source; weighting by the phase space gives a rate ratio.
    def spectrum(s):
        cut=np.zeros_like(s)
        for fourth in weights_m4:
            cut+=fourth*np.sqrt(np.maximum(0.,1-4*sqrt(fourth)/s))
        return cut/(8*pi**2*v*v*(s-mh*mh)**2)
    mass_weight=lambda s: zh_kinematic(s)/zh_kinematic(mh**2)
    return {"source_normalized":integrate(spectrum,QLO**2,QHI**2),
            "rate_relative_to_LO_Zh125":integrate(lambda s: spectrum(s)*mass_weight(s),QLO**2,QHI**2)}


def weak_tree_fb(v,mh,fourth):
    """Use only common running weak parameters in the dimensionful LO cross section."""
    gweak, gy = .65096, sqrt(.6)*.461228
    mz = v*sqrt(gweak**2+gy**2)/2
    sw2 = gy**2/(gweak**2+gy**2)
    sr=ECM**2
    z=mz*mz/sr
    gf=1/(sqrt(2)*v*v)
    prefactor=gf**2*mz**4*((-1+4*sw2)**2+1)/(96*pi*sr*(1-z)**2)
    gev2_to_fb=0.3893793721e12
    def unit_sigma(s):
        lam=(1-z-s/sr)**2-4*z*s/sr
        return gev2_to_fb*prefactor*np.sqrt(lam)*(lam+12*z)
    def production(s):
        cut=fourth*np.sqrt(1-4*sqrt(fourth)/s)
        return unit_sigma(s)*cut/(8*pi**2*v*v*(s-mh**2)**2)
    return {"tree_Z_mass_GeV":mz,"LO_Zh_fb":unit_sigma(mh*mh),
            "single_real_window_floor_fb":integrate(production,QLO**2,QHI**2)}


def compute():
    cmin,c,cmax=known_trace_values()
    d=d_target(V,MH,c)
    m=d**.25
    rows={}
    for name,fractions in (("one_real", [1.]), ("four_real_equal", [.25]*4),
                           ("three_real_unequal",[.15,.3,.55])):
        parts=d*np.array(fractions)
        rows[name]={"masses":(parts**.25).tolist(),
                    "pole_delta":sum(single_delta(h**.25,V,MH) for h in parts),
                    "recoil":recoil(V,MH,parts)}
    lower=rows["one_real"]
    for another in list(rows.values())[1:]:
        assert another["pole_delta"]>lower["pole_delta"]
        for key in lower["recoil"]:
            assert another["recoil"][key]>lower["recoil"][key]
    # Source box endpoints propagated into both the required fourth-moment
    # and the pole denominator at once; no likelihood claim is attached.
    corners=[(v,mh,d_target(v,mh,trace)) for v,mh,trace in
             product([V-.032,V+.032],[MH-.11,MH+.11],[cmin,cmax])]
    values=[single_delta(fourth**.25,v,mh) for v,mh,fourth in corners]
    dmin=min(q[2] for q in corners)
    dmax=max(q[2] for q in corners)
    # Positive cut gives a conservative floor for every point of the full
    # input box, without assuming monotonicity in D near the pair threshold.
    portal_floor=recoil(V+.032,MH-.11,[dmax])
    # rho used the heavier threshold and maximum numerator; replace D by Dmin.
    # For the actual rate ratio take the largest possible unit on-shell rate,
    # attained here at the smallest m_h (Zh phase space decreases with mass).
    recoil_floor={k:y*dmin/dmax for k,y in portal_floor.items()}
    assert 2*dmax**.25<QLO<QHI<ECM-MZ
    assert lower["pole_delta"]>sqrt(d)/(48*pi**2*V**2)
    # Unknown extra widths allow equal visible rates for different residues.
    # Positive widths are not available in the strict closed-portal-only theory.
    same_mu=.94
    width_alias=[]
    for entry in (lower,rows["four_real_equal"]):
        k=1-entry["pole_delta"]
        extra=k*k/same_mu-k
        assert extra>=0 and np.isclose(k*k/(k+extra),same_mu)
        width_alias.append(extra)
    return {"inputs":"same charged MZ values and leading mass match as round 1372",
            "v":V,"mh":MH,"Mstar":m,"single_real_portal":2*m*m/V**2,
            "window_GeV": [QLO,QHI],"ecm_GeV":ECM,
            "sharply_minimized_pole_delta":lower["pole_delta"],
            "zero_momentum_minimum":sqrt(d)/(48*pi**2*V**2),
            "full_box_sharp_minimum_interval":[min(values),max(values)],
            "full_box_conservative_recoil_floor":recoil_floor,
            "common_EW_scheme_recoil":weak_tree_fb(V,MH,d),
            "fixed_fourth_moment_examples":rows,
            "extra_widths_over_SM_for_same_visible_mu_0.94":width_alias}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(i) for i in x]
    return float(format(x,".9g")) if isinstance(x,float) else x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true",help="verify frozen results (default)")
    options=p.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    actual=rounded(compute())
    if options.write:
        with file.open("x",encoding="utf-8") as out:
            json.dump(actual,out,indent=2,ensure_ascii=False)
            out.write("\n")
    else:
        assert actual==json.loads(file.read_text(encoding="utf-8"))
    print("1374 independent spectral, recoil and rate-alias checks passed.")
