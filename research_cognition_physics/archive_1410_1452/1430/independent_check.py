"""Spin traces, a fixed charge assignment, and response transfer with independent quadrature."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def chiral_trace(s,mu,el,gl,gr):
    eye=np.eye(2)
    zero=np.zeros((2,2))
    sigma=np.diag([1.,-1.])
    gamma0=np.block([[zero,eye],[eye,zero]])
    gamma3=np.block([[zero,sigma],[-sigma,zero]])
    momentum=sqrt((s-(mu+el)**2)*(s-(mu-el)**2)/(4*s))
    em=(s+mu**2-el**2)/(2*sqrt(s))
    ee=sqrt(s)-em
    plus=gamma0*em-gamma3*momentum+np.eye(4)*mu
    minus=gamma0*ee+gamma3*momentum-np.eye(4)*el
    g=np.diag([gl,gl,gr,gr])
    return float(np.trace(plus@g@minus@gamma0@g.conj().T@gamma0).real)


def calculation():
    cfg=read("inputs.json")
    fit=read(cfg["angles_masses"])
    atlas=read(cfg["atlas_et"])
    mK,mpi,mp,ep=[cfg[k] for k in ("MK_GeV","Mpi_GeV","mu_pole_GeV","e_pole_GeV")]
    mh,v=cfg["mh_GeV"],cfg["v_GeV"]
    limits=[(mp+ep)**2,(mK-mpi)**2]
    smin=(mp-ep)**2
    lm=fit["yukawa"]["charged_leptons"]
    dm=fit["yukawa"]["down"]
    def convert(vals):
        return np.array(vals)*v/sqrt(2)
    e,m,t=convert(np.array(lm)[:,0])
    d,st=convert(np.array(dm)[:2,0])
    n=atlas["width_multipliers"]
    # Analytic derivative of (Aq+Bq cos2phi)(Al-Bl cos2phi), minimizing derivative
    # with a ratio bound uniform down to the actual unequal-mass pair threshold.
    r_max=(dm[0][0]+n*dm[0][1])/(dm[1][0]-n*dm[1][1])
    r_min=(dm[0][0]-n*dm[0][1])/(dm[1][0]+n*dm[1][1])
    le=(lm[0][0]+n*lm[0][1])/(lm[1][0]-n*lm[1][1])
    le_b=2*le/(1+le**2)
    q_b=2*r_min/(1+r_min**2)
    derivative_margin=q_b/(1+2*q_b)-le_b
    assert derivative_margin>0
    for phase in (0.,.37,pi/2):
        gl=m*np.exp(-1j*phase)
        gr=e*np.exp(1j*phase)
        for s in np.linspace(*limits,6)[1:-1]:
            trace=chiral_trace(s,mp,ep,gl,gr)
            mirror=chiral_trace(s,ep,mp,gr.conjugate(),gl.conjugate())
            expected=(e**2+m**2)*(s-mp**2-ep**2)-4*ep*mp*e*m*np.cos(2*phase)
            assert np.isclose(trace,expected,atol=1e-15) and np.isclose(mirror,expected,atol=1e-15)
    shape=np.array(cfg["ETMC_table_IV_s_f0"])
    # Square-root endpoints are removed before Gaussian quadrature; split the
    # interpolated illustrative form factor at every listed knot.
    knots=np.r_[limits[0],shape[:,0][(shape[:,0]>limits[0])&(shape[:,0]<limits[1])],limits[1]]
    arc=np.arcsin(np.sqrt((knots-limits[0])/(limits[1]-limits[0])))
    xx,ww=np.polynomial.legendre.leggauss(40)
    all_s,all_w=[],[]
    for l,r in zip(arc,arc[1:]):
        a=(l+r)/2+xx*(r-l)/2
        s=limits[0]+(limits[1]-limits[0])*np.sin(a)**2
        ds=(limits[1]-limits[0])*np.sin(2*a)*(r-l)/2
        phase_space=np.sqrt(((mK+mpi)**2-s)*(limits[1]-s)*(s-limits[0])*(s-smin))/s
        all_s.extend(s)
        all_w.extend(ww*ds*phase_space)
    s,w=np.array(all_s),np.array(all_w)
    width=cfg["hwidth_GeV"]["direct"]
    def yield_rate(x,y,cab,mls,ratio,f0,algebraic_envelope,accepted=False):
        el,mul,tau=mls
        # Bound on |S_h/v|², two electric charge states counted by ATLAS.
        br=atlas["atlas"]["MC_marginal95"]["et"]
        factor=1-atlas["leptonic_relative_rate_factor_half_width"] if algebraic_envelope else 1.
        gh2=br*width*8*pi/(mh*(tau*tau+el*el)*(1-x)*x*(1-y)*factor)
        pds=cab*cab*(1-cab**2)
        pem=(1-x)*x*y
        quark=((mK*mK-mpi*mpi)*f0*ratio/2)**2
        if algebraic_envelope:
            effective_s=(limits[1] if accepted else s)
            lept=(mul**2+el**2)*(effective_s-smin)
            den=(mh*mh-limits[1])**2
        else:
            lept=(mul**2+el**2)*(s-mp**2-ep**2)-4*mul*el*mp*ep
            den=(mh*mh-s)**2
        return float(np.sum(w*quark*lept/den)*gh2**2*pds*pem*cfg["tauK_s"]/
                     (256*pi**3*mK**3*cfg["hbar_GeV_s"]))
    def xval(k): return (1-sqrt(1-k))/2
    x0=xval(fit["DayaBay_amplitude"][0])
    y0=fit["NO"]["atmospheric_sine_squared"][0]
    cab0=fit["quark_sines_and_phase"][0][0]
    exact=yield_rate(x0,y0,cab0,[e,m,t],(st+d)/(st-d),np.interp(s,shape[:,0],shape[:,1]),False)
    low_amps=fit["DayaBay_amplitude"]
    x=xval(low_amps[0]-n*low_amps[1])
    y=y0+n*fit["NO"]["atmospheric_sine_squared"][-1]
    cab=cab0+n*fit["quark_sines_and_phase"][0][1]
    emu=convert(np.array(lm)[:2,0]+n*np.array(lm)[:2,1])
    tau=(lm[2][0]-n*lm[2][1])*v/sqrt(2)
    # Separate denominator's electron mass from numerator to reproduce a strict outer box.
    elow=(lm[0][0]-n*lm[0][1])*v/sqrt(2)
    adjust=((tau*tau+emu[0]**2)/(tau*tau+elow**2))**2
    def envelope(accepted):
        return adjust*yield_rate(x,y,cab,[*emu,tau],(1+r_max)/(1-r_max),cfg["shape_abs_upper"],True,accepted)
    return {"unweighted_phase_space_moments":[float(w.sum()),float(w@(s-smin))],
            "positive_common_phase_test_margin":derivative_margin,
            "central_real_width330_MeV_single_charge":exact,
            "quadrature_differential_ceiling_at_330_MeV":envelope(False),
            "quadrature_flat_response_ceiling_at_330_MeV":envelope(True),
            "nominal_Ses_upper_at_10p5_MeV":envelope(True)*(.0105/width)**2/cfg["NA62"]["ses"]}


def matches(a,b):
    return a.keys()==b.keys() and all(np.allclose(a[k],b[k],rtol=1e-10,atol=1e-23) for k in a)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    group=cli.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=cli.parse_args()
    res=calculation()
    if args.write:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as h:
            json.dump(res,h,indent=2)
            h.write("\n")
    else:
        assert matches(res,read("independent_results.json"))
        author=read("results.json")
        for a,(lo,hi) in zip(res["unweighted_phase_space_moments"],author["monotone_integral_bounds_I0_I1"]):
            assert lo<a<hi
        val=author["illustration_central_inputs_real_source_saturating_et"]["direct"]
        assert np.isclose(val,res["central_real_width330_MeV_single_charge"],rtol=3e-6,atol=1e-22)
        assert res["quadrature_flat_response_ceiling_at_330_MeV"] < author["scenarios"]["direct"]["accepted_flat_B_upper"]
        print("1430 independent chiral trace, mirror charge, quadrature and accepted-yield bound verified")
