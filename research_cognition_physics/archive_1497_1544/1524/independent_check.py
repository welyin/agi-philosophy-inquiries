"""Joint top/quartic elbow and a weak-current spectral check on the Higgs vacuum angle."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads((BASE/path).read_text(encoding="utf-8-sig"))


def make_frame(n, seed):
    rng=np.random.default_rng(seed)
    raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    return np.linalg.qr(raw)[0]


def mixed_moment(mw, c):
    # The charged and neutral matrices may each mix with singlets. Construct the
    # same active soft block and arbitrary different allowed mixing entries.
    soft=make_frame(3,3)@np.diag([.82e6,1.e6,1.12e6])@make_frame(3,3).conj().T
    lepton=np.diag([0.,.1**2,1.76**2])  # sample Dirac mass-squared, not an input fit
    def hermitian(diagonal, n, seed):
        rng=np.random.default_rng(seed)
        off=1.e4*(rng.normal(size=(3,n))+1j*rng.normal(size=(3,n)))
        return np.block([[diagonal,off],[off.conj().T,np.diag(np.linspace(.9e6,1.2e6,n))]])
    # Common neutral D shift absorbed into 'soft'.
    ec,uc=np.linalg.eigh(hermitian(soft+lepton+mw**2*c*np.eye(3),3,7))
    en,un=np.linalg.eigh(hermitian(soft,2,11))
    flux=abs(uc[:3].conj().T@un[:3])**2
    value=float(np.sum(flux*(ec[:,None]-en[None,:])))
    assert min(ec)>0 and min(en)>0 and abs(flux.sum()-3)<1e-12
    assert abs(value-np.trace(lepton)-3*mw**2*c)<2e-8
    return dict(sum_weights=float(np.sum(flux)),raw_mass_moment_GeV2=value,
                recovered_c=float((value-np.trace(lepton))/(3*mw**2)))


def certificate():
    quartic=read("../1522/results.json")
    transport=read("../1523/results.json")
    data=read("../1523/inputs.json")
    ew=read("../1394/inputs.json")
    mw,mw_error=ew["measured"]["CMS"]
    mz=ew["calibration"]["m_Z"][0]
    ms=quartic["scale_GeV"]
    ma=ms
    hratio=1+(mz/ma)**2
    def vacuum_c(projection):
        return projection/math.sqrt(hratio**2-(hratio**2-1)*projection**2)
    favorable=transport["windows"]["1518_slow_growth"]
    g1,g2,_=favorable["gauge_at_low"]
    dmax=(g2*g2+.6*g1*g1)/8
    ymax=favorable["yt_squared_4pi"]["quadrature_bracket"][1]
    yte,dyte=data["sm_dr_at_one_TeV"]["yt"]
    ctop=2*((yte-2*dyte)/ymax)**2-1
    cases={}
    for s in (1,3):
        deficit=quartic["bounded_cases"][str(s)]["necessary_other_threshold_whole_input_box"]
        row={}
        for budget in (0,.01,.03):
            absolute=math.sqrt(max(0.,(deficit-budget)/dmax))
            signed=absolute if s==1 else None
            if s==1:
                assert 0<ctop<absolute
            # Charged lepton masses removed, both scheme/loop and mass errors
            # in the spectral moment remain to be bounded by a real experiment.
            splitting=(mw-2*mw_error)**2*absolute
            finite_split=(mw-2*mw_error)**2*vacuum_c(absolute)
            row[str(budget)]=dict(abs_c_min=absolute,positive_c_min_with_shared_top_transport=signed,
                                  min_absolute_mean_normalized_mass_difference_GeV2=splitting,
                                  equivalent_gap_at_mass_sum_2000_GeV=splitting/2000,
                                  finite_tree_mA_vacuum_c_min=vacuum_c(absolute),
                                  finite_tree_mA_gap_at_mass_sum_2000_GeV=finite_split/2000)
        cases[str(s)]=row
    c=.98
    proj=hratio*c/math.sqrt(1+(hratio**2-1)*c*c)
    assert abs(vacuum_c(proj)-c)<1e-13
    central=transport["windows"]["1518_central"]["gauge_at_low"]
    dg=(central[1]**2+.6*central[0]**2)/8
    lam=quartic["quartic"]
    pref=quartic["coefficient"]
    required=(lam-dg*proj*proj)/pref
    assert 0<required<3
    xroots=[6-sign*math.sqrt(36-12*required) for sign in (1,-1)]
    gap=math.sqrt(c*(ma**2+mz**2))
    m1=(math.hypot(gap,2*ms)+gap)/2
    m2=ms**2/m1
    weighted=m1*m1+m2*m2-2*ms**2
    assert abs(weighted-(m1-m2)**2)<1e-8
    assert abs(weighted/(ma**2+mz**2)-c)<1e-12
    tanb=math.sqrt((1+c)/(1-c))
    # Full MSSM leading third-family flow above both right slepton thresholds.
    # Gaugeless lepton Yukawas of e and mu are neglected. Decoupling a lone E^c
    # lowers b_1 by 1/5 until its threshold. The light partner beforehand raised
    # g1 relative to the degenerate spectrum, and exact endpoint logs cancel.
    upper_end=transport["first_heavy_doublet_GeV"]
    g3min2=1/(1/central[2]**2+3*math.log(upper_end/ms)/(8*math.pi**2))
    g12,g22=central[0]**2,central[1]**2
    ty,by,ly=.86,.20,.15
    yuk_bound=[6*ty**2+by**2,ty**2+6*by**2+ly**2,3*by**2+4*ly**2]
    loss=[13*g12/15+3*g22+16*g3min2/3,
          7*g12/15+3*g22+16*g3min2/3,9*g12/5+3*g22]
    assert all(u<g for u,g in zip(yuk_bound,loss))
    ht=yte/math.sqrt((1+proj)/2)
    # Independent finite matching bounds on yb and ytau in the chosen witness.
    cosb=math.sqrt((1-proj)/2)
    assert ht<ty and .016/cosb<by and .012/cosb<ly
    z=mw**2*c
    readout=mixed_moment(mw,c)
    # Same two charged eigenmasses and neutral mass with opposite active weights.
    span=6500.
    p=(1+z/span)/2
    assert 0<p<1
    return dict(c_top_two_width_integral_outer=ctop, d_max=dmax, cases=cases,
                common_scale_witness=dict(vacuum_c=c,projection_c=proj,tan_vacuum_beta=tanb,x_squared=xroots,top_initial_DR=ht,
                    assumed_yb_low_upper=.016,assumed_ytau_low_upper=.012,
                    rhs_yuk_upper=yuk_bound,rhs_gauge_lower=loss, mA_GeV=ma,
                    right_slepton_soft_GeV=[m1,m2],mu_threshold_GeV=ms,
                    right_slepton_delta_A1_at_Q=-.2*math.log(ms/m2)/(2*math.pi),
                    soft_trace_GeV2=weighted,charged_minus_neutral_after_lepton_GeV2=z,
                    weak_weighted_spectrum=readout),
                opposite_c_same_masses=dict(neutral_GeV=ms,
                    charged_GeV=[math.sqrt(ms**2-span),math.sqrt(ms**2+span)],
                    heavy_charged_left_weights=[p,1-p],c_options=[c,-c]),
                missing_current_example=dict(missing_weight_per_active=.001,
                    minimum_unseen_mass_squared_span_to_mask_sign_GeV2=z/.001))


def equal(x,y):
    if isinstance(x,dict):
        return isinstance(y,dict) and x.keys()==y.keys() and all(equal(x[k],y[k]) for k in x)
    if isinstance(x,list):
        return isinstance(y,list) and len(x)==len(y) and all(equal(a,b) for a,b in zip(x,y))
    if x is None:
        return y is None
    return math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-8)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    command=parser.add_mutually_exclusive_group()
    command.add_argument("--check",action="store_true")
    command.add_argument("--write",action="store_true")
    cfg=parser.parse_args()
    evidence=certificate()
    if cfg.write:
        with (BASE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(evidence,f,indent=2,allow_nan=False)
            f.write("\n")
    elif not equal(evidence,read("independent_results.json")):
        raise SystemExit("Mismatch in the independent spectral certificate.")
    print("Complete charged-current moment, two-width bounds and soft eigenstate witnesses verified.")
