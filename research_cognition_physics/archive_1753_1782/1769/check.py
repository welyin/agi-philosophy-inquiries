"""Use a direct charged-doublet mass floor in the one-loop three-gauge projection."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG=read("inputs.json")
LOW=read(CFG["ew_gauge"])
OLD=read(CFG["odd_spectrum"])["shared_witness"]
MC, NEW_WEAK=read(CFG["even_gate"])["pole_mass_GeV"]
COLOR=np.array([*OLD["colored_fermions_GeV"],MC,MC])
MZ=LOW["mu_ref_GeV"]
SOFT=CFG["archived_soft_GeV"]
PI=math.pi
SM=np.array([41/10,-19/6,-7.])
SUSY=np.array([33/5,1.,-3.])
C=np.array([.4,0.,1.])
L=np.array([.6,1.,0.])
ADJ=np.array([0.,2.,3.])
W=np.array([5/7,-12/7,1.])
np.testing.assert_allclose([W@SUSY,W@SM,W@C,W@L],[0.,19/14,9/7,-9/7],atol=1.e-14)


def logmass(m,soft):
    return (2*np.log(m/MZ)+np.log(np.hypot(m,soft)/MZ))/3


def electroweak(sign1=0,sign2=0):
    return np.array([LOW[key][0]+sign*LOW[key][1] for key,sign in
                     (("g1_GUT",sign1),("g2",sign2))])


def gauge(weak,soft=SOFT,gs=None,mt=None):
    if gs is None:
        gs=electroweak()
    if mt is None:
        mt=CFG["top_MS_matching_GeV"][0]
    a12=4*PI/gs**2
    offs=(ADJ/(12*PI)+(SM-SUSY)*math.log(soft/MZ)/(2*PI)
          -(C*sum(logmass(COLOR,soft))+L*sum(logmass(weak,soft)))/(2*PI))
    a3=W@offs-W[:2]@a12
    coupling=math.sqrt(4*PI/a3)
    # b6-b5=2/3. At the matching scale m_t the one-loop decoupling constant has no log.
    a5=a3-math.log(mt/MZ)/(3*PI)
    above=SUSY+len(COLOR)*C+len(weak)*L
    t=(a12[0]-a12[1]-offs[0]+offs[1])/(above[0]-above[1])
    U=MZ*math.exp(2*PI*t)
    inv_U=a12[0]-offs[0]-above[0]*t
    np.testing.assert_allclose(np.r_[a12,a3],inv_U+above*t+offs,rtol=1e-14)
    assert max(COLOR.max(),weak.max(),np.hypot(max(COLOR.max(),weak.max()),soft))<U
    return dict(g3_SM=coupling,a3_SM=a3,a5_QCD=a5,alpha5=1/a5,
                crossing_GeV=U,inverse_alpha_U=inv_U)


def split(m,order):
    # Thomas-Wells Eq. 2-3. Use u^3=1-x to regularize the logarithmic endpoint.
    z,w=np.polynomial.legendre.leggauss(order)
    u=(z+1)/2
    one_minus_x=u**3
    x=1-one_minus_x
    r=(m/MZ)**2
    integral=(w/2)@(3*u**2*(1+one_minus_x)*np.log1p(x/(r*one_minus_x**2)))
    return MZ/CFG["alpha_em_MZ_inverse"]/2*math.sqrt(r)/PI*integral


def doublet(m):
    massdiff=split(m,100)
    assert abs(massdiff-split(m,64))<3e-10
    charged=read(CFG["charged_pion"])
    GF=read(CFG["weak_constants"])["const_GeV"]["GF"]
    fpi=charged["decay_constants_MeV"]["pi"][0]/1000
    mp=charged["meson_masses_GeV"]["pi"]
    gamma_pi=GF**2*fpi**2*CFG["Vud_for_lifetime"]**2/PI*massdiff**3*math.sqrt(1-(mp/massdiff)**2)
    v2=1/(math.sqrt(2)*GF)
    return dict(radiative_splitting_GeV=massdiff,leading_single_pion_width_GeV=gamma_pi,
                approximate_ctau_mm=1.973269804e-13*CFG["pion_branching_template"]/gamma_pi,
                triplet_operator_Lambda_over_abs_c_TeV_to_pion_closure=v2/(4*(massdiff-mp))/1000,
                triplet_operator_Lambda_over_abs_c_TeV_to_one_GeV_gap=v2/(4*(1-massdiff))/1000)


def floor_spectrum(floor,new_floor):
    return np.array([floor,floor,new_floor,new_floor])


def bounds(weak,soft,width):
    mt,dmt=CFG["top_MS_matching_GeV"]
    rows=[gauge(weak,soft,electroweak(i*width,j*width),mt+k*width*dmt)
          for i,j,k in product((-1,1),repeat=3)]
    return {tag:[min(row[tag] for row in rows),max(row[tag] for row in rows)]
            for tag in ("alpha5","g3_SM","a3_SM","crossing_GeV")}


def solve_monotone(fn,lo,hi):
    if not fn(lo)>=0>=fn(hi):
        return None
    for _ in range(60):
        d=(lo+hi)/2
        if fn(d)>0:
            lo=d
        else:
            hi=d
    return (lo+hi)/2


def calculate():
    floor=CFG["step_mass_floor_GeV"]
    light=OLD["weak_fermions_GeV"][0]
    archived=floor_spectrum(light,NEW_WEAK)
    live=floor_spectrum(floor,MC)
    test=gauge(archived)
    target=LOW["g3"][0]-.67*LOW["g3"][1]
    assert abs(test["g3_SM"]-target)<5e-10
    assert np.all(np.sort(live)<=np.sort(COLOR))
    out=dict(original_calibrated=test,at_225_floor_frozen_high_pair=gauge(floor_spectrum(floor,NEW_WEAK)),
             at_225_floor_relax_high_pair=gauge(live),
             independent_endpoint_mass_role_GeV=MC,low_doublet_at_original=doublet(light),
             low_doublet_at_experimental_line=doublet(CFG["pure_doublet_ATLAS_charged_limit_GeV"]),
             central_input_EW_weak_floor_sensitivity={
                 str(m):gauge(floor_spectrum(m,MC))["alpha5"] for m in [100.,light,200.,floor,250.]})
    comparisons={}
    for name,(alpha,err) in CFG["strong_lattice_5f"].items():
        table={}
        for n in (0,1,2):
            observed_low=alpha-n*err
            envelope=bounds(live,SOFT,n)
            high_alpha=envelope["alpha5"][1]
            mt=CFG["top_MS_matching_GeV"][0]+n*CFG["top_MS_matching_GeV"][1]
            favorable=electroweak(-n,n)
            np.testing.assert_allclose(gauge(live,SOFT,favorable,mt)["alpha5"],high_alpha,rtol=2e-15)
            # At the time the collider datum was taken, mass values far below M_Z
            # would require another threshold model. Restrict the inversion to m>=M_Z.
            old_allowed=solve_monotone(
                lambda m:gauge(floor_spectrum(m,MC),SOFT,favorable,mt)["alpha5"]-observed_low,
                MZ,2000.)
            soft_max=solve_monotone(lambda s:gauge(live,s,favorable,mt)["alpha5"]-observed_low,
                                    1000.,100000.)
            table[str(n)]=dict(electroweak_and_top_box=envelope,
                         minimum_lattice_alpha5=observed_low,
                         minimum_extra_projected_inverse_coupling=max(0,1/high_alpha-1/observed_low),
                         max_common_soft_GeV_if_two_pure_light=soft_max,
                         max_common_odd_doublet_step_mass_GeV_no_collider=old_allowed,
                         nominal_tree_step_must_be_below_MZ_at_chosen_inputs=old_allowed is None)
        comparisons[name]=table
    out["independent_strong_tests"]=comparisons
    trial=floor_spectrum(250.,7300.)
    lattice=CFG["strong_lattice_5f"]["ALPHA_2026"][0]
    new_soft=solve_monotone(lambda s:gauge(trial,s)["alpha5"]-lattice,1000.,30000.)
    assert new_soft is not None
    out["gauge_only_soft_recalibration"]=dict(weak_masses_GeV=trial.tolist(),soft_GeV=new_soft,
                                              **gauge(trial,new_soft))
    return out


def clean(v):
    if isinstance(v,dict):
        return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list):
        return [clean(x) for x in v]
    if v is None or isinstance(v,bool):
        return v
    return float(f"{v:.11g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=clean(calculate())
    if args.check:
        if result != read("results.json"):
            raise SystemExit("1769 saved result differs")
        print("1769 OK: two-coupling unification projection, same-spectrum collider template, 5-flavour and soft limits")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

