"""Infer the light family axis from full tree CKM inputs and test a common source."""
import argparse
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location("textures1486", HERE.parent/"1486"/"check.py")
old = module_from_spec(spec)
spec.loader.exec_module(old)
I = old.I
hull = old.hull


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def published(row):
    x = I.point(row[0])
    w = sum((I.point(z) for z in row[1:]),I.point(0))
    return hull([x-2*w,x+2*w])


def cosine_positive_degrees(value):
    """Cos with a rational pi enclosure and upper/lower alternating Taylor sums."""
    pi_min = Fraction("3.14159265358979323846264338327950288")
    pi_max = pi_min+Fraction(1,10**35)
    val = Fraction(str(value))/180
    assert 0 < val*pi_max < Fraction(3,2)
    def partial(arg,last):
        return sum((-1)**k*arg**(2*k)/math.factorial(2*k) for k in range(last+1))
    return hull([I.point(partial(pi_max*val,11)),I.point(partial(pi_min*val,10))])


def tree_axis(us,cb,ratio,cos_gamma):
    s13 = cb*ratio
    c13 = (1-s13.square()).sqrt()
    s12,s23 = us/c13,cb/c13
    c12,c23 = (1-s12.square()).sqrt(),(1-s23.square()).sqrt()
    k = c12*s23*s13/(s12*c23)
    sin2_gamma = 1-cos_gamma.square()
    cd = cos_gamma*(1-k.square()*sin2_gamma).sqrt()-k*sin2_gamma
    z = c23*c13*ratio
    return ((s12.square()+(c12*z).square()-2*s12*c12*z*cd)/(1+z.square())).sqrt()


def central_ckm(us,cb,r,gamma):
    s13=cb*r
    c13=math.sqrt(1-s13**2)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12**2),math.sqrt(1-s23**2)
    delta=gamma+math.asin(c12*s23*s13/(s12*c23)*math.sin(gamma))
    phase=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13*phase.conjugate()],
        [-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13],
        [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13]])


def fit_textures(sd,params):
    masses = read(HERE/"../1480/inputs.json")["lattice"]
    q = old.width(masses["md"])/old.width(masses["ms"])
    osc = read(HERE/"../1484/inputs.json")
    R = (1+old.width(osc["T2K_dm32_NO_eV2"])/old.width(osc["KamLAND_dm21_eV2"])).sqrt()
    p3 = (1-(1-old.width(osc["DayaBay_sin2_2theta13"])).sqrt())/2
    p1 = (1-p3)*(1-old.width(osc["JUNO_sin2_theta12_NO"]))
    out = {}
    for name, factors in params.items():
        boxes = [old.charged(qi,si,factors) for qi in q.parts(24) for si in sd.parts(24)]
        se = hull(pair[0] for pair in boxes)
        spectrum = hull(pair[1] for pair in boxes)
        if name == "small_left":
            # For this Clebsch placement L+1/L decreases with q and increases with sd.
            def sharp_value(q0,s0):
                b2 = s0.square()*(1-q0.square())+q0.square()
                trace_over_det = 12*(q0+1/q0)-143*q0/(12*b2)
                return (trace_over_det+(trace_over_det.square()-4).sqrt())/2
            assert 0 < q.lo < q.hi < sd.lo
            spectrum=hull([sharp_value(I(q.hi,q.hi),I(sd.lo,sd.lo)),
                           sharp_value(I(q.lo,q.lo),I(sd.hi,sd.hi))])
        passed = []
        for branch in ("small","large"):
            bounds = old.p3_bounds(se,R,p1,branch)
            if any(max(row.lo,p3.lo)<=min(row.hi,p3.hi) for row in bounds):
                passed.append(branch)
        out[name] = {"angle_sine":se.values(),"mass_ratio_window":spectrum.values(),
                     "reactor_overlap_possible":passed}
    return out


def solve():
    p = read(HERE/"inputs.json")
    data = read(HERE/p["quark_inputs"])
    bella = read(HERE/p["belle"])["BelleII_Vcb_BGL"]
    texture = read(HERE/p["lepton_model"])["clebsch"]
    us = published(data["K_product"])/published(data["K_form"])
    cb,r = published(bella),published(data["ub_over_cb"])
    angles = p["chargedB"]
    g_lo = angles["central_degrees"]-2*angles["lower_error"]-angles["additional_degrees"]
    g_hi = angles["central_degrees"]+2*angles["upper_error"]+angles["additional_degrees"]
    cg = hull([cosine_positive_degrees(round(g_hi,5)),cosine_positive_degrees(round(g_lo,5))])
    mix = published(p["mixing"]["td_over_ts"])
    mix_axis = hull(x/(1+x.square()).sqrt() for x in (I(mix.lo,mix.lo),I(mix.hi,mix.hi)))
    use_samples = [us.parts(4),cb.parts(2),r.parts(4),cg.parts(6)]
    boxes = []
    for sign in (1,-1):
        sd = hull(tree_axis(a,b,v,sign*c) for a in use_samples[0] for b in use_samples[1]
                  for v in use_samples[2] for c in use_samples[3])
        tangent = sd/(1-sd.square()).sqrt()
        hmix = (mix/tangent).square()
        lo,hi = max(sd.lo,mix_axis.lo),min(sd.hi,mix_axis.hi)
        boxes.append({"sd":sd.values(),"td_over_ts":tangent.values(),
                      "oscillation_ratio_multiplier":hmix.values(),
                      "intersects_mixing":lo<=hi,
                      "common_texture_test":fit_textures(I(lo,hi),texture) if lo<=hi else None})
    assert boxes[0]["intersects_mixing"] and not boxes[1]["intersects_mixing"]
    assert [key for key,row in boxes[0]["common_texture_test"].items()
            if row["reactor_overlap_possible"]] == ["small_left"]
    inputs = read(HERE/"../1480/inputs.json")
    central_ratio=inputs["lattice"]["md"][0]/inputs["lattice"]["ms"][0]
    osc=read(HERE/"../1484/inputs.json")
    R=math.sqrt(1+osc["T2K_dm32_NO_eV2"][0]/osc["KamLAND_dm21_eV2"][0])
    lep=inputs["muonium"]["pole_ratio"][0]
    lep_ms=lep*(1+3*math.log(lep)/(2*math.pi*p["QED_alpha_inverse"]))
    low_mass=boxes[0]["common_texture_test"]["small_left"]["mass_ratio_window"][0]
    upper_lep=lep+2*inputs["muonium"]["pole_ratio"][1]
    upper_ms=upper_lep*(1+3*math.log(upper_lep)/(2*math.pi*p["QED_alpha_inverse"]))
    max_transfer=upper_ms/low_mass
    vector=read(HERE/"../1488/inputs.json")
    s2,err=read(HERE/"../1477/inputs.json")["SLD_sin2_eff_e"]
    low_angle,high_angle=s2-2*err,s2+2*err
    high_t=max_transfer*1.005
    left=max(high_angle,abs(low_angle-high_t**2/2),abs(high_angle-high_t**2/2))
    zmax=(left**2+high_angle**2)/((high_angle-.5)**2+high_angle**2)
    re,ere=vector["lep"]["R_e"]
    rm,erm=vector["lep"]["R_mu"]
    zfloor=(re-2*ere)/(rm+2*erm)
    belle_tau=read(HERE/"../1482/inputs.json")["Belle_II"]
    wfloor=belle_tau["ratio_central"]-2*belle_tau["ratio_width"]
    ef=vector["amplitude_or_rate_transfer"][1]
    assert high_t<wfloor<=high_t*ef and zmax*ef<zfloor
    center=central_ckm(data["K_product"][0]/data["K_form"][0],bella[0],data["ub_over_cb"][0],
                       math.radians(angles["central_degrees"]))
    top=center[2].conj()
    bottom=np.array([0.,0.,1.])
    line=np.array([center[2,1],-center[2,0],0.])
    line/=np.linalg.norm(line)
    np.testing.assert_allclose([top.conj()@line,bottom@line],[0,0],atol=1e-13)
    np.testing.assert_allclose(center.conj().T@center,np.eye(3),atol=1e-13)
    sd=abs(line[1])
    su=abs((center@line)[1])
    assert math.isclose(su,data["ub_over_cb"][0]/math.hypot(1,data["ub_over_cb"][0]))
    witness=old.witness(central_ratio,sd,R,texture["small_left"])
    witness["common_scale_transfer_at_one_loop"]=lep_ms/witness["lepton_mass_ratio"]
    witness["pole_ratio_transfer"]=lep/witness["lepton_mass_ratio"]
    c=abs(top.conj()@bottom)
    pt=np.eye(3)-np.outer(top,top.conj())
    pb=np.eye(3)-np.outer(bottom,bottom)
    pi=np.outer(line,line.conj())
    np.testing.assert_allclose(np.linalg.norm(np.linalg.matrix_power(pt@pb,3)-pi,ord=2),c**5,atol=2e-13)
    goal=p["filter_amplitude_example"]
    n=math.ceil((1+math.log(goal)/math.log(c))/2)
    assert c**(2*n-1)<=goal<c**(2*n-3)
    tag=p["symmetric_tag_amplitude"]
    return {"tree_and_mixing":boxes,"center_light_sines_up_down":[su,sd],
            "conditional_electron_witness":witness,"max_one_loop_transfer":max_transfer,"heavy_overlap":c,
            "filter_cycles_amplitude":n,
            "higgs_only_muon_singlet":{"max_T_with_transport":high_t,"max_W_with_vertex_budget":high_t*ef,
            "Belle_W_lower":wfloor,"max_Z_with_vertex_budget":zmax*ef,"LEP_Z_lower":zfloor},
            "max_axis_error_degrees_at_tag":math.degrees(math.asin(tag*math.sqrt(2/(1-c)))),
            "tag_for_requested_angle":math.sin(math.radians(p["angular_tolerance_degrees"]))*math.sqrt((1-c)/2)}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    out=HERE/"results.json"
    result=solve()
    if args.write_new:
        with out.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert result==read(out)
        print("1490: tree-only heavy-axis reconstruction, mixing holdout, projected textures and filter precision checked.")
