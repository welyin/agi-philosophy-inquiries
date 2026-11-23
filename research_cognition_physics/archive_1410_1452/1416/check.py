"""Rebuild both charged currents, scalar flavour projectors and circular holdouts."""
import argparse
import json
from math import asin, atan2, cos, pi, sin, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def rotation(i,j,angle):
    r=np.eye(3)
    r[i,i]=r[j,j]=cos(angle)
    r[i,j],r[j,i]=sin(angle),-sin(angle)
    return r


def phase_and_j(m):
    s13=abs(m[0,2])
    c13=sqrt(1-s13*s13)
    s12,s23=abs(m[0,1])/c13,abs(m[1,2])/c13
    c12,c23=sqrt(1-s12*s12),sqrt(1-s23*s23)
    j=(m[0,0]*m[1,1]*m[0,1].conjugate()*m[1,0].conjugate()).imag
    sk=j/(s12*c12*s23*c23*s13*c13*c13)
    ck=(s12**2*s23**2+c12**2*c23**2*s13**2-abs(m[2,0])**2)/(2*s12*s23*c12*c23*s13)
    return atan2(sk,ck),j,[s12,s13,s23]


def sector(ang,alpha,swap=False,reverse=False):
    t12,t13,t23=ang
    o_left=(rotation(1,2,t23)@rotation(0,2,-t13 if reverse else t13)).T
    o_right=rotation(0,1,t12)
    if swap:
        # Real weak rows of the mass seeds change relative to a fixed selector.
        change=np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,1.]])
        o_left,o_right=change@o_left,change@o_right
    p=np.diag([np.exp(1j*alpha),1,1])
    curr=o_left.T@p@o_right
    tag=np.diag([1.,0.,0.])
    return curr, o_right.T@tag@o_right


def arc_test(predicted,reference):
    a,b=sorted(predicted)
    intervals=[(reference[0]+2*k,reference[1]+2*k) for k in range(-2,3)]
    distance=min(max(0.,lo-b,a-hi) for lo,hi in intervals)
    contained=any(lo<=a and b<=hi for lo,hi in intervals)
    return {"overlap":distance==0., "whole_inside":contained, "gap_pi":distance}


def run():
    data=read("inputs.json")
    weak=read(data["quark_and_lepton_angle_inputs"])
    db=read(data["theta13_from"])["daya_bay"]
    q=weak["quark_sines_and_phase"]
    width=data["work_half_widths"]
    alpha,sigma=q[3]
    q_angles=[asin(x[0]) for x in q[:3]]
    qmat, qtag=sector(q_angles,alpha)
    np.testing.assert_allclose(phase_and_j(qmat)[0],alpha,rtol=1e-10)
    guard=data["joint_published_rounding_guard_pi"]
    out={}
    for ordering in ("NO","IO"):
        solar=weak[ordering]["solar_sine_squared"]
        atmospheric=weak[ordering]["atmospheric_sine_squared"]
        theta13=asin(sqrt((1-sqrt(1-db["sin_sq_2theta13"]))/2))
        angles=[asin(sqrt(solar[0])),theta13,asin(sqrt(atmospheric[0]))]
        ca=weak["quark_sines_and_phase"][0]
        cab=[ca[0]+j*width*ca[1] for j in (-1,1)]
        sol=[solar[0]+j*width*solar[1] for j in (-1,1)]
        assert 0<sol[0]<sol[1]<.5 and 0<cab[0]<cab[1]<1/sqrt(2)
        ratio=[cab[0]*sqrt(1-cab[0]**2)/sqrt(sol[1]*(1-sol[1])),
               cab[1]*sqrt(1-cab[1]**2)/sqrt(sol[0]*(1-sol[0]))]
        modes={}
        for swap in (False,True):
            for rev in (False,True):
                mat,tag=sector(angles,-alpha,swap,rev)
                delta,j,sines=phase_and_j(mat)
                np.testing.assert_allclose(sines,np.sin(angles))
                analytic=(-1 if not swap else 1)*alpha+(pi if rev else 0.)
                np.testing.assert_allclose(np.exp(1j*delta),np.exp(1j*analytic),atol=1e-13)
                arc=[delta/pi-width*sigma/pi,delta/pi+width*sigma/pi]
                comparison={}
                for name,limits in data["phase_credible_arcs_pi"].items():
                    v=arc_test(arc,limits)
                    v["whole_inside_under_rounding"]=arc_test(arc,[limits[0]+guard,limits[1]-guard])["whole_inside"]
                    v["overlap_under_inward_rounding"]=arc_test(arc,[limits[0]+guard,limits[1]-guard])["overlap"]
                    comparison[name]=v
                modes[f"axis{2 if swap else 1}_13{'minus' if rev else 'plus'}"]={
                    "delta_pi":delta/pi,"quark_two_width_arc_pi":arc,"J":float(j),
                    "nu_scalar_diagonal_tag":tag.diagonal().tolist(),"readouts":comparison}
        out[ordering]={"phase_cases":modes,"scalar_amplitude_ratio_ds_over_nu12":float(abs(qtag[0,1])/sqrt(solar[0]*(1-solar[0]))),
                       "scalar_ratio_two_widths":ratio}
    return out


def compare(got,expected):
    if isinstance(got,dict):
        assert got.keys()==expected.keys()
        for key in got:
            compare(got[key],expected[key])
    elif isinstance(got,(bool,np.bool_)):
        assert bool(got) is expected
    else:
        np.testing.assert_allclose(got,expected,rtol=5e-10,atol=2e-13)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write",action="store_true")
    result=run()
    if cli.parse_args().write:
        with (HERE/"results.json").open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(result,read("results.json"))
    print(json.dumps({k:[v["scalar_amplitude_ratio_ds_over_nu12"],v["scalar_ratio_two_widths"],
                          {m:[t["delta_pi"],t["J"]] for m,t in v["phase_cases"].items()}] for k,v in result.items()}))
