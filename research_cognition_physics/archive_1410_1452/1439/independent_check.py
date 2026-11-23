"""Exact tree-graph messenger case and independent finite-loop check."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def heavy_triangle(h,m):
    # A direct positive Feynman-parameter integral instead of an arcsine formula.
    nodes,weights=np.polynomial.legendre.leggauss(48)
    x=(nodes+1)/2
    weights=weights/2
    y=(1-x[:,None])*x[None,:]
    z=h*h/(4*m*m)
    product=4*x[:,None]*y
    return float(4*np.sum(weights[:,None]*weights[None,:]*(1-x[:,None])*(1-product)/(1-z*product)))


def point(yuk,bare,sm,bottom,h):
    left,m,right=np.linalg.svd(bare+yuk)
    response=np.diag(left.T@yuk@right.T)/m
    step=1e-5
    hi=np.linalg.svd(bare+math.exp(step)*yuk,compute_uv=False)
    lo=np.linalg.svd(bare+math.exp(-step)*yuk,compute_uv=False)
    fd=(np.log(hi)-np.log(lo))/(2*step)
    np.testing.assert_allclose(fd,response,rtol=2e-7,atol=2e-8)
    assert abs(response.sum()-3)<1e-11
    loops=[heavy_triangle(h,heavy) for heavy in m[:2]]
    amplitude=sm+(response[-1]-1)*bottom+np.dot(response[:2],loops)
    return {"masses_ascending_GeV":m[::-1].tolist(),"kappas_ascending":response[::-1].tolist(),
            "left_b_singlet_probability":float(left[1,-1]**2),
            "right_b_doublet_probability":float(right[-1,2]**2),
            "heavy_loop_forms":loops,"pure_gg_ratio":float(abs(amplitude/sm)**2)}


def compute():
    cfg=read(HERE/"inputs.json")
    old=read(HERE/cfg["production"])
    author=read(HERE/"results.json")
    radial=read(HERE/cfg["radial"])
    form,masses=runpy.run_path(str(HERE.parent/"1377"/"check.py"))["inherited_spectrum"]()
    h=cfg["h_GeV"]
    m=masses[-1]
    sm=sum(form(h,mi)[1] for mi in masses)
    bottom=form(h,m)[1]
    a=c=650.
    md=mq=300.
    # det(M^* M-m^2 I)=0 solved explicitly for the middle edge, all units in GeV.
    middle=m*math.sqrt((a*a+md*md-m*m)*(c*c+mq*mq-m*m)/((a*a-m*m)*(c*c-m*m)))
    bare=np.diag([0.,md,mq])
    tree=point(np.array([[0.,a,0.],[0.,0.,middle],[c,0.,0.]]),bare,sm,bottom,h)
    assert np.isclose(tree["masses_ascending_GeV"][0],m)
    assert 0<=tree["kappas_ascending"][-1]<=1 and tree["masses_ascending_GeV"][1]>300
    k=radial["bb"]["kappa_upper"]
    assert tree["kappas_ascending"][0]<k and tree["kappas_ascending"][1]>=2-k
    floor=sm.real+(k-1)*bottom.real+(3-k)*4/3
    rate=(1-old["central_upper_bb_fraction"])*cfg["min_gg_ratio_transfer"]*floor**2/abs(sm)**2
    assert rate>old["marginal_rectangles"]["2"]["upper_ratio"]
    edge=author["admitted_example"]["outer_mass_insertion_GeV"]
    loop=point(np.array([[0.,edge,0.],[0.,0.,700.],[edge,-700.,0.]]),bare,sm,bottom,h)
    np.testing.assert_allclose(loop["pure_gg_ratio"],author["admitted_example"]["pure_gg_rate_over_SM"],rtol=1e-10)
    return {"tree_middle_edge_GeV":middle,"tree_instance":tree,"closed_loop_author_point":loop,
            "tree_finite_mass_production_floor":rate,"T_mass_GeV":mq,"T_Higgs_log_response":0.}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif args.check:
        def equal(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
            return np.allclose(a,b,atol=1e-12,rtol=4e-10)
        assert equal(result,read(path))
        print("1439 independent mass, weak representations and finite coloured responses verified")
    else:
        print(json.dumps(result,indent=2))
