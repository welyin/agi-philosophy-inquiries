"""Independent spin traces and stability of the one-direction cross-family identity."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def gamma():
    eye,zero=np.eye(2),np.zeros((2,2))
    matrices=[np.block([[eye,zero],[zero,-eye]])]
    for p in [np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]:
        matrices.append(np.block([[zero,p],[-p,zero]]))
    return matrices


def spin_trace(l,r,g0,vertex):
    conjugate_vertex=g0@vertex.T.conj()@g0
    return np.trace(l@vertex@r@conjugate_vertex).real


def compute():
    gf=1.1663787e-5
    mt,mh,mu=172.5,125.,.0015
    mB,mm=5.27965,.1056583755
    v2=1/(math.sqrt(2)*gf)
    g0,g1,g2,g3=gamma()
    ident=np.eye(4)
    g5=1j*g0@g1@g2@g3
    p=math.sqrt(mB*mB/4-mm*mm)
    l=mB/2*g0-p*g3+mm*ident
    r=mB/2*g0+p*g3-mm*ident
    scalar=spin_trace(l,r,g0,ident)
    pseudo=spin_trace(l,r,g0,g5)
    assert math.isclose(scalar,2*(mB*mB-4*mm*mm))
    assert math.isclose(pseudo,2*mB*mB)
    for angle in [0,.5,2.2]:
        vertex=.3*ident+np.exp(1j*angle)*g5
        assert math.isclose(spin_trace(l,r,g0,vertex),.3**2*scalar+pseudo)
    eu=(mt**2+mu**2-mh**2)/(2*mt)
    pu=math.sqrt(eu**2-mu**2)
    vertex=(mt*(ident-g5)+mu*(ident+g5))/(2*math.sqrt(v2))
    avg=spin_trace(eu*g0-pu*g3+mu*ident,mt*(g0+ident),g0,vertex)/2
    exact_k=pu*avg/(8*math.pi*mt**2)
    zero_k=mt**3/(32*math.pi*v2)*(1-(mh/mt)**2)**2
    beta=2*p/mB
    norm=(1.512/6.582119569e-13)*mB*beta/(8*math.pi)*(.186*mm/(2*v2))**2
    w=beta*mB**2*.995/mh**2
    bound=math.sqrt(1.9e-10*.00019*1.4/(norm*w*w*zero_k))
    hypothetical_charge_product=.01
    theta=.5*math.acos(math.sqrt(bound/hypothetical_charge_product))
    plus=np.array([math.cos(theta),math.sin(theta)])
    minus=np.array([math.cos(theta),-math.sin(theta)])
    # Scale by 0.1 to make both same-sector Gram blocks contractions.
    up=np.array([.1*plus,.1*minus])
    down=up.copy()
    neutral=abs((up@up.T)[0,1]*(down@down.T)[0,1])
    charged=abs((up@down.T)[0,0]*(up@down.T)[1,1])
    assert np.linalg.norm(up,2)<1 and math.isclose(neutral/charged,math.cos(2*theta)**2)
    return {"Bd_spin_weight_ratio_S_over_P":scalar/pseudo,
            "Kt_massless_unit_GeV":zero_k,"Kt_ratio_at_up_mass_1p5MeV":exact_k/zero_k,
            "CMS_factorized_product_limit":bound,"hypothetical_right_product":.01,
            "minimum_maximum_ray_angle_deg_if_point01":math.degrees(theta),
            "minimum_maximum_orthogonal_fraction_if_point01":math.sin(theta)**2}


if __name__=="__main__":
    cmd=argparse.ArgumentParser()
    mode=cmd.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cmd.parse_args()
    data=compute()
    out=HERE/"independent_results.json"
    if args.write:
        with out.open("x",encoding="utf-8") as f:
            f.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    else:
        saved=json.loads(out.read_text(encoding="utf-8"))
        assert saved.keys()==data.keys() and all(math.isclose(saved[k],v,rel_tol=2e-12) for k,v in data.items())
    print("PASS Dirac spin traces, independently normalized bounds, and second-source angles")
