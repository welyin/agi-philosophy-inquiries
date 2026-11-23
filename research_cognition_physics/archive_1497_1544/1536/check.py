"""Tree-level background determinant, Yukawa assignments and held-out tau spin."""
import argparse
import cmath
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def evaluate():
    cfg=read("inputs.json")
    prev=read(cfg["earlier"])
    types={"I":(2,2),"II":(1,1),"X":(2,1),"Y":(1,2)}
    data={}
    # Two complex neutral backgrounds. The common up matrix always uses the conjugate of 2.
    for label,(down,lepton) in types.items():
        beta=math.pi/4
        amp=.37
        up_odd=-amp/math.tan(beta)
        non_up={1:-amp*math.tan(beta),2:amp/math.tan(beta)}
        z1,z2=1.12*cmath.exp(.17j),.69*cmath.exp(-.21j)
        factor=(z2.conjugate()*(z1 if down==1 else z2))**3
        data[label]=dict(odd_top_equals_tau_at_tanbeta1=abs(up_odd-non_up[lepton])<1e-12,
                          protected=down==2,finite_phase=cmath.phase(factor),
                          mass_phase_slope=3*(up_odd+non_up[down]))
        assert (abs(cmath.phase(factor))<1e-12)==(down==2)
    assert [k for k,v in data.items() if v["protected"] and v["odd_top_equals_tau_at_tanbeta1"]]==["X"]
    b0,_,bw=cfg["bb_over_ZZ_SM_normalized"]
    t0,_,tw=prev["cms_tau_over_z_sm"]
    eta_lo,eta_hi=prev["rate_transfer"]
    mh=prev["higgs_mass_GeV"]
    bt=1-(2*prev["tau_mass_GeV"]/mh)**2
    b_fast=1-(2*min(cfg["bottom_kin_range_GeV"])/mh)**2
    b_slow=1-(2*max(cfg["bottom_kin_range_GeV"])/mh)**2
    boxes={}
    for n in range(3):
        b,t=(b0+n*bw)/eta_lo,(t0+n*tw)/eta_lo
        cap=(math.sqrt(b)+math.sqrt(t))
        row=dict(sqrt_rate_sum_upper=cap,required_transfer_or_lower=
                 ((math.sqrt(b0+n*bw)+math.sqrt(t0+n*tw))/2)**2)
        if cap>=2:
            assert t<4-b
            cos_lower=(4+t-b)/(4*math.sqrt(t))
            tan2=min(bt*(t-1),b_fast*(b-1))
            row.update(common_odd_angle_outer_deg=math.degrees(math.acos(cos_lower)),
                       full_vertex_angle_outer_deg=math.degrees(math.atan(math.sqrt(tan2))),
                       minimum_VV_coupling_squared=1/(1+tan2))
        boxes[str(n)]=row
    witness=cfg["shared_odd_witness_ratio"]
    rates={"bottom":1+witness**2/b_slow,"tau":1+witness**2/bt}
    for name,(mid,lo,hi) in {"bottom":cfg["bb_over_ZZ_SM_normalized"],
                             "tau":prev["cms_tau_over_z_sm"]}.items():
        assert mid-2*lo <= rates[name]*eta_lo <= mid+2*hi
    spin=math.degrees(math.atan(witness))
    direct,err=cfg["holdout_ATLAS_tau_spin_deg"]
    assert abs(spin-direct)<err and boxes["1"]["sqrt_rate_sum_upper"]<2
    # Check the finite mass identity and circle bound on explicit complex projections.
    for u in np.linspace(.85,1.1,15):
        for w in np.linspace(-.3,.3,23):
            t=u*u+w*w/bt
            b=(2-u)**2+w*w/b_slow
            np.testing.assert_allclose(b-t,4*(1-u)+w*w*(1/b_slow-1/bt),atol=1e-14)
            assert math.sqrt(b)+math.sqrt(t)>=2-1e-12
    return dict(types=data,rate_boxes=boxes,full_witness=dict(w=witness,rates=rates,angle_deg=spin),
                spin_holdout_deg=[direct,err])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--check",action="store_true")
    g.add_argument("--write-new",action="store_true")
    opt=p.parse_args()
    result=evaluate()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a:
                    compare(a[k],b[k])
            else:
                np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-15)
        compare(result,read("results.json"))
        print("1536: four assignments, determinant test, two-rate upper envelope and withheld tau spin passed.")
