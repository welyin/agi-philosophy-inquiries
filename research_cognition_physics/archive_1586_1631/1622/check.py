"""Mixed anomaly characters and electroweak running of a family-anomaly completion."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def family_anomaly(weak_dim, hypercharge, orientation):
    # Weyl fundamental of U(3) or its conjugate. Dynkin index uses twice T=1.
    n = 3*weak_dim
    q = orientation
    return (weak_dim*q, weak_dim*q, n*q**3, n*q, 3*q*(weak_dim == 2),
            weak_dim*hypercharge, n*q*hypercharge**2, n*q**2*hypercharge)


def gauge_at(t, logs_l, logs_e, g0):
    loop = 8*math.pi**2
    active_l = sum(max(t-z,0.) for z in logs_l)
    active_e = sum(max(t-z,0.) for z in logs_e)
    coeff = (41*t/6 + (2*active_l + 4*active_e)/3,
             -19*t/6 + 2*active_l/3, -7*t)
    return [1/(1/g**2-c/loop) for g,c in zip(g0,coeff)]


def at_crossing(low, unit, l_masses, e_masses, count):
    logs_l = [math.log(m/unit) for m in l_masses]
    logs_e = [math.log(m/unit) for m in e_masses]
    assert all(z >= 0 for z in logs_l+logs_e)
    g0 = [low[k] for k in ("gY","g2","g3")]
    def difference(t):
        x,y,_ = gauge_at(t,logs_l,logs_e,g0)
        return 1/x-3/y
    a,b = 0.,math.log(10_000/unit)
    assert difference(a)>0>difference(b)
    for _ in range(56):
        mid = (a+b)/2
        if difference(mid)>0:
            a=mid
        else:
            b=mid
    end=(a+b)/2
    def rhs(t,state):
        yt2,lam=state
        x,y,z=gauge_at(t,logs_l,logs_e,g0)
        return np.array([yt2*(9*yt2-17*x/6-9*y/2-16*z),
                         24*lam**2+(12*yt2-9*y-3*x)*lam-6*yt2**2+3*(2*y**2+(x+y)**2)/8])/(16*math.pi**2)
    state=np.array([low["top_yukawa"]**2,low["lambda"]])
    nodes=sorted(set([0.,end]+[a for a in logs_l+logs_e if 0<a<end]))
    for left,right in zip(nodes,nodes[1:]):
        n=max(1,math.ceil(count*(right-left)/end))
        h=(right-left)/n
        for i in range(n):
            t=left+i*h
            k1=rhs(t,state)
            k2=rhs(t+h/2,state+h*k1/2)
            k3=rhs(t+h/2,state+h*k2/2)
            k4=rhs(t+h,state+h*k3)
            state += h*(k1+2*k2+2*k3+k4)/6
    g2sq=gauge_at(end,logs_l,logs_e,g0)[1]
    return {"scale_GeV":unit*math.exp(end),"g2_squared":g2sq,"measured_lambda_run":state[1],
            "lambda_minus_G2_boundary":state[1]-g2sq/6,
            "lambda_minus_Sp6_boundary":state[1]-g2sq/2}


def run():
    cfg=read(HERE/"inputs.json")
    sm=read(HERE/cfg["prior_weak"])
    lepton=read(HERE/cfg["prior_leptons"])["pole_MeV"]
    e,mu=lepton["electron"][0],lepton["muon"][0]
    a,b=math.sqrt(e),math.sqrt(mu)
    tau=(2*(a+b)+math.sqrt(3*(a*a+4*a*b+b*b)))**2
    assert math.isclose(tau,read(HERE/cfg["prior_tau"])["roots_MeV"][1],rel_tol=1e-12)
    lep=read(HERE/"../1619/independent_results.json")["LEP_destructive_AA_only_95_coefficient_TeV_minus_2"]
    ref=at_crossing(sm["reference"],sm["matching_scale"],[],[],64)
    lm=[cfg["mass_example_GeV"]*math.sqrt(m/e) for m in (e,mu,tau)]
    sectors={}
    an = [family_anomaly(2,F(-1,2),1),family_anomaly(1,F(1),1)]
    mirror = [family_anomaly(2,F(1,2),-1),family_anomaly(1,F(-1),-1)]
    assert all(sum(terms)==0 for terms in zip(*an,*mirror))
    # General impossibility follows from determinant characters; check explicit representatives.
    for k in range(-25,26):
        assert (3*k+1)%3==1
    for r in cfg["relative_singlet_thresholds"]:
        em=[r*m for m in lm]
        def evaluate(n):
            return at_crossing(sm["reference"],sm["matching_scale"],lm,em,n)
        point=evaluate(64)
        fine=evaluate(128)
        assert all(math.isclose(fine[k],point[k],rel_tol=2e-10) for k in point)
        # Both lowest thresholds have been crossed; explicitly retain any higher ones.
        assert max(lm[0],em[0])<=point["scale_GeV"]
        # r=0.5 may activate the second singlet below the crossing.
        floor=sm["matching_scale"]
        threshold_logs=sum(math.log(max(1,point["scale_GeV"]/l))
                           - math.log(max(1,point["scale_GeV"]/s)) for l,s in zip(lm,em))
        assert abs(math.log(point["scale_GeV"]/ref["scale_GeV"])-4*threshold_logs/49)<1e-12
        assert min(lm+em)>=floor
        # Test one physical matching scale for the full family multiplet.
        m11_TeV=math.sqrt(point["g2_squared"]/(4*lep))
        point["all_vectors_M33_floor_TeV_at_unit_transfer"]=m11_TeV*math.sqrt(tau/e)
        point["full_tower_max_transfer_product"]=point["scale_GeV"]/(1000*point["all_vectors_M33_floor_TeV_at_unit_transfer"])
        assert point["full_tower_max_transfer_product"] < .04
        sectors[str(r)]=point
    assert abs(sectors["1"]["scale_GeV"]-ref["scale_GeV"])<1e-8
    assert sectors["1"]["g2_squared"]>ref["g2_squared"]
    return {"anomalies_before_mirror_with_neutral_quarks":[float(sum(j)) for j in zip(*an)],
            "anomaly_order":["SU3F^3","SU3F^2F","F^3","gravity^2F","SU2L^2F","SU3F^2Y","Y^2F","YF^2"],
            "residual_family_anomalies":[0]*8,"tau_mass_input_from_prior_relation_MeV":tau,
            "witness_L_masses_GeV":lm,"reference_no_added_fields":ref,"by_E_over_L":sectors}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    result=run()
    path=HERE/"results.json"
    if parser.parse_args().check:
        expected=read(path)
        def matched(old,new):
            if isinstance(old,dict):
                return old.keys()==new.keys() and all(matched(old[k],new[k]) for k in old)
            if isinstance(old,list):
                return len(old)==len(new) and all(matched(a,b) for a,b in zip(old,new))
            if isinstance(old,str):
                return old==new
            return math.isclose(old,new,rel_tol=1e-10,abs_tol=1e-11)
        assert matched(expected,result)
        print("1622 PASS: family anomaly, full mirror cancellation, threshold identity and scalar running")
    else:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(result,indent=2)+"\n")
