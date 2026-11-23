"""Compare protected PMNS columns with published CP intervals on the circle."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PI = math.pi


def read(name):
    return json.loads((HERE/name).read_bytes())


def intersection(periodic_interval, lo, hi):
    a,b = periodic_interval
    overlaps = [(max(lo,a+2*i*PI),min(hi,b+2*i*PI)) for i in range(-2,3)]
    return [[x,y] for x,y in overlaps if x<=y]


def upper_mu(solar, reactor, atm, cosine):
    return (solar*(1-atm)+(1-solar)*reactor*atm
            +2*math.sqrt(solar*(1-solar)*reactor*atm*(1-atm))*cosine)


def decreasing_cut(f, lo, hi):
    assert f(lo)>0 and f(hi)<0
    for _ in range(65):
        mid=(lo+hi)/2
        if f(mid)>0:
            lo=mid
        else:
            hi=mid
    return hi


def compute():
    cfg=read("inputs.json")
    orig=read(cfg["port_results"])["cases"]
    reactor,solar,atmos=read(cfg["angle_windows"])["joint_windows"]
    vlo,vhi=[(1-math.sqrt(1-a))/2 for a in reactor]
    joint=cfg["joint_2025"]
    peak,minus,plus=joint["hpd1_center_minus_plus_in_pi"]
    # Two rounded input digits enter each derived endpoint.
    j1wide=[PI*(peak-minus-.01),PI*(peak+plus+.01)]
    j3core=[PI*(joint["hpd3_endpoints_in_pi"][0]+.005),PI*(joint["hpd3_endpoints_in_pi"][1]-.005)]
    current=cfg["t2k_updated"]
    p,m,q=current["fc1_center_minus_plus_rad"]
    fc1wide=[p-m-.01,p+q+.01]
    prev=cfg["t2k_2023"]
    rnd=prev["printed_endpoint_half_unit_rad"]
    fc90hi=prev["fc90_rad"][1]+rnd
    fc90core=[prev["fc90_rad"][0]+rnd,prev["fc90_rad"][1]-rnd]
    fc2hi=prev["fc2_near_zero_rad"][1]+rnd
    fc3hi=prev["fc3_near_zero_rad"][1]+rnd
    rows={}
    for label,point in orig.items():
        if not point["solar_reactor_compatible"]:
            continue
        probs=point["protected_probability_envelopes"]
        # Outward slack for decimal endpoints of the archived envelopes.
        pmlo=probs[1][0]-2e-9
        xhigh=min(solar[1],1-(probs[0][0]-2e-9)/(1-vlo))
        zlo,zhi=atmos
        coeff=2*math.sqrt(xhigh*(1-xhigh)*vhi*zlo*(1-zlo))
        cosine=(pmlo-xhigh*(1-zlo)-(1-xhigh)*vhi*zlo)/coeff
        # Upper bound on dp_mu/dz for positive cos, all x>solar_lo, v<vhi.
        slope=-solar[0]+vhi+math.sqrt(.25*vhi)*(1-2*zlo)/math.sqrt(zlo*(1-zlo))
        assert slope<0 and 1>cosine>point["necessary_cos_delta_lower"]>0
        half=math.acos(cosine)+1e-9
        assert intersection(j1wide,-half,half)==[]
        assert intersection(j3core,-half,half)==[[-half,half]]
        assert intersection(fc1wide,-half,half)==[]
        assert -half<fc90hi<fc2hi<0<fc3hi<half
        upper_atm=decreasing_cut(lambda z:upper_mu(xhigh,vhi,z,math.cos(fc90hi))-pmlo,zlo,zhi)
        # A real unitary PMNS construction at the archived central column (no CP fit).
        pe,pm,pt=point["p_MZ"]
        lower_v=max(vlo,1-pe/(1-solar[1]))
        upper_v=min(vhi,1-pe/(1-solar[0]))
        v=(upper_v+lower_v)/2
        x=1-pe/(1-v)
        assert vlo<=v<=vhi and solar[0]<=x<=solar[1]
        delta=cfg["witness_delta_rad"]
        z=decreasing_cut(lambda z:upper_mu(x,v,z,math.cos(delta))-pm,zlo,zhi)
        assert fc90core[0]<delta<fc90core[1] and z<=upper_atm
        assert abs((1-x)*(1-v)-pe)<1e-13 and abs(upper_mu(x,v,z,math.cos(delta))-pm)<1e-13
        assert abs(1-pe-pm-pt)<3e-9
        rows[label]=dict(min_cos_from_column=cosine, necessary_abs_delta_deg=half*180/PI,
                        joint_normal_HPD1_intersection_empty=True,joint_normal_HPD3_contains_entire_arc=True,
                        updated_T2K_normal_FC1_intersection_empty=True,
                        t2k_2023_FC90_outer_rad=[-half,fc90hi],
                        t2k_2023_FC2_outer_rad=[-half,fc2hi],
                        t2k_2023_FC3_outer_rad=[-half,fc3hi],
                        t2k_2023_FC90_and_column_s23sq_upper=upper_atm,
                        exact_kinematic_witness_at_delta_minus_point55=[x,v,z])
    assert len(rows)==6
    return dict(cases=rows, number_of_initial_ports=len(orig),
                joint_NO_phase_posterior_mass_upper=1-joint["hpd1_mass"],
                largest_remaining_abs_delta_degrees=max(r["necessary_abs_delta_deg"] for r in rows.values()),
                worst_updated_FC1_endpoint_gap_rad=abs(fc1wide[1])-max(r["necessary_abs_delta_deg"]*PI/180 for r in rows.values()))


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,float) else x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=tidy(compute())
    if p.parse_args().check:
        assert result==read("results.json")
        print("1705: six circular CP sets, overlapping experimental intervals and charged-current witnesses verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
