"""Absolute-mass and Majorana-phase tests of four published normal-ordering points."""
import argparse
import cmath
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
# Dueck & Rodejohann, table 13. (scenario, m0 in meV, m_bb in meV).
POINTS = [("MN_no_RGE",3.03,.35),("MN_RGE",2.40,.49),
          ("MS_no_RGE",2.58,.38),("MS_RGE",6.83,.44)]


def read(path):
    return json.loads((HERE.parent/path).read_text(encoding="utf-8-sig"))


def interval(x, width=2):
    return x[0]-width*x[1], x[0]+width*x[-1]


def small_root(amplitude):
    return (1-math.sqrt(1-amplitude))/2


def phase_witness(a,b,c,t):
    low,high=max(abs(a-b),abs(c-t)),min(a+b,c+t)
    assert low<=high
    r=(low+high)/2
    alpha=math.acos((r*r-a*a-b*b)/(2*a*b))
    two=a+b*cmath.exp(1j*alpha)
    third=cmath.phase(two)+math.acos((t*t-r*r-c*c)/(2*r*c))
    assert abs(abs(two+c*cmath.exp(1j*third))-t)<2e-13
    return [math.degrees(alpha),math.degrees(third)%360]


def solution():
    archived=read("1484/inputs.json")
    juno=read("1446/inputs.json")["juno_207_days"]
    amplitude=archived["DayaBay_sin2_2theta13"]
    x0,x1=interval(archived["JUNO_sin2_theta12_NO"])
    z0,z1=map(small_root,interval(amplitude))
    s0,s1=(q*1e6 for q in interval(juno["solar_21"]))
    d0,d1=(q*1e6 for q in interval(juno["normal_atmospheric_31"]))
    b0=(1-z0)*x0*s0+z0*d0
    b1=(1-z1)*x1*s1+z1*d1
    x=archived["JUNO_sin2_theta12_NO"][0]
    z=small_root(amplitude[0])
    b_mid=(1-z)*x*juno["solar_21"][0]*1e6+z*juno["normal_atmospheric_31"][0]*1e6
    rows={}
    nuc=read("1370/inputs.json")["kamland_zen"]
    for key,m,t in POINTS:
        ml,mh=m-.005,m+.005
        tl,th=t-.005,t+.005
        beta=(math.sqrt(ml*ml+b0),math.sqrt(mh*mh+b1))
        a=((1-x1)*(1-z1)*ml,(1-x0)*(1-z0)*mh)
        b=(x0*(1-z1)*math.sqrt(ml**2+s0),x1*(1-z0)*math.sqrt(mh**2+s1))
        c=(z0*math.sqrt(ml**2+d0),z1*math.sqrt(mh**2+d1))
        min_ab_difference=max(0.,a[0]-b[1],b[0]-a[1])
        cos_upper=min(1.,-1.+((c[1]+th)**2-min_ab_difference**2)/(2*a[0]*b[0]))
        a_mid=(1-x)*(1-z)*m
        bmid=x*(1-z)*math.sqrt(m*m+juno["solar_21"][0]*1e6)
        cmid=z*math.sqrt(m*m+juno["normal_atmospheric_31"][0]*1e6)
        witness=phase_witness(a_mid,bmid,cmid,t)
        # NME-specific sensitivities required at the printed point value.
        half_lives=[nuc["half_life_lower_years"]*(q*1000/t)**2
                    for q in nuc["m_bb_upper_eV_endpoints_by_nme"]]
        rows[key]={"m0_meV":m,"m_bb_meV":t,"m_beta_interval_meV":beta,
                   "m_beta_central_meV":math.sqrt(m*m+b_mid),
                   "abs_alpha21_min_deg":math.degrees(math.acos(cos_upper)),
                   "central_effective_alpha21_alpha31minus2delta_deg":witness,
                   "Xe_half_life_for_printed_value_years":half_lives}
    mb_sorted=sorted(q[2] for q in POINTS)
    ratio=min((high-.005)/(low+.005) for low,high in zip(mb_sorted,mb_sorted[1:]))
    m0_sorted=sorted(q[1] for q in POINTS)
    gap=min((high-.005)**2-(low+.005)**2 for low,high in zip(m0_sorted,m0_sorted[1:]))
    first_three_upper=max(rows[q[0]]["m_beta_interval_meV"][1] for q in POINTS[:3])
    robust=rows["MS_RGE"]["m_beta_interval_meV"][0]-first_three_upper
    assert robust>0 and gap>0
    return {"points":rows,"m_beta_squared_common_offset_interval_meV2":[b0,b1],
            "m_beta_gap_fourth_to_three_meV":robust,
            "all_four_m0_squared_min_separation_meV2":gap,
            "all_four_amplitude_ratio_max_to_min_for_separation":ratio,
            "common_xe_amplitude_symmetric_rel_halfwidth_max":(ratio-1)/(ratio+1)}


def compact(x):
    if isinstance(x,dict):
        return {k:compact(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [compact(v) for v in x]
    if isinstance(x,float):
        return float(f"{x:.9g}")
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compact(solution())
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif args.check:
        assert result==json.loads(path.read_text(encoding="utf-8"))
        print("1555 independent: four absolute-mass and isotope/phase certificates verified")
    else:
        print(json.dumps(result,indent=2))
