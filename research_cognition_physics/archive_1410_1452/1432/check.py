"""Joint image of the intrinsic neutrino 13 element and solar angle for one charged 12 rotation."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(filename):
    return json.loads((HERE/filename).read_text(encoding="utf-8"))


def box(row,n):
    return [row[0]-n*row[1], row[0]+n*(row[2] if len(row)>2 else row[1])]


def ratio(a,b):
    return [a[0]/b[1],a[1]/b[0]]


def compute():
    cfg=load("inputs.json")
    low=load(cfg["charged_low"])
    reactor=load(cfg["reactors"])
    tree=load(cfg["charged_tree"])
    leps=low["yukawa"]["charged_leptons"]
    juno=reactor["juno_primary"]
    db=reactor["daya_bay"]
    out={}
    for n in cfg["widths"]:
        re=ratio(box(leps[0],n),box(leps[1],n))
        q=ratio(box(tree["K_product"],n),box(tree["K_form"],n))
        dn=low["yukawa"]["down"]
        rd=ratio(box(dn[0],n),box(dn[1],n))
        t_down=[math.atan(math.sqrt(v)) for v in rd]
        up=low["yukawa"]["up"]
        ru=ratio(box(up[0],n),box(up[1],n))
        t_up=[math.atan(math.sqrt(v)) for v in ru]
        z=[(1-math.sqrt(1-v))/2 for v in box([db["sin_sq_2theta13"],db["one_sigma"]],n)]
        y=box(low["NO"]["atmospheric_sine_squared"],n)
        solar=box([juno["sin_sq_theta12"],juno["one_sigma"]],n)
        t12=[math.asin(math.sqrt(v)) for v in solar]
        all_points=[]
        for th,x,k in itertools.product([math.atan(math.sqrt(v)) for v in re],z,y):
            ce,se=math.cos(th),math.sin(th)
            a=ce*math.sqrt(x)
            b=se*math.sqrt(k*(1-x))
            assert a>b and k>re[1]*x/(1-x)
            intrinsic_low=a-b
            intrinsic_high=a+b
            solar_radius=math.atan(se*math.sqrt(1-k)/(ce*math.sqrt(1-x)-se*math.sqrt(k*x)))
            needed=math.atan(math.sqrt(x/(k*(1-x))))
            all_points.append((intrinsic_low,intrinsic_high,solar_radius,needed,needed-th))
        radius=max(row[2] for row in all_points)
        arc=[t12[0]-radius,t12[1]+radius]
        r=dict(
            observed_Cabibbo_sine=q,
            down_only_sine=[math.sin(v) for v in t_down],
            two_quark_rotations_sine_envelope=[math.sin(t_down[0]-t_up[1]),math.sin(t_down[1]+t_up[1])],
            charged_rotation_deg=[math.degrees(math.atan(math.sqrt(v))) for v in re],
            intrinsic_s13=[min(r[0] for r in all_points),max(r[1] for r in all_points)],
            intrinsic_theta12_deg=list(map(math.degrees,arc)),largest_solar_angle_change_deg=math.degrees(radius),
            intrinsic_zero13_min_charged_angle_deg=math.degrees(min(v[3] for v in all_points)),
            intrinsic_zero13_extra_charged_angle_deg=math.degrees(min(v[4] for v in all_points)),
            required_mass_ratio_for_zero13_lower=z[0]/(y[1]*(1-z[0])),
            measured_mass_ratio=re,
            bimaximal_gap_deg=45-math.degrees(arc[1]))
        if n==0:
            td,tu=t_down[0],t_up[0]
            a,b=math.cos(tu)*math.sin(td),math.sin(tu)*math.cos(td)
            r["quark_phase_cos_for_cabibbo"]=(a*a+b*b-q[0]**2)/(2*a*b)
            # A joint witness: preserve the exact TB solar angle, permit intrinsic 13.
            te=math.atan(math.sqrt(re[0]))
            ce,se=math.cos(te),math.sin(te)
            c13,s13=math.sqrt(1-z[0]),math.sqrt(z[0])
            c23,s23=math.sqrt(1-y[0]),math.sqrt(y[0])
            c12,s12=math.cos(t12[0]),math.sin(t12[0])
            # Traditional solar sum rule applies on the separately assumed intrinsic-13 zero slice.
            if_zero={}
            for name,tan_sq in (("TB_solar",.5),("bimaximal",1.)):
                num=s12**2*s23**2+c12**2*c23**2*s13**2-tan_sq*(c12**2*s23**2+s12**2*c23**2*s13**2)
                den=2*s12*c12*s23*c23*s13*(1+tan_sq)
                if_zero[name]=num/den
            r["cos_delta_if_zero13_and_mass_angle_released"]=if_zero
            A=ce*c13+se*s23*s13
            B=se*c23
            ph=math.acos(((A*A+B*B)/3 - s12*s12*A*A-c12*c12*B*B)/(2*s12*c12*A*B))
            de=math.pi-ph
            em=cmath.exp(1j*ph)
            ep=cmath.exp(1j*de)
            row=[ce*c12*c13+se*em*(-s12*c23-c12*s23*s13*ep),
                 ce*s12*c13+se*em*(c12*c23-s12*s23*s13*ep),
                 ce*s13*ep.conjugate()+se*em*s23*c13]
            assert abs(sum(abs(v)**2 for v in row)-1)<1e-12
            assert abs(abs(row[1])**2/(1-abs(row[2])**2)-1/3)<1e-12
            assert abs(abs(row[2])-r["intrinsic_s13"][0])<1e-12
            r["TB_solar_only_witness"]={"physical_dirac_deg":math.degrees(de),
                                         "charged_exchange_phase_deg":math.degrees(ph),
                                         "intrinsic_s13":abs(row[2])}
        out[str(n)]=r
    return out


def near(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(near(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(near(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-11,abs_tol=2e-14)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    flags=p.add_mutually_exclusive_group()
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=compute()
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert near(result,load("results.json"))
    print(json.dumps(result))
