"""Integer-certified interval comparison with an arbitrary singlet threshold menu."""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, pi, sqrt
from pathlib import Path


DIR = Path(__file__).resolve().parent
def load(file):
    return json.loads(file.read_text(encoding="utf-8"), parse_float=F)

P = load(DIR / "inputs.json")
Q = int(P["integer_scale"])
N = P["n_cells"]
D = 2*N
PI_L, PI_H = F("3.14159"), F("3.14160")


def up(a,b):
    return -((-a)//b)


def ceil_scaled(f):
    return up(Q*f.numerator, f.denominator)


def root_scaled(f):
    k=isqrt((Q*Q*f.numerator)//f.denominator)
    return k if k*k*f.denominator==Q*Q*f.numerator else k+1


def endpoint(row, sign):
    return row[0]+sign*2*sum(row[1:])


def initial():
    a=load(DIR/P["higgs_and_weak"])
    g=load(DIR/P["gauge"])
    y=load(DIR/P["charged"])["yukawa"]
    B=8*PI_H**2*endpoint(a["higgs_mass_GeV"],1)**2/a["higgs_transfer_lower"]
    C=2*a["eta_trace_upper"]/a["v_over_sqrt2_GeV"]**2
    su=sum(endpoint(i,1)**2 for i in y["up"][:-1])
    sd=sum(endpoint(i,1)**2 for i in y["down"])
    se=sum(endpoint(i,1)**2 for i in y["charged_leptons"])
    return B,C,g,ceil_scaled(endpoint(y["up"][-1],1)**2),ceil_scaled(3*(su+sd)+se)


def exponential_box(z):
    term=total=F(1)
    for i in range(1,151):
        term=term*z/i
        total+=term
    return total, total+term*z/151/(1-z/152)


def verify_time_window(mz):
    target=F(P["limit_GeV"])
    assert mz*exponential_box(8*PI_H**2*F(49,100))[1] < target
    assert mz*exponential_box(8*PI_L**2/2)[0] > target
    assert mz*exponential_box(8*PI_L**2*F(41,100))[0] > 10**16
    # Lower bound on log(2) from the first two positive atanh-series terms.
    assert 8*PI_H**2*F(7,800) < 2*(F(1,3)+F(1,81))


@lru_cache(maxsize=1)
def controls():
    _,_,g,_,_=initial()
    mz=g["mu_ref_GeV"]
    verify_time_window(mz)
    a=[endpoint(g[key],-1)**2 for key in ("g1_GUT","g2","g3")]
    # g_i^2(x)=g_i^2(0)/(1-b_i g_i^2(0) x).
    beta=[F(41,10),F(-19,6),F(-7)]
    result=[]
    for j in range(N):
        lower=[]
        for v,b in zip(a,beta):
            point=j if b>0 else j+1
            num=v.numerator*b.denominator*D
            den=v.denominator*b.denominator*D-b.numerator*v.numerator*point
            assert den>0
            lower.append((Q*num)//den)
        g1,g2,g3=lower
        gt=(17*g1+45*g2+160*g3)//20
        gn=(9*g1+45*g2)//20
        dyadic_slot=up(800*(j+1),7*D)
        result.append((gt,gn,dyadic_slot))
    return result


def end_bound(start, rate):
    # exp(z) <= 1/(1-z) for signed z<1. One exact outward rounding.
    den=D*Q-rate
    if den<=0:
        raise ArithmeticError("cell too large")
    return up(start*D*Q,den)


def cell(u,a,s,gt,gn):
    upper=[2*u,2*a,2*s]
    for _ in range(3):
        uu,aa,ss=upper
        ru=up(9*uu,2)+aa+ss-gt
        ra=max(0,up(5*aa,2)+3*uu+ss-gn)
        rs=3*uu+aa+up(5*ss,2)
        out=[end_bound(u,ru),end_bound(a,ra),end_bound(s,rs)]
        new=[max(u,out[0]),out[1],out[2]]
        # Strict first pass establishes an invariant cell. Subsequent refinements
        # may coincide exactly after rounding.
        if any(v>old for v,old in zip(new,upper)):
            raise ArithmeticError("positive invariant box not closed")
        upper=new
    return out,upper


def envelope(Bscale,Cscale,calendar):
    B,C,g,u,s=initial()
    B*=F(Bscale)
    C*=F(Cscale)
    w=root_scaled(B*C)
    mz=g["mu_ref_GeV"]
    slots={j:min(w,ceil_scaled(C*mz**2*4**j)) for j in range(1,60)}
    e=0
    umax,amax,smax=u,0,s
    near_top=0
    integrated=[0,0,0]
    for j,(gt,gn,slot) in enumerate(controls()):
        cap=slots[slot] if calendar else w
        a=e+cap
        (un,an,sn),(ucell,acell,scell)=cell(u,a,s,gt,gn)
        e+=an-a
        umax,amax,smax=max(umax,ucell),max(amax,acell),max(smax,scell)
        if (j+1)*100 >= 49*D:
            near_top=max(near_top,ucell)
        if 100*j < 41*D:
            for i,v in enumerate((ucell,acell,scell)):
                integrated[i]+=v
        u,s=un,sn
    iu,ia,iss=[F(v,Q*D) for v in integrated]
    y=load(DIR/P["charged"])["yukawa"]
    d0=(endpoint(y["down"][0],-1)*endpoint(y["charged_leptons"][1],-1)/
        (endpoint(y["down"][1],1)*endpoint(y["charged_leptons"][0],1)))
    exp_hi=exponential_box(3*(iu+ia+iss/3)/4)[1]
    # All summed cells are upper bounds through 10^16 GeV, whichever singlets are active.
    high_D=F((10**9*d0/exp_hi).__floor__(),10**9)
    return dict(weight_cap=w/Q,all_top_squared=umax/Q,all_neutral_trace=amax/Q,
                all_other_weighted_squared=smax/Q,top_squared_UV_slab=near_top/Q,
                integrals_dx_to_1e16=[float(iu),float(ia),float(iss)],
                lower_mass_double_ratio_1e16=float(high_D),
                comparison_at_x_half=[u/Q,(w+e)/Q,s/Q])


def matching_example():
    """Two degenerate singlets match two light splittings and the central Higgs term."""
    info=load(DIR/P["higgs_and_weak"])
    nu=load(DIR/"../1370/inputs.json")["NO"]
    msum=(sqrt(nu["solar"][0])+sqrt(nu["atm_abs"][0]))*1.e-9
    u=float(info["v_over_sqrt2_GeV"])
    h=float(info["higgs_mass_GeV"][0])
    M=(8*pi*pi*h*h*u*u/msum)**(1/3)
    strength=M*msum/(u*u)
    assert abs(M*M*strength/(8*pi*pi*h*h)-1)<1.e-12
    data={"M_GeV":M,"sum_gate_weight":strength,"weak_trace":msum/(2*M)}
    return {k:float(f"{v:.10g}") for k,v in data.items()}


def calculate():
    B,C,g,u,s=initial()
    cases={"both_nominal":("1","1"),"higgs_137pct":("1.37","1"),
           "weak_135p5pct":("1","1.355"),"both_116pct":("1.16","1.16")}
    output={key:envelope(x,y,True) for key,(x,y) in cases.items()}
    output["no_calendar_nominal"]=envelope("1","1",False)
    assert all(x["all_top_squared"]<1 and x["all_neutral_trace"]<4*float(PI_L)
               for x in output.values())
    return dict(round=1587, B_outward_GeV2=float(B),C_inverse_GeV2=float(C),
                initial_top_squared=u/Q, initial_other_weighted_squared=s/Q,
                scenarios=output,nonempty_matching_example=matching_example())


def main():
    a=argparse.ArgumentParser()
    group=a.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=a.parse_args()
    out=calculate()
    path=DIR/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as h:
            json.dump(out,h,indent=2)
            h.write("\n")
    elif json.loads(path.read_text(encoding="utf-8"))!=out:
        raise SystemExit("1587 replay mismatch")
    print("1587 certified integer bounds verified")


if __name__=="__main__":
    main()
