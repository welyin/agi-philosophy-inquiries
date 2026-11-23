"""Test one common reciprocal right-handed current using the normalized vector and axial rates."""
import argparse
import itertools
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"independent_results.json"


def read(n,name):
    return json.loads((ROOT.parent/str(n)/name).read_bytes())


LIGHT=read(1613,"inputs.json")
CKM=read(1702,"inputs.json")["direct_semileptonic_moduli"]


def limits(x,t):
    e=math.hypot(*x[1:])
    return (x[0]-t*e,x[0]+t*e)


def rectangle(t,pion):
    d=limits([.9739,.0029] if pion else LIGHT["new_nuclear"],t)
    s=limits(CKM["Vus"],t)
    b=limits(CKM["Vub"],t)
    l=LIGHT["leptonic"]
    raw=limits([l["combination"],l["branching_error"],l["radiative_error"]],t)
    A=(raw[0]/(l["decay_constant_ratio"]+t*l["upper"]),
       raw[1]/(l["decay_constant_ratio"]-t*l["lower"]))
    return d,s,A,b


def free_geometry(d,s,A,B,sign):
    t=math.sqrt(A*d/s)  # |1-lambda|/|1+lambda|
    a=math.hypot(d,s*t,B*(1 if sign<0 else t))
    b=a/t
    return (a+b)/2,a,b,[d/a,s/b,B/(a if sign<0 else b)]


def certify_full_box(t,pion,sign):
    d,s,A,B=rectangle(t,pion)
    # For the light two contributions increasing d and s increases the norm.
    # At fixed d,s, q=sqrt(A*s/d) has q^3<s/d, so the light norm decreases in A.
    assert (A[1]*s[1]/d[0])**1.5<s[0]/d[1]
    q=math.sqrt(A[0]*s[1]/d[1])
    light=math.hypot(d[1]+s[1]/q,s[1]+d[1]*q)/2
    tlo=math.sqrt(A[0]*d[0]/s[1])
    thi=math.sqrt(A[1]*d[1]/s[0])
    bterm=B[1]*(1+(1/tlo if sign<0 else thi))/2
    return math.hypot(light,bterm)


def free_witness(pion,sign):
    def eval_path(k):
        d,s,A,B=rectangle(k,pion)
        return free_geometry(d[1],s[1],A[0],B[1],sign)
    lo,hi=0.,4.
    assert eval_path(lo)[0]<1<eval_path(hi)[0]
    for _ in range(70):
        m=(lo+hi)/2
        if eval_path(m)[0]<1:lo=m
        else:hi=m
    _,a,b,moduli=eval_path(hi)
    lam=(b-a)/2
    d,s,A,B=rectangle(hi,pion)
    npred=observables(moduli,lam,sign)
    assert max(abs(i-j) for i,j in zip(npred,[d[1],s[1],A[0],B[1]]))<1e-12
    return {"reachable_box_scale_on_explicit_path":hi,"real_right_strength":lam,"unitary_row":moduli}


def observables(v,z,sign):
    a,b=abs(1-z),abs(1+z)
    return [v[0]*a,v[1]*b,v[1]/v[0]*a/b,v[2]*(a if sign<0 else b)]


def fixed_strength(t,pion,sign,v):
    d,s,A,B=rectangle(t,pion)
    dlo,dhi=(x/v[0] for x in d)
    slo,shi=(x/v[1] for x in s)
    if sign<0:
        dlo=max(dlo,B[0]/v[2])
        dhi=min(dhi,B[1]/v[2])
    else:
        slo=max(slo,B[0]/v[2])
        shi=min(shi,B[1]/v[2])
    r=v[1]/v[0]
    tlo,thi=A[0]/r,A[1]/r
    if dhi<dlo or shi<slo or dhi<tlo*slo or dlo>thi*shi:
        return None
    a=min(dhi,shi*thi)
    b=min(shi,dhi/tlo)
    if a+b<2 or abs(a-b)>2:
        return None
    re=(b*b-a*a)/4
    im=math.sqrt(max(0,(a*a+b*b)/2-1-re**2))
    z=complex(re,im)
    for y,(lo,hi) in zip(observables(v,z,sign),(d,s,A,B)):
        assert lo-2e-12<=y<=hi+2e-12
    return [re,im]


def fixed_threshold(v,pion,sign):
    lo,hi=0.,7.
    assert fixed_strength(hi,pion,sign,v) is not None
    for _ in range(65):
        m=(lo+hi)/2
        if fixed_strength(m,pion,sign,v) is None:lo=m
        else:hi=m
    return {"minimum_nested_width":hi,"right_strength_at_onset":fixed_strength(hi,pion,sign,v)}


def calculate():
    row=read(1712,"results.json")["cases"]["interior_all_one_width"]
    u=row["quark_neutrino"]["Vus_MZ"]
    b=limits(CKM["Vub"],row["input_widths"]["Vub"])[1]
    v=[math.sqrt(1-u*u-b*b),u,b]
    results={}
    for pion,sign in itertools.product([False,True],[-1,1]):
        key=("pion_beta" if pion else "nuclear_beta")+("_b_minus" if sign<0 else "_b_plus")
        cap=certify_full_box(1.,pion,sign)
        lower,upper=0.,4.
        for _ in range(65):
            m=(lower+upper)/2
            if certify_full_box(m,pion,sign)<1:lower=m
            else:upper=m
        report=free_witness(pion,sign)
        assert lower<=report["reachable_box_scale_on_explicit_path"]+1e-10
        dbox,sbox,abox,bbox=rectangle(1.,pion)
        L,aa,bb,unit_row=free_geometry(dbox[1],sbox[1],abox[0],bbox[1],sign)
        rr=(bb-aa)/2
        assert abs((L-rr)-aa)+abs((L+rr)-bb)<1e-14
        results[key]={"one_width_upper_for_half_sum_distances":cap,
                      "one_width_real_competing_current_witness":{"left_normalization":L,
                          "right_strength":rr,"unitary_row":unit_row},
                      "rigorous_no_solution_below_nested_width":lower,
                      "free_unitary_row":report,"frozen_1712_row":fixed_threshold(v,pion,sign)}
    assert max(results[key]["one_width_upper_for_half_sum_distances"] for key in results if key.startswith("nuclear"))<.999
    d,s,A,B=[e[0] for e in rectangle(0.,False)]
    r=math.sqrt(A*s/d)
    factor=math.sqrt(A*d/s)
    return {"scope":"Single lepton-universal right-vector amplitude, no other distortions, signs d-, s+ and indicated b. Nested marginal boxes, not joint confidence or an EFT truncation.",
            "frozen_down_strange_bottom_moduli":v,"true_Cabibbo_ratio_from_light_central_observables":r,
            "b_axial_over_vector_indicated_by_light_central_ratios":{"sign_minus":1/factor,"sign_plus":factor},
            "BelleII_Btaunu_relative_rate_width":math.hypot(.41,.19)/1.24,
            "branches":results}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.check:
        assert result==json.loads(OUT.read_bytes())
        print("1713 independent reciprocal right-current bounds passed")
    else:
        OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

