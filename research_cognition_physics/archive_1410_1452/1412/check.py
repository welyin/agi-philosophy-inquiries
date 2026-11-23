"""Three-generation FX textures, tree-level flavor readouts and interval bounds."""
import argparse
import cmath
import json
from itertools import product
from math import acos, asin, atan, atan2, cos, degrees, isclose, sin, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def box(measurement,k):
    c,*errors = measurement
    return c-k*sum(errors),c+k*sum(errors)


def bound(f,*intervals):
    values = [f(*corner) for corner in product(*intervals)]
    return min(values),max(values)


def ratio(n,d):
    assert d[0]>0
    return bound(lambda a,b:a/b,n,d)


def times(*xs):
    return bound(lambda *p:float(np.prod(p)),*xs)


def rotation(t):
    c=1/sqrt(1+t*t)
    return np.array([[c,-c*t,0],[c*t,c,0],[0,0,1]])


def trig(t):
    return t/sqrt(1+t*t),1/sqrt(1+t*t)


def model(up,down,u,d,k,w):
    su,cu=trig(u)
    sd,cd=trig(d)
    s=w/cu
    c=sqrt(1-s*s)
    A,B=su*cd*c,cu*sd
    cosine=(A*A+B*B-k*k)/(2*A*B)
    assert abs(cosine)<1
    phi=acos(cosine)
    W=np.array([[cmath.exp(-1j*phi),0,0],[0,c,s],[0,-s,c]])
    Ru,Rd=rotation(u),rotation(d)
    V=Ru.T@W@Rd
    assert np.allclose(V.conj().T@V,np.eye(3))
    assert isclose(abs(V[0,1]),k) and isclose(abs(V[1,2]),w)
    assert isclose(abs(V[0,2]/V[1,2]),u) and isclose(abs(V[2,0]/V[2,1]),d)
    # Both genuine three-flavor Yukawas. These matrices allow comparisons even if u was freed.
    Du=np.diag([-up[0],*up[1:]])
    Dd=np.diag([-down[0],*down[1:]])
    Yu=Ru@Du@Ru.T
    Yd=W@Rd@Dd@Rd.T@W.conj().T
    assert np.allclose(np.linalg.svd(Yu,compute_uv=False),up[::-1])
    assert np.allclose(np.linalg.svd(Yd,compute_uv=False),down[::-1])
    assert np.allclose((W@Rd).conj().T@Yd@(W@Rd),Dd)
    J=(V[0,0]*V[1,1]*V[0,1].conj()*V[1,0].conj()).imag
    assert isclose(J,su*cu*sd*cd*s*s*c*sin(phi),rel_tol=1e-10)
    gamma=degrees(cmath.phase(-V[0,0]*V[0,2].conj()/(V[1,0]*V[1,2].conj())))
    return {"phi_degrees":degrees(phi),"J":J,"gamma_DK_mod180_for_two_CP_branches":[gamma,180-gamma],
            "up_11_over_charm":abs(Yu[0,0])/up[1],"down_11_over_strange_in_named_frame":
            abs((W.conj().T@Yd@W)[0,0])/down[1]}


def cp_bounds(u,d,k,w):
    su=bound(lambda t:trig(t)[0],u)
    cu=bound(lambda t:trig(t)[1],u)
    sd=bound(lambda t:trig(t)[0],d)
    cd=bound(lambda t:trig(t)[1],d)
    s=ratio(w,cu)
    c=bound(lambda x:sqrt(1-x*x),s)
    A,B=times(su,cd,c),times(cu,sd)
    # All four positive input factors are bounded first. The quotient denominator is positive.
    numerator=bound(lambda a,b,k_:a*a+b*b-k_*k_,A,B,k)
    co=ratio(numerator,times((2,2),A,B))
    assert -1<co[0]<co[1]<1
    min_abs_cos=0 if co[0]<=0<=co[1] else min(map(abs,co))
    sf=(sqrt(1-max(map(abs,co))**2),sqrt(1-min_abs_cos**2))
    J=times(su,cu,sd,cd,s,s,c,sf)
    # tan-angle formula for the signed real and imaginary parts of the gamma quartet.
    # In each Cartesian range below all factors are nonnegative except cos(phi).
    real1=bound(lambda u_,d_,c_:u_*(1-(d_*c_)**2),u,d,c)
    real2=times(d,c,bound(lambda x:1-x*x,u),co)
    x=(real1[0]-real2[1],real1[1]-real2[0])
    y=times(d,c,bound(lambda t:1+t*t,u),sf)
    gamma=bound(lambda h,v:degrees(atan2(v,h)),x,y)
    return {"phi_deg":list(bound(lambda x:degrees(acos(x)),co)),
            "abs_J":list(J),"gamma_plus_deg":list(gamma),
            "gamma_conjugate_mod180_deg":[180-gamma[1],180-gamma[0]]}


def departures(r,t):
    assert 0<r[0]<=r[1]<t[0]**2 and t[1]<1
    return {"minimum_up_turn_deg":degrees(atan(t[0])-atan(sqrt(r[1]))),
            "minimum_symmetric_11_over_charm":(t[0]**2-r[1])/(1+t[0]**2),
            "minimum_zero11_forward_over_reverse":(r[1]**2+t[0]**2)/(r[1]*(1+t[0]**2))}


def calculate():
    cfg=read("inputs.json")
    tree=read(cfg["tree"])
    low=read(cfg["masses_low"])
    high=read(cfg["masses_high"])["high_scale"]
    output={}
    for scale,data,angles,radian in [("MZ",low["yukawa"],low["quark_sines_and_phase"],False),
                                     ("high",high,high["angles_theta12_theta13_theta23_delta"],True)]:
        up,down=data["up"],data["down"]
        m=[a for a,_ in up]
        n=[a for a,_ in down]
        a,b=sqrt(m[0]/m[1]),sqrt(n[0]/n[1])
        ck=[x[0] for x in angles]
        s12,s13,s23=(map(sin,ck[:3]) if radian else ck[:3])
        c23=sqrt(1-s23*s23)
        c12=sqrt(1-s12*s12)
        td=abs((s12*s23-c12*c23*s13*cmath.exp(1j*ck[3]))/
               (c12*s23+s12*c23*s13*cmath.exp(1j*ck[3])))
        result={"square_root_ratios":[a,b],"CKM_same_scale_reference":[s13/(s23*sqrt(1-s13*s13)),td]}
        for mult in (1,2):
            r=ratio(box(up[0],mult),box(up[1],mult))
            t13=box(angles[1],mult)
            t23=box(angles[2],mult)
            if radian:
                t13,t23=bound(sin,t13),bound(sin,t23)
            r_ckm=bound(lambda x,y:x/(sqrt(1-x*x)*y),t13,t23)
            analog=departures(r,r_ckm)
            bounds={"sqrt_mu_mc":list(bound(sqrt,r)),"sqrt_md_ms":list(bound(sqrt,ratio(box(down[0],mult),box(down[1],mult)))),
                    "same_scale_fit_ub_cb":list(r_ckm),"same_scale_fit_minimum_cost":analog}
            if scale=="MZ":
                direct=box(tree["ub_over_cb"],mult)
                bounds["tree_ub_cb"]=list(direct)
                bounds["tree_minimum_cost"]=departures(r,direct)
                us=ratio(box(tree["K_product"],mult),box(tree["K_form"],mult))
                w=box(cfg["BelleII_Vcb_BGL"],mult)
                bounds["both_textures_CP"]=cp_bounds(bounds["sqrt_mu_mc"],bounds["sqrt_md_ms"],us,w)
                bounds["down_only_CP"]=cp_bounds(direct,bounds["sqrt_md_ms"],us,w)
            result[str(mult)]=bounds
        if scale=="MZ":
            us=tree["K_product"][0]/tree["K_form"][0]
            w=cfg["BelleII_Vcb_BGL"][0]
            result["conditional_CP_central"]=model(m,n,a,b,us,w)
            result["down_only_CP_central"]=model(m,n,tree["ub_over_cb"][0],b,us,w)
            dm=[box(down[0],1)[0],box(down[1],1)[1],n[2]]
            U=box(tree["ub_over_cb"],1)[1]
            K=ratio(box(tree["K_product"],1),box(tree["K_form"],1))[1]
            witness=model(m,dm,U,sqrt(dm[0]/dm[1]),K,w)
            G=witness["gamma_DK_mod180_for_two_CP_branches"][0]
            assert cfg["gamma_expanded_1_width"][0]<G<cfg["gamma_expanded_1_width"][1]
            result["down_only_one_width_witness"]={"down_Yukawas":dm,"ub_cb":U,"us":K,"cb":w,
                                                   "readings":witness}
        output[scale]=result
    return output


def agrees(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(agrees(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return isclose(a,b,rel_tol=1e-10,abs_tol=1e-13)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    modes=ap.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write-new",action="store_true")
    args=ap.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert agrees(result,read("results.json"))
        print("1412 physical textures, tree ratios, complex quartets and matching costs verified")
