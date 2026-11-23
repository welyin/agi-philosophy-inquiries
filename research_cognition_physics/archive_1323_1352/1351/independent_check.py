"""Inverse bound on the CKM phase from the held-out frequencies, without a gamma fit."""
import argparse
import cmath
import json
import math
from fractions import Fraction as F
from pathlib import Path


def f(value):
    return F(str(value))


def certified_bounds():
    # All enclosures in this part are rational; coarse c12 and h bounds below
    # are certified by squaring positive expressions.
    qlo,qhi=f(.2151)/f(.9736), f(.2163)/f(.9672)
    rho0,rho1=f(.075),f(.091)
    wmax=f(.0495)
    c13sq=1-(rho1*wmax)**2
    alo,ahi=f(.2209),f(.224)
    blo,bhi=f(.974),f(.976)
    hlo,hhi=f(.0748),rho1
    assert alo<=qlo and ahi*ahi*c13sq>=qhi*qhi
    assert blo**2<=1-ahi**2 and bhi**2>=1-alo**2
    assert hlo**2<=rho0**2*(c13sq-wmax**2)
    assert (f(5279.58)-f(.43))/(f(5366.90)+f(.51)) > f(.98)
    assert (f(5279.58)+f(.43))/(f(5366.90)-f(.51)) < f(.985)
    mlo,mhi=f(.98),f(.985)
    xilo,xihi=f(1.182),f(1.230)
    reslo,reshi=f(.99)/f(1.01),f(1.01)/f(.99)
    rlo=f(.5019)/f(17.7766)
    rhi=f(.5081)/f(17.7600)
    # X = |Vtd/Vts|^2, rather than the hadronic frequency ratio.
    xlo,xhi=rlo*xilo**2/(mhi*reshi), rhi*xihi**2/(mlo*reslo)
    Cmin,Cmax=blo**2+alo**2*hlo**2,bhi**2+ahi**2*hhi**2
    Nmin=alo**2+blo**2*hlo**2-xhi*Cmax
    Nmax=ahi**2+bhi**2*hhi**2-xlo*Cmin
    assert Nmin>0
    cmin=Nmin/(2*ahi*bhi*hhi*(1+xhi))
    cmax=Nmax/(2*alo*blo*hlo*(1+xlo))
    mirror_floor=mlo*reslo/xihi**2*qlo**2/(1-qlo**2)
    extra_max=rhi/mirror_floor
    return {
        "data_frequency_ratio": [float(rlo),float(rhi)],
        "inferred_ckm_ratio_squared":[float(xlo),float(xhi)],
        "cos_delta_outer":[float(cmin),float(cmax)],
        "principal_delta_deg_outer":[math.degrees(math.acos(float(cmax))),
                                     math.degrees(math.acos(float(cmin)))],
        "mirror_frequency_floor":float(mirror_floor),
        "max_further_mirror_d_over_s_multiplier":float(extra_max),
    }


def construct():
    q,w,rho=.2157/.9704,.04,.083
    z=w*rho
    c13=math.sqrt(1-z*z)
    a,b=q/c13,math.sqrt(1-(q/c13)**2)
    s23=w/c13
    c23=math.sqrt(1-s23*s23)
    h=rho*c13*c23
    k=b*s23*z/(a*c23)
    target=.5050/17.7683
    mass=5279.58/5366.90
    xi=1.206
    x=target*xi*xi/mass
    cos_delta=(a*a+b*b*h*h-x*(b*b+a*a*h*h))/(2*a*b*h*(1+x))
    delta=math.acos(cos_delta)
    gamma=cmath.phase(cmath.exp(1j*delta)+k)
    # A forward evaluation uses the *complex* entries of the CKM matrix.
    def direct(phase):
        unit=cmath.exp(1j*phase)
        cd=-a*c23-b*s23*z*unit
        td=a*s23-b*c23*z*unit
        ts=-b*s23-a*c23*z*unit
        measured_angle=cmath.phase(-(b*c13)*z*unit/(cd*w))
        return mass*abs(td/ts)**2/xi**2,measured_angle
    result,angle=direct(delta)
    assert abs(result-target)<1e-14 and abs(angle-gamma)<1e-13
    assert 62.6<math.degrees(gamma)<74.9
    mirror=gamma+math.pi
    other=mirror+math.asin(k*math.sin(mirror))
    other_ratio,other_angle=direct(other)
    assert abs(cmath.exp(1j*other_angle)+cmath.exp(1j*gamma))<1e-13
    return {"w":w,"gamma_deg":math.degrees(gamma),"delta_deg":math.degrees(delta),
            "zero_extra_correction_ratio":result,"mirrored_ratio_same_other_inputs":other_ratio,
            "absolute_frequency_not_fixed":True}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=dict(round=1351,status="PASS",bounds=certified_bounds(),existence_witness=construct())
    path=Path(__file__).resolve().with_name("independent_results.json")
    if args.write_new:
        with path.open("x",encoding="utf8",newline="\n") as out:
            out.write(json.dumps(data,indent=2)+"\n")
    else:
        assert json.loads(path.read_text(encoding="utf8"))==data
    print(json.dumps(data))

