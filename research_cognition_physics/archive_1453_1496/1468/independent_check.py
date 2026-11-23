"""Phase read by a real tagged decay mode, including the same Z -> c cbar source."""
import argparse
import json
from math import exp, pi, sin, sqrt
import cmath
from pathlib import Path
import runpy

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def required_ratio(delta,theta):
    if abs(delta)<1.e-15:
        return 0.
    denom=sin(delta/2-theta)
    return min(1.,abs(sin(delta/2)/denom)) if abs(denom)>1.e-15 else 1.


def run():
    g=read("../1369/inputs.json")["g3"][0]
    eta=(g*g/(4*pi))**(6/23)
    m=read("../1404/inputs.json")
    c=m["const_GeV"]
    short=m["short_distance"]
    box=c["GF"]**2*c["MW"]**2*short["S0_published"][0]*short["etaB"][0]/(4*pi*pi)
    tree=c["GF"]*eta/sqrt(2)
    p=read("../1465/inputs.json")
    data=read("../1465/"+p["three_rates"])
    lib=runpy.run_path(str(HERE/"../1465/check.py"))
    build=runpy.run_path(str(HERE/"../1467/check.py"))["build"]
    results={}
    for sign in ("free_0.36","free_-0.36"):
        F=complex(*read("../1467/results.json")[sign]["Fsb"])
        V,*_=build(F,p,data,lib["hermitian_power"],lib["complete_basis"])
        lc=V[1,1].conjugate()*V[1,2]
        lt=V[2,1].conjugate()*V[2,2]
        phi=cmath.phase((box*lt**2+tree*F**2)/lc**2)
        theta=cmath.phase(F/lc)
        vals={}
        for name,mean,sigma in [("CMS13",-.073,.022+.010),("LHCb2015_18_JpsiKK",-.039,.022+.006),
                                ("LHCb_pipi",-.077,.034+.007)]:
            def min_over_interval(low,high):
                if low<=phi<=high:
                    return 0.
                # In these small real intervals there is no pole or stationary
                # point away from zero; verify explicitly.
                a,b=(phi-high)/2-theta,(phi-low)/2-theta
                assert sin(a)*sin(b)>0 and abs(b-a)<pi
                return min(required_ratio(phi-low,theta),required_ratio(phi-high,theta))
            delta=phi-mean
            z=cmath.exp(1j*delta)
            h=(z-1)/(cmath.exp(1j*theta)-z*cmath.exp(-1j*theta))
            assert abs(h.imag)<1e-10
            zcheck=(1+h*cmath.exp(1j*theta))/(1+h*cmath.exp(-1j*theta))
            assert abs(z-zcheck)<1e-12 and abs(abs(h)-required_ratio(delta,theta))<1e-10
            vals[name]={"relative_amplitude_for_central_phase":abs(h),
                        "one_width_floor_for_Z_only_relative_amplitude":min_over_interval(mean-sigma,mean+sigma)}
        results[sign]={"phase_before_decay_pollution":phi,"phase_Z_minus_tree":theta,
                       "weak_ratio_abs_F_over_charm":abs(F/lc),"decay_modes":vals}
    return json.loads(json.dumps(results),parse_float=lambda x:float(f"{float(x):.11g}"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
    else:
        assert result==read("independent_results.json")
        print("1468: mixing phases and the independently matched charmonium decay phase checked")
