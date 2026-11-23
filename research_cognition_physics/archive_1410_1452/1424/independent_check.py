"""Independent verification. A pair of vacuum phase bins eliminating the slow term."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONV = 1.26693268


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def compute():
    cfg=read("inputs.json")
    old=read(cfg["angles"])
    z=(1-math.sqrt(1-old["DayaBay_amplitude"][0]))/2
    alpha=old["quark_sines_and_phase"][-1][0]
    rows={}
    for label,sign in (("NO",1),("IO",-1)):
        if label=="NO":
            x=old[label]["solar_sine_squared"][0]
            solar=old[label]["solar_gap"][0]
        else:
            x=cfg["solar_IO_fallback"]["s12_squared"][0]
            solar=cfg["solar_IO_fallback"]["dm21_eV2"][0]
        y=old[label]["atmospheric_sine_squared"][0]
        m31=sign*cfg["daya_bay_ee_eV2"][0]+x*solar
        m32=m31-solar
        u1=math.pi/(2*CONV*abs(m31+m32))
        ds=CONV*u1*solar
        r=(math.sin(ds)/math.sin(3*ds))**2
        B=2*math.sqrt(x*(1-x)*y*(1-y)*z)
        p3=y*(1-z)
        moment=[]
        comb=[]
        for delta in (-alpha,math.pi+alpha):
            phase=cmath.exp(1j*delta)
            p1=abs(math.sqrt(x*(1-y))+math.sqrt((1-x)*y*z)*phase)**2
            p2=abs(math.sqrt((1-x)*(1-y))-math.sqrt(x*y*z)*phase)**2
            assert abs(p1+p2+p3-1)<2e-15
            def survival(u):
                return abs(p1+p2*cmath.exp(-2j*CONV*solar*u)+p3*cmath.exp(-2j*CONV*m31*u))**2
            comb.append(survival(u1)-r*survival(3*u1))
            moment.append(abs(m31-x*solar)-abs((p1*m31+p2*m32)/(p1+p2)))
        exact=8*p3*B*abs(math.cos(alpha))*(math.sin(ds)+r*math.sin(3*ds))
        assert math.isclose(abs(comb[0]-comb[1]),exact,rel_tol=1e-12)
        # Any two positive normalized response distributions in 5% L/E windows.
        rmin=(math.sin(.95*ds)/math.sin(3*1.05*ds))**2
        rmax=(math.sin(1.05*ds)/math.sin(3*.95*ds))**2
        gap=8*p3*B*abs(math.cos(alpha))*(
            math.cos(.05*math.pi/2)*math.sin(.95*ds)
            + rmin*math.cos(.15*math.pi/2)*math.sin(3*.95*ds))
        rows[label]={"central_moment_differences_eV2":moment,"two_L_over_E_km_per_GeV":[u1,3*u1],
                     "monochromatic_gap_probability":exact,"monochromatic_weight":r,
                     "five_percent_weight_bounds":[rmin,rmax],"five_percent_gap_lower":gap,
                     "per_bin_absolute_error_ceiling_if_angles_fixed":gap/(2*(1+rmax))}
    return rows


def near(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(near(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(near(i,j) for i,j in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-10,abs_tol=1e-15)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    opts=parser.add_mutually_exclusive_group()
    opts.add_argument("--write",action="store_true")
    opts.add_argument("--check",action="store_true")
    args=parser.parse_args()
    values=compute()
    if args.write:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2)
            stream.write("\n")
    else:
        assert near(values,read("independent_results.json"))
    print(json.dumps(values))
