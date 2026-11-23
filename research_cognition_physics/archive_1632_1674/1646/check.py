"""Three-port mass-squared transport, with complete angular and spectral boxes."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def lower(x):
    return np.nextafter(x, -np.inf)


def upper(x):
    return np.nextafter(x, np.inf)


class Box:
    """Outward binary64 intervals with optional broadcast subinterval axes."""
    def __init__(self, left, right=None):
        self.lo = np.asarray(left)
        self.hi = self.lo if right is None else np.asarray(right)

    def __add__(self, other):
        other = wrap(other)
        return Box(lower(self.lo+other.lo), upper(self.hi+other.hi))

    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-wrap(other))

    def __rsub__(self, other):
        return wrap(other) + (-self)

    def __mul__(self, other):
        other = wrap(other)
        terms = [a*b for a,b in itertools.product((self.lo,self.hi),(other.lo,other.hi))]
        return Box(lower(np.minimum.reduce(terms)), upper(np.maximum.reduce(terms)))

    __rmul__ = __mul__

    def sq(self):
        low = np.where((self.lo<=0) & (self.hi>=0), 0., np.minimum(self.lo**2,self.hi**2))
        return Box(np.maximum(0.,lower(low)), upper(np.maximum(self.lo**2,self.hi**2)))

    def sqrt(self):
        return Box(np.maximum(0.,lower(np.sqrt(np.maximum(0.,self.lo)))),
                   upper(np.sqrt(np.maximum(0.,self.hi))))

    def magnitude(self):
        return self.sq().sqrt()


def wrap(x):
    return x if isinstance(x,Box) else Box(x)


def measured(x, n):
    return x[0]-n*x[1], x[0]+n*x[-1]


def parameter_boxes(cfg, source, case, width):
    if case=="quark":
        down=[measured(row,width) for row in source["yukawa"]["down"]]
        d,s,b=down
        r=[(s[0]**2-d[1]**2)/(b[1]**2-d[1]**2),
           (s[1]**2-d[0]**2)/(b[0]**2-d[0]**2)]
        return [r]+[measured(row,width) for row in source["quark_sines_and_phase"][:3]], 1
    data=source[case]
    solar=measured(data["solar_gap"],width)
    atm=measured(data["atmospheric_gap"],width)
    ratios=[solar[0]/atm[1], solar[1]/atm[0]]
    normal=case=="NO"
    r=[v/(1+v) for v in ratios] if normal else ratios
    p=[math.sqrt(v) for v in measured(data["solar_sine_squared"],width)]
    t=[math.sqrt(v) for v in measured(data["atmospheric_sine_squared"],width)]
    u=[math.sqrt((1-math.sqrt(1-v))/2) for v in measured(source["DayaBay_amplitude"],width)]
    return [r,p,u,t], 1 if normal else -1


def split(pair, count, axis):
    # Input functions use elementary floating operations. Inflate by 2e-13,
    # much larger than their endpoint error on the registered parameter boxes.
    left,right=pair
    points=np.linspace(left-2e-13,right+2e-13,count+1)
    dims=[1]*5
    dims[axis]=count
    return Box(lower(points[:-1].reshape(dims)),upper(points[1:].reshape(dims)))


def enclose_slice(cfg, physics, case, width, mirrored, r_index, radial_count):
    params,sign=parameter_boxes(cfg,physics,case,width)
    r_points=np.linspace(*params[0],radial_count+1)
    ranges=[[r_points[r_index],r_points[r_index+1]],*params[1:]]
    r,p,u,t=[split(pair,1 if i==0 else radial_count,i) for i,pair in enumerate(ranges)]
    if case=="quark":
        if width==0:
            ang=read(cfg["charged"])["gamma_Bpm"]["deg"]
            dlo,dhi=ang-.01,ang+.01
        else:
            dlo,dhi=cfg["quark_phase_outer_deg_by_width"][str(width)]
        zrange=sorted(math.cos(math.radians(a+180*mirrored)) for a in (dlo,dhi))
    else:
        zrange=[-1.,1.]
    z=split(zrange,cfg["splits"]["delta"],4)
    c=(1-p.sq()).sqrt(); C=(1-u.sq()).sqrt(); T=(1-t.sq()).sqrt()
    a=r*p*c; b=u*(sign-r*p.sq()); d=1-2*t.sq()
    X=t*T*(b.sq()-a.sq())+a*b*d*z
    Y=t*T*(sign*C.sq()-r*c.sq()+r*p.sq()*u.sq())-a*u*d*z
    re=X*Y-a.sq()*b*u*(1-z.sq())
    im=a*t*T*u*sign*(sign-r)*(1-z.sq()).sqrt()
    re_abs,im_abs=re.magnitude(),im.magnitude()
    alo=math.degrees(float(np.min(np.arctan2(im_abs.lo,re_abs.hi))))-1e-9
    ahi=math.degrees(float(np.max(np.arctan2(im_abs.hi,re_abs.lo))))+1e-9
    # Bounds independent of the phase, for the uniform-gain subcase.
    abs_b=b.magnitude()
    third=t*T*(sign*C.sq()-r*c.sq()+r*p.sq()*u.sq())
    # H23 correction has modulus <= a*u via T^2+t^2=1.
    min23=float(np.min((third.magnitude()-a*u).lo))
    max12=float(np.max((C*(a*T+abs_b*t)).hi))
    min12=float(np.min((C*(abs_b*t-a*T)).lo))
    min13=float(np.min((C*(abs_b*T-a*t)).lo))
    all_min=min(min12,min13,min23)
    spread_lo=min23/max12
    spread_hi=(1+params[0][1]+2e-12)/(2*all_min) if all_min>0 else None
    return {"distance_to_real_axis_deg":[max(0.,alo),min(90.,ahi)],
            "edge_modulus_max_over_min_lower":spread_lo-1e-9,
            "edge_modulus_max_over_min_upper":spread_hi+1e-9 if spread_hi else None,
            "scaled_real_part_range":[float(np.min(re.lo)),float(np.max(re.hi))]}


def enclose(cfg, physics, case, width, mirrored=False):
    count=cfg["splits"]["radial"] if width else 1
    slices=[enclose_slice(cfg,physics,case,width,mirrored,i,count) for i in range(count)]
    def outer(key):
        return [min(s[key][0] for s in slices), max(s[key][1] for s in slices)]
    return {"distance_to_real_axis_deg":outer("distance_to_real_axis_deg"),
            "scaled_real_part_range":outer("scaled_real_part_range"),
            "edge_modulus_max_over_min_lower":min(s["edge_modulus_max_over_min_lower"] for s in slices),
            "edge_modulus_max_over_min_upper":max(s["edge_modulus_max_over_min_upper"] for s in slices)}


def mixing(p,u,t,phase):
    c,C,T=(math.sqrt(1-s*s) for s in (p,u,t))
    e=np.exp(1j*phase)
    return np.array([[c*C,p*C,u/e],[-p*T-c*t*u*e,c*T-p*t*u*e,t*C],
                     [p*t-c*T*u*e,-c*t-p*T*u*e,T*C]])


def loop(h):
    return h[0,1]*h[1,2]*h[2,0]


def verify_closed_expression(cfg, physics):
    for case in ("quark","NO","IO"):
        par,sg=parameter_boxes(cfg,physics,case,0)
        r,p,u,t=[pair[0] for pair in par]
        c,C,T=[math.sqrt(1-v*v) for v in (p,u,t)]
        a,b,d=r*p*c,u*(sg-r*p*p),1-2*t*t
        for phase in np.linspace(0.,2*math.pi,29):
            mat=mixing(p,u,t,phase)
            h=(mat*np.array([0.,r,sg]))@mat.conj().T
            z,y=math.cos(phase),math.sin(phase)
            x=t*T*(b*b-a*a)+a*b*d*z
            v=t*T*(sg*C*C-r*c*c+r*p*p*u*u)-a*u*d*z
            formula=C*C*complex(x*v-a*a*b*u*y*y,a*t*T*u*sg*(sg-r)*y)
            assert abs(formula-loop(h))<=2e-14*abs(loop(h))


def inverse_example():
    angles=(.22,.035,.19,.97)
    v=mixing(*angles)
    h=(v*np.array([.006,.21,1.]))@v.conj().T
    rng=np.array([.4,2.,1.3])*np.exp(1j*np.array([.16,-.42,.75]))
    target=(rng[:,None]*h)*rng.conj()[None,:]+np.diag([.15,.21,.32])
    # Only observed Hermitian kernels, never the original complex gains.
    ratios=np.array([target[i,j]/h[i,j] for i,j in [(0,1),(1,2),(2,0)]])
    a,b,c=abs(ratios)
    gains=np.sqrt([a*c/b,a*b/c,b*c/a])
    phases=[0.,-np.angle(ratios[0]),np.angle(ratios[2])]
    recovered=gains*np.exp(1j*np.array(phases))
    lifted=recovered[:,None]*h*recovered.conj()[None,:]
    lifted+=np.diag(np.diag(target-lifted))
    residual=float(np.max(abs(lifted-target)))
    # One directed edge can make any nonzero loop phase while retaining positivity.
    rotated=target.copy()
    rotated[0,1]*=np.exp(.09j)
    rotated[1,0]=rotated[0,1].conj()
    mismatch=float(np.angle(loop(rotated)/loop(h)))
    assert residual<2e-14 and np.min(np.linalg.eigvalsh(rotated))>0 and abs(mismatch-.09)<2e-14
    # All real diagonal shifts and absolute masses cancel from off-diagonal products.
    return {"recovered_positive_gains":gains.tolist(),"reconstruction_error":residual,
            "single_edge_rotation_result_rad":mismatch,"smallest_modified_eigenvalue":float(np.min(np.linalg.eigvalsh(rotated)))}


def check_phase(cfg, physics):
    g=read(cfg["charged"])["gamma_Bpm"]
    p,u,t=[measured(v,2)[1] for v in physics["quark_sines_and_phase"][:3]]
    shift=math.degrees(math.asin(p*u*t/math.sqrt((1-p*p)*(1-t*t))))
    for n in (1,2):
        old=[g["deg"]-n*g["down"]-g["transfer_cap_deg"]-shift,
             g["deg"]+n*g["up"]+g["transfer_cap_deg"]+shift]
        lo,hi=cfg["quark_phase_outer_deg_by_width"][str(n)]
        assert lo<old[0] and old[1]<hi
    assert shift<.01
    return shift


def run():
    cfg=read("inputs.json"); phys=read(cfg["reference"])
    verify_closed_expression(cfg,phys)
    out={"gamma_s_minus_delta_upper_deg":check_phase(cfg,phys), "constructive_test":inverse_example(),"widths":{}}
    for width in (0,1,2):
        cases={key:enclose(cfg,phys,case,width,mirror) for key,case,mirror in
               (("q_positive","quark",False),("q_shifted","quark",True),("NO","NO",False),("IO","IO",False))}
        for name in ("NO","IO"):
            l_hi=cases[name]["distance_to_real_axis_deg"][1]
            gaps={branch:cases[branch]["distance_to_real_axis_deg"][0]-l_hi for branch in ("q_positive","q_shifted")}
            assert min(gaps.values())>0
            gap=min(gaps.values())
            # Maximal allowed loop-phase improvement at fixed per-edge fractional amplitude.
            cases[name]["relative_to_quark"]={"separation_deg":gaps,
                "one_edge_amplitude_fraction_floor":math.sin(math.radians(gap)),
                "three_edges_max_fraction_floor":math.sin(math.radians(gap/3))}
            assert cases["q_positive"]["edge_modulus_max_over_min_lower"]>cases[name]["edge_modulus_max_over_min_upper"]
        out["widths"][str(width)]=cases
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda s:float(f"{float(s):.10g}"))
    if args.check:
        assert result==read("results.json")
    elif args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as target:
            json.dump(result,target,indent=2)
            target.write("\n")
    else:
        print(json.dumps(result,indent=2))
    print("1646 PASS: continuous quark/lepton phase bounds and constructive port equivalence")
