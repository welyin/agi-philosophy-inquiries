"""Directed interval enclosure for a three-meson rephasing invariant."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE/relative).read_text(encoding="utf-8-sig"))


class Interval:
    def __init__(self, low, high=None):
        self.low=low
        self.high=low if high is None else high

    @staticmethod
    def bound(low,high):
        return Interval(math.nextafter(low,-math.inf),math.nextafter(high,math.inf))

    def __add__(self,other):
        return self.bound(self.low+other.low,self.high+other.high)

    def __neg__(self):
        return Interval(-self.high,-self.low)

    def __sub__(self,other):
        return self+(-other)

    def __mul__(self,other):
        v=[a*b for a,b in product((self.low,self.high),(other.low,other.high))]
        return self.bound(min(v),max(v))

    def __truediv__(self,other):
        assert other.low>0
        return self*self.bound(1/other.high,1/other.low)

    def abs_range(self):
        low=0 if self.low<=0<=self.high else min(abs(self.low),abs(self.high))
        return low,max(abs(self.low),abs(self.high))


def rectangle(row):
    error=2*sum(row[1:])
    return Interval.bound(row[0]-error,row[0]+error)


def cos_from_sin(s):
    return Interval.bound(math.sqrt(1-s.high**2),math.sqrt(1-s.low**2))


def split(s, count):
    points=[s.low+(s.high-s.low)*k/count for k in range(count)]+[s.high]
    return [Interval.bound(x,y) for x,y in zip(points,points[1:])]


ZERO=Interval(0)
def mul(x,y):
    a,b=x
    c,d=y
    return (a*c-b*d,a*d+b*c)


def sc(x,a):
    return (a*x[0],a*x[1])


def add(x,y):
    return (x[0]+y[0],x[1]+y[1])


def co(x):
    return (x[0],-x[1])


def mixing(us,cb,ubcb,gamma,branch):
    # Reconstruct the CKM phase from gamma with a rigorous wider angular interval.
    s13=cb*ubcb
    c13=cos_from_sin(s13)
    s12,s23=us/c13,cb/c13
    c12,c23=cos_from_sin(s12),cos_from_sin(s23)
    k=c12*s23*s13/(s12*c23)
    d0=math.radians(gamma.low)-math.asin(k.high)-1e-13
    d1=math.radians(gamma.high)+math.asin(k.high)+1e-13
    assert 0<d0<d1<math.pi/2
    z=sc((Interval.bound(math.cos(d1),math.cos(d0)),
          Interval.bound(math.sin(d0),math.sin(d1))),Interval(branch))
    c=[add((-s12*c23,ZERO),sc(z,-c12*s23*s13)),
       add((c12*c23,ZERO),sc(z,-s12*s23*s13)),(s23*c13,ZERO)]
    t=[add((s12*s23,ZERO),sc(z,-c12*c23*s13)),
       add((-c12*s23,ZERO),sc(z,-s12*c23*s13)),(c23*c13,ZERO)]
    return c,t


def run():
    cfg=read("inputs.json")
    semi=read(cfg["semileptonic"])
    us=rectangle(semi["K_product"])/rectangle(semi["K_form"])
    cb=rectangle(read(cfg["Belle_BGL"])["BelleII_Vcb_BGL"])
    ratio=rectangle(semi["ub_over_cb"])
    g=read(cfg["gamma"])["chargedB"]
    gam=Interval(g["central_degrees"]-2*g["lower_error"]-g["additional_degrees"],
                 g["central_degrees"]+2*g["upper_error"]+g["additional_degrees"])
    u,c,t=[rectangle(row) for row in read(cfg["ordinary_yukawa_MZ"])["yukawa"]["up"]]
    lowest=[1.,1.,1.]
    calls=0
    for params in product(split(us,2),split(cb,4),split(ratio,4),split(gam,8),[-1,1]):
        charm,top=mixing(*params)
        for su,sc_sign in product([-1,1],repeat=2):
            # A global sign on all three eigenvalues changes P's sign, not its square.
            a=c*Interval(sc_sign)-u*Interval(su)
            b=t-u*Interval(su)
            edges=[add(sc(mul(co(charm[i]),charm[j]),a),sc(mul(co(top[i]),top[j]),b))
                   for i,j in ((0,1),(1,2),(2,0))]
            p=mul(mul(edges[0],edges[1]),edges[2])
            (rl,rh),(il,ih)=[z.abs_range() for z in p]
            assert il>0
            # Direct lower bound for |sin(2 arg P)|. If Re crosses zero,
            # the box alone cannot discriminate the modulo-pi phase.
            bound=2*rl*il/(rh*rh+ih*ih)
            lowest=[min(lowest[k],x) for k,x in enumerate((rl,il,bound))]
            calls+=1
    assert lowest[0]>0 and lowest[2]>.005
    return {"interval_boxes_including_discrete_branches":calls,
            "abs_real_product_floor":lowest[0],"abs_imaginary_product_floor":lowest[1],
            "abs_sine_of_three_scalar_mixing_phases_floor":lowest[2],
            "distance_from_real_axis_mod_pi_floor_degrees":math.degrees(math.asin(lowest[2]))}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    out=HERE/"independent_results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
    else:
        assert json.loads(out.read_text())==result
    print(json.dumps(result))
