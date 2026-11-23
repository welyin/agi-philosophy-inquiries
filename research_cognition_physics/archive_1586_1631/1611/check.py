"""Gauge meeting versus a co-located Higgs double zero, with a rational ODE subsolution."""
import argparse
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCALE = 10**16


class Interval:
    """Closed interval with integer endpoints in units of 10^-16."""
    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi

    def __add__(self, rhs):
        rhs = box(rhs)
        return Interval(self.lo + rhs.lo, self.hi + rhs.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, rhs):
        return self + (-box(rhs))

    def __rsub__(self, lhs):
        return box(lhs) + (-self)

    def __mul__(self, rhs):
        rhs = box(rhs)
        terms = [x*y for x in (self.lo,self.hi) for y in (rhs.lo,rhs.hi)]
        return Interval(min(terms)//SCALE,-((-max(terms))//SCALE))

    __rmul__ = __mul__

    def __truediv__(self, rhs):
        rhs = box(rhs)
        assert rhs.lo > 0
        return self * Interval(SCALE**2//rhs.hi,-((-SCALE**2)//rhs.lo))

    def __rtruediv__(self, lhs):
        return box(lhs)/self

    def sq(self):
        high = max(self.lo**2,self.hi**2)
        low = 0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2)
        return Interval(low//SCALE,-((-high)//SCALE))

    def sqrt(self):
        assert self.lo>=0
        return Interval(math.isqrt(self.lo*SCALE),math.isqrt(self.hi*SCALE)+1)

    def positive(self):
        return Interval(max(0,self.lo),max(0,self.hi))

    def values(self):
        return [self.lo/SCALE,self.hi/SCALE]


def box(lo, hi=None):
    if isinstance(lo,Interval):
        return lo
    l = Fraction(str(lo))*SCALE
    h = l if hi is None else Fraction(str(hi))*SCALE
    return Interval(math.floor(l),math.ceil(h))


def point(i):
    return Interval(i,i)


def hull(boxes):
    return Interval(min(b.lo for b in boxes),max(b.hi for b in boxes))


def log_ratio(m, mu):
    # Decimal.ln is correctly rounded in ROUND_HALF_EVEN. Ratios here terminate.
    with localcontext() as ctx:
        ctx.prec = 50
        t = (Decimal(str(m))/Decimal(str(mu))).ln()
        val = box(t)
    return Interval(val.lo-1,val.hi+1)


PI = box("3.1415926535897932","3.1415926535897933")
LOOP = 16*PI.sq()
BY = box(41)/6
BW = -box(19)/6
BD_Y = box(4)/9
BD_C = box(2)/3
BT_W = box(4)/3


def setup():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    old = json.loads((HERE/cfg["reference"]).read_text(encoding="utf-8"))
    mu = cfg["matching_scale_GeV"]
    ref, cal, widths = old["reference"], old["matching"], old["marginal_width"]
    sstep = cfg["widths"]*box(cal["g3_per_alpha"])*box(widths["alpha_s"])
    gs = box(ref["g3"])*(1+hull([-sstep,sstep]))
    gx,gw = box(ref["gY"]),box(ref["g2"])
    lstep = cfg["widths"]*box(cal["lambda_per_relative_GeV"])*box(widths["mh"])
    measured = box(ref["lambda"])*(1+hull([-lstep,lstep]))
    mr,m2 = cfg["two_B_masses_GeV"]
    d1 = hull([log_ratio(m,mu) for m in mr])
    d2 = log_ratio(m2,mu)
    U,S,A = [],[],[]
    def matching(strong, logs):
        i1,i2,i3 = box(3)/(5*gx.sq()),1/gw.sq(),1/strong.sq()
        t = (LOOP/2*(i1-i3)-box(2)/5*logs)/(box(103)/10)
        w = ((box(39)/5)*t-box(4)/15*logs-LOOP/2*(i1-i2))/BT_W
        return t,w,t-w/2
    for strong in [point(gs.lo),point(gs.hi)]:
        for gate in [point(d1.lo),point(d1.hi)]:
            t,w,a = matching(strong,gate+d2)
            U.append(t); S.append(w); A.append(a)
    whole = (hull(U),hull(S),hull(A))
    central = matching(box(ref["g3"]),log_ratio(mr[0],mu)+d2)
    return cfg, (gx,gw,gs), (d1,d2), whole, central, measured


def gauge_squared(t, g, d, a):
    dlog = (t-d[0]).positive()+(t-d[1]).positive()
    # For a fixed mean log a, hinge sums range from coincident to maximally split.
    minimum = 2*(t-a).positive()
    maximum = t+(t-2*a).positive()
    wlog = Interval(minimum.lo,maximum.hi)
    return (1/(1/g[0].sq()-2*(BY*t+BD_Y*dlog)/LOOP),
            1/(1/g[1].sq()-2*(BW*t+BT_W*wlog)/LOOP),
            1/(1/g[2].sq()-2*(-7*t+BD_C*dlog)/LOOP))


def quartic_source(x, y):
    return (2*y.sq()+(x+y).sq())*3/8


def subsolution(count, g, d, t_end, a):
    x,y,_ = gauge_squared(t_end,g,d,a)
    q = point((quartic_source(x,y)/6).sqrt().lo)
    lam = box(0)
    step = box(1)/count
    radius = (point(t_end.hi)*step).hi
    pad = Interval(-radius,radius)
    for i in range(count):
        t = t_end*(box(count-i-1,count-i)/count)
        x,y,z = gauge_squared(t,g,d,a)
        qa, la = q+pad, lam+pad
        assert qa.lo>max(la.hi,0)
        dq = t_end*qa*(17*x/6+9*y/2+16*z-9*qa)/LOOP
        dl = t_end*(-24*la.sq()+(9*y+3*x-12*qa)*la+
                    6*qa.sq()-quartic_source(x,y))/LOOP
        assert abs(dq.lo)<=t_end.hi and abs(dl.lo)<=t_end.hi
        q = point((q+step*point(dq.lo)).lo)
        lam = point((lam+step*point(dl.lo)).lo)
    return q.lo/SCALE, lam.lo/SCALE


def solve():
    cfg,g,d,(t,w,a),central,measured = setup()
    qmin,lmin = subsolution(cfg["steps"],g,d,t,a)
    # w/t lower bound deliberately decorrelates numerator and denominator.
    assert (w/t).lo>SCALE and (w/t).hi<2*SCALE and a.lo>0 and (2*a).hi<t.lo
    assert lmin > 0.129 and measured.hi/SCALE<0.12425
    mu = cfg["matching_scale_GeV"]
    return {"gauge_log_meeting":t.values(),"summed_triplet_log":w.values(),
            "active_triplets_equivalent":(w/t).values(),
            "triplet_mean_log_mass":a.values(),
            "central_meeting_GeV":mu*math.exp(central[0].lo/SCALE),
            "central_triplet_geometric_mean_GeV":mu*math.exp(central[2].lo/SCALE),
            "whole_triplet_geometric_mean_GeV":[mu*math.exp(x) for x in a.values()],
            "low_quartic_calibration":measured.values(),
            "uniform_low_quartic_lower_bound":lmin,"uniform_low_top_squared_lower_bound":qmin,
            "certified_gap_lower_rounded_down":math.floor((lmin-measured.hi/SCALE)*1e6)/1e6,
            "uniform_slices":cfg["steps"]}


if __name__=="__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args = p.parse_args()
    result = solve()
    target = HERE/"results.json"
    if args.check:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    else:
        with target.open("x",encoding="utf-8") as output:
            output.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("1611: gauge matching and continuum-of-splits Higgs comparison verified")
