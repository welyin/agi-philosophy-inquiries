"""Exact outward interval bounds for two Clebsch placements and the joint electron row."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SCALE = 10**12


def ceil_div(a,b):
    return -(-a//b)


class I:
    """Closed interval, with endpoints integer multiples of 10**-12."""
    def __init__(self, low, high):
        assert low <= high
        self.lo, self.hi = low, high

    @classmethod
    def point(cls,x):
        if isinstance(x,I):
            return x
        value = Fraction(str(x))*SCALE
        return cls(math.floor(value), math.ceil(value))

    def __add__(self, other):
        b = I.point(other)
        return I(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi,-self.lo)

    def __sub__(self,other):
        return self+-I.point(other)

    def __rsub__(self,other):
        return I.point(other)+-self

    def __mul__(self,other):
        b = I.point(other)
        edge = [x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(edge)//SCALE,ceil_div(max(edge),SCALE))

    __rmul__ = __mul__

    def __truediv__(self,other):
        b = I.point(other)
        assert b.lo > 0
        return self*I(SCALE**2//b.hi,ceil_div(SCALE**2,b.lo))

    def __rtruediv__(self,other):
        return I.point(other)/self

    def sqrt(self):
        assert self.lo >= 0
        return I(math.isqrt(self.lo*SCALE),math.isqrt(self.hi*SCALE)+1)

    def square(self):
        a,b = self.lo,self.hi
        low = 0 if a<=0<=b else min(a*a,b*b)
        return I(low//SCALE,ceil_div(max(a*a,b*b),SCALE))

    def parts(self,n):
        split = [self.lo+(self.hi-self.lo)*i//n for i in range(n+1)]
        return [I(a,b) for a,b in zip(split,split[1:])]

    def values(self):
        return [self.lo/SCALE,self.hi/SCALE]


def hull(boxes):
    boxes = list(boxes)
    return I(min(a.lo for a in boxes),max(a.hi for a in boxes))


def width(data):
    middle, low = map(I.point,data[:2])
    high = I.point(data[-1])
    return hull((middle-2*low,middle+2*high))


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def charged(q,s,clebsch):
    ca,cb,cc = clebsch
    b = (q.square()*(1-s.square())+s.square()).sqrt()
    a,c = q/b,(1-q.square())*s*(1-s.square()).sqrt()/b
    x,y,z = ca*a,cb*b,cc*c
    diag = y.square()+z.square()-x.square()
    tan2 = 2*x*z/diag
    se = ((1-1/(1+tan2.square()).sqrt())/2).sqrt()
    trace = x.square()+y.square()+z.square()
    det = x*y
    mass_ratio = (trace+(trace.square()-4*det.square()).sqrt())/(2*det)
    return se,mass_ratio


def projection(R,branch):
    B = 6*(R+1/R)+16
    h = B+(B.square()-1936).sqrt()
    r = 8/h if branch=="small" else h/242
    lam = 6*r/R
    gap = 6*r*(R-1/R)
    he = (11*r.square()-lam/3)/gap
    hm = (2+99*r.square()-12*r-5*lam/6)/gap
    re = (33*r.square()-2*r-lam/3)/gap
    im = 2*I.point(3).sqrt()*r/gap
    return he,hm,re,im


def p3_bounds(se, R, p1, branch):
    bands = [[],[]]
    for angle in se.parts(12):
        ce = (1-angle.square()).sqrt()
        sc = angle*ce
        cosphi = (4-3*angle.square()-6*p1)/(4*sc)
        if cosphi.lo > SCALE or cosphi.hi < -SCALE:
            continue
        cosphi = I(max(-SCALE,cosphi.lo),min(SCALE,cosphi.hi))
        sinphi = (1-cosphi.square()).sqrt()
        for ratio in R.parts(72):
            he,hm,re,im = projection(ratio,branch)
            mid = ce.square()*he+angle.square()*hm+2*sc*re*cosphi
            offset = 2*sc*im*sinphi
            for values,box in zip(bands,(mid-offset,mid+offset)):
                values.append(box)
    return [hull(row) for row in bands]


def witness(q,s,R,clebsch):
    b = math.sqrt(q*q*(1-s*s)+s*s)
    a,c = q/b,(1-q*q)*s*math.sqrt(1-s*s)/b
    ca,cb,cc = clebsch
    E = np.array([[0,ca*a],[cb*b,cc*c]])
    eig, left = np.linalg.eigh(E@E.T)
    se = abs(left[0,1])
    ce = math.sqrt(1-se*se)
    B = 6*(R+1/R)+16
    r = 8/(B+math.sqrt(B**2-1936))
    a0,b0 = np.array([0,1,1]),np.array([1,3,1])
    M = np.outer(a0,a0)+r*np.exp(2j*np.pi/3)*np.outer(b0,b0)
    _,U = np.linalg.eigh(M@M.conj().T)
    measured = abs(ce*U[0]-se*U[1])**2
    atmospheric = 1-abs(U[2,2])**2/(1-measured[2])
    return {"lepton_mass_ratio":math.sqrt(eig[1]/eig[0]),"rotation_degrees":math.degrees(math.asin(se)),
            "reactor_amplitude":4*measured[2]*(1-measured[2]),"solar":measured[1]/(1-measured[2]),
            "atmospheric":atmospheric,"tau_first":abs(U[2,0])**2}


def solve():
    data = read(ROOT/"inputs.json")
    masses = read(ROOT/data["quark_lepton"])
    lattice = masses["lattice"]
    q = width(lattice["md"])/width(lattice["ms"])
    lepton = width(masses["muonium"]["pole_ratio"])
    charm = read(ROOT/data["charm"])
    _,v,errors = next(row for row in charm["edges"] if row[0]=="cd")
    error = sum(map(I.point,errors),I.point(0))
    s = hull(((I.point(v)-2*error)/I.point("1.02"),(I.point(v)+2*error)/I.point(".98")))
    osc = read(ROOT/data["oscillations"])
    ds,atm = osc["KamLAND_dm21_eV2"],osc["T2K_dm32_NO_eV2"]
    R = (1+width(atm)/width(ds)).sqrt()
    p3 = (1-(1-width(osc["DayaBay_sin2_2theta13"])).sqrt())/2
    p1 = (1-p3)*(1-width(osc["JUNO_sin2_theta12_NO"]))
    out = {}
    for name,clebsch in data["clebsch"].items():
        boxes = [charged(qi,si,clebsch) for qi in q.parts(16) for si in s.parts(16)]
        se = hull(row[0] for row in boxes)
        mass = hull(row[1] for row in boxes)
        branches = {}
        for branch in ("small","large"):
            bands = p3_bounds(se,R,p1,branch)
            branches[branch] = {"reactor_weights": [a.values() for a in bands],
                                "overlap_possible":any(max(a.lo,p3.lo)<=min(a.hi,p3.hi) for a in bands)}
        ex = witness(lattice["md"][0]/lattice["ms"][0],v,math.sqrt(1+atm[0]/ds[0]),clebsch)
        ex["lepton_ratio_transfer"] = masses["muonium"]["pole_ratio"][0]/ex["lepton_mass_ratio"]
        out[name]={"sine_rotation":se.values(),"lepton_ratio":mass.values(),"transfer":(lepton/mass).values(),
                   "branches":branches,"real_minus_witness":ex}
    passing = [(name, branch) for name,row in out.items() for branch,test in row["branches"].items()
               if test["overlap_possible"]]
    assert passing == [("small_left", "small")]
    ex=out["small_left"]["real_minus_witness"]
    assert .0803<ex["reactor_amplitude"]<.0899 and .2908<ex["solar"]<.3164
    assert data["mass_matching_window"][0]<ex["lepton_ratio_transfer"]<data["mass_matching_window"][1]
    return {"requested_electron_weights":{"first":p1.values(),"third":p3.values()},"textures":out}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    result=solve()
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==read(path)
        print("1486: outward weak-basis, mass matching and both phased reactor branches verified.")
