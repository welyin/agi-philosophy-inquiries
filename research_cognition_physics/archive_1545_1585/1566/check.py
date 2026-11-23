"""Running signed down mass at the physical lightest up-flavor port."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def current(angles):
    a,b,c,delta = angles
    sa,sb,sc = np.sin([a,b,c])
    ca,cb,cc = np.cos([a,b,c])
    z = np.exp(1j*delta)
    return np.array([[ca*cb,sa*cb,sb/z],[-sa*cc-ca*sc*sb*z,ca*cc-sa*sc*sb*z,sc*cb],
                     [sa*sc-ca*cc*sb*z,-ca*sc-sa*cc*sb*z,cc*cb]])


def at(row,n):
    return row[0]+n*row[1]


def initial(cfg,phase=None):
    high = read(HERE/cfg["high_quarks"])["high_scale"]
    positions = cfg["width_coordinate"]
    a = high["angles_theta12_theta13_theta23_delta"]
    angles = [at(a[0],positions["theta12"]),at(a[1],positions["theta13"]),
              a[2][0],a[3][0] if phase is None else phase]
    v = current(angles)
    p = abs(v[0])**2
    down = high["down"]
    s = at(down[1],positions["strange"])
    b = at(down[2],positions["bottom"])
    eps = np.array(cfg["signs"])
    d = -(p[1]*eps[1]*s+p[2]*eps[2]*b)/(p[0]*eps[0])
    assert down[0][0]-2*down[0][1] <= d <= down[0][0]+2*down[0][1]
    up = np.diag([x[0]**2 for x in high["up"]]).astype(complex)
    dl = np.array([d,s,b])
    hd = (v*dl**2)@v.conj().T
    le = read(HERE/cfg["high_leptons"])["scales"]["1e16"][3:]
    he = np.diag([x[0]**2 for x in le]).astype(complex)
    return np.array([up,hd,he]),v,eps*dl


def flow(t,state,weak,cfg):
    # t=log(mu/Mhigh). Every derivative is with respect to *increasing* log mu.
    a,b,e = state
    l = math.log(cfg["high_GeV"]/weak["mu_ref_GeV"])+t
    g0 = np.array([weak[k][0] for k in ("g1_GUT","g2","g3")])
    g1,g2,g3 = 1/(1/g0**2-2*l*np.array(cfg["gauge_b"])/LOOP)
    trace = np.trace(3*(a+b)+e).real
    au = trace-(.85*g1+2.25*g2+8*g3)
    ad = trace-(.25*g1+2.25*g2+8*g3)
    ae = trace-2.25*(g1+g2)
    cross = a@b+b@a
    return np.array([2*au*a+3*a@a-1.5*cross,
                     2*ad*b+3*b@b-1.5*cross,2*ae*e+3*e@e])/LOOP


def measured(state,signs):
    a,b,_ = state
    _,q = np.linalg.eigh(a)
    e,v = np.linalg.eigh(b)
    assert min(e)>0
    signed = signs*np.sqrt(e)
    prob = abs(q[:,0].conj()@v)**2
    f = float(prob@signed)
    return dict(F=f,F_over_s=f/math.sqrt(e[1]), F_over_d=f/math.sqrt(e[0]))


def linear(cfg,weak,phase=None):
    state,v,d = initial(cfg,phase)
    adot,bdot,_ = flow(0.,state,weak,cfg)
    a = np.diag(state[0]).real
    root = (v*d)@v.conj().T
    moving_port = sum(2*(adot[j,0]*root[0,j]).real/(a[0]-a[j]) for j in (1,2))
    ds = (v.conj().T@bdot@v)/(d[:,None]+d[None,:])
    root_change = (v@ds@v.conj().T)[0,0].real
    assert np.linalg.norm(root@(v@ds@v.conj().T)+(v@ds@v.conj().T)@root-bdot) < 1e-17
    f1,f3,f5 = [float((abs(v[0])**2)@d**power) for power in (1,3,5)]
    delta = .01
    finite = (measured(state+delta*flow(0,state,weak,cfg),np.sign(d))["F"] -
              measured(state-delta*flow(0,state,weak,cfg),np.sign(d))["F"])/(2*delta)
    total = float(moving_port+root_change)
    assert math.isclose(finite,total,rel_tol=3e-4,abs_tol=2e-15)
    return dict(signed_high=d.tolist(), F=f1,F_cubed=f3,F_fifth=f5,
                moving_up_port=float(moving_port),signed_root=float(root_change),
                derivative_total=total, derivative_over_s=total/d[1],finite_difference=finite)


def run(cfg,weak,steps):
    state,_,d = initial(cfg)
    t = 0.
    out = {}
    goals = [math.log(x/cfg["high_GeV"]) for x in cfg["check_scales_GeV"]]
    grid = sorted(set(np.linspace(0,goals[-1],steps+1).tolist()+goals),reverse=True)
    for next_t in grid[1:]:
        delta = next_t-t
        k1 = flow(t,state,weak,cfg)
        k2 = flow(t+delta/2,state+delta*k1/2,weak,cfg)
        k3 = flow(t+delta/2,state+delta*k2/2,weak,cfg)
        k4 = flow(next_t,state+delta*k3,weak,cfg)
        state += delta*(k1+2*k2+2*k3+k4)/6
        state = (state+state.conj().transpose(0,2,1))/2
        t = next_t
        if next_t in goals:
            out[str(cfg["check_scales_GeV"][goals.index(t)])] = measured(state,np.sign(d))
    return out


def signed_interval(cfg, sign):
    """Outward rounded interval for 16 pi² dF/dlog(mu) at a real CKM endpoint."""
    from fractions import Fraction
    mesh = 10**34

    class Interval:
        def __init__(self, x, y=None):
            if isinstance(x,Interval):
                self.lo, self.hi = x.lo, x.hi
            else:
                x = Fraction(str(x)) if not isinstance(x,Fraction) else x
                y = x if y is None else y
                self.lo = Fraction(math.floor(x*mesh),mesh)
                self.hi = Fraction(math.ceil(y*mesh),mesh)

        def __add__(self, other):
            k=Interval(other)
            return Interval(self.lo+k.lo, self.hi+k.hi)
        __radd__ = __add__

        def __neg__(self):
            return Interval(-self.hi,-self.lo)

        def __sub__(self,other):
            return self+(-Interval(other))

        def __rsub__(self,other):
            return -self+other

        def __mul__(self,other):
            k=Interval(other)
            p=[a*b for a in (self.lo,self.hi) for b in (k.lo,k.hi)]
            return Interval(min(p),max(p))
        __rmul__ = __mul__

        def __truediv__(self,other):
            k=Interval(other)
            assert not k.lo<=0<=k.hi
            return self*Interval(1/k.hi,1/k.lo)

    def trig(x, odd):
        # Alternating series through power 16/17, with a signed next term.
        x=Fraction(x)
        terms=[(-1)**k*x**(2*k+odd)/math.factorial(2*k+odd) for k in range(10)]
        a,b=sum(terms[:-1]),sum(terms)
        return Interval(min(a,b),max(a,b))

    def value(row,delta=0):
        return Fraction(str(row[0]))+delta*Fraction(str(row[1]))

    h=read(HERE/cfg["high_quarks"])["high_scale"]
    work=cfg["width_coordinate"]
    q=h["angles_theta12_theta13_theta23_delta"]
    ang=[value(q[0],work["theta12"]),value(q[1],work["theta13"]),value(q[2])]
    sa,sb,sc=[trig(x,1) for x in ang]
    ca,cb,cc=[trig(x,0) for x in ang]
    v=[[ca*cb,sa*cb,sign*sb],
       [-sa*cc-sign*ca*sc*sb,ca*cc-sign*sa*sc*sb,sc*cb],
       [sa*sc-sign*ca*cc*sb,-ca*sc-sign*sa*cc*sb,cc*cb]]
    p=[x*x for x in v[0]]
    s=Interval(value(h["down"][1],work["strange"]))
    b=Interval(value(h["down"][2],work["bottom"]))
    signs=cfg["signs"]
    d=[-(p[1]*signs[1]*s+p[2]*signs[2]*b)/p[0],signs[1]*s,signs[2]*b]
    a=[Interval(value(row)**2) for row in h["up"]]
    port=Interval(0)
    for j in (1,2):
        hj=sum(v[j][i]*v[0][i]*d[i]*d[i] for i in range(3))
        dj=sum(v[j][i]*v[0][i]*d[i] for i in range(3))
        port+=-3*(a[j]+a[0])/(a[0]-a[j])*hj*dj
    root=sum(Fraction(3,2)*p[i]*d[i]*d[i]*d[i] for i in range(3))
    for i in range(3):
        for j in range(3):
            g=sum(v[m][i]*a[m]*v[m][j] for m in range(3))
            root-=Fraction(3,2)*v[0][i]*v[0][j]*g*(d[i]*d[i]+d[j]*d[j])/(d[i]+d[j])
    pi=Interval(Fraction(314159,100000),Fraction(314160,100000))
    bound=(port+root)/(16*pi*pi)
    assert bound.hi < 0
    return [float(bound.lo),float(bound.hi)]


def compute():
    cfg = read(HERE/"inputs.json")
    weak = read(HERE/cfg["weak_gauge"])
    main = linear(cfg,weak)
    e0,ep = [linear(cfg,weak,x)["derivative_total"] for x in (0.,math.pi)]
    mean,cosine = (e0+ep)/2,(e0-ep)/2
    for phase in (.71,2.1,4.4):
        p = linear(cfg,weak,phase)
        assert math.isclose(p["derivative_total"],mean+cosine*math.cos(phase),abs_tol=2e-17)
    coarse, fine = [run(cfg,weak,steps) for steps in cfg["rk4_steps"]]
    error = max(abs(fine[k]["F_over_s"]-coarse[k]["F_over_s"]) for k in fine)
    assert error<1e-7 and abs(main["F"])<1e-17
    return dict(initial=main,phase_harmonic=dict(mean=mean,cosine=cosine,
             full_phase_range=[mean-abs(cosine),mean+abs(cosine)], certified_endpoint_minmax=[signed_interval(cfg,z) for z in (1,-1)]),
             scale_transport=fine,max_step_refinement_F_over_s=error)


def compare(new,old):
    if isinstance(new,dict):
        assert old.keys()==new.keys()
        for k,v in new.items():
            compare(v,old[k])
    else:
        np.testing.assert_allclose(new,old,rtol=1e-7,atol=3e-15)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(result,read(target))
        print("1566 signed root, moving port and three finite-scale transports verified")
