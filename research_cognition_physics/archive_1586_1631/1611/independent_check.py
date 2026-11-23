"""Two triplet thresholds, neutral Yukawa area and a distinct critical scale."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
LOW = json.loads((HERE.parent/"1608"/"inputs.json").read_text(encoding="utf-8"))
OSC = json.loads((HERE.parent/"1370"/"inputs.json").read_text(encoding="utf-8"))
BASE = 200.
V = 246.22
PI2 = math.pi**2
B = np.array([4.1, -19/6, -7.])
DB = np.array([4/15, 0., 2/3])
DS = np.array([0., 4/3, 0.])


def gauge_parameters(md, shift):
    ref = LOW["reference"]
    g = np.array([math.sqrt(5/3)*ref["gY"], ref["g2"], ref["g3"]*(1+shift)])
    inv = 8*PI2/g**2
    logs = np.log(np.array([md, 2100.])/BASE)
    # Three linear rows: common inverse coupling, full log and summed triplet exposure.
    sol = np.linalg.solve(np.column_stack((np.ones(3), B+2*DB, DS)),
                          inv+sum(logs)*DB)
    common, meeting, exposure = sol
    avg_log = meeting-exposure/2
    assert 0 < avg_log < meeting and common>0
    return inv, logs, meeting, avg_log


def gauge(t, model):
    inv, logs, triplets = model
    strength = inv - B*t
    for start in logs:
        strength = strength-DB*max(0.,t-start)
    for start in triplets:
        strength = strength-DS*max(0.,t-start)
    return 8*PI2/strength


def terminal(t, model, step=.16):
    """Integrate y itself with analytic affine 1/g² at each threshold."""
    nodes=sorted({0., t, *[float(n) for n in (*model[1],*model[2]) if 0<n<t]},reverse=True)
    g1, g2, _ = gauge(t, model)
    gy = .6*g1
    q = math.sqrt(2*g2*g2+(g2+gy)**2)/4
    yt, lam = math.sqrt(q), 0.
    for high,low in zip(nodes,nodes[1:]):
        hi = gauge(high,model)
        lo = gauge(low,model)
        inverse = 1/hi
        slope = (1/lo-inverse)/(low-high)
        n = int(math.ceil((high-low)/step))
        dt = (low-high)/n
        slope = tuple(slope)
        inverse = tuple(inverse)
        def rhs(rel,y,x):
            y1, w, col = [1/(a+b*rel) for a,b in zip(inverse,slope)]
            hy=.6*y1
            dy=y*(4.5*y*y-17*hy/12-9*w/4-8*col)
            dx=24*x*x+(12*y*y-9*w-3*hy)*x-6*y**4+3*(2*w*w+(hy+w)**2)/8
            return dy/(16*PI2),dx/(16*PI2)
        for j in range(n):
            u = dt*j
            a,b=rhs(u,yt,lam)
            c,d=rhs(u+dt/2,yt+dt*a/2,lam+dt*b/2)
            e,f=rhs(u+dt/2,yt+dt*c/2,lam+dt*d/2)
            p,h=rhs(u+dt,yt+dt*e,lam+dt*f)
            yt+=dt*(a+2*c+2*e+p)/6
            lam+=dt*(b+2*d+2*f+h)/6
    return yt,lam


def separate_root(model):
    left,right=math.log(1e8/BASE),math.log(1e16/BASE)
    target=LOW["reference"]["lambda"]
    assert terminal(left,model)[1]<target<terminal(right,model)[1]
    for _ in range(40):
        mid=(left+right)/2
        if terminal(mid,model)[1]>target:
            right=mid
        else:
            left=mid
    t=(left+right)/2
    y,l=terminal(t,model,.08)
    assert abs(l-target)<5e-9
    dy = 1e-5
    def at(x,yt):
        g1,g2,_=gauge(x,model)
        gy=.6*g1
        return (-6*yt**4+3*(2*g2**2+(gy+g2)**2)/8)/(16*PI2)
    g1,g2,gs=gauge(t,model)
    hy=.6*g1
    yh=(math.sqrt(2*g2*g2+(g2+hy)**2)/4)**.5
    yp=yh*(4.5*yh*yh-17*hy/12-9*g2/4-8*gs)/(16*PI2)
    curvature=(at(t+dy,yh+dy*yp)-at(t-dy,yh-dy*yp))/(2*dy)
    assert curvature>0
    return {"scale_GeV":BASE*math.exp(t),"yt_200":y,"yt_at_double_zero":yh,
            "beta_lambda_prime":curvature}


def product(m1,m2,mean):
    return 2e-9*mean*math.sqrt(m1*m2)/V**2


def masses(order, solar, atm):
    if order=="NO":
        return math.sqrt(solar),math.sqrt(atm)
    return math.sqrt(atm-solar),math.sqrt(atm)


def compute():
    md0=1100
    shift=2*LOW["marginal_width"]["alpha_s"]*LOW["matching"]["g3_per_alpha"]
    rows=[]
    for md,delta in itertools.product([md0,1600],[-shift,0.,shift]):
        inv,logs,t,avg=gauge_parameters(md,delta)
        for frac in (0.,.5,1.):
            d=frac*min(avg,t-avg)
            model=(inv,logs,(avg-d,avg+d))
            assert np.ptp(gauge(t,model))<1e-13
            top,lam=terminal(t,model)
            rows.append({"B1_GeV":md,"g3_relative_shift":delta,"split_fraction":frac,
                         "scale_GeV":BASE*math.exp(t),"mean_GeV":BASE*math.exp(avg),
                         "yt_200":top,"lambda_200":lam})
    inv,logs,t,avg=gauge_parameters(md0,0.)
    mean=BASE*math.exp(avg)
    amins,amaxs=[f(r["mean_GeV"] for r in rows) for f in (min,max)]
    spectrum={}
    for ordering in ("NO","IO"):
        o=OSC[ordering]
        light=masses(ordering,o["solar"][0],o["atm_abs"][0])
        prods=[product(*masses(ordering,so,at),m) for so,at,m in itertools.product(
            o["solar"][1:],o["atm_abs"][1:],(amins,amaxs))]
        y_area=product(*light,mean)
        samples={}
        # Explicit Takagi pair, a complex unitary 3x2 left isometry and unequal masses.
        left=np.fft.fft(np.eye(3))[:,:2]/math.sqrt(3)
        heavy=np.array([mean/4,4*mean])
        for z in (0.,.7+1.6j):
            rot=np.array([[np.cos(z),np.sin(z)],[-np.sin(z),np.cos(z)]])
            yy=1j*math.sqrt(2)/V*left@np.diag(np.sqrt(1e-9*np.array(light)))@rot@np.diag(np.sqrt(heavy))
            sing=np.linalg.svd(yy,compute_uv=False)
            eff=-(V**2/2)*yy@np.diag(1/heavy)@yy.T
            assert np.allclose(eff,left@np.diag(1e-9*np.array(light))@left.T,rtol=1e-12,atol=1e-24)
            assert math.isclose(math.prod(sing),y_area,rel_tol=1e-12)
            samples[str(z)]=sing.tolist()
        spectrum[ordering]={"central_masses_eV":light,"central_product":y_area,
                            "full_product_range":[min(prods),max(prods)],
                            "central_lower_bound_smax":math.sqrt(y_area),
                            "central_mixing_singular_product":math.sqrt(math.prod(light))*1e-9/mean,
                            "illustrative_singular_values":samples}
    assert min(r["lambda_200"] for r in rows)>0.130124285
    return {"threshold_tests":rows,"neutrino":spectrum,
            "distinct_central_double_zero":separate_root((inv,logs,(avg,avg)))}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [rounded(v) for v in x]
    return float(f"{x:.8g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=rounded(compute())
    if args.write:
        with OUT.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert result==json.loads(OUT.read_text("utf-8"))
    print("1611 independent threshold, Yukawa-area and double-zero checks passed")
