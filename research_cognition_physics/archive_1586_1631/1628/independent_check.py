"""Independent charged-vacuum column audit in the positive two-column cone."""
import argparse
import json
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OUT = HERE / "independent_results.json"
# Only the outward integer interval arithmetic is reused from the earlier audit.
arithmetic = runpy.run_path(str(STAGE / "1625" / "independent_check.py"))
Box, ONE, SCALE, scaled = (arithmetic[k] for k in ["Box", "ONE", "SCALE", "scaled"])


def point(value):
    return Box(scaled(value), scaled(value))


def enclose_v():
    source = json.loads((STAGE / "1619/inputs.json").read_bytes())["pole_MeV"]
    roots = []
    for flavor in ["electron", "muon"]:
        center, sigma = map(scaled, source[flavor])
        roots.append(Box(center-2*sigma, center+2*sigma).sqrt())
    e, mu = roots
    tau = point(2)*(e+mu) + (point(3)*(e*e+point(4)*e*mu+mu*mu)).sqrt()
    norm = (e*e + mu*mu + tau*tau).sqrt()
    return e/norm, mu/norm, tau/norm


def residual(order, x12, x13, x23, v):
    e, mu, tau = v
    s12, c12 = x12.sqrt(), (ONE-x12).sqrt()
    s13, c13 = x13.sqrt(), (ONE-x13).sqrt()
    s23, c23 = x23.sqrt(), (ONE-x23).sqrt()
    if order == "NO":
        return (tau*(s12*s23 - c12*c23*s13) - mu*(s12*c23 + c12*s23*s13)
                - e*c12*c13)
    return c13*(tau*c23-mu*s23)-e*s13


def angular_threshold(order, x12, x13, v):
    lo, hi = 0, SCALE
    while hi-lo > 8:
        mid = (lo+hi)//2
        test = residual(order,x12,x13,Box(mid,mid),v)
        if test.lo <= 0 <= test.hi:
            break
        if (test.lo > 0) == (order == "NO"):
            hi = mid
        else:
            lo = mid
    # lo,hi bracket the zero for every mass point in the outward box.
    return [lo/SCALE, hi/SCALE]


def abstract_extremizers():
    # The background is unrestricted in this sharp two-level inequality.
    for t in [.01,.1,.366,.55,.70,2**(-.5),.8,.99]:
        c = np.sqrt(1-t*t)
        n = np.array([c,t])
        d = max(0.,t*c-t*t)
        h = np.diag([d,1+d])
        hn = h@n
        background = np.outer(hn,hn)/(n@hn)
        target = 2*t*c if t<=2**(-.5) else 1.
        assert abs(np.ptp(np.linalg.eigvalsh(background))-target) < 3e-14
        assert np.linalg.norm((h-background)@n) < 2e-15
        assert np.linalg.eigvalsh(h-background)[0] > -3e-15


def fit_numerical_witnesses():
    """Find a charged-column eigenplane through a basis-invariant Krylov determinant."""
    masses = json.loads((STAGE/"1619/inputs.json").read_bytes())["pole_MeV"]
    e, mu = np.sqrt([masses["electron"][0],masses["muon"][0]])
    tau = 2*(e+mu)+np.sqrt(3*(e*e+4*e*mu+mu*mu))
    amps = np.array([e,mu,tau])/np.linalg.norm([e,mu,tau])
    q = np.diag((amps/amps[2])**2)
    x,y,z = (.30,.022,.50)
    sx,sy,sz = np.sqrt([x,y,z])
    cx,cy,cz = np.sqrt([1-x,1-y,1-z])
    U = np.array([[cx*cy,sx*cy,sy],
                  [-sx*cz-cx*sz*sy,cx*cz-sx*sz*sy,sz*cy],
                  [sx*sz-cx*cz*sy,-cx*sz-sx*cz*sy,cz*cy]])
    results = {}
    for order,bracket in [("NO",(-.0006,-.0005)),("IO",(-.0055,-.0054))]:
        levels = np.array([0,.000075, .0025 if order=="NO" else -.0025])
        levels -= levels.min()
        target = (U*levels)@U.T
        v = amps*np.array([1,1 if order=="NO" else -1,1])

        def krylov(b):
            block = target-b*q
            first = block@v
            return np.linalg.det(np.column_stack((v,first,block@first)))

        lo,hi = bracket
        assert krylov(lo)*krylov(hi)<0
        for _ in range(60):
            mid = (lo+hi)/2
            if krylov(lo)*krylov(mid)<=0:
                hi=mid
            else:
                lo=mid
        b=(hi+lo)/2
        eig,frame=np.linalg.eigh(target-b*q)
        d=eig[0]
        assert abs(frame[:,0]@v)<3e-12
        source=target-b*q-d*np.eye(3)
        inverse=(frame[:,1:]*(1/(eig[1:]-d)))@frame[:,1:].T
        a=1/(v@inverse@v)
        ew,rotation=np.linalg.eigh(source-a*np.outer(v,v))
        assert np.max(abs(ew[:2]))<1e-16 and a>0 and ew[2]>0
        w=np.sqrt(ew[2])*rotation[:,2]
        background=d*np.eye(3)+b*q
        shift=-np.min(np.diag(background))
        assert shift>0
        rebuilt=background+shift*np.eye(3)+a*np.outer(v,v)+np.outer(w,w)
        error=np.linalg.norm(rebuilt-(target+shift*np.eye(3)),ord=2)
        assert error<1e-16
        results[order]={"slope_eV2":round(b,12),"column_strength_eV2":round(a,12),
                        "lightest_mass_squared_eV2":round(shift,12),
                        "reconstruction_operator_error_upper_eV2":1e-16}
    return results


def calculate():
    v = enclose_v()
    box = json.loads((STAGE / "1625/inputs.json").read_bytes())["box"]
    x12 = point(box["sin2_theta12"][0])
    x13 = point(box["sin2_theta13"][1])
    assert (v[2]*x12.sqrt()-v[1]*(ONE-x12).sqrt()*x13.sqrt()).lo > 0
    result = {"family_source":"1619 heavy Koide root; independent two-width charged inputs",
              "orders":{}}
    for order, index in (("NO",0),("IO",1)):
        x23 = point(box["sin2_theta23"][index])
        r = residual(order,x12,x13,x23,v)
        assert r.lo > 0 and r.hi < SCALE/2
        rlo = Box(r.lo,r.lo)
        fraction = point(2)*rlo*(ONE-rlo*rlo).sqrt()
        gap = point(box["abs31_eV2"][0])
        if order == "NO":
            gap = gap*point(box["solar_over_abs31"][0])
        result["orders"][order] = {
            "orthogonality_boundary_sin2theta23_enclosure":angular_threshold(order,x12,x13,v),
            "overlap_lower":rlo.lo/SCALE, "light_mass_squared_over_a_lower":(rlo*rlo).lo/SCALE,
            "background_spread_over_gap_lower":fraction.lo/SCALE,
            "background_spread_eV2_lower":(fraction*gap).lo/SCALE}
    abstract_extremizers()
    result['fit_numerical_witnesses']=fit_numerical_witnesses()
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = calculate()
    if args.check:
        assert results == json.loads(OUT.read_bytes())
    else:
        with OUT.open("x",encoding="utf-8") as out:
            json.dump(results,out,indent=2)
            out.write("\n")
    print("1628 independent PASS: charged error box, phase triangle and spectral spread")
