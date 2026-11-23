"""Carry physical column phases through the same Higgs selector and two real Majorana sources."""
import argparse
import itertools
import json
from math import cos, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rot(i, j, sine):
    q = np.eye(3)
    q[i,i] = q[j,j] = sqrt(1-sine**2)
    q[i,j], q[j,i] = sine, -sine
    return q


def amplitude_box(a, b, cross):
    """Nonnegative lengths; extrema of the quadratic on each colored rectangle."""
    def vals(t):
        candidates = list(itertools.product(a,b))
        candidates += [(x, float(np.clip(-t*x,*b))) for x in a]
        candidates += [(float(np.clip(-t*y,*a)), y) for y in b]
        return [sqrt(max(0.,x*x+y*y+2*x*y*t)) for x,y in candidates]
    return [min(vals(min(cross))),max(vals(max(cross)))]


def calculate():
    get = lambda rel: json.loads((HERE/rel).read_text(encoding="utf-8"))
    cfg=get("inputs.json")
    osc=get(cfg["oscillation"])
    alpha, err=get(cfg["quark_phase"])["quark_sines_and_phase"][3]
    cos_box=sorted(cos(2*(alpha+j*cfg["CKM_error_multiple"]*err)) for j in (-1,1))
    output={"round":1426,"absolute_interference_axis1":sorted(abs(t) for t in cos_box),
            "minimum_X_abs_separation_at_registered_CKM":1-max(abs(t) for t in cos_box),"orderings":{}}
    for ordering in ("NO","IO"):
        par=osc[ordering]
        x, xl, xh=par["s12"]
        y, yl, yh=par["s13"]
        ds, dsl, dsh=par["solar"]
        da, dal, dah=par["atm_abs"]
        ids = [1,2] if ordering=="NO" else [0,1]
        masses = np.array([0.,sqrt(ds),sqrt(da)]) if ordering=="NO" else np.array([sqrt(da-ds),sqrt(da),0.])
        heavy = np.array(cfg["heavy_GeV"])*1e9
        rows={}
        for swap,reverse,parity in itertools.product([False,True],[False,True],[1,-1]):
            axes = np.eye(3)[[1,0,2]] if swap else np.eye(3)
            charged = axes @ (rot(1,2,sqrt(cfg["representative_s23sq"])) @ rot(0,2,(-1 if reverse else 1)*sqrt(y))).T
            neutral = axes @ rot(0,1,sqrt(x))
            charged_left = np.diag([np.exp(.5j*alpha),1,1]) @ charged
            dirac = np.diag([np.exp(-.5j*alpha),1,1]) @ neutral[:,ids] @ np.diag(np.sqrt(heavy*masses[ids]))
            signed_heavy = np.diag(heavy*[1,parity])
            mnu = -dirac @ np.linalg.solve(signed_heavy,dirac.T)
            takagi = np.ones(3,complex)
            takagi[ids] = [1j, 1j if parity==1 else 1.]
            un = np.diag([np.exp(-.5j*alpha),1,1]) @ neutral @ np.diag(takagi)
            pmns = charged_left.conj().T @ un
            assert np.allclose(un.conj().T@mnu@un.conj(),np.diag(masses),atol=5e-16)
            mass_basis = charged_left.conj().T@mnu@charged_left.conj()
            amp=abs((pmns[0]**2)@masses)
            assert np.isclose(amp,abs(mass_basis[0,0]))
            beta2=(abs(pmns[0])**2)@(masses**2)
            assert np.isclose(beta2,np.sum(abs(mass_basis[0])**2))
            if ordering=="NO":
                a=(1-y)*x*masses[1]
                b=y*masses[2]
                ratio=(amp**2-a*a-b*b)/(2*a*b)
                assert np.isclose(ratio,parity*(1 if swap else cos(2*alpha)))
            else:
                assert np.isclose(amp,(1-y)*abs((1-x)*masses[0]+parity*x*masses[1]))
            key=("axis2" if swap else "axis1")+("_same" if parity==1 else "_opposite")
            if key in rows:
                assert np.isclose(rows[key],amp*1000)
            else:
                rows[key]=float(amp*1000)
        if ordering=="NO":
            a_box=[1000*(1-yh)*xl*sqrt(dsl),1000*(1-yl)*xh*sqrt(dsh)]
            b_box=[1000*yl*sqrt(dal),1000*yh*sqrt(dah)]
            bands={key:amplitude_box(a_box,b_box,
                    ([sign] if "axis2" in key else [sign*t for t in cos_box]))
                   for key in rows for sign in [1 if "same" in key else -1]}
        else:
            # Deliberately outer-bound dependence on the shared solar angle and mass gap.
            a_box=[1000*(1-yh)*(1-xh)*sqrt(dal-dsh),1000*(1-yl)*(1-xl)*sqrt(dah-dsl)]
            b_box=[1000*(1-yh)*xl*sqrt(dal),1000*(1-yl)*xh*sqrt(dah)]
            bands={key:amplitude_box(a_box,b_box,[1 if "same" in key else -1]) for key in rows}
        assert all(bands[key][0]<=v<=bands[key][1] for key,v in rows.items())
        output["orderings"][ordering]={"central_meV":rows,"outer_marginal_boxes_meV":bands}
    return output


def matches(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(matches(a[k],b[k]) for k in a)
    return bool(np.allclose(a,b,rtol=1e-10,atol=1e-12))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=calculate()
    if args.check:
        assert matches(result,json.loads((HERE/"results.json").read_text(encoding="utf-8")))
        print("1426 Takagi phases, physical selector and independent mass bands verified")
    else:
        print(json.dumps(result,indent=2))
