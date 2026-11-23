"""Check a charged-scalar holdout and a blind operator using full weak doublets."""
import argparse
import json
from math import cos,sin,atan
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def operators(h1,h2):
    x=np.vdot(h1,h2)
    return h2*x,h2*x*np.vdot(h1,h1),h1*x*x


def run():
    rng=np.random.default_rng(1437)
    for _ in range(50):
        h1,h2=np.zeros(2,complex),np.zeros(2,complex)
        h1[1],h2[1]=rng.normal(size=2)+1j*rng.normal(size=2)
        _,a,b=operators(h1,h2)
        assert np.allclose(a,b)
    rows=[]
    for t in [.5,1.,2.]:
        angle=atan(t)
        c,s=cos(angle),sin(angle)
        v1,v2=np.array([0.,c]),np.array([0.,s])
        base,f1,f2=operators(v1,v2)
        step=1e-6
        def evaluate(x, direction):
            a,b=direction
            return operators(v1+np.array([a*x,0]),v2+np.array([b*x,0]))
        def diff(direction):
            return [(plus[0]-minus[0])/(2*step)
                    for plus,minus in zip(evaluate(step,direction),evaluate(-step,direction))]
        deriv=diff((-s,c))
        gold=diff((c,s))
        assert abs((gold[1]-gold[2])/f1[1])<1e-9
        shift=(deriv[1]-deriv[2])/f1[1]
        assert np.isclose(deriv[0]/base[1],1/t) and np.isclose(shift,t+1/t)
        # Add the identically neutral-blind operator to cancel the physical charged muon vertex.
        deform=-1/(1+t*t)
        assert abs(deriv[0]/base[1]+deform*shift)<1e-9
        rows.append({"tan_beta":t,"charged_over_unit_original":round(1/t,9),
                     "neutral_A_mu_over_tau":round(2+t*t,9),
                     "difference_increment_per_unit":round(shift,9),
                     "relative_deformation_for_charged_zero":round(deform,9)})
    return {"neutral_complex_tests":50,"cases":rows}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group()
    g.add_argument("--write-new",action="store_true")
    g.add_argument("--check",action="store_true")
    out=run()
    if p.parse_args().write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert out==json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        print("1437 independent full-doublet charged response and exact blind operator verified")
