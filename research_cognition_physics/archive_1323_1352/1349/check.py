"""Exact reduced dynamics of the two-site gauge-invariant, finite-bank sector."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def branch(detuning, mass, time):
    # |x>, |y> with fixed chirality-X sign. Spin-L has equal X components.
    def propagate(shift):
        h = np.array([[0.,1.],[1.,shift]])
        e, v = np.linalg.eigh(h)
        return v @ (np.exp(-1j*time*e)*v[0].conj())
    a, b = propagate(detuning+mass), propagate(detuning-mass)
    left, right = (a+b)/2, (a-b)/2
    return np.array([abs(left[0])**2,abs(right[0])**2,abs(left[1])**2,abs(right[1])**2])


def record(group, matched, mu, time):
    charge, counts = ([1.,2.5],[3,6]) if group=="SU3" else ([0.,1.,3.],[1,3,5])
    symmetric = 1. if matched else .5
    out = np.zeros(4)
    for z,k in zip(charge,counts):
        for sign,weight in [(1,symmetric),(-1,1-symmetric)]:
            out += k/9*weight*branch(1.3*(z-1),.7*(mu+sign),time)
    assert abs(np.sum(out)-1)<1.e-13
    return out


def run():
    rows = []
    for mu,time in itertools.product([1.,2.],[1.4,2.2]):
        grid = np.array([[record(group,matched,mu,time) for matched in (True,False)]
                         for group in ("SU3","SO3")])
        det = np.linalg.det(grid[:,:,3])
        if mu==1:
            assert abs(det)<1.e-15
            np.testing.assert_allclose(grid[:,0,3]/grid[:,1,3],2)
        else:
            assert abs(det)>1.e-7
        p = grid[0,0]
        rows.append(dict(mu=mu,time=time,complete_probabilities=grid.tolist(),
                         group_by_bank_determinant=float(det),
                         joint_minus_marginal_product=float(p[3]-(p[2]+p[3])*(p[1]+p[3]))))
    return dict(round=1349,groups=["SU3","SO3"],banks=["aligned","orthogonal"],
                records=["x L","x R","y L","y R"],rows=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    result=run()
    path=Path(__file__).with_name("results.json")
    if parser.parse_args().write:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        saved=json.loads(path.read_text(encoding="utf8"))
        for k in ("round","groups","banks","records"):
            assert saved[k]==result[k]
        for x,y in zip(result["rows"],saved["rows"],strict=True):
            assert x.keys()==y.keys()
            for k in x:
                np.testing.assert_allclose(x[k],y[k],atol=1e-13,rtol=1e-13)
    print([(r["mu"],r["time"],r["group_by_bank_determinant"],
            np.array(r["complete_probabilities"])[:,:,3].tolist()) for r in result["rows"]])
