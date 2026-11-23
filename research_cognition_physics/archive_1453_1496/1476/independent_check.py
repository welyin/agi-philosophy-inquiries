"""Resolve common left-current Gram matrices from three virtualities in a finite pole model."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DIRECT = {
    "ud": (.97373, .00031), "cd": (.2238, .0011+.0015+.0004+.0002+.0022),
    "cs": (.9589, .0023+.0040+.0015+.0005+.0095),
    "cb": (.04057, .00031+.00095+.00058), "tb": (.988, .051),
}


def seed():
    def rot(i,j,s,delta=0):
        c=(1-s*s)**.5
        arr=np.eye(3,dtype=complex)
        arr[i,i]=arr[j,j]=c
        arr[i,j]=s*np.exp(-1j*delta)
        arr[j,i]=-arr[i,j].conjugate()
        return arr
    return rot(1,2,.040) @ rot(0,2,.0035,1.1) @ rot(0,1,.224)


def solve_graph():
    order = "uctdsb"
    lap=np.zeros((6,6))
    intervals = {edge: [(v-2*s)/1.02, (v+2*s)/.98] for edge,(v,s) in DIRECT.items()}
    for (i,j), (mag, _) in intervals.items():
        x=np.eye(6)[order.index(i)] - np.eye(6)[order.index(j)]
        lap += mag**2 * np.outer(x,x)
    resist=np.linalg.pinv(lap)
    distances={}
    for i,j in itertools.combinations(range(6),2):
        vec=np.eye(6)[i]-np.eye(6)[j]
        distances[order[i]+order[j]]=(vec @ resist @ vec)**.5
    return intervals, distances


def response(a,x,t):
    return (a*(t/(t+x))) @ a.conj().T


def coefficients(x, t):
    return np.array([np.prod(s+x)/(s*np.prod([s-y for y in t if s!=y])) for s in t])


def replay():
    windows, graph=solve_graph()
    # Illustrative constant tree pole spectrum in GeV, without RG or threshold inference.
    mu=np.array([.0022,1.27,172.5])
    md=np.array([.0047,.093,4.18])
    u=seed()
    assert np.allclose(u.conj().T@u,np.eye(3))
    momenta=np.array([100.,250.,600.])
    ts=momenta**2
    t=ts[-1]
    a=(np.sqrt(1+mu**2/t)[:,None]*u)*np.sqrt(1+md**2/t)
    edges={key:abs(a["uct".index(key[0]),"dsb".index(key[1])]) for key in DIRECT}
    assert all(windows[k][0] <= val <= windows[k][1] for k,val in edges.items())
    results = {}
    for name,mat,mass,external in [("up",a,md,mu),("down",a.conj().T,mu,md)]:
        target=mat@mat.conj().T
        rows=np.array([response(mat,mass**2,s) for s in ts])
        assert np.allclose(rows[-1],np.diag(1+external**2/t))
        weight=coefficients(mass**2,ts)
        restored=np.tensordot(weight,rows,axes=(0,0))
        assert np.linalg.norm(restored-target,2)<5e-14
        assert np.allclose(weight @ (ts[:,None]/(ts[:,None]+mass[None,:]**2)),np.ones(3))
        raw_off=target-np.diag(target.diagonal())
        extra_mass=1000.
        unseen=float(weight@(ts/(ts+extra_mass**2)))
        predicted=1-extra_mass**6*np.prod([(extra_mass**2-x**2)/extra_mass**2 for x in mass])/np.prod(ts+extra_mass**2)
        assert abs(unseen-predicted)<1e-13
        results[name]=dict(coefficients=weight.tolist(),noise_l1=sum(abs(weight)),
                           one_point_max_offdiagonal=float(np.max(abs(rows[-1]-np.diag(rows[-1].diagonal())))),
                           true_max_offdiagonal=float(np.max(abs(raw_off))),
                           unseen_TeV_pole_response=unseen)
    sing=np.linalg.svd(a,compute_uv=False)
    results.update(direct_moduli=edges, tree_diameter=max(graph.values()),
                   singular_values=sing.tolist(), maximum_norm=max(sing),
                   squared_strength_spread=sing[0]**2-sing[-1]**2)
    return results


def tidy(v):
    if isinstance(v,dict): return {k:tidy(a) for k,a in v.items()}
    if isinstance(v,list): return [tidy(a) for a in v]
    return round(float(v),11)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    opts=p.add_mutually_exclusive_group()
    opts.add_argument("--check",action="store_true")
    opts.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    val=tidy(replay())
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(val,indent=2)+"\n")
    elif args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == val
        print("1476 independent check passed")
    else:
        print(json.dumps(val,indent=2))
