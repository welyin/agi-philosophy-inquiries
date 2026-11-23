"""Compare the full colour-space anticommutator with independently overlaid singlets."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

PAIRS = list(itertools.combinations(range(4), 2))


def generators(name):
    e = np.eye(3, dtype=complex)
    if name == "SO3":
        return [1j*(np.outer(e[a],e[b])-np.outer(e[b],e[a]))
                for a,b in ((0,1),(1,2),(2,0))]
    t = []
    for a,b in itertools.combinations(range(3),2):
        m = np.outer(e[a],e[b])
        t.extend([(m+m.T)/2, 1j*(m-m.T)/2])
    return t + [np.diag([1.,-1.,0.])/2,np.diag([1.,1.,-2.])/(2*np.sqrt(3))]


def pair_action(t, i, j, state):
    # One generator on each specified tensor factor.
    out = np.zeros_like(state)
    for a in t:
        v = state
        for k in (i,j):
            op = -a.T if k % 2 else a
            v = np.moveaxis(np.tensordot(op,v,axes=([1],[k])),0,k)
        out += v
    return out


def test(name):
    t = generators(name)
    cr = float(np.trace(sum(a@a for a in t)).real/3)
    dim = len(t)
    singlet = np.einsum("ij,kl->ijkl",np.eye(3),np.eye(3)).astype(complex)/3
    actions = np.stack([pair_action(t,i,j,singlet).reshape(-1) for i,j in PAIRS])
    psi = singlet.reshape(-1)
    assert np.isclose(np.linalg.norm(psi),1)
    is_cross = np.array([i//2 != j//2 for i,j in PAIRS])
    assert np.max(abs(actions[is_cross] @ psi.conj())) < 1e-13
    assert np.isclose(np.linalg.norm(actions[1])**2, cr**2/dim)
    hard = np.array([[1,0,0],[-1,0,0],[.5,np.sqrt(.75),0],[-.5,-np.sqrt(.75),0]])
    def weights(direction):
        return np.array([(1-hard[i]@hard[j])/((1-hard[i]@direction)*(1-hard[j]@direction))
                         for i,j in PAIRS])
    def compare(w,v):
        k, l = -2*w@actions, -2*v@actions
        independent = (2*cr)**2*(w[0]+w[-1])*(v[0]+v[-1])
        signal = (np.vdot(k,l).real-independent)/(4*cr**2*w[0]*v[-1])
        def contrast(s):
            return s[1]-s[2]-s[3]+s[4]
        expected = contrast(w)*contrast(v)/(dim*w[0]*v[-1])
        np.testing.assert_allclose(signal,expected,rtol=1e-10,atol=1e-11)
        return float(signal)
    rng = np.random.default_rng(1347)
    for _ in range(40):
        n,m = rng.normal(size=(2,3))
        compare(weights(n/np.linalg.norm(n)),weights(m/np.linalg.norm(m)))
    values = {}
    for angle in (20.,30.):
        a = np.array([1-hard[i]@hard[j] for i,j in PAIRS])
        integrated = np.log1p(a*np.tan(np.radians(angle))**2/2)
        values[str(int(angle))] = compare(integrated,integrated)
    return dict(normalized_cross=values,normal_direction=compare(weights([0,0,1]),weights([0,0,-1])),
                fourth_pointer_abs=(4+4/dim)/3)


def run():
    u,s = test("SU3"),test("SO3")
    for region in ("20","30"):
        assert np.isclose(s["normalized_cross"][region]/u["normalized_cross"][region],8/3)
    return dict(round=1347, su3=u,so3=s, randomized_direction_pairs_per_group=40)


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    result=run()
    out=Path(__file__).with_name("results.json")
    if p.parse_args().write:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        old=json.loads(out.read_text(encoding="utf8"))
        for group in ("su3","so3"):
            for a in ("20","30"):
                np.testing.assert_allclose(old[group]["normalized_cross"][a],
                                           result[group]["normalized_cross"][a],atol=1e-12)
            for k in ("normal_direction","fourth_pointer_abs"):
                assert abs(old[group][k]-result[group][k])<1e-12
        assert all(old[k]==result[k] for k in ("round","randomized_direction_pairs_per_group"))
    print(json.dumps(result))
