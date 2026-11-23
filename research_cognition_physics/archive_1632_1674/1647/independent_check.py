"""Independent axis and helicity-row obstruction for real bipartite exchanges."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def interval(row):
    return np.array([row[0]-2*sum(row[1:]),row[0]+2*sum(row[1:])])


def modulus_box():
    first = read("../1381/inputs.json")
    tree = read("../1445/inputs.json")
    us = interval(first["K_product"])/interval(first["K_form"])[::-1]
    ub = interval(first["ub_over_cb"])*interval(tree["BelleII_Vcb"])
    return np.array([interval(first["ud_beta"]),us,ub])


def first_row_distance():
    box = modulus_box()
    smin,bmax = box[1,0],box[2,1]
    # For the nearest helicity row, y is its smaller endpoint and sqrt(2 y (1-y))
    # is its middle component. Their gaps to the s lower bound and b upper bound agree.
    T = smin+bmax
    y = T*T/(1+T+math.sqrt((1+T)**2-3*T*T))
    hel = np.array([1-y,math.sqrt(2*y*(1-y)),y])
    nearest = np.clip(hel,box[:,0],box[:,1])
    delta = y-bmax
    assert 0<y<1/3 and 0<delta<.03
    np.testing.assert_allclose(abs((hel-nearest)[1:]),delta,atol=2e-14)
    assert abs(hel[0]-nearest[0])<delta
    # For larger y, the smallest endpoint already enforces this gap. If y>=1/3
    # no pairing has a sufficiently large first component. Any repeated pair stays apart.
    zeq = min(max(a[0]-b[1],b[0]-a[1],0)/2 for a,b in itertools.combinations(box,2))
    assert zeq>delta and box[0,0]-max(2/3,1/math.sqrt(2))>delta
    return {"row_intervals":box.tolist(),"sharp_Linfinity_distance_to_box":delta,
            "axis_row_lower_bound":zeq,"optimizing_helicity_row":hel.tolist(),
            "nearest_rectangular_point":nearest.tolist(),"axis_cos":1-2*y}


def offdiag(i,j):
    A = np.zeros((3,3))
    A[i,j] = A[j,i] = 1.
    return A


def pulse(H,a):
    d,W = np.linalg.eigh(H)
    return (W*np.exp(-1j*a*d))@W.T


def pulse_test():
    P = np.diag([1.,1j,1.])
    A,B = offdiag(0,1),offdiag(1,2)
    actual = []
    vectors = []
    for idx in range(2):
        U = np.eye(3,dtype=complex)
        scalar = 0.
        for n in range(17+idx):
            H = math.sin((idx+1.2)*(n+1)) * A + math.cos((idx+.35)*(n+2)) * B
            t = math.sin(n+.4)*.55
            d = math.cos(n+idx)*.19
            U = pulse(H+d*np.eye(3),t)@U
            scalar += t*d
        R = np.exp(1j*scalar)*P.conj().T@U@P
        np.testing.assert_allclose(R.imag,0,atol=2e-14)
        np.testing.assert_allclose(R.T@R,np.eye(3),atol=2e-14)
        assert abs(np.linalg.det(R)-1)<5e-14
        ev,V = np.linalg.eig(R)
        # Order eigenlines as +, axis, -, use arbitrary eigenvector phases.
        order = np.argsort(ev.imag)[::-1]
        assert abs(ev[order[0]].imag)>.01 and abs(ev[order[1]].imag)<1e-13
        vectors.append(V[:,order])
        # Axis extraction and its independently normalized *real* orientation.
        W = R.real
        n = np.array([W[2,1]-W[1,2],W[0,2]-W[2,0],W[1,0]-W[0,1]])
        actual.append(n/np.linalg.norm(n))
    V = vectors[0].conj().T@vectors[1]
    np.testing.assert_allclose(V.conj().T@V,np.eye(3),atol=4e-13)
    c = float(np.dot(*actual))
    t = math.sqrt((1-c*c)/2)
    d = np.array([[(1+c)/2,t,(1-c)/2],[t,abs(c),t],[(1-c)/2,t,(1+c)/2]])
    np.testing.assert_allclose(abs(V),d,atol=3e-13)
    J = (V[0,0]*V[1,1]*V[0,1].conj()*V[1,0].conj()).imag
    assert abs(J)<1e-13
    return {"axis_dot":c,"absolute_weak_current":np.round(abs(V),10).tolist(),
            "quartet_magnitude_below":1e-13}


def closure_test():
    X = [-1j*offdiag(*pair) for pair in [(0,1),(1,2),(2,0)]]
    def comm(a,b):
        return a@b-b@a
    K = [comm(X[i],X[(i+1)%3]) for i in range(3)]
    # K[0] is skew on 0,2; K[1] on 0,1; K[2] on 1,2.
    diag = [comm(X[0],K[1]),comm(X[1],K[2])]
    def rank(basis):
        return int(np.linalg.matrix_rank([np.r_[x.real.ravel(),x.imag.ravel()] for x in basis]))
    assert rank(X[:2]+[K[0]])==3
    assert rank(X+K+diag)==8
    return {"bipartite_generator_dimension":3,"independent_triangle_dimension":8}


def solve():
    return {"sharp_first_row_test":first_row_distance(),"timed_pulses":pulse_test(),
            "generator_audit":closure_test()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    path = HERE/"independent_results.json"
    result = solve()
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1647 independent PASS: scalar drift, spectral lines, sharp row distance, generated Lie algebra")
    else:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")

