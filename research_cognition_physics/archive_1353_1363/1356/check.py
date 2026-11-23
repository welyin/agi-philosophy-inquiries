"""Finite charged matter processes and indistinguishable center lifts."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np


def complex_response(n, t, u, z):
    return u**n * z**t


def record(transport, angle):
    path = np.array([1, np.exp(1j*angle)*transport])/np.sqrt(2)
    ports = (np.array([[1, 1], [1, -1]])/np.sqrt(2)) @ path
    values = np.abs(ports)**2
    np.testing.assert_allclose(values.sum(), 1, atol=2e-14)
    return float(values[1])


def bounds(value):
    # |delta p| <= |delta angle|/2 + (1-v)/2 + readout error.
    epsilon = .02/2 + .04/2 + .005
    return [max(0, value-epsilon), min(1, value+epsilon)]


def calculate():
    root = np.exp(2j*np.pi/3)
    sources = {"unit": (3,0), "upper": (2,1), "lower": (-1,1),
               "pair_ud": (1,2), "external_singlet": (1,0)}
    loops = [(1,root), (root.conjugate(),1)]
    outputs = {}
    for name,(n,t) in sources.items():
        rows = []
        for u,z in loops:
            transport = complex_response(n,t,u,z)
            probabilities = [record(transport,angle) for angle in (0,np.pi/2)]
            rows.append({"minus": probabilities,"bounds": [bounds(p) for p in probabilities]})
        outputs[name] = rows
        if name!="external_singlet":
            np.testing.assert_allclose(rows[0]["minus"],rows[1]["minus"],atol=2e-14)
    assert all(a[1] < b[0] or b[1] < a[0] for a,b in
               zip(outputs["external_singlet"][0]["bounds"],outputs["external_singlet"][1]["bounds"]))

    # Arbitrary intertwiners between same-N input and output multiplicities.
    rng = np.random.default_rng(1356)
    labels_in = np.array([-2, 1, 1, 3, 4, 4])
    labels_out = np.array([-2,-2, 1, 3, 3, 4, 4])
    free = rng.normal(size=(7,6))+1j*rng.normal(size=(7,6))
    block = free * (labels_out[:,None] == labels_in)
    block /= max(1, np.linalg.norm(block,ord=2))
    residual = np.linalg.norm((root**(-labels_out))[:,None]*block - block*root**(-labels_in))
    assert residual < 1e-13
    # Coarse center-preserving alone admits n and n+3 even at equal dimension.
    assert abs(root**4 - root) < 1e-13 and 4 != 1
    # Exact equal charge, distinct defects, perfect SO(3)-scalar conversion.
    delta = np.eye(3).reshape(9) / np.sqrt(3)
    bridge = delta.reshape(1,9)
    center_defect = float(np.linalg.norm(bridge - root**2*bridge,ord=2))
    assert abs(center_defect - np.sqrt(3)) < 1e-13
    phase_checks = 0
    for word in product(list(sources)[:-1], repeat=3):
        n=sum(sources[s][0] for s in word)
        t=sum(sources[s][1] for s in word)
        for u,z in loops:
            np.testing.assert_allclose(complex_response(n,t,u,z),(u/z)**n,atol=1e-13)
            phase_checks += 1
    return {"round":1356, "two_paths":outputs, "finite_error_bound":.035,
            "three_factor_tests":phase_checks,"vertex_norm_for_unlocked_sink":center_defect}


def verify(actual, expected):
    if isinstance(actual,dict):
        assert actual.keys()==expected.keys()
        for key in actual:
            verify(actual[key],expected[key])
    elif isinstance(actual,list):
        assert len(actual)==len(expected)
        for x,y in zip(actual,expected):
            verify(x,y)
    else:
        np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    data=calculate()
    path=Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x",encoding="utf8") as out:
            json.dump(data,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        verify(data,json.loads(path.read_bytes()))
    print({k:[r["minus"] for r in v] for k,v in data["two_paths"].items()})
