"""Independent finite reference replay. All reported fidelities use squared overlap."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def face_optimum(p, rate):
    """Enumerate saturated box faces, maximizing over the unpinned simplex by Cauchy."""
    caps = np.roll(p, 1)
    candidates = []
    for mask in itertools.product((False, True), repeat=len(p)):
        fixed = np.array(mask)
        w = np.zeros_like(p)
        w[fixed] = caps[fixed]
        leftover = rate - w.sum()
        if leftover < -1e-13:
            continue
        if np.any(~fixed):
            w[~fixed] = leftover * p[~fixed] / sum(p[~fixed])
        if abs(w.sum() - rate) < 1e-13 and np.all(w <= caps + 1e-13):
            amplitude = sum(np.sqrt(np.maximum(0, w) * p))
            candidates.append(amplitude**2/rate)
    return min(1.0, max(candidates))


def root_psd(a):
    eig, vec = np.linalg.eigh(a)
    return (vec * np.sqrt(np.maximum(eig, 0))) @ vec.conj().T


def mixed_witness():
    rng = np.random.default_rng(1354)
    p = np.array([1/2, 1/3, 1/6])
    # Each conditional ket is a 2x2 coefficient matrix (accessible R, inaccessible E).
    blocks = []
    for _ in p:
        a = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
        blocks.append(a/np.linalg.norm(a))
    actual = 0j
    upper = 0.
    no_partner_bound = 0.
    unitary = np.zeros((6,6), complex)
    for j, a in enumerate(blocks):
        target = (j+1)%len(p)
        b = blocks[target]
        left, singular, vh = np.linalg.svd(a @ b.conj().T)
        polar = vh.conj().T @ left.conj().T
        unitary[2*target:2*target+2,2*j:2*j+2] = polar
        coefficient = np.sqrt(p[j]*p[target])
        actual += coefficient*np.vdot(b,polar @ a)
        no_partner_bound += coefficient
        sqrta = root_psd(a.conj().T @ a)
        sqrtb = root_psd(b.conj().T @ b)
        root_fidelity = np.linalg.svd(sqrta @ sqrtb,compute_uv=False).sum()
        np.testing.assert_allclose(singular.sum(),root_fidelity,atol=5e-14)
        upper += coefficient*root_fidelity
    np.testing.assert_allclose(unitary.conj().T @ unitary,np.eye(6),atol=5e-14)
    np.testing.assert_allclose(actual,upper,atol=5e-14)
    assert 0 < upper < no_partner_bound
    return {"unrestricted_charge_bound":float(no_partner_bound),"joint_polar_attained":float(upper)}


def calculate():
    p = np.array([1/2,1/3,1/6])
    common = sum(np.sqrt(p*np.roll(p,1)))
    rates = [1/3,1/2,3/4,1.]
    fidelities = [face_optimum(p,s) for s in rates]
    np.testing.assert_allclose(fidelities[-1],common**2,atol=1e-14)
    probs = [s*(1-f)/2 for s,f in zip(rates,fidelities)]
    assert probs[1] < 2*.005 < probs[2]
    # A nonzero shift yields orthogonal joint states for a perfectly charge-tagged E.
    psi=np.diag(np.sqrt(p))
    joint=[abs(np.vdot(psi,np.roll(psi,t,axis=0)))**2 for t in [1,2,3]]
    np.testing.assert_allclose(joint,[0,0,1],atol=2e-15)
    return {"round":1354,"face_success_values":rates,"face_squared_fidelities":fidelities,
            "unconditional_swap_minus":probs,"tagged_joint_repeated":joint,
            "conditional_partner_witness":mixed_witness()}


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for k in actual:
            compare(actual[k],saved[k])
    else:
        np.testing.assert_allclose(actual,saved,rtol=1e-12,atol=2e-14)


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=ap.parse_args()
    data=calculate()
    target=Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with target.open("x",encoding="utf8") as f:
            json.dump(data,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(data,json.loads(target.read_bytes()))
    print(json.dumps(data))

