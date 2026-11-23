"""Canonical three particles among six modes read out by one flavour-symmetric exchange."""
import argparse
import itertools
import json
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def move(mask, origin, target):
    if not mask & (1 << origin):
        return None
    sign = (-1)**((mask & ((1 << origin)-1)).bit_count())
    mask ^= 1 << origin
    if mask & (1 << target):
        return None
    sign *= (-1)**((mask & ((1 << target)-1)).bit_count())
    return mask | (1 << target), sign


def solve():
    api = runpy.run_path(str(HERE.parent/"1748"/"check.py"))
    bg = api["OLD"]["mother"](.125)
    L, y, _ = np.linalg.svd(bg["yW"][2, :, :3])
    L, y = L[:, ::-1], y[::-1]
    h = y**2 / np.sum(y**2)
    h -= h[0]
    levels = h/h[1]  # Units of first energy gap; absolute Hamiltonian scale remains free.
    masks = [m for m in range(64) if m.bit_count() == 3]
    index = {m:i for i,m in enumerate(masks)}
    n = np.array([[(m >> i & 1) + (m >> (i+3) & 1) for i in range(3)] for m in masks])
    energies = n @ levels
    energies -= energies.min()
    ground = np.flatnonzero(energies == 0)
    np.testing.assert_array_equal(n[ground], [[2,1,0],[2,1,0]])
    signs = np.array([1, -1, 1])
    event = np.all((-1)**n == signs, axis=1)
    assert event.sum() == 4

    def probs(temperature):
        weight = np.exp(-energies/temperature)
        return (weight @ event)/weight.sum()

    # Positive temperature examples at the specified threshold, not a monotonicity claim.
    benchmarks = {}
    for target in (.99, .999):
        lo, hi = .01, 1.
        assert probs(lo) > target > probs(hi)
        for _ in range(60):
            mid = (lo+hi)/2
            if probs(mid) > target:
                lo = mid
            else:
                hi = mid
        benchmarks[str(target)] = (lo+hi)/2

    # V = sum_ij |i><j| ⊗ sum_a c^dagger_(aj) c_(ai).
    dim = len(masks)
    V = np.zeros((3*dim,3*dim))
    for m, p in index.items():
        for j,i,a in itertools.product(range(3),range(3),range(2)):
            hop = move(m,3*a+i,3*a+j)
            if hop:
                target, sign = hop
                V[i*dim+index[target],j*dim+p] += sign
    np.testing.assert_array_equal(V,V.T)
    gind = np.array([i*dim+j for i in range(3) for j in ground])
    qind = np.setdiff1d(np.arange(3*dim),gind)
    np.testing.assert_array_equal(V[np.ix_(gind,gind)],np.diag(np.repeat([2,1,0],2)))
    goal = np.eye(3*dim)[:,gind]*np.repeat(signs,2)
    dynamics = []
    for g in (.1,.01,.001):
        H = np.diag(np.tile(energies,3))+g*V
        eig,U = np.linalg.eigh(H)
        evolution = (U*np.exp(-1j*np.pi*eig/g)) @ U[gind].conj().T
        dynamics.append(dict(exchange_over_first_gap=g,
                             full_isometry_error=np.linalg.norm(evolution-goal,2),
                             maximum_reference_leak_probability=np.linalg.norm(evolution[qind],2)**2))
    assert dynamics[-1]["full_isometry_error"] < .005

    # Equilibrium of independent copies with one chemical potential:
    # odd probabilities a,b,c cannot exceed 1/2 and acceptance is multilinear.
    maximum = max(b*(1-a)*(1-c) + (1-b)*a*c
                  for a,b,c in itertools.product((0.,.5),repeat=3))
    assert maximum == .5
    # A common energy splitting between the two copies leaves the eigenlines unchanged.
    def ground_at_split(split):
        energies_two = np.concatenate([levels-split/2,levels+split/2])
        lowest = np.argsort(energies_two)[:3]
        return np.bincount(lowest%3,minlength=3)
    np.testing.assert_array_equal(ground_at_split(levels[-1]*.999),[2,1,0])
    np.testing.assert_array_equal(ground_at_split(levels[-1]*1.001),[1,1,1])
    # Exact charged-basis parity read back into the actual heavy rows.
    Xbase,Fbase,_,_ = api["terms"](bg,0,(1,1,1))
    Xwant,Fwant,_,_ = api["terms"](bg,0,(1,-1,1))
    P = L @ np.diag(signs) @ L.conj().T
    np.testing.assert_allclose(Xbase@P,Xwant,atol=2e-13)
    match = np.linalg.norm(P@Fbase@P-Fwant,2)
    assert match < 4e-13
    _,actual = api["daughter"](bg,P@Fbase@P)
    obs = actual["neutrinos"]
    return dict(reference_level_spacing_over_coefficient=h.tolist(),
                levels_in_first_gap_units=levels.tolist(),
                canonical_ground_degeneracy=len(ground),
                critical_copy_splitting_in_first_gap_units=levels[-1],
                canonical_infinite_temperature_probability=event.mean(),
                chemical_potential_only_upper_bound_up_to_common_sign=maximum,
                temperature_in_first_gap_units_at_registered_probabilities=benchmarks,
                probability_at_half_a_gap=probs(.5),
                exchange_gates=dynamics,
                physical_matching_operator_residual=match,
                strong_gap_limit_neutrinos=dict(angles=obs["angle_readings"],
                                                 dm31=obs["dm31_eV2"],dm32=obs["dm32_eV2"]))


def simple(obj):
    if isinstance(obj,dict):
        return {k:simple(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple,np.ndarray)):
        return [simple(x) for x in obj]
    if isinstance(obj,(float,np.floating)):
        return float(f"{obj:.8g}") if abs(obj)>1e-13 else 0.
    return obj


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("--check",action="store_true")
    args = parse.parse_args()
    result = simple(solve())
    output = HERE/"independent_results.json"
    if args.check:
        assert result == json.loads(output.read_bytes())
        print("1750 canonical flavour reference, exchange, finite gap, thermal and endpoint checks passed.")
    else:
        with output.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")

