"""Construct the two readout bases from labeled permutation operations."""
import argparse
import itertools
import json
from collections import deque
from math import asin, degrees, isclose, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def permutation(p):
    return np.eye(len(p))[:, p]


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    members = {identity}
    pending = deque([identity])
    while pending:
        current = pending.popleft()
        for g in generators:
            product = tuple(g[i] for i in current)
            if product not in members:
                members.add(product)
                pending.append(product)
    return members


def verify_port_graphs():
    transitive = {}
    for size in range(4, 9):
        cycle = (1, 2, 0, *range(3, size))
        valid = 0
        for i, j in itertools.combinations(range(size), 2):
            swap = list(range(size))
            swap[i], swap[j] = j, i
            group = generated_group((cycle, swap))
            if len({g[0] for g in group}) == size:
                assert len(group) == 24
                valid += 1
        transitive[str(size)] = valid
    counter = generated_group(((1, 2, 0, 3, 4), (3, 4, 2, 0, 1)))
    assert len({g[0] for g in counter}) == 5
    return {"transitive_single_swaps_by_size": transitive,
            "two_swaps_five_ports_group_order": len(counter)}


def charged_spectrum(T, contrast, modes, setup, even):
    c, b = setup["offset"], setup["oriented"]
    H = c*contrast + even*contrast@(T+T.T) + 1j*b*contrast@(T-T.T)
    energies = np.real(np.diag(modes.conj().T@H@modes))
    assert np.allclose(H@modes, modes*energies)
    assert min(energies) > 0
    return H, energies


def quark_constraint(k, measurement):
    center, width = measurement
    amin = center-k*width
    high, low = amin**2, 1-amin**2
    small_gap = asin(sqrt(1/6))-asin(sqrt(low))
    large_gap = asin(amin)-asin(sqrt(2/3))
    assert small_gap > 0 and large_gap > 0
    return {"Vud_min": amin, "large_probability_floor": high, "other_probability_ceiling": low,
            "small_entry_gap_degrees": degrees(small_gap),
            "large_entry_gap_degrees": degrees(large_gap),
            "least_total_ray_rotation_degrees": degrees(min(small_gap,large_gap))}


def calculate():
    cfg = read("inputs.json")
    T, P = permutation(cfg["cycle"]), permutation(cfg["swap"])
    eye = np.eye(4)
    common = np.ones(4)/2
    contrast = eye-np.outer(common,common)
    anchor = contrast@eye[:,3]
    e = anchor/np.linalg.norm(anchor)
    z = (eye[:,3]-eye[:,0])/sqrt(2)
    omega = np.exp(2j*pi/3)
    mu = np.array([1,omega,omega**2,0])/sqrt(3)
    modes = np.column_stack([e,mu,mu.conj()])
    assert np.allclose(modes.conj().T@modes, np.eye(3))
    assert np.allclose(T@modes, modes*np.array([1,omega**2,omega]))
    assert np.allclose(P@z,-z)
    eigenweights = abs(modes.conj().T@z)**2
    assert np.allclose(eigenweights,[2/3,1/6,1/6])
    odd, even = (eye-P)/2,(eye+P)/2
    t_zero = (eye+T+T@T)/3
    from_ref = odd@eye[:,3]
    assert np.allclose(from_ref,z/sqrt(2))
    full_prob = np.linalg.norm(t_zero@from_ref)**2
    assert isclose(full_prob,1/3)
    assert np.allclose(t_zero@from_ref,e/sqrt(3))
    paths = [even@T@z,even@T@T@z]
    assert all(isclose(np.vdot(x,x),.75) for x in paths)
    echo = np.vdot(z,T@z)
    assert isclose(echo,.5)
    setup = cfg["charged_square"]
    charged, charged_sq = charged_spectrum(T, contrast, modes, setup, setup["even"])
    _, alt = charged_spectrum(T, contrast, modes, setup, cfg["alternate_even_coefficient"])
    light = int(np.argmin(alt))
    assert light != 0 and np.argmin(charged_sq) == 0
    charged_dirac = (modes*np.sqrt(charged_sq))@modes.conj().T
    assert np.allclose(charged_dirac@charged_dirac.conj().T, charged)

    # Columns are real contrast differences accessible to a P-even filter.
    v = np.array([0,1,-1,0])/sqrt(2)
    w = np.array([1,-1,-1,1])/2
    assert np.allclose((paths[0]+paths[1])/np.linalg.norm(paths[0]+paths[1]),w)
    assert np.allclose((paths[1]-paths[0])/np.linalg.norm(paths[1]-paths[0]),v)
    sine = cfg["neutral"]["sine_of_even_plane_rotation"]
    cosine = sqrt(1-sine**2)
    directions = np.column_stack([z,cosine*w+sine*v,-sine*w+cosine*v])
    mass = np.array(cfg["neutral"]["masses"])
    assert np.allclose(directions@directions.T,contrast)
    neutral_dirac = directions[:,1:]*mass[1:]
    mass_sq = neutral_dirac@neutral_dirac.T
    assert np.allclose(P@neutral_dirac, neutral_dirac)
    assert np.allclose(mass_sq@z,0)
    assert not np.allclose(charged@P,P@charged)
    assert not np.allclose(mass_sq@T,T@mass_sq)
    # Between charged and neutral spaces the weak current is the common identity on contrast.
    cc = modes.conj().T@contrast@directions
    assert np.allclose(cc.conj().T@cc,np.eye(3))
    probabilities = []
    for time in cfg["neutral"]["propagation_times_in_pi"]:
        prop = (directions*np.exp(-1j*pi*time*mass**2))@directions.T
        survival = abs(np.vdot(e,prop@e))**2
        sumrule = sum(abs(modes.conj().T@prop@e)**2)
        assert isclose(sumrule,1)
        probabilities.append(survival)
    assert np.allclose(probabilities,[(31**2+15**2)/48**2,9/64])

    # Finite group average removes contrast; averaging probabilities would not do so.
    gs = generated_group((cfg["cycle"],cfg["swap"]))
    assert len(gs)==24
    assert np.allclose(sum(permutation(g) for g in gs)/len(gs),eye-contrast)
    quark = read(cfg["quark_input"])["ud_beta"]
    return {"port_graph":verify_port_graphs(),"reference_contrast_success":float(np.vdot(anchor,anchor)),
            "reference_to_tfixed_via_odd_success":full_prob,"one_cycle_echo_amplitude":echo,
            "even_success_after_rotating_odd":float(np.vdot(paths[0],paths[0])),
            "null_flavor_weights":eigenweights.tolist(),"charged_square_in_mode_order":charged_sq.tolist(),
            "alternate_electron_zero_weight":eigenweights[light],
            "neutral_e_weights":(abs(cc[0])**2).tolist(),
            "two_electron_survival_readings":probabilities,
            "quark_transfer":{str(k):quark_constraint(k,quark) for k in cfg["rectangular_widths"]}}


def agrees(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(agrees(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return isclose(a,b,rel_tol=1e-11,abs_tol=1e-14)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as output:
            json.dump(result,output,indent=2,allow_nan=False)
            output.write("\n")
    else:
        assert agrees(result,read("results.json"))
        print("1410 port action, two residual masses, operational probes and cross-sector bounds verified")
