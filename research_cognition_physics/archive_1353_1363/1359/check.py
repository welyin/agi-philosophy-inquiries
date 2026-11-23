"""Construct the two physical vertices on all 236 one-path Fock states."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def parity(order):
    inversions = sum(order[i] > order[j] for i in range(3) for j in range(i+1, 3))
    return (-1)**inversions


def vertices():
    triples = list(itertools.product(range(3), range(2), range(3), range(2),
                                     range(3), range(2)))
    middle = list(itertools.product(range(3), range(3), range(2)))
    first = np.zeros((18, 216))
    second = np.zeros((2, 18))
    for r, (c, k, spin) in enumerate(middle):
        if c == k:
            second[spin, r] = 1/math.sqrt(3)
        for col, (i, alpha, j, beta, l, sigma) in enumerate(triples):
            if (l, sigma) != (k, spin) or len({i, j, c}) < 3 or alpha == beta:
                continue
            first[r, col] = parity([i,j,c]) * (1 if alpha == 0 else -1) / 2
    assert np.allclose(first @ first.T, np.eye(18))
    assert np.allclose(second @ second.T, np.eye(2))
    return first, second, triples


def make_parent(first, second, g=1., h=1.):
    parent = np.zeros((236,236), dtype=complex)
    parent[216:234,:216] = g*first
    parent[234:,216:234] = h*second
    return parent + parent.conj().T


def check_roles(first, second):
    rng = np.random.default_rng(1359)
    su, _ = np.linalg.qr(rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)))
    su[:, 0] /= np.linalg.det(su)
    so, _ = np.linalg.qr(rng.normal(size=(3,3)))
    so[:, 0] *= np.linalg.det(so)
    worst = 0.
    for m in [su, so, np.exp(2j*math.pi/3)*so]:
        one = np.kron(m, np.eye(2))
        triple = np.kron(np.kron(one, one), one)
        middle = np.kron(m.conj(), one)
        worst = max(worst, float(np.linalg.norm(first@triple-middle@first)),
                    float(np.linalg.norm(second@middle-second)))
    assert worst < 1.e-12
    return worst


def check_preparation():
    basis = list(itertools.product(range(3), repeat=3))
    positions = {v:i for i,v in enumerate(basis)}
    antisym = np.zeros((27,27))
    for p in itertools.permutations(range(3)):
        for i, v in enumerate(basis):
            j = positions[tuple(v[k] for k in p)]
            antisym[j,i] += parity(p)/6
    assert np.allclose(antisym@antisym, antisym)
    prepared = antisym @ (np.eye(27)/math.sqrt(27))
    lift = float(np.linalg.norm(prepared)**2)
    assert np.isclose(lift, 1/27)
    return lift


def calculate():
    first, second, triples = vertices()
    parent = make_parent(first, second)
    values, basis = np.linalg.eigh(parent)
    logical = (second@first).T
    assert np.allclose(logical.T@logical, np.eye(2))
    dark = np.zeros(216)
    dark[triples.index((0,0,0,1,0,0))] = 1
    assert np.linalg.norm(first@dark) == 0
    settings = []
    for angle, phase in [(math.pi, math.pi/3), (math.pi/2, 2*math.pi/3)]:
        t = angle/math.sqrt(2)
        u = (basis*np.exp(-1j*t*values))@basis.conj().T
        evolved = u[:,:216] @ logical
        amp = (math.cos(angle)-1)/2
        assert np.linalg.norm(u[234:,:216]-amp*logical.T) < 1.e-12
        assert np.linalg.norm(u[:216,234:]-amp*logical) < 1.e-12
        # Matrix columns encode the two spin amplitudes entangled with a passive partner.
        wave = evolved/math.sqrt(2)
        occ = [float(np.linalg.norm(wave[s])**2) for s in [slice(216),slice(216,234),slice(234,None)]]
        assert abs(sum(occ)-1) < 1.e-12
        rec = wave[234:]
        left, right = rec/math.sqrt(2), np.exp(1j*phase)*rec/math.sqrt(2)
        minus = float(np.linalg.norm((left-right)/math.sqrt(2))**2)
        plus = float(np.linalg.norm((left+right)/math.sqrt(2))**2)
        assert abs(minus - amp**2*(1-math.cos(phase))/2) < 1.e-12
        assert abs(plus+minus-occ[-1]) < 1.e-12
        recovered = u[:216,234:] @ rec
        assert np.linalg.norm(recovered - amp**2*logical/math.sqrt(2)) < 1.e-12
        leak = math.sqrt(.88)*logical[:,0] + math.sqrt(.12)*dark
        leak_p = float(np.linalg.norm(u[234:,:216] @ leak)**2)
        assert np.isclose(leak_p, .88*occ[-1])
        settings.append({"omega_t":angle,"receiver_phase":phase,
                         "source_middle_receiver":occ,"receiver_plus_minus":[plus,minus],
                         "unconditioned_from_three_independent_pair_seeds_minus":minus/27,
                         "orthogonal_path_log_minus":occ[-1]/2,
                         "12_percent_orthogonal_leak_receiver":leak_p,
                         "roundtrip_double_success":float(np.linalg.norm(recovered)**2)})
    eigen = {label:int(np.count_nonzero(np.isclose(values,v))) for label,v in
             [("-sqrt2",-math.sqrt(2)),("-1",-1.),("0",0.),("1",1.),("sqrt2",math.sqrt(2))]}
    assert list(eigen.values()) == [2,16,200,16,2]
    return {"dimensions":[216,18,2],"single_path_spectrum_multiplicities":eigen,
            "settings":settings,"prepare_colour_probability":check_preparation(),
            "colour_intertwining_max_residual":check_roles(first,second),
            "receiver_readout_total_error_budget":.02}


def compare(x, y):
    if isinstance(x, dict):
        assert x.keys() == y.keys()
        for k in x:
            compare(x[k], y[k])
    elif isinstance(x, list):
        assert len(x) == len(y)
        for a,b in zip(x,y):
            compare(a,b)
    else:
        assert np.isclose(x,y,atol=1.e-12,rtol=1.e-10)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args=parser.parse_args()
    record=calculate()
    target=Path(__file__).with_name("results.json")
    if args.write:
        with target.open("x",encoding="utf8") as out:
            json.dump(record,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(record,json.loads(target.read_text(encoding="utf8")))
    print(json.dumps(record,indent=2))
