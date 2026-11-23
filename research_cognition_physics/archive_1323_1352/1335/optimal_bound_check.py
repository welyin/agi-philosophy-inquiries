"""Round 1335 supplemental audit, independent from the author and first checker."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

I = np.eye(3)
UP = np.array([1., 0.])
DOWN = np.array([0., 1.])
OPT = (np.sqrt(2.) - 1.) / np.sqrt(3.)


def cross(v):
    x, y, z = v
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=complex)


# Real Cartesian flavor basis. i cross(e_k) obeys the spin-1 commutation rule.
J = [1j * cross(e) for e in I]
PAULI = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.diag([1, -1])]
X = sum(np.kron(p, j) for p, j in zip(PAULI, J))
E, V = np.linalg.eigh(X)


def column(s):
    return np.kron(s[:, None], I)


def player(a, b):
    return np.array([np.cos(a), np.exp(1j*b) * np.sin(a)])


def formation(theta):
    u = (V * np.exp(-1j * theta * E)) @ V.conj().T
    return u @ column(UP), u @ column(DOWN)


def lower_mass(s, cosine_complement, t, heavy):
    total = t*t + heavy*heavy
    disc = np.hypot(t*t-heavy*heavy, 2*t*heavy*cosine_complement) / total
    # Stable at s=0 and small t/M.
    return np.sqrt(2.) * t * heavy * s / (total * (1. + disc))**.5


def parent(w, complement, x, t, heavy):
    access = column(x)
    k = access.conj().T @ w
    sigma = np.linalg.svd(k, compute_uv=False)
    full = np.column_stack([t * access, heavy * complement])
    all_masses = np.linalg.svd(full, compute_uv=False)[::-1]
    comp = np.linalg.svd(access.conj().T @ complement, compute_uv=False)[::-1]
    lower = lower_mass(sigma, comp, t, heavy)[::-1]
    upper = np.sqrt(np.maximum(0., t*t+heavy*heavy-lower*lower))
    error = np.max(np.abs(np.sort(np.r_[lower, upper]) - all_masses))
    scalar = np.trace(k)/3
    skew_error = np.max(np.abs((k + k.T) / 2 - scalar*I))
    return k, sigma, all_masses, error, skew_error


def run():
    rng = np.random.default_rng(133502)
    assert np.allclose(E, [-2,-2,1,1,1,1])
    rows = 0
    max_error = max_skew = max_isometry = 0.
    for theta in (0., -0.51, .63, np.pi/3, 1.5, 2.4):
        w, c = formation(theta)
        full_basis=np.column_stack([w, c])
        max_isometry=max(max_isometry,np.max(np.abs(full_basis.conj().T@full_basis-np.eye(6))))
        for a,b in ((0.,0.),(np.pi/2,1.1),(.47,.23),(.88,1.31),(.83,-1.98)):
            for t,heavy in ((.08,2.),(1.,1.),(1.,1.2),(4.,.33)):
                _, s, mass, error, skew_error = parent(w,c,player(a,b),t,heavy)
                assert s[1] >= s[0]/3-1e-13
                assert mass[1] >= OPT*mass[2]-2e-13
                if error > 3e-13*max(t, heavy):
                    raise ValueError((theta, a, t, heavy, error))
                max_error, max_skew = max(max_error,error), max(max_skew,skew_error)
                rows += 1
    assert max_skew < 1e-13 and max_isometry < 1e-13

    # Test the general complex-transpose lemma, including nonnormal matrices.
    for _ in range(90):
        z = rng.normal(size=3)+1j*rng.normal(size=3)
        a = rng.normal()+1j*rng.normal()
        s=np.linalg.svd(cross(z)+a*I,compute_uv=False)
        assert s[1] >= s[0]/3-1e-12

    w,c=formation(np.pi/3)
    sharp_cases=[]
    for heavy in (1.,1.1,1.01,1.001,1.0001):
        _,s,mass,_,_=parent(w,c,UP,1.,heavy)
        sharp_cases.append(dict(heavy=heavy,signal_ratio=float(s[1]/s[0]),
                                full_ratio=float(mass[1]/mass[2]),middle_gap=float(mass[3]-mass[2])))
    assert abs(sharp_cases[0]['signal_ratio']-1/3)<1e-13
    assert abs(sharp_cases[0]['full_ratio']-OPT)<2e-13
    assert sharp_cases[0]['middle_gap']<1e-12
    assert all(z['middle_gap']>0 for z in sharp_cases[1:])

    w,c=formation(.63)
    ku=column(player(.47,.23)).conj().T@w
    hu=ku.conj().T@ku
    probes=[]
    for b in (1.31,1.31+np.pi):
        kd=column(player(.88,b)).conj().T@w
        hd=kd.conj().T@kd
        comm=hu@hd-hd@hu
        probes.append(dict(down_port_phase=float(b),
                           mass_shape=np.linalg.eigvalsh(hd).tolist(),
                           joint_trace=float(np.trace(hu@hd).real),
                           cp_imag_trace=float(np.trace(np.linalg.matrix_power(comm,3)).imag)))
    assert np.allclose(probes[0]['mass_shape'],probes[1]['mass_shape'])
    assert abs(probes[0]['joint_trace']-probes[1]['joint_trace'])>1e-4
    assert abs(probes[0]['cp_imag_trace'])>1e-8
    return dict(parent_cases=rows,complex_skew_cases=90,
                max_paired_spectrum_error=float(max_error),max_cartesian_symmetry_error=float(max_skew),
                max_complete_basis_error=float(max_isometry),best_uniform_constant=float(OPT),
                sharp_cases=sharp_cases,phase_probe=probes)


def compare(real, expected):
    if isinstance(real,dict):
        assert real.keys()==expected.keys()
        for key in real:
            compare(real[key],expected[key])
    elif isinstance(real,list):
        assert len(real)==len(expected)
        for a,b in zip(real,expected):
            compare(a,b)
    else:
        assert np.isclose(real,expected,rtol=1e-7,atol=1e-8),(real,expected)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    out=Path(__file__).with_name('optimal_bound_results.json')
    if args.check:
        compare(result,json.loads(out.read_text(encoding='utf8')))
    else:
        with out.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
    print(json.dumps(result,allow_nan=False))
