"""Independent audit of static transfer, finite endpoint moment bounds and flavour readout."""
import argparse
import json
from pathlib import Path

import numpy as np


def propagator(h, t):
    energies, modes = np.linalg.eigh(h)
    return (modes * np.exp(-1j * energies * t)) @ modes.conj().T


def chain(a, middle, end, left, right):
    return np.array([[a, left.conjugate(), 0], [left, middle, right.conjugate()],
                     [0, right, end]], complex)


def occupations(h, t):
    return (abs(propagator(h, t)[:, 0]) ** 2).tolist()


def finite_bounds():
    # Different moduli, phases, and energies in a candidate that violates both exact endpoint conditions.
    h = chain(0.25, 0.42, -0.12, np.exp(.3j), .71*np.exp(-.7j))
    ev, modes = np.linalg.eigh(h)
    c = ((h[0, 0] + h[2, 2])/2).real
    ediff = abs(h[0, 0] - h[2, 2])
    gdiff = abs(abs(h[0, 1])**2 - abs(h[1, 2])**2)
    energy_deficit = (ediff / np.ptp(ev)) ** 2
    strength_deficit = (gdiff / np.ptp((ev-c)**2)) ** 2
    # Simultaneous measurement of the spectral projectors is stronger than either moment.
    frequency_envelope = float(np.sum(abs(modes[2] * modes[0].conj())) ** 2)
    rates = np.array([occupations(h, t)[2] for t in np.linspace(0, 27, 1001)])
    assert rates.max() <= frequency_envelope + 2.e-12
    assert frequency_envelope <= 1-max(energy_deficit, strength_deficit) + 2.e-12
    return {"population_ceiling_first_moment": float(1-energy_deficit),
            "population_ceiling_centered_second_moment": float(1-strength_deficit),
            "population_ceiling_spectral_weights": frequency_envelope,
            "largest_on_registered_grid": float(rates.max())}


def flavor():
    # Free flavour reference, and an independently fixed basis in which conversion is read.
    angle = np.pi/6
    rotate = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    masses = np.diag([.375, 1.625])
    other = rotate @ masses @ rotate.T
    site = chain(0., 0., 0., 1., 1.)
    bell = np.eye(2)/np.sqrt(2)
    initial = np.zeros((6, 2), dtype=complex)
    initial[:2] = bell
    response = {}
    for phase, key in [(np.pi/2, "halfpi"), (np.pi, "pi")]:
        t = phase/(masses[1, 1]-masses[0, 0])
        receiver_site = occupations(site, t)[2]
        for name, m in [("no_mixing", masses), ("rotated", other)]:
            h = np.kron(site, np.eye(2)) + np.kron(np.eye(3), m)
            joint = propagator(h, t) @ initial
            h_state = joint[4:]
            received = float(np.linalg.norm(h_state)**2)
            assert np.isclose(received, receiver_site)
            free_pair = propagator(m, t) @ bell
            fidelity = float(abs(np.vdot(free_pair, h_state))**2 / received)
            assert np.isclose(fidelity, 1.)
            # Fresh input is one independently identified flavour; the numerator includes receipt.
            physical_one = propagator(h, t)[4:, 0]
            trans = float(abs(physical_one[1])**2 / received)
            expected = (0 if name == "no_mixing" else (3/4)*np.sin(phase/2)**2)
            assert np.isclose(trans, expected)
            response.setdefault(key, {})[name] = {"receiver": received, "free_reference_fidelity": fidelity,
                                                 "changed_flavor_given_receiver": trans}
    return response


def run():
    delta = np.sqrt(8/3)
    omega = np.sqrt(delta**2+8)
    t_star = 4*np.pi / omega
    cases = {}
    for name, detuning in [("resonant", 0.), ("detuned", delta)]:
        h = chain(0., detuning, 0., 1., 1.)
        cases[name] = {"half": occupations(h, t_star/2), "full": occupations(h, t_star)}
    assert np.allclose(cases["detuned"]["full"], [0,0,1])
    assert np.allclose(cases["detuned"]["half"], [.5,0,.5])
    # Author's separate flavour demonstration uses this detuned site and unit flavour gap.
    # Check its two absolute occupations using the closed even-branch evolution.
    author_times = np.array([np.pi / 2, np.pi])
    even_amplitude = np.exp(-1j * delta * author_times/2) * (
        np.cos(omega * author_times/2) + 1j*delta/omega*np.sin(omega * author_times/2))
    author_absolute = abs(even_amplitude-1)**2/4
    assert np.allclose(author_absolute,[.43280398066661696,.6380103306968355], atol=1.e-13)
    # Ratio r=1/3 gives an all-time upper bound, not just failure at a registered time.
    odd = chain(0, 1, 0, 1, 1)
    for t in np.linspace(0, 2*np.pi, 45):
        c = np.cos(t)
        assert np.isclose(occupations(odd,t)[2], (1-c)**2*(5+4*c)/9)
    assert np.isclose(occupations(odd,2*np.pi/3)[2], .75)
    # Irrational delta/omega. Search a single explicit time with a recorded accuracy.
    irrational_delta = np.sqrt(2)
    w = np.sqrt(10)
    candidates = np.arange(1, 20001)
    rates = np.sin(np.pi*candidates*(1-irrational_delta/w)/2)**2
    i = np.flatnonzero(rates > 1-1.e-7)[0]
    k = int(candidates[i])
    witnessed = occupations(chain(0, irrational_delta, 0, 1, 1), 2*np.pi*k/w)[2]
    assert np.isclose(witnessed,rates[i]) and witnessed > 1-1.e-7
    return {"registered_site_occupations": cases, "odd_rational_all_time_max_probability": .75,
            "irrational_explicit_k": k, "irrational_witness_probability": witnessed,
            "unequal_endpoints": finite_bounds(), "flavour": flavor()}


def compare(got, saved):
    if isinstance(got,dict):
        assert got.keys() == saved.keys()
        for label, values in got.items():
            compare(values,saved[label])
    else:
        assert np.allclose(got,saved,rtol=3.e-10,atol=1.e-12)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    value=run()
    path=Path(__file__).with_name("independent_results.json")
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(value,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        compare(value,json.loads(path.read_text(encoding="utf-8")))
    print("PASS",value)
