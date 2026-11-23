"""Compare two charge-summed rare B decays for any common neutral-scalar amplitudes."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
MASSES = np.array([0.000510999, 0.1056584, 1.77686])
CHANNELS = {"emu": (0, 1), "mutau": (1, 2)}
MESONS = {"Bs": (5.367, [6.3e-9, 7.2e-9], 4.2e-5),
          "Bd": (5.280, [1.3e-9], 1.4e-5)}


def gram(texture, pair, meson_mass):
    i, j = pair
    x, y = MASSES[[i, j]] / meson_mass
    phase = np.sqrt((1-(x+y)**2)*(1-(x-y)**2))
    p, q = texture[i, j], texture[j, i]
    diag = (1-x*x-y*y)*(abs(p)**2+abs(q)**2)
    off = -4*x*y*np.conj(p*q)
    return phase*np.array([[diag, off], [off.conjugate(), diag]])


def spin_trace(texture, pair, mass, amplitudes):
    total = 0.
    for i, j in (pair, pair[::-1]):
        x, y = MASSES[[i, j]] / mass
        r = amplitudes[0]*texture[i, j]
        l = amplitudes[1]*texture[j, i].conjugate()
        total += ((1-x*x-y*y)*(abs(r)**2+abs(l)**2)
                  -4*x*y*(r*l.conjugate()).real)
    return total*np.sqrt((1-(x+y)**2)*(1-(x-y)**2))


def evaluate(texture, mass, limits, direct_limit):
    light, heavy = [gram(texture, p, mass) for p in CHANNELS.values()]
    # The eigenvalues of the Hermitian pencil bound every complex two-port amplitude.
    values, vectors = np.linalg.eigh(heavy)
    assert min(values) > 0
    inverse_root = (vectors/np.sqrt(values)) @ vectors.conj().T
    bounds, eigenvectors = np.linalg.eigh(inverse_root @ light @ inverse_root)
    for ratio, eigvec in zip(bounds, eigenvectors.T):
        a = inverse_root @ eigvec
        rates = [spin_trace(texture, p, mass, a) for p in CHANNELS.values()]
        np.testing.assert_allclose(rates[0]/rates[1], ratio, rtol=1e-12)
    for a in ([1, 0], [0, 1], [1+2j, -3j]):
        for p in CHANNELS.values():
            rate = spin_trace(texture, p, mass, a)
            np.testing.assert_allclose(np.vdot(a, gram(texture, p, mass)@a), rate,
                                       rtol=1e-12, atol=1e-18)
    caps = [limit / bounds[0] for limit in limits]
    return {"charge_summed_emu_over_mutau": bounds.tolist(),
            "emu_limits_in_primary_order": limits,
            "conditional_mutau_caps": caps,
            "published_direct_single_species_mutau_limit": direct_limit,
            "direct_limit_over_weakest_conditional_cap": direct_limit/max(caps)}


def calculate():
    source = runpy.run_path(str(HERE.parent/'1793'/'independent_check.py'))
    result = {}
    for step in (.1, .05):
        textures = source['matrices'](step)[0]
        result[str(step)] = {
            matching: {species: evaluate(texture, *spec) for species, spec in MESONS.items()}
            for matching, (texture, _) in textures.items()}
    for matching in result['0.1']:
        for species in MESONS:
            field = 'charge_summed_emu_over_mutau'
            np.testing.assert_allclose(result['0.1'][matching][species][field],
                                       result['0.05'][matching][species][field], rtol=5e-6)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    raw = json.dumps(calculate(), indent=2, ensure_ascii=False) + '\n'
    path = HERE/'independent_results.json'
    if args.check:
        assert json.loads(path.read_bytes()) == json.loads(raw)
    else:
        path.write_bytes(raw.encode('utf-8'))
    print('1795 independent', hashlib.sha256(raw.encode()).hexdigest()[:12])
