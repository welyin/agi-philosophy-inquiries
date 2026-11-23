"""Compare tree-level lepton-number-conserving matching with equal light Majorana masses."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def weight(d, m):
    f = np.linalg.solve(m.T, d.T).T
    return f @ f.conj().T


def operator(d, m):
    return -d @ np.linalg.solve(m, d.T)


def compare(d, t, c):
    m = np.block([[np.zeros_like(t), t], [t.T, c]])
    k = -t @ np.linalg.solve(c, t.T)
    endpoints = np.column_stack([d, np.zeros_like(d)])
    light = operator(endpoints, m)
    np.testing.assert_allclose(light, operator(d, k), rtol=2e-12, atol=1e-22)
    w, w0 = weight(endpoints, m), weight(d, k)
    delta = w-w0
    f = np.linalg.solve(t, d.T).T
    np.testing.assert_allclose(delta, f @ f.conj().T, rtol=2e-10, atol=1e-27)
    sig = np.linalg.svd(m, compute_uv=False)
    eig = np.linalg.eigvalsh(delta)
    assert abs(eig[0]) < eig[-1]*1e-12 and eig[1] > 0
    normalized = eig[1]*eig[2]*np.prod(sig)/np.linalg.det(d.conj().T @ d).real
    np.testing.assert_allclose(normalized, 1., atol=1e-12)
    floor = np.linalg.eigvalsh(delta-d@d.conj().T/sig[0]**2)
    assert floor[0] > -eig[-1]*1e-12
    # An unconstrained two-singlet competitor matches both Wilson matrices.
    ge, gu = np.linalg.eigh(w)
    root = gu[:, -2:]*np.sqrt(ge[-2:])
    pinv = np.linalg.pinv(root)
    fake_m = -pinv @ light @ pinv.T
    fake_d = root @ fake_m
    np.testing.assert_allclose(operator(fake_d, fake_m), light, rtol=3e-12, atol=1e-22)
    np.testing.assert_allclose(weight(fake_d, fake_m), w, rtol=3e-12, atol=2e-27)
    return dict(heavy_masses_TeV=(sig/1000).round(7).tolist(),
                trace_current_weight=float(f'{np.trace(w).real:.9e}'),
                excess_eigenvalues=[float(f'{v:.9e}') for v in eig[1:]],
                trace_ratio=float(f'{np.trace(w).real/np.trace(w0).real:.9f}'),
                pseudo_determinant_ratio=round(float(normalized), 11),
                unconstrained_two_state_masses_TeV=(np.linalg.svd(fake_m, compute_uv=False)/1000).round(7).tolist())


def compute():
    cfg = read('../1685/inputs.json')
    old = read('../1688/results.json')['rank_two_witness']
    witness = read('../1678/independent_results.json')['joint_witness']
    eps, x = cfg['scalar_magnitude_over_cutoff'], cfg['real_X_over_cutoff']
    s = eps*np.exp(1j*np.pi*cfg['scalar_angle_over_pi'])
    kappa = read('../1688/inputs.json')['majorana_cross_coefficient']
    t = np.array([[1., 0.], [kappa*s*s, x]])
    c = -np.eye(2)
    r = old['small_root_relative_amplitude']
    a, b = np.array([0., 1., 1.]), np.array([1., 3., 1.])
    legs = (s/eps)**np.array([2, 0, 1])
    db = math.sqrt(r)*b
    da = a*x/eps**2
    d = legs[:,None] * np.column_stack([db, s*s*(da+kappa*db)])
    base_mass = cfg['cutoff_TeV']*1000
    reduced = operator(d, -t@np.linalg.solve(c,t.T))
    svals = np.linalg.svd(reduced, compute_uv=False)
    np.testing.assert_allclose(svals[0]/svals[1], witness['R'], rtol=1e-10)
    source = math.sqrt(math.sqrt(witness['d21'])*1e-9 * base_mass/svals[1])
    d = source*d
    first = compare(d,base_mass*t,base_mass*c)
    second = compare(d,2*base_mass*t,4*base_mass*c)
    np.testing.assert_allclose(np.array(first['excess_eigenvalues'])/4,second['excess_eigenvalues'],rtol=1e-8,atol=0)
    return dict(benchmark_dirac_scale_GeV=round(source, 11), initial=first,
                scaled_same_low_mass=second, mass_gap_solar_eV2=witness['d21'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    result = compute()
    path = HERE/'independent_results.json'
    if ap.parse_args().check:
        assert result == json.loads(path.read_bytes())
        print('PASS: current matching, heavy-spectrum product, physical two-singlet competitor')
    else:
        with path.open('x', encoding='utf8') as fp:
            fp.write(json.dumps(result, indent=2)+'\n')
