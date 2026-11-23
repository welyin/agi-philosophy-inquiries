"""Sharp positive Higgs moment bound in tree matching when the cubic LNV current moment vanishes."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def invariant(port, heavy):
    f = np.linalg.solve(heavy.T, port.T).T
    mass = f @ port.T
    trace_w = np.linalg.norm(f)**2
    moment = np.linalg.norm(port @ heavy.conj().T)**2
    masses = np.linalg.svd(mass, compute_uv=False)
    return moment * trace_w**3 / sum(masses)**4


def unitary(rng, n):
    return np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]


def compute():
    # Saturation for an arbitrary pair of nonzero light masses and arbitrarily rotated flavors.
    cap = 27/4
    rng = np.random.default_rng(1695)
    y = unitary(rng, 3)[:,:2] @ np.diag(np.sqrt([.009,.05]))
    t = np.eye(2)
    c = math.sqrt(2.)*np.eye(2)
    heavy = np.block([[0*t, t],[t,c]])
    endpoint = np.column_stack([y, 0*y])
    ratio = invariant(endpoint,heavy)
    np.testing.assert_allclose(ratio,cap)
    np.testing.assert_allclose(invariant(y,t),1)
    for _ in range(32):
        z = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        c = z+z.T
        t = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        m = np.block([[0*t,t],[t.T,c]])
        v = unitary(rng,4)
        port = endpoint @ v
        m = v.T @ m @ v
        np.testing.assert_allclose(port @ m.conj() @ port.T,0,atol=5e-15)
        assert invariant(port,m)>=cap-1e-11
    # Frozen 1689 mass graph and light-flavor directions. The ratio is invariant under common units.
    cfg = json.loads((HERE.parent/'1689/inputs.json').read_bytes())
    prior = json.loads((HERE.parent/'1685/inputs.json').read_bytes())
    s = prior['scalar_magnitude_over_cutoff']*np.exp(1j*np.pi*prior['scalar_angle_over_pi'])
    x = prior['real_X_over_cutoff']
    m0, g, hx = [cfg[k] for k in ['m_N0H0','coupling_to_s2','coupling_to_x']]
    t = np.array([[m0,0],[g*s**2,hx*x]])
    c = np.diag(cfg['bare_H0_H6_masses_in_sterile_units'])
    spec = json.loads((HERE.parent/'1688/results.json').read_bytes())
    r = spec['rank_two_witness']['small_root_relative_amplitude']
    a, b = np.array([0,1,1]), np.array([1,3,1])
    legs = (s/abs(s))**np.array([2,0,1])
    db = np.sqrt(r)*b
    d = legs[:,None]*np.column_stack([db,s**2*(a/abs(s**4/x**2)**.5+g/m0*db)])
    heavy = np.block([[0*t,t],[t.T,c]])
    frozen = invariant(np.column_stack([d,0*d]), heavy)
    # Reconstruct the additional weak-current weight using exactly the same ports.
    p = t.conj() @ t.T
    extra = d @ np.linalg.solve(p, d.conj().T)
    bmat = d @ d.conj().T
    rebuilt = np.trace(bmat @ np.linalg.pinv(extra,hermitian=True) @ bmat).real
    np.testing.assert_allclose(rebuilt,np.trace(d@p@d.conj().T).real)
    lightdata = json.loads((HERE.parent/'1678/independent_results.json').read_bytes())['joint_witness']
    total = 1e-9*math.sqrt(lightdata['d21'])*(1+lightdata['R'])
    uu = json.loads((HERE.parent/'1558/inputs.json').read_bytes())['original_u_example_GeV']
    target = 125.11**2
    wbound = (cap*total**4/(8*math.pi**2*uu**2*target))**(1/3)
    return dict(sharp_zero_cubic_constant=ratio,sharp_unrestricted_constant=1.,
                frozen_1689_moment_ratio=round(float(frozen),8),
                minimum_trace_W_at_calibrated_Higgs_equal_positive_weights=float(f'{wbound:.10e}'),
                weighted_W_floor_if_wmin_half=float(f'{wbound*2**(-1/3):.10e}'),
                fixed_ports_excess_reconstructs_positive_moment=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    output = compute()
    output['sharp_zero_cubic_constant'] = round(float(output['sharp_zero_cubic_constant']), 10)
    dest = HERE/'independent_results.json'
    if parser.parse_args().check:
        assert output == json.loads(dest.read_bytes())
        print('PASS: sharp 27/4 trace inequality, two-heavy control, frozen graph and weak-weight inverse')
    else:
        with dest.open('x',encoding='utf8') as fp:
            fp.write(json.dumps(output,indent=2)+'\n')
