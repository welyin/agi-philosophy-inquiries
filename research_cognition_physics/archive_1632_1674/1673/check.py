"""Compare static source matching, spacelike propagation and finite-endpoint poles."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def mass(z, M, g):
    return np.array([[M, g*z],[-g*z.conjugate(), M]],dtype=complex)


def transfer(z, M, g, Q):
    block = mass(z, M, g)
    zero = np.zeros((2,2),dtype=complex)
    hamiltonian = np.block([[zero, block],[block.conj().T,zero]])
    return -np.linalg.inv(hamiltonian-1j*Q*np.eye(4))[2,1]


def poles(z, M, g, t):
    ext = np.zeros((3,3),dtype=complex)
    ext[0,1] = ext[2,0] = t
    ext[1:,1:] = mass(z,M,g)
    return np.sort(np.linalg.svd(ext,compute_uv=False))


def stable(obj):
    if isinstance(obj,dict):
        return {k:stable(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [stable(v) for v in obj]
    if isinstance(obj,(float,np.floating)):
        return float(f'{obj:.12g}')
    return obj


def compute():
    cfg = read('inputs.json')
    g = cfg['g']
    phi = -np.pi/6
    z = (3+2*np.exp(4j*phi))/5
    q = abs(z)
    M = g*np.sqrt(q)
    H1, Hq = [np.linalg.svd(mass(v,M,g),compute_uv=False) for v in (1.,z)]
    np.testing.assert_allclose(Hq/H1, np.sqrt(q))
    point = read(cfg['frozen_point'])['branches'][cfg['family']]['witness']
    old = runpy.run_path(str(HERE.parent/'1668/check.py'))
    E,_ = old['OLD']['charged_basis'](old['read'](old['INPUT']['charged']))
    matrix = old['OLD']['rank_two'](point['d21'],point['d32'],2*np.pi/3,'small')[0]
    frame = np.diag(np.exp(1j*np.array([0.,2*phi,3*phi])))
    f0 = transfer(1.,M,g,0)
    bridge = frame @ np.diag([f0, f0, transfer(z,M,g,0)])/f0
    physical = bridge @ matrix @ bridge.T
    masses,modes = old['OLD']['takagi'](physical)
    U = E.conj().T@modes
    angles = old['OLD']['angles'](U)
    delta = old['phase_delta'](U)
    np.testing.assert_allclose(np.diff(masses**2),[point['d21'],point['d32']],rtol=1e-11)
    np.testing.assert_allclose(angles,point['three_angles'],atol=1e-10)
    assert abs(delta-point['delta_over_pi'])<2e-10
    momenta = []
    for u in cfg['euclidean_Q2_over_Hq2']:
        Q = np.sqrt(u)*Hq[0]
        ratio = transfer(z,M,g,Q)/transfer(1.,M,g,Q)
        np.testing.assert_allclose(ratio, (z/q)*(1+q*u)/(1+u))
        momenta.append([u,abs(ratio)])
    endpoint = []
    for t_over_g in cfg['endpoint_t_over_g']:
        t = t_over_g*g
        target = t**2/(g*(1+q))
        M2 = (g-target)*(q*g-target)
        assert 0 < target < q*g and M2 > 0
        small = [poles(v,M,g,t)[0] for v in (1.,z)]
        exact = [poles(v,np.sqrt(M2),g,t) for v in (1.,z)]
        np.testing.assert_allclose([d[0] for d in exact],target,rtol=2e-13)
        assert min(d[1] for d in exact)>target and small[1]<small[0]
        endpoint.append({'t_over_g':t_over_g,'fixed_static_light_ratio':small[1]/small[0],
                         'exact_equal_pole_M2_over_g2':M2/g**2,'equal_light_pole_over_g':target/g,
                         'nearest_heavy_over_g':min(d[1] for d in exact)/g})
    return stable({'q':q,'Hq_over_H1':Hq[0]/H1[0],'spacelike_response':momenta,
                   'restored_leading_mass':{'three_angles':angles.tolist(),'delta_over_pi':delta,
                      'd21':masses[1]**2-masses[0]**2,'d32':masses[2]**2-masses[1]**2},
                   'finite_endpoint_comparison':endpoint})


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    out = compute()
    if args.check:
        assert out==read('results.json')
        print('PASS: heavy Dirac resolvent, mass texture, spacelike ratio and finite physical poles')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
