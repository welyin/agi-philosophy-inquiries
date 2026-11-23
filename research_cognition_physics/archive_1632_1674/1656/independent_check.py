"""Recompute scalar conversion with the inverse mass matrix in the flavour basis."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def compute():
    cfg = read('inputs.json')
    portal, charged = read(cfg['portal']), read(cfg['texture'])
    v, f, s = [portal[k] for k in ('v_GeV','f_GeV','sin_theta')]
    c = math.sqrt(1-s*s)
    m = read(cfg['texture_inputs'])['mh_GeV']
    mmu = read(cfg['muon'])['mu_mass_GeV']
    q, sine = charged['quark_q_s']
    b = math.hypot(q*math.sqrt(1-sine**2), sine)
    yf = np.array([[0, b],[q/b, (1-q*q)*sine*math.sqrt(1-sine**2)/b]])
    _, right = np.linalg.eigh(yf.T@yf)
    powers = 2*right[0]**2
    mat = cfg['quark_fractions']
    light = np.array([[mat[f'up_{i}'][0], mat[f'down_{i}'][0], mat['strange_common'][0]]
                      for i in ('p','n')])
    masses = np.array(cfg['nucleon_masses_GeV'])
    anomaly = (1-light.sum(axis=1))*2/27
    higgs = masses/v*(light.sum(axis=1)+3*anomaly)
    singlet = masses/f*(light[:,1:]@powers)
    source = np.array([higgs, singlet])
    extra = np.column_stack([masses/f*light[:,0], masses/f*anomaly])
    targets = [cfg['targets'][k] for k in cfg['target_order']]
    norm = np.array([t['overlaps'] for t in targets])
    cap = np.array([t['capture_per_s']*cfg['hbar_GeV_s'] for t in targets])
    limits = np.array([t['limit90'] for t in targets])
    rotation = np.array([[c,s],[-s,c]])
    main = read('results.json')
    results = {}
    for M in cfg['mH_tests_GeV']:
        inverse = np.linalg.inv(rotation.T@np.diag([m*m,M*M])@rotation)
        kernel = inverse[1]@source
        group = norm@kernel
        # Eq. (1) coefficient C = -(G_F/sqrt(2))*g_S, eq. (14) gives 4C^2.
        coeff = 4*mmu**5*(group/s)**2 / cap
        widths = {}
        for key, lepton in charged['branches'].items():
            chirals = lepton['h_e_mu_and_mu_e_abs_right']
            prediction = coeff*np.dot(chirals,chirals)
            widths[key] = prediction.tolist()
            np.testing.assert_allclose(prediction, main['reference'][str(M)]['rates_at_ref'][key],
                                       rtol=2e-8,atol=1e-27)
        caps = np.sqrt(limits/coeff)
        np.testing.assert_allclose(caps, main['reference'][str(M)]['Yh_norm_caps'],rtol=2e-8)
        free = inverse[1,1]*extra
        cancel = np.linalg.solve(free,-kernel)
        np.testing.assert_allclose(kernel+free@cancel,0,atol=5e-23)
        # With one real heavy-quark coefficient, the minimax of the two
        # absolute normalized amplitudes is attained at opposite signs.
        w = np.sqrt(4*mmu**5/(limits*cap))
        a, b_heavy = w*group/s, w*(norm@free[:,1])/s
        best = -(a[0]+a[1])/(b_heavy[0]+b_heavy[1])
        floors = a+b_heavy*best
        np.testing.assert_allclose(floors[0],-floors[1],rtol=1e-10)
        assert b_heavy.min()>0 and abs(floors[0])>0
        results[str(M)] = {'norm_h_coupling_caps':caps.tolist(),
                           'simultaneous_n_u_and_sum_n_heavy_zero':cancel.tolist(),
                           'one_heavy_profile_h_norm_upper':1/abs(floors[0]),
                           'rates':widths}
    assert np.all(powers>0) and np.linalg.det(norm)!=0 and np.linalg.det(extra)!=0
    return {'down_quark_radial_powers':powers.tolist(), 'nuclear_response_det':np.linalg.det(norm),
            'by_mass_GeV':results}


def tidy(o):
    if isinstance(o, dict):
        return {k:tidy(v) for k,v in o.items()}
    if isinstance(o, list):
        return [tidy(i) for i in o]
    return float(f'{o:.10g}') if isinstance(o,(float,np.floating)) else o


if __name__ == '__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    data=tidy(compute())
    if args.check:
        assert data==read('independent_results.json')
        print('PASS: KKO normalization, non-diagonal propagator, rates, two-nucleus profile and null')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(data,indent=2)+'\n')
