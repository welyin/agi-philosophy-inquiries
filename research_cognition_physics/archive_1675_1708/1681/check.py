"""Strong-angle sign, discrete color test, and exact sector budget on the common flavor point."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def spread(values):
    return np.max(values) / np.min(values)


def best_budget(sectors):
    # Log width is the maximum of all intradepartmental affine differences.
    # Parallel lines keep only the largest intercept. A finite nonconstant
    # minimum on the positive half-line occurs at its origin or a crossing.
    lines = {}
    for amplitudes, powers in sectors:
        logs = np.log(amplitudes)
        for i,j in itertools.permutations(range(len(logs)),2):
            m = int(powers[i] - powers[j])
            lines[m] = max(lines.get(m,-np.inf),float(logs[i]-logs[j]))
    candidates = [0.0]
    for (m,a),(n,b) in itertools.combinations(lines.items(),2):
        x = (a-b)/(n-m)
        if x>0:
            candidates.append(x)
    value,x = min((max(b+m*x for m,b in lines.items()),x) for x in candidates)
    return float(np.exp(value)),float(np.exp(-x))


def compute():
    cfg = read('inputs.json')
    uk = read(cfg['up_yukawas_at_1e16'])['high_scale']
    up,du = np.array(uk['up']).T
    obs = read(cfg['quark_joint_witness'])['joint_witness']
    p = np.array(read(cfg['down_powers'])['powers_A_B_C_E_F_R_D'])
    phi = -np.pi/6
    omega = np.exp(4j*phi)
    dc,de = np.array(uk['down']+uk['angles_theta12_theta13_theta23_delta'][:3]).T
    old = runpy.run_path(str(HERE.parent/'1678/independent_check.py'))
    short,long,_,_ = old['get_matrices'](dc+de*np.array(obs['quark_pulls']),obs['quark_gamma_high'])
    d = short+long
    d[1,1] /= 3
    z = d[1,2]
    imag_path = z.imag/omega.imag
    real_path = z.real-imag_path*omega.real
    edges = np.array([d[0,1].real,d[1,0].real,d[1,1].real,d[2,1].real,
                      d[2,2].real,real_path,imag_path])
    kd = int(p[0]+p[1]+p[4])
    assert kd==4 and np.all(edges!=0)
    # Rebuild BEFORE weak-basis rephasing; reduced Z alone cannot give the QCD angle.
    a,b,c,e,f,r,s = edges*np.exp(1j*p*phi)
    cp_basis = np.array([[0,a,0],[b,c,r+s],[0,e,f]])
    np.testing.assert_allclose(np.linalg.det(cp_basis),-np.prod(edges[[0,1,4]])*np.exp(1j*kd*phi))
    h = cp_basis@cp_basis.conj().T
    assert abs((h[0,1]*h[1,2]*h[2,0]).imag)>1e-12
    group = []
    for powers in itertools.product(range(4),repeat=3):
        n = kd+sum(powers)
        if n%4==0 and n%6==0:
            group.append(powers)
            y = up*np.exp(1j*np.array(powers)*phi)
            assert abs(np.sin(np.angle(np.prod(y)*np.linalg.det(cp_basis))))<1e-12
    assert set(group)==set(itertools.permutations([2,3,3]))
    leading = min(best_budget([(up,np.array(row))])[0] for row in group)
    lead_joint = min(best_budget([(up,np.array(row)),(abs(edges),p)])[0] for row in group)
    lifted = []
    for i in range(9):
        for j in range(9-i):
            q = np.array([i,j,8-i-j])
            for_up = best_budget([(up,q)])
            joint = best_budget([(up,q),(abs(edges),p)])
            lifted.append((q.tolist(),for_up,joint))
    free_up = min(lifted,key=lambda a:a[1][0])
    free_joint = min(lifted,key=lambda a:a[2][0])
    eps = .845
    all_order_exponents = np.array([45,8,3])
    assert (kd+all_order_exponents.sum())%12==0
    floor = abs(edges[4]/edges[1])
    trial = max(spread(up/eps**all_order_exponents), spread(abs(edges)/eps**p))
    np.testing.assert_allclose(trial,floor,rtol=1e-12)
    common = read('../1653/inputs.json')
    ps = (2*np.array(common['left']))%4
    ns = kd+ps.sum()
    assert (ns,ns%4)==(6,2)
    ex = common['epsilon']
    # q=c2/(c6*eps^4), total determinant proportional to 1+q*exp(-4i phi).
    theta = read(cfg['neutron_edm_convention'])['conditional_edm']['theta_limit_rad']
    tan = math.tan(theta)
    lo,hi = -2*tan/(math.sqrt(3)-tan),2*tan/(math.sqrt(3)+tan)
    for q in (lo,hi):
        np.testing.assert_allclose(abs(np.angle(1+q/omega)),theta,rtol=1e-12)
    return {'allowed_leading_powers':group,'central_up_only_min_spread':leading,
            'up_two_halfwidth_box_lower_bound':(up[1]-2*du[1])/(up[0]+2*du[0]),
            'central_leading_at_epsilon_half_min':min(spread(up/ex**np.array(q)) for q in group),
            'joint_leading_min_spread':lead_joint,'quark_path_amplitudes':edges.tolist(),
            'total_power_12_best_up':free_up,'total_power_12_best_joint':free_joint,
            'unrestricted_global_joint_floor_at_frozen_point':floor,
            'saturating_power_epsilon_and_two_spreads':[all_order_exponents.tolist(),eps,
                    spread(up/eps**all_order_exponents),spread(abs(edges)/eps**p)],
            'SU5_light_exponents':ps.tolist(),'light_index':int(ns),'single_heavy_min_power':6,
            'heavy_c2_over_c6_conditional_range':[ex**4*lo,ex**4*hi]}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [rounded(v) for v in x]
    if isinstance(x,(np.floating,float)):
        return float(format(x,'.12g'))
    return x


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    expected = rounded(compute())
    if parser.parse_args().check:
        assert expected==read('results.json')
        print('PASS: joint weak and strong flavor charge residues, coefficient bounds and spectator quality')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(expected,indent=2)+'\n')
