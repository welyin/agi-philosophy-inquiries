"""One-loop rho from full weak currents and rotation-independent outer bounds."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def rotation(angles):
    a, b, d, phase = angles
    A, B, D = np.sqrt(1-np.square([a,b,d]))
    z = np.exp(1j*phase)
    return np.array([[A*D,a*D,d/z],[-a*B-A*b*d*z,A*B-a*b*d*z,b*D],
                     [a*b-A*B*d*z,-A*b-a*B*d*z,B*D]])


def kernels(m, n):
    x, y = m*m, n*n
    if x == y:
        return 0., 0.
    if min(x,y) == 0:
        return x+y, 0.
    t = (x-y)/(x+y)
    if abs(t)<.005:
        q = t*t
        plus = 2*q*(1/3+q*(1/15+q*(1/35+q/63)))
        minus = q*(1/3+q*(1/5+q*(1/7+q/9)))
        return (x+y)*plus, 4*m*n*minus
    l = math.log(x/y)
    return x+y-2*x*y*l/(x-y), 2*m*n*((x+y)*l/(x-y)-2)


def current(m,n,left,right):
    p,q = kernels(m,n)
    return (abs(left)**2+abs(right)**2)*p + 2*(left*right.conjugate()).real*q


def x_box(cfg, widths):
    def value(sign):
        w = cfg['Mw_GeV'][0]+sign*widths*cfg['Mw_GeV'][1]
        w0 = cfg['Mw_ref'][0]-sign*widths*cfg['Mw_ref'][1]
        d = cfg['s_eff'][0]-cfg['s_eff_ref'][0]
        return (w/w0)**2-1 + 2*(d+sign*widths*(cfg['s_eff'][1]+cfg['s_eff_ref'][1]))
    return [value(-1),value(1)]


def analyze(cfg, packed, y_mass, original):
    mat = np.asarray(packed)
    mass_matrix = mat[:,:,0]+1j*mat[:,:,1]
    ul, masses, vh = np.linalg.svd(mass_matrix)
    ix = np.argsort(masses)
    masses, ul, vr = masses[ix], ul[:,ix], vh.conj().T[:,ix]
    weak = ul[:3]
    shape = original.conj().T @ weak
    ref = rotation(cfg['reference_rotation'])
    norm = 3/(16*math.pi**2*cfg['vev_GeV']**2)
    z_left = ul.conj().T @ np.diag([1,1,1,0,-1]) @ ul
    z_right = -np.outer(vr[4].conj(), vr[4])
    charged_exotic = sum(current(y_mass,m,ul[4,j],vr[4,j]) for j,m in enumerate(masses))
    neutral = sum(current(masses[i],masses[j],z_left[i,j],z_right[i,j])
                  for i in range(5) for j in range(i+1,5))
    theta = np.array([[kernels(up,m)[0] for m in masses] for up in cfg['up_GeV']])
    charged = np.sum(abs(weak)**2*theta)
    fixed_sm = np.sum(abs(ref)**2*theta[:,:3])
    moving_sm = np.sum(abs(original)**2*theta[:,:3])
    fixed = norm*(charged+charged_exotic-neutral-fixed_sm)
    same = norm*(charged+charged_exotic-neutral-moving_sm)
    def bound(moving):
        ops = [shape @ np.diag(row) @ shape.conj().T for row in theta]
        common = charged_exotic-neutral
        if moving:
            ops = [h-np.diag(row[:3]) for row,h in zip(theta,ops)]
        else:
            common -= fixed_sm
        endpoints = []
        for i in range(3):
            eig = np.array([np.linalg.eigvalsh(h-ops[i]) for h in ops])
            base = common+np.trace(ops[i]).real
            endpoints.append(norm*(base+np.array([eig[:,0].sum(),eig[:,-1].sum()])))
        return [max(x[0] for x in endpoints),min(x[1] for x in endpoints)]
    fixed_range, same_range = bound(False),bound(True)
    assert fixed_range[0]-2e-11<=fixed<=fixed_range[1]+2e-11
    assert same_range[0]-2e-11<=same<=same_range[1]+2e-11
    assert max(abs(shape@shape.conj().T-np.eye(3)).flat)<1e-9
    c2 = 1-cfg['s2']
    remainder = [[x_box(cfg,w)[0]-c2*fixed_range[1], x_box(cfg,w)[1]-c2*fixed_range[0]] for w in [1,2]]
    return {'M_Y_GeV':float(y_mass), 'sum_terms_GeV2': [charged,charged_exotic,neutral,fixed_sm],
            'rho_fixed_SM':fixed,'rho_same_U_SM':same, 'T_same_U_SM':same*cfg['alpha_inverse'],
            'any_rotation_fixed_SM':fixed_range,'any_rotation_same_U_SM':same_range,
            'remainder_any_rotation_one_two_widths':remainder}


def calculate():
    cfg = read('inputs.json')
    old = read(cfg['old'])
    new = read(cfg['released'])
    points = {}
    for name in ['common_witness','high_mass_witness']:
        row = old[name]
        points['1635_'+name] = analyze(cfg,row['weak_basis_mass_GeV'],row['M_Y_GeV'],rotation(cfg['reference_rotation']))
    for row in new['points']:
        points['1640_'+str(row['M_GeV'])] = analyze(cfg,row['weak_basis_mass_GeV'],row['M_GeV'],rotation(cfg['released_rotation']))
    diagonal = np.diag([.005,.1,4.8,2100.,1600.]).astype(complex)
    diagonal[:3,:3] = rotation(cfg['reference_rotation']) @ diagonal[:3,:3]
    blank = analyze(cfg,np.stack((diagonal.real,diagonal.imag),-1),1600.,rotation(cfg['reference_rotation']))
    assert abs(blank['rho_fixed_SM']) < 2e-12
    assert abs(points['1635_common_witness']['rho_fixed_SM']-.00068124)<2e-8
    assert points['1635_high_mass_witness']['remainder_any_rotation_one_two_widths'][1][1]<-.00017
    assert all(kernels(m,m)==(0.,0.) for m in [0,3,1e3])
    return {'X_central':x_box(cfg,0)[0], 'X_one_two_widths':[x_box(cfg,1),x_box(cfg,2)],
            'decoupled_rho':blank['rho_fixed_SM'], 'points':points}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):
        return [clean(v) for v in x]
    return float(format(x,'.10g'))


if __name__ == '__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=clean(calculate())
    target=HERE/'results.json'
    if args.check:
        assert json.loads(target.read_bytes())==result
        print('PASS: complete one-loop currents, five spectra, arbitrary U(3) envelopes and EW windows')
    else:
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
