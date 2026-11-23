"""Weak harmonics, tagged scalar vertices and the QCD phase in a CP-real matching basis."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def run():
    down = read('../1681/results.json')['quark_path_amplitudes']
    powers = np.array(read('../1653/results.json')['powers_A_B_C_E_F_R_D'])
    phi = -np.pi/6
    angles = np.arange(64)*2*np.pi/64
    def weak(ang):
        a,b,c,e,f,r,d = down*np.exp(1j*powers*ang)
        y = np.array([[0,a,0],[b,c,r+d],[0,e,f]])
        h = y@y.conj().T
        return (h[0,1]*h[1,2]*h[2,0]).imag
    vals = np.array([weak(q) for q in angles])
    expected = down[0]**2*down[2]*down[3]*down[4]*down[6]
    np.testing.assert_allclose(vals, expected*np.sin(4*angles),rtol=1e-12,atol=1e-18)
    # At a protected sign of the full tree-level determinant the extra strong phase
    # equals Arg of the mass relative to its S^6 seed, in this fixed UV convention.
    eps = .5
    ref = eps*np.exp(1j*phi)
    points = []
    for ratio in [-.04,-.005,0.,.012,.1]:
        def mass(y):
            return y**6+ratio*y**2
        m = mass(ref)
        h = (6*ref**6+2*ratio*ref**2)/m
        finite = (mass(ref*math.exp(1e-5))-mass(ref*math.exp(-1e-5)))/(2e-5*m)
        np.testing.assert_allclose(finite,h,rtol=2.e-9)
        theta = np.angle(m/ref**6)
        assert abs(np.angle(np.exp(1j*(theta+np.angle(h-2)))))<2.e-14
        points.append([ratio,theta,h.real,h.imag])
    # Add a neutral real polynomial which has any chosen phase, but zero derivative.
    z0 = np.exp(4j*phi)
    delta = .17
    b = -4*math.sin(delta)/math.sqrt(3)
    a = math.cos(delta)+5*b/12
    def alt(y):
        z = (y/eps)**4
        return y**6*(a+b*(z+z*z/2+z**3/3))
    phase = alt(ref)/ref**6
    z_deriv = 1+z0+z0**2
    np.testing.assert_allclose(phase,np.exp(1j*delta),atol=2.e-15)
    assert abs(z_deriv)<1.e-14
    h = (alt(ref*math.exp(1e-5))-alt(ref*math.exp(-1e-5)))/(2e-5*alt(ref))
    np.testing.assert_allclose(h,6.,atol=2.e-8,rtol=0)
    return {'weak_harmonic_normalized_quark_witness':float(expected),
            'two_path_points_ratio_theta_realH_imagH':np.round(points,11).tolist(),
            'neutral_polynomial_competitor':{'phase_rad':delta,'a':a,'b':b,'physical_radial_realH':6.0,
                                             'physical_radial_imagH':0.0}}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    out = run()
    target = HERE/'independent_results.json'
    if p.parse_args().check:
        assert out == json.loads(target.read_bytes())
        print('PASS: physical weak harmonic, scalar vertex and hidden QCD phase at fixed masses')
    else:
        with target.open('x',encoding='utf-8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
