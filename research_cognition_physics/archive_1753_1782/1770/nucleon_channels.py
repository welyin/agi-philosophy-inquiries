"""Recover the two-flavour rate Gram matrix seen by three nucleon-decay channels."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def read(rel):
    return json.loads((STAGE/rel).read_bytes())


def observables(gram, phase, a, b, c):
    d, s = gram[0,0].real, gram[1,1].real
    cross = (np.exp(1j*phase)*gram[0,1]).real
    p = np.array([b,a*np.exp(1j*phase)])
    n = np.array([-b,c*np.exp(1j*phase)])
    readings = np.array([d,(p.conj()@gram[:2,:2]@p).real,(n.conj()@gram[:2,:2]@n).real])
    measurement = np.array([[1,0,0],[b*b,a*a,2*a*b],[b*b,c*c,-2*b*c]])
    recovered = np.linalg.solve(measurement,readings)
    np.testing.assert_allclose(recovered,[d,s,cross],atol=1e-14)
    np.testing.assert_allclose(np.linalg.det(measurement),-2*a*b*c*(a+c))
    return readings


def fixed_phase_rescaling(d, f, heavy):
    """Vary the physical heavy frame; preserve the entire complex d,s block."""
    generators = np.array([[[1,0],[0,0]], [[0,0],[0,1]],
                           [[0,1],[1,0]], [[0,-1j],[1j,0]]])
    def pack(matrix):
        return np.array([matrix[0,0].real, matrix[1,1].real,
                         matrix[0,1].real, matrix[0,1].imag])
    def value_and_jacobian(u):
        heavy_in_frame = u @ heavy
        loss = d - f[:,None]*heavy_in_frame
        value = np.eye(2) - loss.conj().T@loss
        columns = []
        for generator in generators:
            dl = -1j*f[:,None]*(generator@heavy_in_frame)
            columns.append(pack(-dl.conj().T@loss-loss.conj().T@dl))
        return value, np.column_stack(columns)
    u = np.eye(2, dtype=complex)
    reference, jac0 = value_and_jacobian(u)
    sigma = np.linalg.svd(jac0, compute_uv=False)[-1]
    # For exp(i t.G), ||G(z)||_2 <= sqrt(2)||z||. The twice differentiated
    # exponential and Gram matrix yield this bound for the packed Hessian.
    h = np.linalg.norm(heavy, 2)
    hessian_bound = 4*math.sqrt(2)*h*(1+h)
    jacobian_floor = math.floor(sigma*10000)/10000
    guaranteed = jacobian_floor**2 / (4*hessian_bound*np.linalg.norm(pack(reference)))
    assert guaranteed > 0
    scale = 1.0001
    for _ in range(10):
        value, jac = value_and_jacobian(u)
        difference = pack(scale*reference-value)
        if np.linalg.norm(difference)<2e-14:
            break
        step = np.linalg.solve(jac, difference)
        eig, axes = np.linalg.eigh(np.einsum('k,kij->ij', step, generators))
        u = ((axes*np.exp(1j*eig))@axes.conj().T)@u
    np.testing.assert_allclose(u.conj().T@u, np.eye(2), rtol=0, atol=2e-14)
    np.testing.assert_allclose(value_and_jacobian(u)[0], scale*reference, rtol=0, atol=2e-14)
    return dict(jacobian_determinant=np.linalg.det(jac0),
                jacobian_minimum_singular_value=sigma,
                guaranteed_two_sided_relative_gram_scale_radius=guaranteed,
                example_gram_scale=scale, relative_squared_coefficient=1/scale,
                example_heavy_unitary_real=u.real.tolist(),
                example_heavy_unitary_imag=u.imag.tolist())


def run():
    old = runpy.run_path(str(STAGE/'1768/check.py'))
    data = read('1766/independent_results.json')['refits']
    fixed = read('1768/results.json')
    had = read('1483/inputs.json')
    a,b0,c = [abs(had['W0_GeV2'][key][0]) for key in ('us_d','ud_s','ds_u')]
    phase_space = runpy.run_path(str(STAGE/'1483/check.py'))['phase']
    m = had['mass_GeV']
    phase_ratio = phase_space(m['n'],m['K0'])/phase_space(m['p'],m['K+'])
    out = {}
    for name,row in data.items():
        q,g,f,w,ckm = old['tagged_frame'](row,name)
        d = np.sqrt(1-f*f)[:,None]*(w.conj().T@q)*g
        j = np.flatnonzero(g<1)
        heavy = np.zeros((2,3))
        heavy[np.arange(2),j] = np.sqrt(1-g[j]**2)
        def gram(u):
            loss = d-f[:,None]*(u@heavy)
            return np.eye(3)-loss.conj().T@loss
        b = b0*ckm[1]/ckm[0]
        h = gram(np.eye(2))
        example = observables(h,.7,a,b,c)
        specfloor = fixed[name]['any_coherent_flavour_floor']
        max_n_over_p = phase_ratio*(c*c+b*b)/(a*a+b*b)/specfloor
        direct_ratio = had['proton']['partial_lifetime_lower_yr']/had['neutron']['partial_lifetime_lower_yr']
        assert max_n_over_p < direct_ratio
        blind = name.index('1')
        record = dict(any_frame_n_over_p_upper=max_n_over_p,
                      neutron_limit_from_p_years=had['proton']['partial_lifetime_lower_yr']/max_n_over_p,
                      current_limit_switch_ratio=direct_ratio)
        if blind in (0,1):
            norm = h[blind,blind].real
            for phases in ((2.7,1.9),(.6,-2.)):
                u = np.array([[0,np.exp(1j*phases[0])],[np.exp(1j*phases[1]),0]])
                other = gram(u)
                np.testing.assert_allclose(other[blind,blind],norm,atol=1e-14)
                observables(other,-.2,a,b,c)
            record.update(fixed_channel='ds'[blind], multiplier=norm)
        else:
            # At diagonal relative U the two independently adjustable phases let both
            # diagonal responses be rescaled. The unmeasured common source phase cancels Re Hds.
            centers = 1-np.sum(abs(d[:,:2])**2,axis=0)-f**2*(1-g[:2]**2)
            shifts = 2*abs(np.diag(d[:,:2]))*f*np.sqrt(1-g[:2]**2)
            eps = .75*min(shifts/centers)
            assert eps>0
            duplicate = []
            for scale in (1.,1.+eps):
                cos = (scale-1)*centers/shifts
                angles = np.angle(np.diag(d[:,:2]))+np.arccos(cos)
                other = gram(np.diag(np.exp(1j*angles)))/scale
                duplicate.append(observables(other,math.pi/2-np.angle(other[0,1]),a,b,c))
            np.testing.assert_allclose(duplicate[0],duplicate[1],rtol=2e-12)
            record.update(relative_squared_coefficient=1/(1+eps),
                          equal_reduced_rates=duplicate[0].tolist())
            record['fixed_source_phase'] = fixed_phase_rescaling(d[:,:2], f, heavy[:,:2])
        out[name]=record
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda v:float(f'{float(v):.11g}'))
    path=HERE/'nucleon_channels_results.json'
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1770 nucleon rates: three-channel inversion, four anchors and two full complex Gram rescalings")
    else:
        with path.open('x',encoding='utf8') as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write('\n')
