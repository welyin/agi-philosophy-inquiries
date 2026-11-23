"""Check the exact residual groups and the two singlet completions of the same inputs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def center():
    # a mod 5, b mod 10: z_5^a and exp(2*pi*i*b/10)
    spin = [(2, -1), (-1, 3), (0, -5)]
    kernel = [(a, b) for a in range(5) for b in range(10)
              if all((2*a*n+b*q) % 10 == 0 for n, q in spin)]
    assert set(kernel) == {(i, 4*i % 10) for i in range(5)}
    charges = [10, 5]
    components = {str(q): sum((q*b) % 10 == 0 for b in range(10)) // 5
                  for q in charges}
    assert components == {"10": 2, "5": 1}
    assert ((0, 5) not in kernel and
            all((2*n*0 + 5*q) % 10 == 5 for n, q in spin))
    sm_kernel = [(a,b,k) for a in range(3) for b in range(2) for k in range(6)
                 if (2*a-2*k) % 6 == 0 and (3*b+3*k) % 6 == 0]
    assert set(sm_kernel) == {(k % 3, k % 2, k) for k in range(6)}
    return dict(spin_kernel_mod5_mod10=kernel, sm_kernel_mod3_mod2_mod6=sm_kernel,
                residual_components_by_absolute_singlet_charge=components)


def norm(a):
    return float(np.linalg.norm(a, 2))


def rel(a, b):
    return norm(a-b)/norm(b)


def witness(x, data, api, cfg, previous):
    _, _, _, f, _, ue, _, _ = api['matrices'](x, data)
    _, _, _, _, y, tensor_i = api['symmetric_sources'](x, data)
    v = cfg['effective_dirac_vev_GeV']
    z = complex(*x[7:9])
    target = f/np.linalg.norm(f) + z*tensor_i/np.linalg.norm(tensor_i)
    lab = ue@data[4]
    scale = api['neutral'](lab.conj().T@target@lab.conj())[1]['common_light_scale_eV']*1e-9
    left = scale*f/np.linalg.norm(f)
    dirac = v*y/math.sqrt(2)
    right = -v*v*np.linalg.norm(tensor_i)*f/(2*scale*z)
    poles = np.sort(np.linalg.svd(right, compute_uv=False))
    np.testing.assert_allclose(poles, previous['heavy_singular_masses_GeV'], rtol=2e-12)
    tilt = np.linalg.solve(right, dirac.T)
    light = left-dirac@tilt
    kinetic = tilt.conj().T@tilt
    assert rel(light, scale*target) < 3e-12

    cases = {}
    for k in cfg['k_values']:
        heavy = np.kron([[0., math.sqrt(k)], [math.sqrt(k), -k]], right)
        both = np.hstack([dirac, np.zeros_like(dirac)])
        wide_tilt = np.linalg.solve(heavy, both.T)
        new_light = left-both@wide_tilt
        new_kinetic = wide_tilt.conj().T@wide_tilt
        assert rel(new_light, light) < 1e-10
        assert rel(new_kinetic, (1+1/k)*kinetic) < 1e-11
        minus = 2/(math.sqrt(1+4/k)+1)
        plus = k/minus
        spectrum = np.sort(np.linalg.svd(heavy, compute_uv=False))
        np.testing.assert_allclose(spectrum, np.sort(np.r_[minus*poles, plus*poles]), rtol=7e-12)
        cases[str(k)] = dict(heavy_neutral_poles_GeV=spectrum.tolist(),
                             relative_schur_mismatch=rel(new_light, light),
                             all_flavor_kinetic_ratio=norm(new_kinetic)/norm(kinetic),
                             norm_extra_kinetic=norm(new_kinetic-kinetic),
                             spinor_singlet_yukawa_norm=math.sqrt(k)*poles[-1]/cfg['spinor_singlet_vev_GeV'],
                             finite_light_mass_relative_error_bound=math.expm1(2*math.asinh(norm(wide_tilt)/2)))
    return dict(even_heavy_poles_GeV=poles.tolist(), even_dimless_kinetic_norm=norm(kinetic),
                odd=cases)


def compute():
    cfg=json.loads((HERE/'inputs.json').read_bytes())
    old_file=HERE/cfg['source_points']
    api=runpy.run_path(str(old_file.with_name('check.py')))
    old=json.loads(old_file.read_bytes())
    saved=json.loads((HERE/cfg['original_finite_spectrum']).read_bytes())['cases']
    data=api['backend']()
    return dict(group=center(), cases={key:witness(np.array(x),data,api,cfg,saved[key])
                                      for key,x in old['mixed_type_I_II'].items()})


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(compute(),indent=2,allow_nan=False)+'\n'
    out=HERE/'results.json'
    if args.check:
        assert json.loads(raw)==json.loads(out.read_bytes())
    else:
        out.write_bytes(raw.encode())
    print('1803',hashlib.sha256(raw.encode()).hexdigest()[:12])
