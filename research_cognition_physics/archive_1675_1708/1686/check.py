"""Bound the full determinant phase and test a tagged two-path counterexample."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_bytes())


def odd_necklaces(length):
    def canonical(word):
        return min(word[i:]+word[:i] for i in range(length))
    return sorted({canonical(w) for t in itertools.product('AB',repeat=length)
                   if canonical(w:=''.join(t))!=canonical(w[::-1])})


def compute():
    cfg = load('inputs.json')
    high = load(cfg['yukawa_data'])['high_scale']
    a,b = [(row[-1][0]+2*row[-1][1])**2 for row in [high['up'],high['down']]]
    h = max(a,b)
    assert h<1
    # For m>=7, 1/m<=1/7 and sum_{m=7}^\infty h^(m-6)=h/(1-h).
    harmonic = sum(1/j for j in range(1,7))
    per = 3*(a*b)**3*(harmonic+h/(7*(1-h)))
    direct = (a*b)**3
    all_budget = 2*per+direct
    assert all(not odd_necklaces(n) for n in range(1,6))
    assert len(odd_necklaces(6))==2
    joint = load(cfg['quark_witness'])['joint_witness']
    old = runpy.run_path(str(HERE.parent/'1678/independent_check.py'))
    qc,qe = np.array(high['down']+high['angles_theta12_theta13_theta23_delta'][:3]).T
    raw = qc+qe*np.array(joint['quark_pulls'])
    sh,long,v,_ = old['get_matrices'](raw,joint['quark_gamma_high'])
    sh[1,1] /= 3
    y = raw[2]*(sh+long)
    up = np.array(high['up'])[:,0]
    yuk_u = np.diag(up)
    hu,hd = yuk_u@yuk_u.T,y@y.conj().T
    c = hu@hd-hd@hu
    invariant = np.linalg.det(c).imag/2
    def discriminant(spectrum):
        return math.prod(j*j-i*i for i,j in itertools.combinations(spectrum,2))
    J = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    from_masses = J*discriminant(up)*discriminant(raw[:3])
    np.testing.assert_allclose(abs(invariant)/abs(from_masses),1,rtol=1.e-9)
    # Verify the nonzero degree-six trace directly. Use cyclic commutator form
    # above as the robust reference at the actual mass hierarchy.
    trace_word = np.trace(hu@hu@hd@hd@hu@hd)
    np.testing.assert_allclose(abs(trace_word.imag/invariant),1,rtol=2.e-6)
    reduced = load(cfg['real_path_shapes'])['quark_path_amplitudes']
    R,D = [raw[2]*p for p in reduced[-2:]]
    omega = np.exp(-2j*np.pi/3)
    tagged = float((R*D*omega).imag)
    np.testing.assert_allclose(y[1,2],R+D*omega,rtol=3e-12)
    rotated = y@np.diag([1,1,np.exp(2j*tagged)])
    np.testing.assert_allclose(rotated@rotated.conj().T,hd,atol=2e-20)
    shift = np.angle(np.linalg.det(rotated)/np.linalg.det(y))
    np.testing.assert_allclose(shift,2*tagged,atol=1e-16)
    theta = load(cfg['edm_calibration'])['conditional_edm']['theta_limit_rad']
    return {'upper_yukawa_squared':[a,b],'nonreal_cyclic_words_at_length_six':odd_necklaces(6),
            'uniform_bound_per_determinant':per,'uniform_bound_direct_theta':direct,
            'combined_induced_phase_bound':all_budget,'edm_working_limit':theta,
            'safety_factor_for_combined_bound':theta/all_budget,'J_weak_witness':float(J),
            'J_times_all_six_gaps_abs':abs(from_masses),'raw_Yukawa_path_overlap_imag':tagged,
            'right_unitary_theta_shift_for_coefficient_one':float(shift),
            'necessary_coefficient_bound_near_zero_single_source':theta/abs(2*tagged)}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    out=compute()
    target=HERE/'results.json'
    if p.parse_args().check:
        assert out==json.loads(target.read_bytes())
        print('PASS: all-order determinant bound, actual six-Yukawa invariant and exact tagged rotation')
    else:
        with target.open('x',encoding='utf-8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
