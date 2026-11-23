"""Dirac protection on the frozen two-flavor charge table, plus charged-basis oscillation witnesses."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
CFG = json.loads((HERE/'inputs.json').read_bytes())
OLD = runpy.run_path(str(ARCH/'1688/independent_check.py'))
CHARGED = runpy.run_path(str(ARCH/'1691/check.py'))
SIGNED = (-1,0,1)


def lepton_triality():
    cfg = json.loads((ARCH/'1688/inputs.json').read_bytes())
    # Only terms with nonzero L enter the three sums below. Charges ordered (F,chi,L).
    fixed = [(q,-1,-1) for q in cfg['su5_generation_charges_10']]  # e^c
    fixed += [(q,0,1) for q in cfg['su5_generation_charges_5bar'] for _ in range(2)]
    fixed += [(0,1,-1),(2,1,-1),(0,-1,1),(6,-1,1)]
    pure = sum(l**3 for _,_,l in fixed)
    fll = sum(q*l*l for q,_,l in fixed)
    fxl = sum(q*x*l for q,x,l in fixed)
    assert (pure%9, fll%3, fxl%3)==(3,0,2)
    pure_repairs = [ls for ls in itertools.combinations_with_replacement(SIGNED,4)
                    if (pure+sum(l**3 for l in ls))%9==0]
    assert pure_repairs == [(-1,-1,-1,0)]
    allowed_pairs = set()
    for fs in itertools.product(SIGNED, repeat=4):
        additional = list(zip(fs, [1]*4, pure_repairs[0]))
        allowed_pairs.add(((fll+sum(q*l*l for q,_,l in additional))%3,
                           (fxl+sum(q*x*l for q,x,l in additional))%3))
    assert (0,0) not in allowed_pairs
    return dict(pure_repair=list(pure_repairs[0]), FLL_and_FchiL_residue_pairs=sorted(map(list,allowed_pairs)))


def one_generation():
    # Multiplicity, 6Y, 3(B-L), color 2T times weak multiplicity, weak 2T times color multiplicity.
    rows = [(6,1,1,2,3),(3,-4,-1,1,0),(3,2,-1,1,0),(2,-3,-3,0,1),(1,6,3,0,0),(1,0,3,0,0)]
    sums = [sum(n*p(y,k) for n,y,k,_,_ in rows)
            for p in [lambda y,k:k, lambda y,k:k**3, lambda y,k:k*k*y, lambda y,k:k*y*y]]
    sums += [sum(c*k for _,_,k,c,_ in rows),sum(w*k for _,_,k,_,w in rows)]
    assert sums == [0]*6
    assert 2*3 % 5 != 0 and (3-3)%5==0
    return sums


def fit(w):
    phi = -math.pi/6
    s, x = .5*np.exp(1j*phi), .5
    vevs = np.array([s,x,s**2,s**3])
    unit = CHARGED['PORT']['UNIT'].copy()
    unit[2] *= w['neutral_modulus']
    a, b = unit*[0,1,1],unit*[1,3,1]
    m = np.outer(a,a)+w['source_r']*np.exp(2j*math.pi/3)*np.outer(b,b)
    h = m.conj().T @ m
    invariant = float((h[0,1]*h[1,2]*h[2,0]).imag/h.trace().real**3)
    assert abs(invariant)>1e-6
    fitted = np.zeros_like(m)
    costs = []
    max_coeff = 0.
    for i, q in enumerate([2,0,1]):
        budget, powers, mat = OLD['interpolation_basis'](q,vevs)
        costs.append(budget)
        for j, wanted in enumerate(m[i]):
            coeff = np.linalg.solve(mat,[wanted.real,wanted.imag])
            fitted[i,j] = sum(c*np.prod(vevs**p) for c,p in zip(coeff,powers))
            max_coeff = max(max_coeff, max(abs(coeff)))
    assert costs==[3,3,3]
    np.testing.assert_allclose(fitted.conj() @ fitted.T,h,atol=2e-14)
    charged_a, charged_b = CHARGED['charged_sources'](w)
    mc = np.outer(charged_a,charged_a)+w['source_r']*np.exp(2j*math.pi/3)*np.outer(charged_b,charged_b)
    vals, us = np.linalg.eigh(mc.conj().T@mc)
    scale = math.sqrt(w['d21_eV2']/vals[1])
    p = abs(us)**2
    pe = p[0,2]
    obs = [4*pe*(1-pe),p[0,1]/(1-pe),p[1,2]/(1-pe)]
    np.testing.assert_allclose(obs,w['three_angles'],atol=5e-12)
    return dict(R=w['R'],weak_frame_CP_invariant=round(invariant,10),
                insertions_needed_in_each_row_for_full_interpolation=costs,
                largest_real_Yukawa_coefficient_at_u174=float(f'{max_coeff*scale*1e-9/174:.9e}'),
                same_three_oscillation_readings=np.round(obs,9).tolist(),
                m_beta_meV=round(float(1000*scale*np.linalg.norm(mc[:,0])),7),
                reference_light_Majorana_mbb_meV=round(float(1000*scale*abs(mc[0,0])),7),
                dirac_delta_B_zero_delta_L_two_amplitude=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    ans = dict(triality=lepton_triality(), one_generation_mixed_anomaly_sums=one_generation(),
               charged_basis_examples=[fit(w) for w in json.loads(
                   (HERE/CFG['reference_neutrinos']).read_bytes())['witnesses']])
    dest = HERE/'results.json'
    if parser.parse_args().check:
        assert ans==json.loads(dest.read_bytes())
        print('PASS: mixed triality obstruction, order-five charge table and both actual charged-frame Dirac matrices')
    else:
        with dest.open('x',encoding='utf8') as fp:
            fp.write(json.dumps(ans,indent=2)+'\n')
