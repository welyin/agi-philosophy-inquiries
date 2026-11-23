"""Round 1629. Direct diagonalization of complete neutral-fermion mass matrices."""
import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def charge_minors(h):
    """Order by the omitted row, e, mu, tau."""
    return np.array([np.linalg.det(h[np.ix_(ix, ix)]).real
                     for ix in ((1, 2), (2, 0), (0, 1))])


def max_relative(a, b):
    return float(np.max(np.abs(np.array(a) - b)) / np.max(np.abs(b)))


def rotate(rng, n):
    x = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    return np.linalg.qr(x)[0]


def seesaw(rng, n_heavy):
    active = rotate(rng, 3)
    sterile = rotate(rng, n_heavy)
    y = rng.uniform(.17, .64)*(rng.normal(size=(2, n_heavy)) +
                              1j*rng.normal(size=(2, n_heavy)))
    d = active[:, :2] @ y
    r = sterile @ np.diag(np.linspace(3, 8, n_heavy)) @ sterile.T
    mass = np.block([[np.zeros((3, 3)), d], [d.T, r]])
    # Use right singular vectors, without solving the seesaw or sampling
    # a mixing parametrization. Nonzero Takagi values are nondegenerate here.
    _, sing, vh = np.linalg.svd(mass)
    mix = vh[::-1].conj().T
    diag = np.diag(mix.T @ mass @ mix)
    mix = mix * np.exp(-.5j*np.angle(diag))
    sing = sing[::-1]
    assert max_relative(mix.T @ mass @ mix, np.diag(sing)) < 1e-12
    assert sing[0] < 2e-14 and sing[1] > 1e-8
    assert sing[2] < sing[3]

    n, s = mix[:3, :3], mix[:3, 3:]
    h = s @ s.conj().T
    h_light = np.eye(3) - n @ n.conj().T
    assert np.linalg.norm(h-h_light) < 5e-14
    assert abs(np.linalg.norm(n[:, 0])**2-1) < 5e-13
    assert np.linalg.norm(s.conj().T @ n[:, 0]) < 3e-13
    p = charge_minors(h)
    assert max_relative(p / p.sum(), abs(n[:, 0])**2) < 1e-10
    assert max_relative(4*charge_minors(h/2), p) < 1e-12

    eff = sing[1]*sing[2]*np.linalg.det(np.eye(3)-h).real/p.sum()
    det_correction = np.linalg.det(np.eye(3)-h).real
    if n_heavy == 2:
        assert abs(eff / (sing[3]*sing[4]) - 1) < 1e-10
        # Recover both physical heavy masses using full complex light amplitude
        # and deficit. This check requires more than the three principal minors.
        lam, frame = np.linalg.eigh(h)
        lplus = frame[:, 1:].conj().T / np.sqrt(lam[1:])[:, None]
        light_amp = (n[:, 1:] * sing[1:3]) @ n[:, 1:].T
        t = -lplus @ light_amp @ lplus.T
        assert max_relative(np.sort(np.linalg.svd(t, compute_uv=False)), sing[3:]) < 1e-10
    all_products = [a*b for a, b in itertools.combinations(sing[3:], 2)]
    smallest_pair_error = min(abs(eff / prod - 1) for prod in all_products)
    return det_correction, smallest_pair_error


def dirac_pair_example(d):
    # Ordering active (e, mu, tau), then two Majorana sterile states,
    # then the two partners of the extra Dirac pair.
    mat = np.zeros((7, 7))
    for active, heavy, light_mass, heavy_mass in [(1,3,.1,2.), (2,4,.2,4.)]:
        mat[active, heavy] = mat[heavy, active] = (light_mass*heavy_mass)**.5
        mat[heavy, heavy] = heavy_mass-light_mass
    mat[0,5] = mat[5,0] = d
    mat[5,6] = mat[6,5] = 3
    eig, vectors = np.linalg.eigh(mat)
    order = np.argsort(abs(eig))
    sing = abs(eig[order])
    n, s = vectors[:3, order[:3]], vectors[:3, order[3:]]
    assert np.allclose(sing[:3], [0,.1,.2])
    h = s@s.T
    frac = charge_minors(h) / np.sum(charge_minors(h))
    left, _, right = np.linalg.svd(n, full_matrices=False)
    expected = np.abs((left@right)[:,0])**2
    assert np.allclose(expected, [1,0,0])
    expected_deficit = d*d/(9+d*d)
    assert np.allclose(h, np.diag([expected_deficit, 1/21, 1/21]))
    return [round(float(x),9) for x in frac], int(np.linalg.matrix_rank(h,tol=1e-13))


def check():
    rng = np.random.default_rng(1629)
    n_two = n_three = 64
    corr = [seesaw(rng,2)[0] for _ in range(n_two)]
    pair_errors = [seesaw(rng,3)[1] for _ in range(n_three)]
    assert min(corr) < .9
    assert max(pair_errors) > .05
    fraction0, rank0 = dirac_pair_example(0)
    fraction1, rank1 = dirac_pair_example(1.5)
    assert rank0 == 2 and rank1 == 3 and fraction1[0] < .2
    box = json.loads((HERE.parent/'1625'/'inputs.json').read_bytes(),parse_float=Fraction)['box']
    no = sorted((1-x)*(1-y) for x in box['sin2_theta12'] for y in box['sin2_theta13'])
    io = box['sin2_theta13']
    return {"two_singlet_full_complex_matrices_tested":n_two,
            "three_singlet_rank_two_full_complex_matrices_tested":n_three,
            "two_singlet_identity_and_recovery_passed":True,
            "three_singlet_all_three_minor_fractions_passed":True,
            "three_singlet_examples_with_effective_product_over_five_percent_from_every_pair":
                sum(1 for gap in pair_errors if gap > .05),
            "four_singlet_same_light_spectrum_polar_angles":
                {"before_fraction":fraction0, "after_fraction":fraction1, "deficit_ranks":[rank0,rank1]},
            "normal_electron_fraction_exact":[str(no[0]), str(no[-1])],
            "inverted_electron_fraction_exact":[str(x) for x in io]}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args=parser.parse_args()
    result=check()
    out=HERE/'independent_results.json'
    if args.check:
        assert result == json.loads(out.read_bytes())
    else:
        encoded=json.dumps(result,indent=2)+'\n'
        with out.open('x',encoding='utf-8') as fp:
            fp.write(encoded)
    print("1629 independent exact neutral mass matrix, cofactor and extra-pair checks passed")
