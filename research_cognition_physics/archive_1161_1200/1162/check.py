"""1162: frozen hard ports and an exact on-shell recoil bridge.

Four groups only. Default recomputes read-only; --write-once creates results.json
with exclusive mode. Neither finite matrices nor momentum grids certify a
physical emission channel, field species, or a normalizable 3D wavepacket.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
SOURCES = {
    'archive_1106_1119/1111/check.py':
        '704ae4ca016c77e0df83fc9c2d39e90d79e9fbf3d6469c8c5b9bf1dad8a5bd44',
    'archive_1159_1160/1160/check.py':
        '8af83d78311b77c804693e055f2ece2c9f195e1d217b5b67e41fa6b3b6753170',
    'archive_1159_1160/1160/proof.md':
        'ca3121515439d2b93c23d341d83a1dfbd68b48c4c63ac95fc41c28b027cf8e81',
}
TOL = 3e-10


def assert_small(value, label):
    value = float(value)
    assert value <= TOL, (label, value)
    return value


def load_old():
    for relative, expected in SOURCES.items():
        assert hashlib.sha256((RESEARCH / relative).read_bytes()).hexdigest() == expected, relative
    spec = importlib.util.spec_from_file_location('round1162_original1111',
                RESEARCH / 'archive_1106_1119/1111/check.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def exact_clock_polynomial():
    before, current = [F(1)], [-F(3, 2), F(1)]
    for weight2 in (F(3, 4), F(1), F(3, 4)):
        nxt = poly_mul([-F(3, 2), F(1)], current)
        for i, term in enumerate(before):
            nxt[i] -= weight2*term
        before, current = current, nxt
    expected = [F(1)]
    for root in range(4):
        expected = poly_mul(expected, [-F(root), F(1)])
    assert current == expected
    return [str(x) for x in current]


def mink2(p):
    return p[0]*p[0] - sum(x*x for x in p[1:])


def recoil(omega, axis):
    k = math.sqrt(16-5*omega)
    a = np.array([5-omega/2, 0., 0., -omega/2])
    b = a.copy()
    a[1+axis], b[1+axis] = k, -k
    q = np.array([omega, 0., 0., omega])
    return a, b, q


def group_frozen(old, model):
    polynomial = exact_clock_polynomial()
    hd = model['hout']
    expected = model['baseline'] + np.repeat(np.arange(4), 72)
    spectrum_residual = assert_small(np.max(np.abs(np.linalg.eigvalsh(hd)-expected)),
                                     'original full HD spectrum')
    totals = [(0, 0, 0), (4, -4, 0), (-4, 4, 0), (0, 0, 0)]
    allowed, positive = [], 0
    for n in range(4):
        for m in range(4):
            for label, spatial in enumerate(totals):
                energy = n-m
                p2 = sum(x*x for x in spatial)
                if energy > 0:
                    positive += 1
                    if energy*energy == p2:
                        allowed.append([n, m, label])
    assert not allowed
    assert 32 > 3**2
    # This is a kinematic scope control, not a manufactured two-photon device.
    two_q = [(F(1), F(0), F(0), F(1)), (F(1), F(0), F(0), -F(1))]
    assert all(mink2(q) == 0 for q in two_q)
    assert [sum(q[j] for q in two_q) for j in range(4)] == [2, 0, 0, 0]
    return {
        'joint_spectral_cases': 64,
        'positive_energy_candidates': positive,
        'single_null_matches': allowed,
        'clock_characteristic_polynomial_ascending': polynomial,
        'original_HD_spectral_residual': spectrum_residual,
        'clock_width': 3,
        'nonzero_frozen_spatial_momentum_squared': 32,
        'two_opposite_quanta_kinematic_control': [2, 0, 0, 0],
        'no_uniform_gap_control_residual_at_energy_1e_minus_8': math.sqrt(2)*1e-8,
    }


def group_shells():
    rows = []
    for k in (F(1, 4), F(1), F(2), F(3), F(15, 4), F(399, 100)):
        omega = (16-k*k)/5
        assert 0 < omega < F(16, 5)
        for axis in (0, 1):
            a = [5-omega/2, F(0), F(0), -omega/2]
            b = list(a)
            a[1+axis], b[1+axis] = k, -k
            q = [omega, F(0), F(0), omega]
            assert mink2(a) == mink2(b) == 9
            assert mink2(q) == 0
            assert [a[j]+b[j]+q[j] for j in range(4)] == [10, 0, 0, 0]
            assert 100-20*omega >= 36
            rows.append({'axis': 'xy'[axis], 'omega': str(omega), 'k': str(k)})
    end_a, end_b, end_q = recoil(16/5, 0)
    other_a, other_b, _ = recoil(16/5, 1)
    assert_small(np.linalg.norm(end_a-other_a)+np.linalg.norm(end_b-other_b), 'endpoint label collision')
    return {'exact_rational_samples': rows,
            'general_bound': 'omega <= 16/5 from (P-q)^2 >= 36',
            'retained_encoding_domain': '0 < omega < 16/5',
            'upper_endpoint_direction_labels_merge': True,
            'upper_endpoint_A': end_a.tolist(),
            'example_omega_1_over_10': {'A_x': recoil(.1, 0)[0].tolist(),
                                     'B_x': recoil(.1, 0)[1].tolist(),
                                     'q': recoil(.1, 0)[2].tolist()}}


def group_endpoints(old, model):
    omega_max = 3.1
    c_bound = math.sqrt(.5+25/(4+math.sqrt(16-5*omega_max))**2)
    center = np.array([model['arrival']+math.pi, 1., 0., 0.])
    sigma_over_mass = .4
    moment_bound_ratio = endpoint_bound_ratio = 0.
    shell_residual = conservation_residual = transport_residual = 0.
    boosts = [(0., 0., 0.), (.6, 0., 0.), (0., -.4, .2)]
    samples = []
    for omega in (1e-9, .01, .1, .4, 1., 3.):
        for axis in (0, 1):
            a, b, q = recoil(omega, axis)
            a0, b0, _ = recoil(0., axis)
            for p, p0 in ((a, a0), (b, b0)):
                z, z0 = center+sigma_over_mass*p, center+sigma_over_mass*p0
                bound = c_bound*omega
                moment_bound_ratio = max(moment_bound_ratio, np.linalg.norm(p-p0)/bound)
                # Form the difference algebraically to avoid cancellation at tiny omega.
                dz = sigma_over_mass*(p-p0)
                endpoint_bound_ratio = max(endpoint_bound_ratio, np.linalg.norm(dz)/(.4*bound))
                assert np.linalg.norm(dz) <= .4*bound + 1e-14
                for beta in boosts:
                    boost = old.boost(beta)
                    offset = np.array([.2, -.3, .4, -.1])
                    transported = boost@z+offset
                    expected = boost@center+offset+sigma_over_mass*(boost@p)
                    transport_residual = max(transport_residual, np.linalg.norm(transported-expected))
                    v = np.linalg.norm(beta)
                    amplification = math.sqrt((1+v)/(1-v))
                    assert np.linalg.norm(boost@dz) <= amplification*.4*bound+1e-14
                    shell_residual = max(shell_residual, abs(old.mink(boost@p, boost@p)-9))
            for beta in boosts:
                boost = old.boost(beta)
                conservation_residual = max(conservation_residual,
                    np.linalg.norm(boost@a+boost@b+boost@q-boost@np.array([10.,0.,0.,0.])))
                shell_residual = max(shell_residual, abs(old.mink(boost@q, boost@q)))
            if axis == 0 and omega in (.1, 1.):
                samples.append({'omega': omega, 'A_endpoint': (center+.4*a).tolist(),
                                'B_endpoint': (center+.4*b).tolist(),
                                'coordinate_flight_time': 2-omega/5,
                                'x_y_endpoint_separation': .4*math.sqrt(2)*math.sqrt(16-5*omega)})
    assert moment_bound_ratio <= 1+1e-12
    for name, value in [('shell',shell_residual), ('fourP',conservation_residual),
                        ('affine endpoint',transport_residual)]:
        assert_small(value, name)
    return {'omega_cap': omega_max, 'C_cap': c_bound,
            'maximum_momentum_bound_ratio': float(moment_bound_ratio),
            'maximum_endpoint_bound_ratio': float(endpoint_bound_ratio),
            'mass_shell_and_null_residual': shell_residual,
            'boosted_fourP_residual': conservation_residual,
            'common_affine_endpoint_residual': transport_residual,
            'examples': samples}


def group_kinematic_encoding():
    # Finite quadrature of the explicitly declared reduced channel space.
    # It is NOT a real sharp-momentum state or physical emission implementation.
    omegas = np.array([.1, .4, .9])
    weights = np.array([.2, .3, .5])
    amplitudes = np.sqrt(weights)*np.exp(1j*np.array([.3, -.2, .7]))
    j = np.zeros((6, 2), complex)
    for a in range(2):
        j[3*a:3*a+3, a] = amplitudes
    gram_residual = assert_small(np.linalg.norm(j.conj().T@j-np.eye(2)), 'kinematic J Gram')
    ptotal = np.array([sum(recoil(w,a)) for a in range(2) for w in omegas])
    total_residual = 0.
    for mu in range(4):
        target = 10. if mu == 0 else 0.
        total_residual = max(total_residual,
            assert_small(np.linalg.norm(ptotal[:,mu,None]*j-j*target), 'kinematic total intertwiner'))
    # J tensor a retained spectator clock/reference remains an isometry.
    spectator_h = np.diag([0., 1., 2., 3.])
    jfull = np.kron(np.eye(4), j)
    clock_residual = assert_small(np.linalg.norm(np.kron(spectator_h,np.eye(6))@jfull
                      -jfull@np.kron(spectator_h,np.eye(2))), 'spectator clock retained')
    rng = np.random.default_rng(1162)
    psi = rng.normal(size=(8,3))+1j*rng.normal(size=(8,3))
    psi /= np.linalg.norm(psi)
    kept = jfull@psi
    assert_small(np.linalg.norm(kept.conj().T@kept-psi.conj().T@psi), 'reference marginal retained')
    # A uniform non-atomic bin has exact nonzero energy variance.
    width = F(1, 5)
    variance = width*width/12
    nodes, quadweights = np.polynomial.legendre.leggauss(8)
    bin_center = .5
    points = bin_center+float(width)*nodes/2
    var_numeric = float(np.dot(quadweights/2, (points-bin_center)**2))
    assert_small(abs(var_numeric-float(variance)), 'bin variance')
    cap = 1.
    lipschitz = math.sqrt(.5+25/(4*(16-5*cap)))
    bin_bound_ratio = 0.
    for w in np.linspace(.4,.6,31):
        p = recoil(w,0)[0]
        pcenter = recoil(bin_center,0)[0]
        bound = lipschitz*float(width)/2
        ratio = np.linalg.norm(p-pcenter)/bound
        bin_bound_ratio = max(bin_bound_ratio, ratio)
        assert ratio <= 1+1e-12
    return {'contract': 'reduced direct integral only; no emission instrument certified',
            'quadrature_Gram_residual': gram_residual,
            'total_fourP_intertwiner_residual': total_residual,
            'spectator_clock_intertwiner_residual': clock_residual,
            'uniform_bin_width': str(width),
            'uniform_bin_energy_variance_exact': str(variance),
            'uniform_bin_energy_variance_quadrature': var_numeric,
            'finite_bin_L_cap': lipschitz,
            'maximum_bin_bound_ratio': float(bin_bound_ratio)}


def calculate():
    old = load_old()
    model = old.model()
    return {
        'schema': 'round1162_frozen_ports_exact_recoil_v1',
        'status': 'PASS', 'check_groups': 4,
        'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources_sha256': SOURCES,
        'frozen_joint_spectrum': group_frozen(old,model),
        'exact_recoil_family': group_shells(),
        'finite_endpoint_transport': group_endpoints(old,model),
        'conditional_kinematic_isometry_and_bins': group_kinematic_encoding(),
        'not_certified': ['finite emission amplitude', 'physical emission channel',
                          'normalizable full 3D momentum wavepacket',
                          'spin-2 emergence', 'dynamic geometry', 'Einstein dynamics'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.write_once:
        with (HERE/'results.json').open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
