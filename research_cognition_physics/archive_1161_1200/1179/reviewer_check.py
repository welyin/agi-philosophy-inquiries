"""Independent 1179 review: no author code imports.

Written before reading author check.py. Rational elementary-function bounds
certify the witness; explicit Fourier and Clifford matrices cross-check the
old physical port dictionary. Default reads and compares, --write is exclusive.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
# Author declared these six files frozen before their final binding.
CORE = {
    'research_note_1179.md': '8c8cc42f1b2c7b2d39579459687a8d082e774e18f2051813152ae87b9d71bbfd',
    '1179/admission.md': '6e892f034895081576748a5f46966c7a1b1e4cb9d01273b8a70e24a6b5f188a8',
    '1179/proof.md': '6154136dc607b8dbba865deb26408133910cd663664bc165d92c7ebf7c66620d',
    '1179/check.py': '31e861284bef9e998f934fe680bbd7ebcd50bd0a20c9375308fc17ea28d15a39',
    '1179/results.json': '8d605e9607c5919fce0d5b17cac1962da7bee74a07c42a233853677ea6d60c9d',
    '1179/NEXT.md': 'e810558d464ffb6ca3946bedb749806ce81166dd75483e8e56da62fa57c8ce6e',
}
SOURCE = {
    '1178/proof.md': 'cd44222a8fedc2d45c497b2e0e02120e18408c3385d5a628e3f0066895f127cf',
    '1178/results.json': 'e6255dd8f47a4bc57e1cf8a1b8cd602d9ce3ac87098013bc0e3e26442ba13b98',
    '1178/check.py': '3c5e06e1a048ca2758c33f1b74f6c67b75e4c19ab3654b516a08e3e115ae6037',
    '1178/review.md': 'c7b76b46abe364ba5c2d031b732766d776937717e29ec6ce6a69cccde641d888',
    '1177/proof.md': '6c3229a8bd9caad35e7d96bc53ced14e867feb89a8aeaff01be8d7d40339c447',
    '1176/proof.md': '966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9',
    '1175/proof.md': '581427de3c8e03176a1145773d6bd0754930d849524110a22971c0e3f8e861f4',
    '_admission/after1174_direction/mean_field_bound.md': '11ac64ce824d3448dec97548efe3dda768bcda30f2e48e1cd88bfe511fcd6b3f',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def iv(x):
    x = Q(x)
    return x, x


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[1], -x[0]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    z = [a*b for a in x for b in y]
    return min(z), max(z)


def scale(x, c):
    return mul(x, iv(c))


def rounded(x, digits=70):
    den = 10**digits
    return Q((x[0]*den).__floor__(), den), Q((x[1]*den).__ceil__(), den)


def atan_point(x, n=70):
    s = sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(n)), Q(0))
    z = s + (-1)**n*x**(2*n+1)/Q(2*n+1)
    return min(s, z), max(s, z)


def pi_interval():
    return rounded(sub(scale(atan_point(Q(1, 5)), 16),
                       scale(atan_point(Q(1, 239)), 4)))


def trig_point(x, cosine=False, n=20):
    offset = 0 if cosine else 1
    s = sum(((-1)**j*x**(2*j+offset)/math.factorial(2*j+offset)
             for j in range(n)), Q(0))
    power = 2*n + offset
    remainder = abs(x)**power/math.factorial(power)
    return rounded((s-remainder, s+remainder))


def sin_interval(x):
    # This explicit monotone branch is independently checked by each caller.
    assert 0 <= x[0] <= x[1] < Q(3, 2)
    return trig_point(x[0])[0], trig_point(x[1])[1]


def cos_interval(x):
    assert 0 <= x[0] <= x[1] < Q(3, 2)
    return trig_point(x[1], True)[0], trig_point(x[0], True)[1]


def sqrt_interval(q):
    den = 10**70
    n = math.isqrt((q*den*den).numerator//(q*den*den).denominator)
    return Q(n, den), Q(n+1, den)


def exp_upper(x, n=80):
    assert 0 <= x < 2
    s = sum((x**j/math.factorial(j) for j in range(n+1)), Q(0))
    tail = x**(n+1)/math.factorial(n+1)/(1-x/Q(n+2))
    return s+tail


def down(x):
    y = float(x)
    while Q.from_float(y) > x:
        y = math.nextafter(y, -math.inf)
    return y


def up(x):
    y = float(x)
    while Q.from_float(y) < x:
        y = math.nextafter(y, math.inf)
    return y


def display(x):
    return {'lower': down(x[0]), 'upper': up(x[1])}


def interval_group():
    pi = pi_interval()
    assert Q(314159, 100000) < pi[0] < pi[1] < Q(314160, 100000)
    time = scale(pi, Q(1, 2))
    k = scale(sin_interval(scale(pi, Q(2, 65))), Q(65, 16))
    s = mul(k, time)
    a = scale(k, Q(1, 5))
    upper_phase = scale(s, Q(31, 20))
    lower_phase = scale(s, Q(9, 20))
    assert 0 < lower_phase[0] < upper_phase[1] < pi[0]/2
    ss, cs, sa, ca = sin_interval(s), cos_interval(s), sin_interval(a), cos_interval(a)
    target_witness = scale(add(mul(ss, sub(iv(1), ca)), scale(mul(cs, sa), 3)), Q(1, 4))
    settings_lipschitz = scale(mul(a, cos_interval(upper_phase)), Q(1, 8))
    linear_tv = scale(target_witness, Q(1, 6))
    actual_linear_tv = sub(linear_tv, iv(Q(16, 100000)))
    assert actual_linear_tv[0] > Q(79, 10000)
    # The pre-existing finite-K proof is uniform on every known pure cluster.
    a_s, a_t = Q(1389, 320), Q(13, 8)
    b2 = 7*a_s*a_s+5*a_t*a_t
    assert b2 == Q(14857247, 102400)
    factor = sqrt_interval(2*b2)[1]+a_s+a_t
    eta_upper = 4*time[1]*exp_upper(Q(9, 10)*time[1])*factor/Q(10**7)
    eta_upper = rounded(iv(eta_upper))[1]
    assert eta_upper < Q(6, 100000)
    return {
        'pi': display(pi), 'native_k': display(k),
        'phase_window': [down(lower_phase[0]), up(upper_phase[1])],
        'target_linear_witness': display(target_witness),
        'random_three_setting_TV_bound': display(settings_lipschitz),
        'convex_preparation_TV_bound': display(linear_tv),
        'actual_convex_TV_bound': display(actual_linear_tv),
        'actual_convex_lower_exact': str(actual_linear_tv[0]),
        'eta_PK_upper': up(eta_upper),
        'eta_PK_upper_exact': str(eta_upper),
        'each_side_probability_budget': '1/12500',
    }


def matrix_group():
    eye2 = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    alpha = [np.kron(x, a) for a in (x, y, z)]
    beta = np.kron(z, eye2)
    eye4 = np.eye(4)
    spin = (np.eye(4)-1j*np.kron(eye2, x+y+z))/2
    # Active cyclic rotation: old x -> y, y -> z, z -> x.
    mode = np.zeros((4, 4), complex)
    mode[0, 0] = 1
    for old, new in ((1, 2), (2, 3), (3, 1)):
        mode[new, old] = 1
    code_rotation = np.kron(mode, spin)
    clifford_error = max(np.linalg.norm(alpha[i]@alpha[j]+alpha[j]@alpha[i]
                                       - (2*eye4 if i == j else 0))
                         for i in range(3) for j in range(3))
    rotation_error = max(np.linalg.norm(spin@alpha[i]@spin.conj().T-alpha[(i+1)%3])
                         for i in range(3))
    assert np.linalg.norm(spin@beta@spin.conj().T-beta) < 1e-14
    k = math.sin(2*math.pi/65)/(16/65)
    t, g = math.pi/2, .2
    original = k*np.array([t+g, t-g, t])
    target = original[[2, 0, 1]]

    def unitary(phases):
        result = np.eye(16, dtype=complex)
        for i in range(3):
            result[4*(i+1):4*(i+2), 4*(i+1):4*(i+2)] = (
                math.cos(phases[i])*eye4-1j*math.sin(phases[i])*alpha[i])
        return result

    intertwining = np.linalg.norm(code_rotation@unitary(original)@code_rotation.conj().T
                                 - unitary(target))
    spinors = [np.kron(np.array([1, 1])/math.sqrt(2), a)
               for a in (np.array([1, 1])/math.sqrt(2),
                         np.array([1, 1j])/math.sqrt(2), np.array([1, 0]))]
    spatial_phase = 2*math.pi*np.arange(65)/65
    effect = (1+np.sin(spatial_phase))/2
    readouts = []
    direct_errors = []
    for i in range(3):
        chi = spinors[i]
        assert np.linalg.norm(alpha[i]@chi-chi) < 1e-14
        psi = ((1+np.exp(1j*spatial_phase-1j*target[i]))/
               math.sqrt(130))[:, None]*chi[None, :]
        p = float(np.sum(effect[:, None]*np.abs(psi)**2))
        predicted = .5+math.sin(target[i])/4
        readouts.append(p)
        direct_errors.append(abs(p-predicted))
        assert abs(float(np.sum(np.abs(psi)**2))-1) < 1e-14
    assert max(clifford_error, rotation_error, intertwining, *direct_errors) < 1e-13
    # Exact minimax certificate, not a numerical optimizer or a new result.
    aa = Q(1, 5)
    wanted = [Q(0), aa, -aa]
    uu, ap = -aa/4, -aa/2
    candidate = [uu+ap, uu-ap, uu]
    residual = [x-y for x, y in zip(wanted, candidate)]
    assert max(map(abs, residual)) == 3*aa/4
    assert sum(w*x for w, x in zip((1, 1, -2), residual)) == 3*aa
    return {
        'clifford_residual': float(clifford_error),
        'proper_cyclic_rotation_residual': float(rotation_error),
        'transported_native_unitary_residual': float(intertwining),
        'target_original_port_probabilities': readouts,
        'maximum_direct_Fourier_formula_error': max(direct_errors),
        'same_prepared_massless_subcode_dimension': 16,
        'full_original_prepared_code_dimension': 48,
        'instantaneous_minimax_for_a_one_fifth_exact': str(3*aa/4),
    }


def convex_group():
    k = math.sin(2*math.pi/65)/(16/65)
    t, g, nu = math.pi/2, .2, .35
    s0, a = k*t, k*g
    target = .5+np.sin([s0, s0+a, s0-a])/4
    weight = np.array([1, 1, -2.])
    rng = np.random.default_rng(117900)
    # This superset of reachable (S,D) gives a kinematic cross-check only.
    s = k*t*(1+nu*rng.uniform(-1, 1, 1000))
    d = k*g*t*rng.uniform(-1, 1, 1000)
    probabilities = .5+np.sin(np.stack([s+d, s-d, s], axis=-1))/4
    witness = probabilities@weight
    closed = .5*np.sin(s)*(np.cos(d)-1)
    assert np.max(witness) <= 1e-14
    assert np.max(np.abs(witness-closed)) < 1e-14
    mixtures = rng.dirichlet(np.ones(1000), size=20)@probabilities
    assert np.max(mixtures@weight) <= 1e-14
    tv = np.sum(np.abs(mixtures-target), axis=1)/3
    theoretical = (target@weight)/6
    assert np.min(tv) >= theoretical
    # Explicit randomized event keeps the three-setting record. For settings
    # x/y accept with probabilities 1/4 or 3/4; for z accept failure bits.
    event = (1.5+probabilities[:, 0]/2+probabilities[:, 1]/2-probabilities[:, 2])/3
    target_event = (1.5+target[0]/2+target[1]/2-target[2])/3
    assert np.max(event) <= 1/2+1e-14
    assert abs(target_event-(1/2+theoretical)) < 1e-14
    return {
        'sampled_identity_residual': float(np.max(np.abs(witness-closed))),
        'sampled_preparation_witness_maximum': float(np.max(witness)),
        'sampled_mixture_witness_maximum': float(np.max(mixtures@weight)),
        'sampled_mixture_TV_minimum': float(np.min(tv)),
        'fixed_event_target_probability': float(target_event),
        'fixed_event_candidate_upper_exact': '1/2',
        'diagnostics_are_not_universal_proof': True,
        'setting_must_be_independent_of_preparation_and_mixture_label': True,
    }


def run():
    assert len(CORE) == 6, 'Wait for author six-core freeze before binding.'
    for table in (CORE, SOURCE):
        for path, digest in table.items():
            assert sha(STAGE/path) == digest, path
    intervals = interval_group()
    archived = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    author_bound = archived['groups']['3_fixed_record_event_separates_the_cycle_target']
    author_eta = archived['groups']['2_all_preparations_uniform_full_output_bound']
    lower = Q(author_bound['formal_fixed_event_gap_lower_exact'])
    upper = Q(author_eta['eta_uniform_prepared']['upper_exact'])
    # An independent, tighter enclosure certifies the author's outward claims.
    assert lower <= Q(intervals['actual_convex_lower_exact'])
    assert upper >= Q(intervals['eta_PK_upper_exact'])
    assert lower > Q(79, 10000) and upper < Q(6, 100000)
    return {
        'status': 'PASS_WITH_SCOPE', 'groups': 3, 'blockers': 0,
        'algorithm': 'independent rational enclosure, Clifford/Fourier, convex record witness',
        'core_sha256': CORE, 'source_sha256': SOURCE,
        'reviewer_script_sha256': sha(Path(__file__)),
        'author_outward_bounds_confirmed_by_independent_tighter_enclosure': True,
        'groups_data': {'strict_intervals': intervals,
                        'original_port_and_rotation': matrix_group(),
                        'convex_preparation_record': convex_group()},
    }


def match(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            match(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            match(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-13), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'review_checks.json'
    if args.write:
        with target.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('WRITE_ONCE PASS_WITH_SCOPE')
    else:
        match(result, json.loads(target.read_text(encoding='utf-8')))
        print('READ_ONLY_MATCH PASS_WITH_SCOPE')
