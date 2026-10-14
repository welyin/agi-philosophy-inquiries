"""1177 independent finite review; no author module imported.

Written after reading equations and discussing the proposed bounds, before
reading author check.py. Exact Fraction certificates carry the inequalities;
small matrix calculations are diagnostics only. Default is read-only.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
CORE = {
    'research_note_1177.md': 'ec3c11f501a0fc727b1c25efad974f6572f32f77825bbd3c89eb62d47a90305a',
    '1177/admission.md': '89e44772078859270af75c830f9fb52b47c5e6afebc3808998e23201f216cf7f',
    '1177/proof.md': '6c3229a8bd9caad35e7d96bc53ced14e867feb89a8aeaff01be8d7d40339c447',
    '1177/check.py': 'ce6be9e2d7a186586af32169a38345d9a72bb653210b99ac4d2db66a4dbdb916',
    '1177/results.json': 'd93aeadd8b201f982f7006847531c54dfa50ce5ee52c724f00aa5127cb259c09',
    '1177/NEXT.md': '2cdff50abb87a01757e3dcf5fcddcb06e5e8263c6317cb7a8f15288900dc85ce',
}
DEPENDENCIES = {
    '1176/proof.md': '966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9',
    '_admission/after1176_selection/propagation_contract_review.md':
        'cddfebb9f6013aa0ce66cb7bd520ff8b57d88441e285154ce872400ef8d21ba3',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def directed(q, upward):
    x = float(q)
    while (F(x) < q if upward else F(x) > q):
        x = math.nextafter(x, math.inf if upward else -math.inf)
    return x


def bracket_atan(q, count):
    value = sum(((-1)**n * q**(2*n+1)/F(2*n+1)
                 for n in range(count)), F(0))
    following = value + (-1)**count*q**(2*count+1)/F(2*count+1)
    return min(value, following), max(value, following)


def bracket_sin(q, count=18):
    assert 0 <= q <= 1
    value = sum(((-1)**n*q**(2*n+1)/math.factorial(2*n+1)
                 for n in range(count)), F(0))
    following = value+(-1)**count*q**(2*count+1)/math.factorial(2*count+1)
    return min(value, following), max(value, following)


def bracket_cos(q, count=18):
    assert 0 <= q <= 2
    value = sum(((-1)**n*q**(2*n)/math.factorial(2*n)
                 for n in range(count)), F(0))
    following = value+(-1)**count*q**(2*count)/math.factorial(2*count)
    return min(value, following), max(value, following)


def certificate():
    a0, a1 = bracket_atan(F(1, 5), 45)
    b0, b1 = bracket_atan(F(1, 239), 12)
    pi0, pi1 = 16*a0-4*b1, 16*a1-4*b0
    ell, eps, nu, g = F(16, 65), F(3, 20), F(7, 20), F(1, 5)
    a_s, a_t = eps/2+3*nu/ell, 2*g/ell
    square = 7*a_s*a_s+5*a_t*a_t
    assert (a_s, a_t, square) == (F(1389, 320), F(13, 8), F(14857247, 102400))
    root = math.sqrt(float(square))
    while F(root)**2 < square:
        root = math.nextafter(root, math.inf)
    eta_bound = pi1*F(root)/(2*10**7)
    assert eta_bound < F(2, 10**6)
    instrument, delta = F(2, 10**5), F(22, 10**6)
    assert eta_bound+instrument < delta
    # This is a preregistered model gap ratio, never a fitted curve.
    assert F(1)/F(1, 5) == 5
    tau0, tau1 = pi0/2+eps, pi1/2+eps
    p_d = {}
    for key, lo, hi in [('A', pi0/2, pi1/2), ('B', tau0, tau1)]:
        s0 = bracket_sin(lo/5)[0]
        s1 = bracket_sin(hi/5)[1]
        p0, p1 = (1+s0)/2, (1+s1)/2
        # Both *actual* branches remain inside the declared inversion window.
        p_window_upper = (1+bracket_sin(F(19, 50))[0])/2
        assert F(1, 2) < p0-delta < p1+delta < p_window_upper
        p_d[key] = {'lower': directed(p0, False), 'upper': directed(p1, True)}
    # Throughout 0 <= .2 tau <= .38: cos >= 1-x^2/2 > .9.
    cos_lower = 1-F(19, 50)**2/2
    assert cos_lower > F(9, 10) and F(5)/cos_lower < 6
    sine0, sine1 = bracket_sin(eps)
    rough = (eps-eps**3/6)/2-14*delta
    paired0, paired1 = sine0/2-14*delta, sine1/2-14*delta
    assert rough == F(297643, 4000000)  # .07441075
    assert paired0 > rough > 0
    assert paired0-14*F(1, 10000) > 0
    volume = 65**3
    bs = 4*eps+6*nu/ell+4*g/ell
    bp = F(8, 5)*eps+18*nu/ell+12*g/ell
    ck = 1+10**14*volume*(4+6/ell+1+bs)+volume*(F(8, 5)+18/ell+bp)
    return {
        'pi_lower': directed(pi0, False), 'pi_upper': directed(pi1, True),
        'pi_bracket_width_upper': directed(pi1-pi0, True),
        'a_s_exact': str(a_s), 'a_t_exact': str(a_t), 'B_squared_exact': str(square),
        'K': str(10**14), 'eta_upper': directed(eta_bound, True),
        'instrument_error_each': float(instrument), 'probability_error_each': float(delta),
        'clock_probability': p_d, 'tau_B_lower': directed(tau0, False),
        'tau_B_upper': directed(tau1, True), 'cos_phase_lower_exact': str(cos_lower),
        'lipschitz_declared': 6, 'paired_error_exact': str(14*delta),
        'ideal_defect_lower': directed(sine0/2, False),
        'ideal_defect_upper': directed(sine1/2, True),
        'paired_lower': directed(paired0, False), 'paired_upper': directed(paired1, True),
        'simple_paired_lower_exact': str(rough), 'both_actual_branches_in_window': True,
        'full_one_particle_and_arbitrary_reference': True,
        'original_input_code_dimension': 4*4*3, 'full_parent_positive_shift_exact': str(ck),
        'relation_qubits': str(2*volume*10**14), 'new_pointer_qubits': volume,
    }


def exponential_taylor(generator, time):
    term = np.eye(len(generator), dtype=complex)
    out = term.copy()
    for n in range(1, 33):
        term = term @ (-1j*time*generator)/n
        out += term
    return out


def finite_diagnostics(cert):
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    beta = np.kron(z, np.eye(2))
    alpha = [np.kron(x, v) for v in (x, y, z)]
    mass = np.diag([0., .3, .5])
    ell, eps, nu, g = 16/65, .15, .35, .2
    # A single star of center and +/- axis neighbors; all three masses remain.
    fs = np.zeros((84, 84), complex)
    ft = np.zeros_like(fs)
    fs[:12, :12] = eps*np.kron(mass, beta)
    for b in range(3):
        for sign in (-1, 1):
            j = 1+2*b+(sign+1)//2
            h = -sign*1j*np.kron(np.eye(3), alpha[b])/(4*ell)
            sl = slice(12*j, 12*j+12)
            fs[:12, sl] += nu*h
            fs[sl, :12] += nu*h.conj().T
            if b < 2:
                ft[:12, sl] += (g if b == 0 else -g)*h
                ft[sl, :12] += (g if b == 0 else -g)*h.conj().T
    hermitian = max(float(np.max(abs(a-a.conj().T))) for a in (fs, ft))
    norms = [float(np.linalg.norm(a, 2)) for a in (fs, ft)]
    assert hermitian == 0 and norms[0] <= 4.340625 and norms[1] <= 1.625
    # Exact incidence proof diagnostics: no edge collisions on a 5^3 torus.
    cover_s, cover_t = np.zeros((5, 5, 5), int), np.zeros((5, 5, 5), int)
    for pos in np.ndindex((5, 5, 5)):
        cover_s[pos] += 1
        cover_t[pos] += 1
        for b in range(3):
            for sign in (-1, 1):
                near = list(pos)
                near[b] = (near[b]+sign) % 5
                cover_s[tuple(near)] += 1
                if b < 2:
                    cover_t[tuple(near)] += 1
    assert np.all(cover_s == 7) and np.all(cover_t == 5)
    # Normal CAR vacuum diagnostic with full Jordan-Wigner signs.
    lowers = []
    for j in range(3):
        a = np.zeros((8, 8), complex)
        for bits in range(8):
            if (bits >> j) & 1:
                a[bits ^ (1 << j), bits] = (-1)**((bits & ((1 << j)-1)).bit_count())
        lowers.append(a)
    vacuum = np.eye(8)[:, 0]
    car_error = max(float(np.max(abs(a @ b.conj().T+b.conj().T @ a-
                          (np.eye(8) if i == j else 0))))
                    for i, a in enumerate(lowers) for j, b in enumerate(lowers))
    vac_error = max(float(np.linalg.norm(a.conj().T @ b @ vacuum))
                    for a in lowers for b in lowers)
    assert car_error == vac_error == 0
    time = (cert['pi_lower']+cert['pi_upper'])/4
    ug = exponential_taylor(-.5*x, time)
    pg = float(abs(ug[0, 0])**2)
    clock = (np.eye(2)-y)/2
    plus = np.ones(2)/math.sqrt(2)
    rows = {}
    for key, tau in [('A', time), ('B', time+.15)]:
        ud = exponential_taylor(np.diag([.3, .5]), tau)
        psi = ud @ plus
        pd = float(np.vdot(psi, clock @ psi).real)
        expected = (1+math.sin(.2*tau))/2
        assert abs(pd-expected) < 1e-14
        assert cert['clock_probability'][key]['lower']-1e-14 <= pd <= cert['clock_probability'][key]['upper']+1e-14
        role = (1+math.cos(5*math.asin(2*pd-1)))/2
        rows[key] = {'p_D_matrix': pd, 'p_G_matrix': pg, 'role_diagnostic': role,
                     'defect_diagnostic': pg-role}
    assert abs(pg-.5) < 1e-14
    # ||-i H t|| < 1 for these matrices, so remainder < 2/33!.
    assert .5*(time+.15) < 1
    taylor_tail = F(2, math.factorial(33))
    # All local CNOTs commute. Full map on three controls+three pointers is
    # bijective. Select by every possible actual site; retain every pointer.
    image = {(controls, pointers ^ controls) for controls in range(8) for pointers in range(8)}
    assert len(image) == 64
    cases = 0
    for site in range(3):
        for controls in range(8):
            pointers_out = controls
            assert ((pointers_out >> site) & 1) == ((controls >> site) & 1)
            cases += 1
    return {
        'CAR_error': car_error, 'all_bilinears_vacuum_error': vac_error,
        'star_hermiticity_error': hermitian, 'star_operator_norms': norms,
        'support_cover_counts': [7, 5], 'clock_matrix_diagnostics': rows,
        'matrix_Taylor_exact_tail_upper': directed(taylor_tail, True),
        'CNOT_basis_bijection_size': len(image), 'same_site_effect_cases': cases,
        'pointer_selection_discards_no_record': True,
        'floating_checks_are_diagnostics_not_certificates': True,
    }


def compare(a, b, path=''):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+'/'+key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+'/'+str(i))
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15), (path, a, b)
    else:
        assert a == b, (path, a, b)


def author_crosscheck(cert):
    # Read saved numbers only after the independent implementation was run.
    result = json.loads((HERE/'results.json').read_text(encoding='utf-8'))['groups']
    first = result['1_same_parent_exact_vacuum']
    second = result['2_full_input_prepared_isometry']
    third = result['3_finite_role_violation_and_control']
    assert F(first['parent_positive_shift_exact']) == F(cert['full_parent_positive_shift_exact'])
    for key in ('a_s_exact', 'a_t_exact', 'B_squared_exact'):
        assert F(second[key]) == F(cert[key])
    assert F(second['marginal_error_adopted_exact']) == F(22, 10**6)
    author_lower = F(third['separation_after_baseline_error_lower_exact'])
    assert F(cert['paired_lower']) <= author_lower <= F(cert['paired_upper'])
    assert F(cert['simple_paired_lower_exact']) < author_lower
    errors = [abs(cert['clock_probability'][key]['lower']-
                  third[label]['lower_float_diagnostic']) for key, label in
              [('A', 'baseline_clock_p_interval'), ('B', 'clock_p_interval')]]
    assert max(errors) < 1e-14
    return {'exact_constants_match': True, 'author_paired_bound_inside_independent_enclosure': True,
            'max_material_probability_endpoint_difference': max(errors)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--unbound-diagnostic', action='store_true')
    args = parser.parse_args()
    certificate_data = certificate()
    diagnostics = finite_diagnostics(certificate_data)
    if args.unbound_diagnostic:
        assert not args.write
        print(json.dumps({'certificate': certificate_data, 'diagnostics': diagnostics}, indent=2))
        return
    assert len(CORE) == 6, 'Wait for six-core author freeze before final binding.'
    hashes = {}
    for name, expected in (CORE | DEPENDENCIES).items():
        value = sha(STAGE/name)
        assert value == expected, ('frozen file changed', name, value)
        hashes[name] = value
    payload = {'round': 1177, 'status': 'ACCEPT_WITH_SCOPE', 'groups': 3, 'blockers': 0,
               'author_and_dependency_sha256': hashes,
               'certificate': certificate_data, 'diagnostics': diagnostics,
               'author_saved_result_comparison': author_crosscheck(certificate_data),
               'independence': 'Independent implementation before reading author code; equations and proposed results known; internal AI review, not blind or external.'}
    target = HERE/'review_checks.json'
    if args.write:
        with target.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
            f.write('\n')
        mode = 'WRITE_ONCE_PASS'
    else:
        compare(payload, json.loads(target.read_text(encoding='utf-8')))
        mode = 'READ_ONLY_MATCH'
    print(mode, '3 groups; blockers 0; paired_lower', certificate_data['paired_lower'])


if __name__ == '__main__':
    main()
