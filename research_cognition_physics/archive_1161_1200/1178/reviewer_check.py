"""1178 independent review. Equations discussed; author code never imported.

The independent analytic implementation was written before reading 1178's
author check.py. Fraction bounds certify inequalities; matrices diagnose the
same generators, reciprocal contraction, energy and original readouts.
Default is read-only; --write exclusively creates review_checks.json.
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
    'research_note_1178.md': '8a3744e80effff3001b5c9b15f5448ba4f4411d99e6b82074a5f8ed23a1d58bf',
    '1178/admission.md': '7820c20c15427c8cdaea24b949ea22cc02bcf448b573ccf806210ef888c96c23',
    '1178/proof.md': 'cd44222a8fedc2d45c497b2e0e02120e18408c3385d5a628e3f0066895f127cf',
    '1178/check.py': '3c5e06e1a048ca2758c33f1b74f6c67b75e4c19ab3654b516a08e3e115ae6037',
    '1178/results.json': 'e6255dd8f47a4bc57e1cf8a1b8cd602d9ce3ac87098013bc0e3e26442ba13b98',
    '1178/NEXT.md': '53fc81c569c250664da5a7c9e8924703c4f69d74a6a43a66004401b6740003b6',
}
DEPENDENCIES = {
    '1176/proof.md': '966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9',
    '1177/proof.md': '6c3229a8bd9caad35e7d96bc53ced14e867feb89a8aeaff01be8d7d40339c447',
    '1177/reviewer_check.py': '78b4fd8915f439f1f77a0cca31927d4a9c271f6eeb778ed6c4f5b56807867788',
    '1170/proof.md': 'c55f5abdc1503801b87fa03336d0710f48d723ecbcd4c133a3fe044a48e1b3af',
    '1171/proof.md': 'ecf48b81b9b42a21bcf7007b386bae16132cb9849e3fc939f02c398860e8d4ea',
    '_admission/after1174_direction/mean_field_bound.md':
        '11ac64ce824d3448dec97548efe3dda768bcda30f2e48e1cd88bfe511fcd6b3f',
    '_admission/after1177_reciprocal_clock/independent_selection.md':
        '1054b07a1579ec899b99e9f17d3487ac542ea27110c0c7bf961765495548e7cb',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(x):
    return (F(x), F(x))


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def scale(a, q):
    return mul(a, exact(q))


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def reciprocal(a):
    assert a[0] > 0
    return 1/a[1], 1/a[0]


def rounded(a, places=65):
    den = 10**places
    return F((a[0]*den).__floor__(), den), F((a[1]*den).__ceil__(), den)


def sqrt_interval(q):
    assert q >= 0
    den = 10**65
    n = math.isqrt((q.numerator*den*den)//q.denominator)
    lo, hi = F(n, den), F(n+1, den)
    assert lo*lo <= q <= hi*hi
    return lo, hi


def atan_interval(q, n):
    part = sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(n)), F(0))
    adjacent = part+(-1)**n*q**(2*n+1)/F(2*n+1)
    return min(part, adjacent), max(part, adjacent)


def trig_interval(a, cosine=False):
    assert 0 <= a[0] <= a[1] <= 2
    count, offset = 20, (0 if cosine else 1)
    out, power = exact(0), exact(1)
    for degree in range(2*count+offset):
        if degree % 2 == offset and degree >= offset:
            j = (degree-offset)//2
            out = add(out, scale(power, F((-1)**j, math.factorial(degree))))
        power = mul(power, a)
    order = 2*count+offset
    remainder = a[1]**order/math.factorial(order)
    return rounded((out[0]-remainder, out[1]+remainder))


def exp_interval(a):
    assert 0 <= a[0] <= a[1] < 4
    count = 70
    lo = sum((a[0]**j/math.factorial(j) for j in range(count)), F(0))
    hi = sum((a[1]**j/math.factorial(j) for j in range(count)), F(0))
    tail = a[1]**count/math.factorial(count)/(1-a[1]/(count+1))
    return rounded((lo, hi+tail))


def outward(q, upward):
    x = float(q)
    while (F(x) < q if upward else F(x) > q):
        x = math.nextafter(x, math.inf if upward else -math.inf)
    return x


def encoded(a):
    lo, hi = rounded(a, 35)
    return {'lower_exact': str(lo), 'upper_exact': str(hi),
            'lower': outward(lo, False), 'upper': outward(hi, True)}


def analytic_certificate():
    pi = rounded(sub(scale(atan_interval(F(1, 5), 48), 16),
                     scale(atan_interval(F(1, 239), 15), 4)))
    time = scale(pi, F(1, 2))
    eps, omega, lam = F(3, 20), F(1, 2), F(-3, 40)
    ell, nu, g, k = F(16, 65), F(7, 20), F(1, 5), 10**14
    assert lam == -eps*omega
    a_s, a_t = eps/2+3*nu/ell, 2*g/ell
    b_square = 7*a_s*a_s+5*a_t*a_t
    assert b_square == F(14857247, 102400)
    w = 2*abs(lam)
    root2b = sqrt_interval(2*b_square)
    exp6wt = exp_interval(scale(time, 6*w))
    # Candidate's stated coarse exact-BG bound, without tightening it to pass.
    eta = scale(mul(mul(time, exp6wt), add(root2b, exact(a_s+a_t))), F(4, 10**7))
    eta_adopted, local_adopted, instrument = F(6, 10**5), F(4, 10**7), F(2, 10**5)
    assert eta[1] < eta_adopted
    alpha_upper = (exp_interval(scale(time, 12*w))[1]-1)/k
    assert alpha_upper < local_adopted**2
    local_upper = sqrt_interval(alpha_upper)[1]
    e_material, e_relation = eta_adopted+instrument, eta_adopted+instrument+local_adopted
    rows = {}
    raw_tau = {}
    for label, q in [('plus', F(1, 2)), ('minus', F(-1, 2))]:
        b = lam*q
        r2 = omega**2+b*b
        r = sqrt_interval(r2)
        angle = mul(scale(r, 2), time)
        zs = scale(sub(exact(1), trig_interval(angle, True)), -omega*b/r2)
        integral = scale(sub(time, mul(trig_interval(angle), reciprocal(scale(r, 2)))), -omega*b/r2)
        tau = rounded(add(time, scale(integral, eps)))
        pd = scale(add(exact(1), trig_interval(scale(tau, F(1, 5)))), F(1, 2))
        visibility = sqrt_interval(1-q*q)
        pg = scale(add(exact(1), mul(visibility, trig_interval(tau, True))), F(1, 2))
        assert tau[0] > 0 and tau[1] < F(19, 10)
        rows[label] = {'q_exact': str(q), 'b_exact': str(b),
                       'zs': encoded(zs), 'integral_zs': encoded(integral),
                       'tau': encoded(tau), 'p_D': encoded(pd), 'p_G': encoded(pg),
                       'Hartree_initial_total_energy_exact': str(-omega-omega*q)}
        raw_tau[label] = tau
    tau_c = add(time, exact(eps))
    pd_c = scale(add(exact(1), trig_interval(scale(tau_c, F(1, 5)))), F(1, 2))
    pg_c = scale(sub(exact(1), trig_interval(exact(eps))), F(1, 2))
    rows['C'] = {'q_exact': '0', 'tau': encoded(tau_c), 'p_D': encoded(pd_c),
                 'p_G': encoded(pg_c), 'Hartree_initial_total_energy_exact': '0'}
    # Use raw certified endpoints, not midpoint subtraction.
    pd_plus = scale(add(exact(1), trig_interval(scale(raw_tau['plus'], F(1, 5)))), F(1, 2))
    pd_minus = scale(add(exact(1), trig_interval(scale(raw_tau['minus'], F(1, 5)))), F(1, 2))
    gap = sub(pd_plus, pd_minus)
    observed = sub(gap, exact(2*e_material))
    assert observed[0] > F(106, 10**5)
    cos_guard = 1-F(19, 50)**2/2
    assert cos_guard > F(9, 10) and 5/cos_guard < 6
    domain_upper = scale(add(exact(1), trig_interval(exact(F(19, 50)))), F(1, 2))[0]
    for row in rows.values():
        lo, hi = F(row['p_D']['lower_exact']), F(row['p_D']['upper_exact'])
        assert F(1, 2) < lo-e_material and hi+e_material < domain_upper
    role_error = e_relation+6*e_material
    # 1177's old single-branch failure certificate remains separated from the
    # new C branch's worst permitted closure error. This is a scope check only.
    old_failure_simple = F(3, 20)/2-F(3, 20)**3/12-7*F(22, 10**6)
    assert old_failure_simple-role_error > F(74, 1000)
    v = 65**3
    added_positive_shift = (k-1)*v*abs(lam)
    return {
        'lambda_exact': str(lam), 'E0_witness_exact': '0', 'K': str(k),
        'time': encoded(time), 'B_squared_exact': str(b_square),
        'w_site_exact': str(w), 'eta_exact_BG_upper': outward(eta[1], True),
        'eta_adopted_exact': str(eta_adopted),
        'local_Gt_counting_upper': outward(local_upper, True),
        'local_Gt_adopted_exact': str(local_adopted),
        'instrument_each_exact': str(instrument), 'p_D_each_error_exact': str(e_material),
        'p_G_each_error_exact': str(e_relation), 'clock_role_curve_error_exact': str(role_error),
        'added_common_positive_shift_exact': str(added_positive_shift),
        'branches': rows, 'ideal_original_material_gap': encoded(gap),
        'actual_original_material_gap_lower': outward(observed[0], False),
        'all_actual_material_branches_in_inverse_window': True,
        'original_code_dimension': 48, 'arbitrary_reference': True,
        'exact_BG_retained': True, 'whole_BG_Hartree_product_not_claimed': True,
    }


def taylor_u(h, time):
    term = np.eye(len(h), dtype=complex)
    result = term.copy()
    for n in range(1, 37):
        term = term @ (-1j*time*h)/n
        result += term
    return result


def matrix_diagnostics(cert):
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    ident = np.eye(2)
    sx, sz, tx = np.kron(x, ident), np.kron(z, ident), np.kron(ident, x)
    ht, bt = -.5*x, -.075*x
    central = .125*ident  # Algebraic phase-kernel check, not a new E0 experiment.
    difference = bt+central-.15*ht
    central_residual = max(float(np.max(abs(difference @ a-a @ difference))) for a in (x, y, z))
    assert central_residual < 1e-15
    columns = []
    for a in (ident, x, y, z):
        columns.append(np.concatenate([(a @ b-b @ a).reshape(-1) for b in (x, y, z)]))
    kernel_rank = int(np.linalg.matrix_rank(np.array(columns).T))
    assert kernel_rank == 3
    bcluster = -.075*tx
    pair = np.kron(sz, bcluster)+np.kron(bcluster, sz)
    swap = np.zeros((16, 16), complex)
    for a in range(4):
        for b in range(4):
            swap[4*b+a, 4*a+b] = 1
    swap_error = float(np.max(abs(swap @ pair @ swap-pair)))
    assert swap_error == 0
    wnorm = float(np.linalg.norm(pair, 2))
    assert abs(wnorm-.15) < 1e-14
    # A generic, known cluster checks BOTH contractions of the same W.
    s0 = np.array([math.sqrt(.8), math.sqrt(.2)], complex)  # z_s=.6
    q0, r0 = .5, math.sqrt(.75)
    t0 = np.array([math.sqrt((1+r0)/2), q0/math.sqrt(2*(1+r0))], complex)
    phi = np.kron(s0, t0)
    contracted = np.einsum('b,abcd,d->ac', phi.conj(), pair.reshape(4, 4, 4, 4), phi)
    wanted = (-.075*q0)*sz+.6*bcluster
    contraction_error = float(np.max(abs(contracted-wanted)))
    assert contraction_error < 1e-14
    clock = (ident-y)/2
    plus, zplus = np.ones(2)/math.sqrt(2), np.array([1., 0.])
    time = (cert['time']['lower']+cert['time']['upper'])/2
    output = {}
    for label in ('plus', 'minus', 'C'):
        row = cert['branches'][label]
        q = float(F(row['q_exact']))
        tau = (row['tau']['lower']+row['tau']['upper'])/2
        initial_s = zplus if label == 'C' else plus
        r = math.sqrt(1-q*q)
        initial_t = np.array([math.sqrt((1+r)/2), q/math.sqrt(2*(1+r))], complex)
        hs = -.5*x-.075*q*z
        state_s = taylor_u(hs, time) @ initial_s
        state_t = taylor_u(ht, tau) @ initial_t
        state_d = taylor_u(np.diag([.3, .5]), tau) @ plus
        xs = float(np.vdot(state_s, x @ state_s).real)
        zs = float(np.vdot(state_s, z @ state_s).real)
        qt = float(np.vdot(state_t, x @ state_t).real)
        pg = float(abs(state_t[0])**2)
        pd = float(np.vdot(state_d, clock @ state_d).real)
        energy = -.5*xs-.5*qt-.075*zs*qt
        energy_error = abs(energy-float(F(row['Hartree_initial_total_energy_exact'])))
        assert energy_error < 1e-14 and abs(qt-q) < 1e-14
        for actual, key in [(pg, 'p_G'), (pd, 'p_D')]:
            assert row[key]['lower']-1e-14 <= actual <= row[key]['upper']+1e-14
        if label != 'C':
            assert row['zs']['lower']-1e-14 <= zs <= row['zs']['upper']+1e-14
        role = (1+r*math.cos(5*math.asin(2*pd-1)))/2
        role_residual = pg-role
        assert abs(role_residual) < 1e-14
        output[label] = {'p_D_matrix': pd, 'p_G_matrix': pg, 'zs_matrix': zs,
                         'q_conservation_error': abs(qt-q), 'energy_error': energy_error,
                         'role_curve_residual': role_residual}
    return {'central_commutator_error': central_residual, 'traceless_response_rank': kernel_rank,
            'pair_exchange_error': swap_error, 'pair_operator_norm': wnorm,
            'same_pair_both_Hartree_contractions_error': contraction_error,
            'branches': output,
            'matrices_are_diagnostics_only': True,
            'original_CNOT_effect_reused_for_arbitrary_BG_state': True}


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


def saved_author_comparison(cert):
    data = json.loads((HERE/'results.json').read_text(encoding='utf-8'))['groups']
    source = data['2_exact_BG_and_unknown_probe_finite_K']['certificate']
    assert F(data['1_selected_response_and_complete_parent']['positive_shift_extra_exact']) == F(cert['added_common_positive_shift_exact'])
    pairs = [('lambda_exact', 'lambda_selected_exact'), ('B_squared_exact', 'B_squared_exact'),
             ('w_site_exact', 'w_site_exact'), ('eta_adopted_exact', 'eta_adopted_exact'),
             ('local_Gt_adopted_exact', 'local_Gt_adopted_exact'),
             ('p_D_each_error_exact', 'D_error_each_exact'),
             ('p_G_each_error_exact', 'G_local_error_each_exact'),
             ('clock_role_curve_error_exact', 'C_role_residual_upper_exact')]
    for mine, theirs in pairs:
        assert F(cert[mine]) == F(source[theirs])
    checked = 0
    for label, fields in [('plus', ['tau_plus', 'p_D_plus', 'p_G_plus_Hartree_local']),
                          ('minus', ['tau_minus', 'p_D_minus', 'p_G_minus_Hartree_local']),
                          ('C', ['C_tau', 'C_p_D', 'C_p_G'])]:
        for mine, theirs in zip(('tau', 'p_D', 'p_G'), fields):
            own = cert['branches'][label][mine]
            theirs = source[theirs]
            assert F(own['lower_exact']) <= F(theirs['lower_exact']) <= F(theirs['upper_exact']) <= F(own['upper_exact'])
            checked += 1
    author_lower = F(source['finite_clock_gap_lower_exact'])
    assert F(cert['actual_original_material_gap_lower']) <= author_lower
    assert author_lower <= F(cert['ideal_original_material_gap']['upper_exact'])-F(1, 6250)
    assert F(source['eta_prepared_interval']['upper_exact']) <= F(cert['eta_exact_BG_upper'])
    assert F(source['sqrt_alpha_interval']['upper_exact']) <= F(cert['local_Gt_counting_upper'])
    return {'exact_constants_match': True, 'author_intervals_within_independent_outward_enclosures': checked,
            'finite_material_gap_lower_verified': True, 'both_background_error_bounds_verified': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--unbound-diagnostic', action='store_true')
    args = parser.parse_args()
    cert = analytic_certificate()
    diagnostics = matrix_diagnostics(cert)
    if args.unbound_diagnostic:
        assert not args.write
        print(json.dumps({'certificate': cert, 'diagnostics': diagnostics}, indent=2))
        return
    assert len(CORE) == 6, 'Author six-core freeze required before formal binding.'
    hashes = {}
    for name, expected in (CORE | DEPENDENCIES).items():
        actual = sha(STAGE/name)
        assert actual == expected, ('frozen file changed', name, actual)
        hashes[name] = actual
    payload = {'round': 1178, 'status': 'ACCEPT_WITH_SCOPE', 'groups': 3, 'blockers': 0,
               'author_and_dependency_sha256': hashes, 'certificate': cert,
               'diagnostics': diagnostics,
               'saved_author_result_comparison': saved_author_comparison(cert),
               'independence': 'Analytic implementation before author code; equations and expected scale discussed; internal AI review, not blind or external.'}
    target = HERE/'review_checks.json'
    if args.write:
        with target.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
            f.write('\n')
        mode = 'WRITE_ONCE_PASS'
    else:
        compare(payload, json.loads(target.read_text(encoding='utf-8')))
        mode = 'READ_ONLY_MATCH'
    print(mode, '3 groups; blockers 0; p_D gap lower', cert['actual_original_material_gap_lower'])


if __name__ == '__main__':
    main()
