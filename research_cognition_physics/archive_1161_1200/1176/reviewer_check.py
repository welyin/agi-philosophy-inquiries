"""1176 independent review: default read-only finite audit of frozen evidence.

The separately implemented Taylor dynamics and post-hoc polar replay were run
before this formal wrapper. --recompute-taylor / --recompute-polar reproduce
them without rewriting their files. No author integration code is imported by
the default calculation. --write creates review_checks.json exclusively.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
CAND = STAGE / '_admission' / 'after1175_selection' / 'tensor_candidate'
CORE = {
    'research_note_1176.md': 'a4edca5555f2cfe8aa64dbc809a38d17558c00f6cc04e87923c2ca763e460d63',
    '1176/admission.md': '4ae5e20a9c35fefaee2bcd6009b9470736d1c1f6b200ff8e4d227854d2349f49',
    '1176/proof.md': '966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9',
    '1176/check.py': '0face41bf7b67a5b3b7d974ae3b025af7241c0bfd3411bec4adc259193d27709',
    '1176/results.json': 'd9b815dec37cd1597202b61a683d768e4d5eb692170992fcaedfb0452f877873',
    '1176/NEXT.md': '8d662827c152698090e9552b044709ef20b0342738244f89aab279f83155410d',
}
EVIDENCE = {
    'candidate.md': '4d550d06b52295d7ff2404089f8e1153aa57a102febc07db37d470f3ce957fef',
    'reviewer_preview.py': 'eaf339ed21acc0535860c0e8143bec9d1bd3cfbbab40965454f737fb02e43722',
    'reviewer_preview.json': '1089c6b33a16f0a0d8d10e95ffcfa513d2dce48498876d1ed64b76e5c36b246f',
    'reviewer_scope.md': '19254e5f9407129b544f7e70517f63c5cc95460d50219ffdb3668aae15c532c1',
    'polar_endpoint_guard.py': '76b3dd71a7b86807ccd123953cdc537b67cb7dc66cd5c1abbb658cc95a273b9b',
    'polar_endpoint_guard.json': '4283338b2058ccb7e0bb07fc0558d79568375a69fedd54b2391920c55e77d093',
    'polar_endpoint_guard.md': '57d8ed0e6f3bb70a86289c603b42b0dfc62894b54b9db1910f5abc55a2babb95',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def compare(a, b, path=''):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            compare(a[k], b[k], path + '/' + k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path + '/' + str(i))
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15), (path, a, b)
    else:
        assert a == b, (path, a, b)


def up(q):
    value = float(q)
    while F(value) < q:
        value = float(np.nextafter(value, np.inf))
    return value


def down(q):
    value = float(q)
    while F(value) > q:
        value = float(np.nextafter(value, -np.inf))
    return value


def algebra():
    # Independently construct the Clifford and parity representation.
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    i2 = np.eye(2)
    beta = np.kron(z, i2)
    alpha = [np.kron(x, a) for a in (x, y, z)]
    gam = [beta] + alpha
    clifford = max(float(np.max(abs(a @ b + b @ a - (2*np.eye(4) if i == j else 0))))
                   for i, a in enumerate(gam) for j, b in enumerate(gam))
    ud = beta @ np.kron(i2, -1j*(x-y)/math.sqrt(2))
    perm = [1, 0, 2]
    parity = max(float(np.max(abs(ud @ alpha[b] @ ud.conj().T - alpha[perm[b]])))
                 for b in range(3))
    pperm = np.eye(4)[[0, 2, 1, 3]]
    r48 = np.kron(pperm, np.kron(np.eye(3), ud))
    p = math.pi/8
    h_a = np.zeros((48, 48), complex)
    h_b = h_a.copy()
    # Arbitrary asymmetric coefficients test the exact structural identity.
    lapse = .97
    velocities = [.82, 1.07, .94]
    for mode in range(4):
        for mass_index, mass in enumerate((0., .3, .5)):
            sl = slice(12*mode+4*mass_index, 12*mode+4*mass_index+4)
            h_a[sl, sl] = lapse*mass*beta
            h_b[sl, sl] = lapse*mass*beta
            if mode:
                axis = mode-1
                h_a[sl, sl] += p*velocities[axis]*alpha[axis]
                h_b[sl, sl] += p*velocities[perm[axis]]*alpha[axis]
    intertwiner = float(np.max(abs(r48 @ h_a @ r48.conj().T - h_b)))
    # Original Fourier menu: neither sampled sine nor special functions needed.
    clock = np.zeros((3, 3), complex)
    clock[1:, 1:] = (np.eye(2)-y)/2
    menu_error = 0.
    for axis in range(3):
        f = np.eye(4, dtype=complex)/2
        f[0, axis+1] = .25j
        f[axis+1, 0] = -.25j
        target = np.eye(4, dtype=complex)/2
        target[0, perm[axis]+1] = .25j
        target[perm[axis]+1, 0] = -.25j
        for internal in (np.eye(3), clock, np.diag([1., 0., 0.])):
            effect = np.kron(f, np.kron(internal, np.eye(4)))
            expected = np.kron(target, np.kron(internal, np.eye(4)))
            menu_error = max(menu_error, float(np.max(abs(r48 @ effect @ r48.conj().T-expected))))
    assert max(clifford, parity, intertwiner, menu_error) < 1e-13
    return {'Clifford_anticommutator_max': clifford, 'Dirac_parity_max': parity,
            'full48_generator_intertwiner_max': intertwiner,
            'all_analytic_menu_transport_max': menu_error}


def calculate():
    for name, expected in CORE.items():
        assert sha(STAGE/name) == expected, ('core changed', name)
    for name, expected in EVIDENCE.items():
        assert sha(CAND/name) == expected, ('review evidence changed', name)
    author = read(HERE/'results.json')
    for name, expected in author['sources'].items():
        assert sha(STAGE/name) == expected, ('author dependency changed', name)
    assert len(author['sources']) == 17
    independent = read(CAND/'reviewer_preview.json')
    numerical = read(CAND/'numerical_bound_checks.json')
    polar = read(CAND/'polar_endpoint_guard.json')
    for name, expected in independent['inputs_sha256'].items():
        if name == 'candidate.md':
            # The preliminary note received subsequent scope clarifications.
            # Preserve its original provenance; the current note and all formal
            # equations have separate frozen bindings, checked above.
            assert expected == 'aa493991057253fa9154b58347132bd01af938527488316bda3b830955eb2140'
        else:
            assert sha(CAND/name) == expected, ('Taylor input changed', name)

    # Exact rational Fock constants and integer source weights.
    ell, eps, nu, gt = F(16, 65), F(3, 20), F(7, 20), F(1, 5)
    volume = 65**3
    bs = 4*eps + (6*nu+4*gt)/ell
    bp = 2*eps*F(4, 5) + (18*nu+12*gt)/ell
    w, lp = 2*volume*bs, volume*bp
    mult = np.array([1]+[2]*32, dtype=np.int64)
    weights = mult[:, None, None]*mult[None, :, None]*mult[None, None, :]
    assert int(weights.sum()) == volume
    assert int(weights[0].sum()) == 65**2
    assert str(bs) == author['group1_parent_source_kernel']['bS_exact']
    assert str(bp) == author['group1_parent_source_kernel']['bP_exact']
    with localcontext() as ctx:
        ctx.prec = 70
        ld = Decimal(lp.numerator)/Decimal(lp.denominator)
        wd = Decimal(w.numerator)/Decimal(w.denominator)
        required = 2*((8*ld/Decimal('0.00001')).ln()/Decimal(10).ln()
                      + 12*wd/Decimal(10).ln())
        assert required+Decimal('1e-50') < 70881071

    # Exact signed-gap budget, independently rounded by rational comparisons.
    p_up = F('0.392699081698725')
    d_up = F('0.000611563643')
    j_up = F('0.045662192326534')
    v_up = F('1.868888431192047')
    assert d_up > ell**2*p_up**3/6
    assert j_up > F(numerical['J_integral_abs_gzt_upper'])
    assert v_up > F(numerical['V_integral_vbar_plus_abs_gzt_upper'])
    numeric_pair = F(9, 100000)
    assert F(numerical['numerical_error_pair']) < numeric_pair
    lattice_q = 2*d_up*j_up+4*p_up*d_up*j_up*v_up
    budget_q = lattice_q+F(2, 100000)+F(4, 100000)+numeric_pair
    lattice, budget = up(lattice_q), up(budget_q)
    menu = numerical['normalized_original_menu']
    y_gap = F(menu['y_clock_minusY_joint']['source_A'])-F(menu['y_clock_minusY_joint']['source_B'])
    x_gap = F(menu['x_massless_position']['source_A'])-F(menu['x_massless_position']['source_B'])
    y_lower = down(abs(y_gap)-budget_q)
    x_lower = down(abs(x_gap)-budget_q)
    separation = down(abs(y_gap)-budget_q-F(6, 100000))
    ab = author['group3_fixed_original_menu_and_budget']
    assert abs(budget-ab['total_pair_probability_budget']) < 2e-18
    assert abs(y_lower-ab['fixed_y_clock_site_lower']) < 2e-18
    assert abs(x_lower-ab['fixed_massless_x_lower']) < 2e-18
    assert abs(separation-ab['new_vs_old_scalar_menu_separation_margin_y']) < 2e-18
    assert y_lower > .00087483 and x_lower > .0030885 and separation > .00081483

    # Independently implemented integration is stored, with frozen provenance.
    max_menu = max(abs(value-independent['active_original_menu'][name][field])
                   for name, item in menu.items() for field, value in item.items())
    assert max_menu < 2e-13
    assert max(independent['convergence'].values()) < 3e-14
    provenance = independent['implementation_provenance']
    assert not provenance['author_preview_imported']
    assert not provenance['author_source_implementation_read_before_writing']
    assert provenance['author_reported_results_seen_before_run']
    records = independent['common_transport_and_reference']
    assert records['fine_outcomes'] == 3*volume
    assert records['record_R_error'] < 1e-13
    assert records['mass0_mass03_coherence_before_record'] > .18
    assert all(item['signed_gap'] == 0 for item in independent['scalar_original_menu'].values())

    # Recheck the post-hoc certificate with rational arithmetic: no eigensolver
    # error assumption is imported into the final probability claim.
    actual_delta = F(polar['exact_stored_minus_raw_entry_l1_rational'])
    assert F(polar['stored_minus_raw_operator_upper']) > actual_delta
    delta = F(polar['stored_to_exact_polar_upper'])
    assert delta >= actual_delta+F(polar['raw_to_exact_polar_upper'])
    assert F(polar['polar_probability_allowance']) >= (2+delta)*delta
    unit_roundoff = F(1, 2**53)
    gamma = 4096*unit_roundoff/(1-4096*unit_roundoff)
    assert F(polar['analytic_Fourier_probability_arithmetic_allowance']) >= 64*gamma
    polar_sum = sum(F(polar[name]) for name in (
        'polar_probability_allowance', 'witness_probability_allowance',
        'analytic_Fourier_probability_arithmetic_allowance'))
    assert polar_sum < F(1, 10**10)
    assert not polar['eigh_accuracy_assumed']
    assert polar['normalized_menu_bitwise_matched_float_leaves'] == 36

    # Weyl invariant is a comparator property. Positive Q suffices analytically.
    q = 2/16/math.tan(math.pi/130)
    zs2 = -2*(.15*4/65+.35*q)
    zt2 = -2*.2*q
    base = (.15-.35)*zs2
    a2 = np.array([base-.2*zt2, base+.2*zt2, base])
    electric = -.5*(a2-np.mean(a2))
    weyl = 8*float(electric@electric)
    assert abs(weyl-4*(.2*zt2)**2) < 1e-14
    assert weyl > .68
    assert abs(weyl-author['group2_autonomous_common_transport']['initial_comparator_Weyl_squared']) < 1e-11
    return {
        'round': 1176, 'date': '2026-10-13',
        'status': 'ACCEPT_WITH_SCOPE', 'blockers': 0, 'groups': 4,
        'author_core_sha256': CORE,
        'author_external_dependency_sha256': author['sources'],
        'independent_evidence_sha256': EVIDENCE,
        'historical_candidate_sha256_at_Taylor_write': independent['inputs_sha256']['candidate.md'],
        'candidate_revision_scope': 'Later prose clarifications; original Taylor evidence retained. Formal equations reviewed and bound through author six-core SHA.',
        'review_text_sha256': sha(HERE/'review.md'),
        'review_code_sha256': sha(HERE/'reviewer_check.py'),
        'group1_complete_parent_and_source': {
            'bS_exact': str(bs), 'bP_exact': str(bp), 'w_exact': str(w), 'L_exact': str(lp),
            'sufficient_log10_K_high_precision': str(required), 'K_power10': 70881071,
            'integer_source_weight_sum': int(weights.sum()),
            'zero_x_plane_weight': int(weights[0].sum()), 'source_particles': 2*volume,
            'full_parent_simulated': False, 'each_exact_BG_retained': True,
            'unknown_prepared_dimension': 48, 'arbitrary_R_by_operator_norm': True,
        },
        'group2_independent_dynamics_and_geometry': {
            'provenance': provenance, 'method': independent['method'],
            'convergence_diagnostic': independent['convergence'],
            'certified2048_menu_max_difference': max_menu,
            'independent_original_menu': independent['active_original_menu'],
            'invariants_diagnostic': independent['active_invariants'],
            'initial_Q_diagnostic': q, 'initial_a_second_diagnostic': a2.tolist(),
            'initial_Weyl_squared_diagnostic': weyl,
            'Weyl_measured_by_full_curvature_instrument': False,
            'traceless_is_velocity_direction_and_initial_metric_tangent': True,
        },
        'group3_original_menu_and_total_budget': {
            'p_upper_exact_decimal': str(p_up), 'd_upper_exact_decimal': str(d_up),
            'J_upper_exact_decimal': str(j_up), 'V_upper_exact_decimal': str(v_up),
            'correlated_lattice_upper': lattice, 'finite_K_pair': 2e-5,
            'instrument_pair_INPUT': 4e-5, 'numeric_pair_adopted': float(numeric_pair),
            'numeric_pair_certificate': numerical['numerical_error_pair'],
            'total_pair_upper': budget,
            'y_clock_site_lower': y_lower, 'massless_x_lower': x_lower,
            'old_scalar_finite_instrument_upper': 6e-5,
            'separation_from_old_scalar_y': separation,
            'posthoc_polar_probability_guard': up(polar_sum),
            'posthoc_eigh_accuracy_assumed': False,
            'mixed_bound_only_controls_A_B_difference': True,
            'all_unknown_inputs_have_nonzero_signal_claimed': False,
        },
        'group4_common_reflection_and_complete_records': {
            'independently_recomputed_algebra': algebra(),
            'record_and_reference_diagnostics': records,
            'scalar_original_menu': independent['scalar_original_menu'],
            'same_full_BG_energy_spectrum_not_same_fixed_probe_total_experiment': True,
            'scalar_exact_zero_only_Hartree_probe_channel': True,
        },
        'execution_scope': {
            'default': 'Read-only SHA binding, exact rational budget and finite algebra; uses already recomputed independent Taylor and polar evidence.',
            'author_integration_imported_in_default': False,
            'independent_taylor_command': 'python -B -X utf8 reviewer_check.py --recompute-taylor',
            'polar_replay_command': 'python -B -X utf8 reviewer_check.py --recompute-polar',
            'numerical_scope': 'Declared gamma_n arithmetic and mathematical residual bounds, not hardware IEEE formal verification.',
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--recompute-taylor', action='store_true')
    parser.add_argument('--recompute-polar', action='store_true')
    args = parser.parse_args()
    if args.recompute_taylor:
        m = module('round1176_independent_taylor', CAND/'reviewer_preview.py')
        fresh, saved = m.run(), read(CAND/'reviewer_preview.json')
        assert fresh['inputs_sha256']['candidate.md'] == EVIDENCE['candidate.md']
        assert saved['inputs_sha256']['candidate.md'] == 'aa493991057253fa9154b58347132bd01af938527488316bda3b830955eb2140'
        # Compare every computed value; only the explicitly recorded prose
        # revision is compared against its own binding, not the historical one.
        fresh['inputs_sha256'].pop('candidate.md')
        saved['inputs_sha256'].pop('candidate.md')
        m.compare(fresh, saved)
    if args.recompute_polar:
        m = module('round1176_polar_guard', CAND/'polar_endpoint_guard.py')
        compare(m.run(), read(CAND/'polar_endpoint_guard.json'))
    result = calculate()
    destination = HERE/'review_checks.json'
    if args.write:
        with destination.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write('\n')
        mode = 'WRITE_ONCE_PASS'
    else:
        compare(result, read(destination))
        mode = 'READ_ONLY_MATCH'
    print(json.dumps({'mode': mode, 'status': result['status'], 'blockers': 0, 'groups': 4,
                      'budget': result['group3_original_menu_and_total_budget']['total_pair_upper'],
                      'y_lower': result['group3_original_menu_and_total_budget']['y_clock_site_lower']}, indent=2))
