"""790: same-Wick first generator anomaly and its anchored path repair.

Exact finite graded calibrations of the analytic map in research_note_790.md.
The full nonminimal quartet and a fermion pair are kept in these diagnostics;
the original continuum state/anomaly is not replaced by these finite matrices.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv790', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
# Extend only this fresh module instance; the frozen 777 engine is unchanged.
bv.EVEN = bv.EVEN+('w', 'psi_star', 'barpsi_star', 'tau')
bv.ODD = bv.ODD+('w_star', 'psi', 'barpsi', 'chi')
bv.NAMES = bv.EVEN+bv.ODD
bv.PAIRS = bv.PAIRS+(('w', 'w_star'), ('psi', 'psi_star'), ('barpsi', 'barpsi_star'))
bv.ZERO = (0,)*len(bv.NAMES)
bv.ONE = {bv.ZERO: Q(1)}
V = {name: bv.var(name) for name in bv.NAMES}
FIELDS = ('q', 'r', 'w', 'b', 'c', 'h', 'psi', 'barpsi')
ANTIS = dict(bv.PAIRS)
project = lambda f: bv.substitute(f, {anti: {} for _, anti in bv.PAIRS})


def gamma(poly, kernel):
    return bv.add(*(bv.scale(bv.diff(bv.diff(poly, a), b), value/2)
                    for (a, b), value in kernel.items()))


def kernel_entry(kernel, a, b, value):
    kernel[(a, b)] = value
    kernel[(b, a)] = -value if a in bv.ODD and b in bv.ODD else value


def setup():
    q, r, w, b, c, h, p, a = [V[x] for x in ('q', 'r', 'w', 'b', 'c', 'h', 'p', 'a')]
    s2 = bv.add(bv.scale(bv.mul(r, r), Q(3, 4)), bv.mul(p, c), bv.mul(b, q),
                bv.scale(bv.mul(b, b), Q(-2, 3)), bv.scale(bv.mul(h, c), -1),
                bv.scale(bv.mul(a, b), -1), bv.scale(bv.mul(V['barpsi'], V['psi']), Q(7, 5)))
    assert not bv.bracket(s2, s2)
    s0 = lambda f: bv.bracket(s2, f)
    vectors = {
        'q': bv.add(bv.mul(q, r), bv.scale(bv.mul(w, w), Q(1, 2))),
        'r': bv.add(bv.mul(q, r), bv.mul(r, w), bv.scale(bv.mul(w, w), Q(1, 2))),
        'w': bv.add(bv.mul(w, r), bv.scale(bv.mul(r, r), Q(1, 2))),
        'b': bv.add(bv.mul(q, b), bv.mul(r, b)),
        'c': bv.mul(q, c), 'h': bv.mul(q, h),
        'psi': bv.mul(bv.add(q, w), V['psi']),
        'barpsi': bv.mul(q, V['barpsi'])}
    generator = bv.add(*(bv.mul(V[ANTIS[name]], vector) for name, vector in vectors.items()))
    # Check the actual canonical vector, including Grassmann signs.
    for name, vector in vectors.items():
        assert bv.bracket(V[name], generator) == vector
    h_kernel = {}
    for a, b, value in (('q', 'q', Q(3, 5)), ('q', 'r', Q(2, 9)), ('q', 'w', Q(1, 6)),
                        ('q', 'b', Q(3, 7)), ('r', 'r', Q(2, 7)), ('r', 'w', Q(4, 11)),
                        ('w', 'w', Q(1, 3)), ('h', 'c', Q(3, 7)),
                        ('barpsi', 'psi', Q(2, 5))):
        kernel_entry(h_kernel, a, b, value)
    return s2, s0, generator, vectors, h_kernel


def full_graded_identity():
    s2, s0, generator, vectors, h_kernel = setup()
    rows, anomalies, contracted = [], [], []
    gamma_field_tests = 0
    # H preserves the FIELD BRST identity; it need not solve the EOM.
    for degree in (1, 2, 3):
        for names in combinations_with_replacement(FIELDS, degree):
            monomial = bv.prod(*(V[x] for x in names))
            if not monomial:
                continue
            assert project(s0(gamma(monomial, h_kernel))) == gamma(project(s0(monomial)), h_kernel)
            gamma_field_tests += 1
    interaction = bv.scale(bv.prod(V['w'], V['w'], V['w']), Q(5, 42))
    ell = bv.scale(V['w'], Q(7, 13))
    assert not s0(interaction) and not s0(ell)
    for sigma in (Q(5, 4), Q(7, 6)):
        w_kernel = {('w', 'w'): sigma}
        c_kernel = {key: -value for key, value in h_kernel.items()}
        c_kernel[('w', 'w')] += sigma
        contract = lambda f: gamma(f, c_kernel)
        for field, anti in bv.PAIRS:
            witness = bv.mul(V[anti], bv.mul(V['w'], V['w']))
            assert s0(gamma(witness, w_kernel)) == gamma(s0(witness), w_kernel)
        contracted_g = contract(generator)
        anomaly = bv.add(s0(contracted_g), bv.scale(contract(s0(generator)), -1))
        assert anomaly == project(anomaly) and anomaly
        k = {name: contract(vec) for name, vec in vectors.items()}
        # Independent variation of the antifield-zero quadratic density.
        body = project(s2)
        delta_body = bv.add(*(bv.mul(bv.diff(body, name, right=True), vectors[name]) for name in FIELDS))
        dk = bv.add(*(bv.mul(bv.diff(body, name, right=True), k[name]) for name in FIELDS))
        assert delta_body == project(s0(generator))
        rho = bv.add(contract(delta_body), bv.scale(dk, -1))
        assert anomaly == bv.scale(rho, -1)
        old_raw = project(contract(interaction))
        new_raw = contract(bv.add(project(interaction), bv.scale(delta_body, -1)))
        new_ell = bv.add(ell, rho)
        assert bv.add(new_raw, new_ell, dk) == bv.add(old_raw, ell)
        no_ghost_h = {key: value for key, value in c_kernel.items() if not set(key) & {'c', 'h'}}
        bad_anomaly = bv.add(s0(gamma(generator, no_ghost_h)),
                             bv.scale(gamma(s0(generator), no_ghost_h), -1))
        assert bad_anomaly != anomaly
        no_fermion_h = {key: value for key, value in c_kernel.items() if not set(key) & {'psi', 'barpsi'}}
        bad_fermion = bv.add(s0(gamma(generator, no_fermion_h)),
                             bv.scale(gamma(s0(generator), no_fermion_h), -1))
        assert bad_fermion != anomaly
        wrong_w = dict(w_kernel)
        wrong_w[('r', 'r')] = Q(1, 5)
        bad_state = bv.add(s0(gamma(generator, wrong_w)),
                           bv.scale(gamma(s0(generator), wrong_w), -1))
        assert bad_state
        rows.append(dict(sigma=str(sigma), contracted_generator=bv.display(contracted_g),
                         source_contact=bv.display(rho), first_anomaly=bv.display(anomaly),
                         source_identity_residual_terms=0,
                         omitted_ghost_error=bv.display(bv.add(bad_anomaly, bv.scale(anomaly, -1))),
                         omitted_fermion_error=bv.display(bv.add(bad_fermion, bv.scale(anomaly, -1))),
                         wrong_nonbisolution_state_error=bv.display(bad_state)))
        anomalies.append(anomaly)
        contracted.append(contracted_g)
    assert anomalies[0] == anomalies[1] and contracted[0] != contracted[1]
    return dict(cases=rows, field_BRST_covariance_monomials=gamma_field_tests,
                complete_nonminimal_and_fermion_sectors=True,
                same_anomaly_different_state_dependent_mean=True,
                arithmetic='Fraction', continuum_state_simulated=False)


def relative_path_signs():
    s2, s0, generator, vectors, h_kernel = setup()
    c_kernel = {key: -value for key, value in h_kernel.items()}
    c_kernel[('w', 'w')] += Q(5, 4)
    contract = lambda f: gamma(f, c_kernel)
    chi, tau = V['chi'], V['tau']
    delta = lambda f: bv.add(s0(f), bv.mul(chi, bv.diff(f, 'tau')))
    first_anomaly = lambda f: bv.add(delta(contract(f)), bv.scale(contract(delta(f)), -1))
    s3 = bv.scale(bv.prod(V['w'], V['w'], V['w']), Q(5, 42))
    ell = bv.scale(V['w'], Q(7, 13))
    a = first_anomaly(generator)
    rho = bv.scale(a, -1)
    extended = bv.add(s3, bv.scale(bv.mul(tau, s0(generator)), -1),
                      bv.scale(bv.mul(chi, generator), -1))
    assert not delta(extended)
    anomaly = first_anomaly(extended)
    expected = bv.add(bv.mul(tau, s0(a)), bv.mul(chi, a))
    assert anomaly == expected
    repair = bv.add(ell, bv.mul(tau, rho))
    assert not bv.add(anomaly, delta(repair))
    assert bv.substitute(repair, {'tau': {}}) == ell
    assert bv.substitute(repair, {'tau': bv.ONE}) == bv.add(ell, rho)
    assert not bv.diff(repair, 'chi')
    wrong_sign = bv.add(ell, bv.scale(bv.mul(tau, rho), -1))
    wrong_defect = bv.add(anomaly, delta(wrong_sign))
    assert wrong_defect
    return dict(classical_extended_identity_residual_terms=0,
                joint_one_loop_first_jet_residual_terms=0,
                endpoint_counterterm=bv.display(bv.add(ell, rho)),
                extended_anomaly=bv.display(anomaly),
                wrong_counterterm_sign_defect=bv.display(wrong_defect),
                no_chi_counterterm_at_first_jet=True,
                original_finite_linear_source='(7/13)*w', anchor_preserved=True,
                scope='Leading one-loop/one-field extended identity only; not higher-source matching.')


def run():
    return dict(round=790, fresh_test_groups=2,
                graded_generator_source=full_graded_identity(),
                anchored_relative_first_jet=relative_path_signs(),
                all_checks_passed=True,
                original_first_generator_anomaly_identified=True,
                original_fixed_gauge_coordinate_first_source_matching_proven=True,
                all_original_anomaly_coefficients_computed=False,
                all_finite_source_contact_coefficients_identified=False,
                auxiliary_gauge_endpoint_matching_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('anchored_source_contact_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
