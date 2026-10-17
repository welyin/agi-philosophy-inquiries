"""1251: finite repeated terminal statistics -> the 1248 weak certificate.

Only exact integer/Fraction certificates are computed.  No quantum trajectory,
detector, preparation, or repeated physical experiment is simulated.  The old
1248 hat-weight implementation is reused read-only; all evidence is hash bound.
Default is a read-only replay.  --write creates results.json exclusively.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FROZEN = {
    '1248/check.py': '69446a470fb0f426f8778bf6813121f936dcd51bdc5a3a6fc608cf88bca96a79',
    '1248/results.json': 'cf61f5d9efc1208541dca163071b106489750fa7cb142ca165d81561f3ac3943',
    '1248/proof.md': '649561632e352c3aa07bf1e75f57edb754b06453e3c53a4b39c564cf5b51f559',
    '1249/canonical.md': '6f7e4417dc96ee8062111cb4778e8e093257c7248916fedbce5db39eeb3cd989',
    '1249/quantum_bridge.md': '037e26ade64a67538b7a359fccb672ded5486795380cc3539a4db66e48a07d5f',
    '1249/quantum_bounds.py': '72db28046211fab6187078914c4034cf6afa70506f28e51312e5fcb01441acda',
    '1249/quantum_bounds_results.json': '77e125812877903dea9625a1fd0c0db571eb70829b28d6bcfa5a44e393ebcaa6',
}
ETA = F(1, 10**8)
TAU = F(5, 10**9)
SYSTEM = F(4, 10**9)
BIN = F(1, 10**10)
DELTA = F(1, 2**40)
K_EFFECT_BIAS = F(1, 10**4)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(x: F) -> dict:
    return {'rational': str(x), 'decimal': float(x)}


def frozen_inputs():
    for rel, expected in FROZEN.items():
        assert sha(ROOT / rel) == expected, f'frozen input changed: {rel}'
    q = json.loads((ROOT / '1249/quantum_bounds_results.json').read_text('utf-8'))
    w = json.loads((ROOT / '1248/results.json').read_text('utf-8'))
    assert q['status'] == w['status'] == 'PASS'
    assert q['science_groups'] == 1 and w['science_groups'] == 2
    return q, w


def quantum_statistics(q):
    old = q['group_constructive_finite_resource_bound']['ledger']
    cbeta = factorial(81) // factorial(40)**2
    assert cbeta * factorial(40)**2 == factorial(81)
    B = 324 * cbeta * 2**40
    Lcut = 56 * B
    M = 2**98 * (1 + Lcut)**2  # exactly old A_s, not a fitted new range
    Ls = 4 * (Lcut + 14)
    assert str(B) == old['B_decimal'] and str(Lcut) == old['Lcut_decimal']
    # ||D^r(chi exp Q)|| <= 27 r! (Lcut+1)^r, 0 <= r <= 20.
    # This is dominated by the original source class (A_s,L_s).
    assert M > 27 and Ls >= Lcut + 1
    assert all(27 * factorial(r) * (Lcut+1)**r <=
               M * factorial(r) * Ls**r for r in range(21))
    assert 48 + 1 <= 64
    EF = int(old['Fourier_transport_exponent_decimal'])
    EH = old['Fourier_H_exponent']
    N = int(old['integer_n_power_exponent_decimal'])
    assert 3*N-EF > 42 and 6*N-EH-EF > 42

    R = 26 * M*M * (2*10**8)**2
    assert R * TAU**2 / (2*M*M) == 13
    batches, tests = 2*129, 2*129*3
    # e > 1+1+1/2+1/6 = 8/3, so this rational upper bounds
    # 2*774 exp(-13), without evaluating a logarithm/exponential.
    statistics_failure_upper = 2*tests*F(3, 8)**13
    # Two additional old 1249 K-classification tasks, not a precision K read.
    ideal_failure_upper = statistics_failure_upper + 2*DELTA
    failure_upper = ideal_failure_upper + 2*K_EFFECT_BIAS
    assert failure_upper < F(1, 100)
    total_error = DELTA + BIN + SYSTEM + TAU
    assert total_error < ETA

    # Finite grid contains both endpoints exactly. Rounding down cannot leave
    # [-M,M]; clipping is retained as a documented defensive convention.
    G = 2*M*10**10+1
    def quantize(x):
        assert -M <= x <= M
        y = (x // BIN) * BIN
        return min(F(M), max(-F(M), y))
    test_points = [-F(M), -F(M)+BIN/3, -BIN/3, F(0), BIN/3,
                   F(M)-BIN/3, F(M)]
    for x in test_points:
        y = quantize(x)
        assert -M <= y <= M and 0 <= x-y < BIN
        assert (y/BIN).denominator == 1
    assert quantize(-F(M)) == -M and quantize(F(M)) == M
    outcome_bound = G**3
    bits = (outcome_bound-1).bit_length()
    assert 2**bits >= outcome_bound and 2**(bits-1) < outcome_bound
    statistical_runs = batches * R
    runs = statistical_runs + 2
    # A transition-region dictionary check: psi=96 gives chi=1/2 by the
    # symmetric original cutoff; other coordinates are in their core.
    # a=2, pi_Q=-8, pi_vartheta=.1, all other material coordinates zero.
    # This off-trajectory algebraic point is not a new physical preparation.
    cut, aval, piQ, piv = F(1, 2), F(2), F(-8), F(1, 10)
    hval = -piQ*piQ/(4*aval**3) + F(3, 20)*aval + F(72, 100) + piv/10
    assert hval == -F(97, 100)
    Jc = -cut*hval
    Mc = F(71, 100)+Jc
    rho_mat = F(3, 20)/aval**2 + (Mc+F(1, 100)+piv/10)/aval**3
    Hval = -piQ/(2*aval**3)
    dictionary_difference = cut*rho_mat-cut*Hval**2
    assert dictionary_difference == cut*(1-cut)*hval/aval**3 == -F(97, 3200)
    return {
        'status': 'PASS',
        'new_symbol': 'a_c=chi*exp(Q); original material rho_c,P_c retained',
        'new_symbol_derivative_class': {
            'orders_checked': [0, 20], 'own_amplitude_upper': '27',
            'own_derivative_scale': str(Lcut+1),
            'old_A_s_M': str(M), 'old_L_s': str(Ls),
            'declared_menu_size': 49, 'old_cap': 64,
            'old_resource_n_unchanged': True,
        },
        'time_menu': {'branches': 2, 'nodes_per_branch': 129,
                      'observables_per_node': 3, 'batches': batches,
                      'scalar_tests': tests, 'window': ['2', '9/4']},
        'finite_POVM': {
            'construction': 'same coherent POVM coarsened by three bounded symbols on a common grid',
            'bin_width': str(BIN), 'range': ['-M', 'M'],
            'rounding': 'floor then clip; M is an integer and both endpoints are grid points',
            'grid_points_per_coordinate': str(G),
            'outcome_count_upper': str(outcome_bound),
            'bits_per_raw_outcome_upper': bits,
            'same_effects_for_both_branches_and_all_nodes': True,
        },
        'statistics': {
            'R_formula': '26*M^2*(200000000)^2', 'R_decimal': str(R),
            'hoeffding_exponent': '13',
            'independence': 'independent repetitions within each batch; outputs within one trial may correlate',
            'joint_failure_rational_upper': exact(failure_upper),
            'ideal_effect_joint_failure_rational_upper': exact(ideal_failure_upper),
            'statistics_failure_rational_upper': exact(statistics_failure_upper),
            'two_old_K_classification_failure_upper': str(2*DELTA),
            'K_effect_probability_bias_contract_per_branch': str(K_EFFECT_BIAS),
            'two_K_effect_probability_bias_upper': str(2*K_EFFECT_BIAS),
            'joint_success_lower': '>99/100',
            'quantum_bias_delta': str(DELTA), 'bin_bias_upper': str(BIN),
            'physical_system_bias_contract': str(SYSTEM),
            'sampling_deviation': str(TAU),
            'sum_of_errors': exact(total_error),
            'accepted_menu_cap': str(ETA),
            'unused_menu_margin': exact(ETA-total_error),
        },
        'material_dictionary_transition_check': {
            'scope': 'off-trajectory symbol identity only; not a new preparation',
            'psi': '96', 'chi': str(cut), 'a': str(aval),
            'pi_Q': str(piQ), 'pi_vartheta': str(piv),
            'h': str(hval), 'Jc': str(Jc),
            'chi_rho_material_minus_chi_H_squared': exact(dictionary_difference),
            'nonzero': True,
        },
        'resources': {
            'statistical_runs_decimal': str(statistical_runs),
            'extra_old_K_classification_runs': 2,
            'runs_decimal': str(runs), 'runs_formula': '258*R+2',
            'raw_result_bits_upper_decimal': str(statistical_runs*bits+4),
            'raw_result_bits_scope': 'result storage only; schedules, devices and preparation records additionally accounted',
            'n_formula': 'n=m*2^N; m=ceil(max(1,hbar)), original 1249 N',
            'N_decimal': str(N),
            'volume_per_run': '1000*n^3',
            'classical_initial_material_energy_per_run': '1110*n^3',
            'volume_inventory_coefficient': str(1000*runs),
            'classical_initial_energy_inventory_coefficient': str(1110*runs),
            'common_n_cubed_factor': 'm^3*2^(3*N)',
            'inventory_not_dissipated_work': True,
            'preparation_ordering_tail_device_sources_included_by_separate_contract': True,
            'no_free_reset_or_simultaneous_noninteracting_universes_claimed': True,
        },
    }


def weak_baseline_certificate():
    spec = importlib.util.spec_from_file_location('frozen_1248_weights', ROOT/'1248/check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    wa, wb, old_weights = module.weights()
    assert len(wa) == len(wb) == 129
    assert sum(wa) == 0 and sum(wb) == 1 and min(wb) >= 0
    absolute_weight = sum(map(abs, wa))
    assert absolute_weight < 384
    h = F(1, 512)
    eta_b = F(123, 20)*ETA + 2*ETA**2
    E = 384*ETA + eta_b + h*h*F(384+20, 8)

    # NEW deletion: only m0=1/100. The K kinetic contribution and all true
    # dynamics are retained. It therefore needs no 1e-10 K estimator.
    m0 = F(1, 100)
    a_lower_during_read_error = F(8, 3)-ETA
    assert a_lower_during_read_error > 2
    derivative_r = m0 / 2**3  # |d [m0/(2 a^2)]/da| at a>=2
    assert derivative_r < F(1, 500)
    deletion_read_error = ETA/F(500)
    deletion_second = m0 * (F(11, 100)**2 + F(3, 8)*F(3, 25))
    assert deletion_second < F(1, 1000)
    assert 20 + deletion_second < 21
    Eminus = 384*ETA + eta_b + deletion_read_error + h*h*F(384+21, 8)
    Dlower = m0*F(1, 3)**2/2
    assert Dlower == F(1, 1800)
    margin = Dlower-E-Eminus
    assert margin > F(16, 100000)
    # Mean tuples remain in the same read domain used by the old product bound.
    assert ETA < F(8, 3)-2
    return {
        'status': 'PASS',
        'weights_reused_from': '1248/check.py:weights',
        'cells': 128, 'nodes': 129, 'h': str(h),
        'sum_a_weights': str(sum(wa)), 'sum_b_weights': str(sum(wb)),
        'sum_abs_a_weights': str(absolute_weight),
        'b_weights_nonnegative': True,
        'a_weights': old_weights['a_weights'], 'b_weights': old_weights['b_weights'],
        'old_same_history_derivative_bounds': {'B_a': 1, 'B_b': 20},
        'new_deletion': 'only m0=1/100; keep K^2/2, J, geometry and all actual dynamics',
        'new_estimator_subtraction': '1/(200*mean_a_j^2)',
        'eta_b': exact(eta_b), 'baseline_read_error_upper': exact(deletion_read_error),
        'baseline_second_derivative_upper': exact(deletion_second),
        'B_b_minus_baseline': 21,
        'full_budget_E': exact(E), 'baseline_omission_budget_Eminus0': exact(Eminus),
        'strict_deletion_gap_lower': exact(Dlower),
        'strict_disjoint_interval_margin': exact(margin),
        'confidence_inherited_from_group1': '>99/100',
        'no_K_precision_contract_added': True,
        'not_old_full_memory_deletion': True,
        'classical_weak_statistic_not_quantum_operator_product': True,
    }


def run():
    q, _ = frozen_inputs()
    return {
        'round': 1251, 'status': 'PASS', 'science_groups': 2,
        'arithmetic': 'exact integers and Fraction; no sampled physical data or quantum simulation',
        'group_joint_terminal_statistics': quantum_statistics(q),
        'group_weak_baseline_separation': weak_baseline_certificate(),
        'inputs_sha256': FROZEN,
        'author_code_sha256': sha(Path(__file__).resolve()),
        'scope': [
            'CT remains an unadopted conjecture; its effective branch is conditional.',
            'The same completed seven-pair quantum process is repeated for known preparations.',
            'Declared material source symbols are not replaced by H^2.',
            'Independent repetitions do not provide a monitored single-run history.',
            'Systematic preparation, timing, measurement, calibration and complete source matching remain explicit conditions.',
            'No apparatus, all-field model, actual precision or empirical GR is claimed.',
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'results.json'
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        print('WROTE_EXCLUSIVELY: 1251/results.json; 2 groups PASS')
    else:
        saved = json.loads(path.read_text('utf-8'))
        assert result == saved, 'saved result mismatch'
        print('READ_ONLY_MATCH: 1251; 2 groups PASS')


if __name__ == '__main__':
    main()
