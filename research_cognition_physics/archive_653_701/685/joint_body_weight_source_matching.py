"""685: original K bulk weight on demonstrably nonzero physical support.

The 615 one-site target has the complete S9 integral. Finite-regulator S9
convergence is proved in the note, not replaced here with fixed E samples.
This does not evaluate the full gauge/thermal integral or establish RP.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_subgroup_measure_source as old
import joint_local_source_lift as body

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_body_weight_source_matching_results.json'
spec = importlib.util.spec_from_file_location('entry685', HERE/'round685_drafts/body_weight_limit_probe.py')
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)


def scalar_body(theta, a, layers):
    """Exact one-site determinant and its theta derivative; no large K."""
    q = old.Q
    d = 4 + 4*a*np.cos(q*theta) + a*a
    dd = -4*a*q*np.sin(q*theta)
    h = a/np.sqrt(d)
    dh = -h*dd/(2*d)
    plus = layers*np.log1p(h)
    minus = layers*np.log1p(-h)
    w = 1/(1+np.exp(minus-plus))
    value = np.sum(2*(layers+1)*np.log(d) + 4*np.logaddexp(plus, minus))
    derivative = np.sum(2*(layers+1)*dd/d +
                        4*layers*dh*(w/(1+h)-(1-w)/(1-h)))
    return float(value), float(derivative)


def body_identity():
    saved = json.loads((HERE/'round685_drafts/body_weight_limit_probe_results.json').read_text('utf8'))
    reproduced = entry.run()
    assert reproduced == saved
    return dict(entry_reproduced=True, inherited_entry_not_counted_as_new_theorem=True,
                no_self_loop_trace_X_constant=True,
                two_by_two_trace_norm_limit=saved['analytic_trace_absolute_H_difference']/2,
                original_joint_limit_requires_aL_diverge_not_La_squared_vanish=True,
                original_direct_determinant_checks=saved['direct_formula_checks'])


def physical_support():
    g5 = np.kron(np.eye(16), old.spin.G5)
    direct = []
    rows = []
    assert int(np.sum(old.Q**2)) == 120
    for theta in (.07, .17, .31):
        _, _, _, H = old.frames(theta)
        X = g5@H
        assert np.linalg.norm(H@H-np.eye(64), 2) < 2e-14
        assert abs(np.trace(X).real-4*np.sum(np.cos(old.Q*theta))) < 1e-13
        # Actual inherited Weyl determinant and full S9 analytic average.
        W = old.physical(theta)
        M, dM = old.measure(theta)
        target = W*M
        assert target > 0 and abs(old.determinant(theta)-W) < 1e-13
        for layers in (1, 3):
            K, _, _ = body.blocks(X, g5, .2, layers)
            phase, value = np.linalg.slogdet(K)
            scalar, deriv = scalar_body(theta, .2, layers)
            generic = entry.determinant(H, g5, .2, layers)
            error = max(abs(value-scalar), abs(generic-scalar), abs(phase-1))
            assert error < 2e-11
            dt = 1e-5
            numeric_deriv = (scalar_body(theta+dt, .2, layers)[0]-scalar_body(theta-dt, .2, layers)[0])/(2*dt)
            assert abs(numeric_deriv-deriv) < 2e-7
            direct.append(dict(theta=theta, L=layers, determinant_residual=float(error),
                               derivative_residual=float(abs(numeric_deriv-deriv))))
        limit = float(2*np.sum(np.cos(old.Q*theta)-1))
        response_limit = float(-2*np.sum(old.Q*np.sin(old.Q*theta)))
        assert limit < 0
        sequence = []
        for a in (.1, .05, .025, .0125, .00625):
            L = int(np.ceil(4/a**2))
            value, derivative = scalar_body(theta, a, L)
            reference = scalar_body(0, a, L)[0]
            scaled = (value-reference)/(a*L)
            response = derivative/(a*L)
            sequence.append(dict(a=a, L=L, aL=a*L,
                scaled_log_body_ratio=scaled, scaled_body_response=response,
                weight_limit_error=abs(scaled-limit), response_limit_error=abs(response-response_limit)))
        assert sequence[-1]['weight_limit_error'] < .01*abs(limit)
        assert sequence[-1]['response_limit_error'] < .01*abs(response_limit)
        assert sequence[-1]['response_limit_error'] < sequence[0]['response_limit_error']
        assert sequence[-1]['weight_limit_error'] < sequence[0]['weight_limit_error']
        target_response = float(-np.sum(old.Q*np.tan(old.Q*theta/2))+dM/M)
        rows.append(dict(theta=theta, exact_full_S9_target_weight=target,
            exact_target_response=target_response, body_ratio_limit=limit,
            body_response_limit=response_limit, sequence=sequence))
    # Closed holonomy, not curvature: the one-site temporal loop is nontrivial.
    loop = np.diag(np.exp(1j*old.Q*.17))
    assert np.linalg.norm(loop-np.eye(16), 2) > .1
    assert all(np.linalg.norm(np.eye(16)@loop@np.eye(16)@loop.conj().T-np.eye(16), 2) < 1e-13 for _ in range(3))
    return dict(one_site_four_directions_are_inherited_inputs=True,
        all_S9_target_integrals_reused_from_615=True,
        finite_soft_S9_not_numerically_integrated=True,
        trace_abs_H_constant=64, trace_X_varies_due_to_wrapped_self_links=True,
        sum_charge_squares=120, scaled_bulk_log_ratio_second_derivative_at_zero=-240,
        nonzero_support_window='abs(theta)<pi/6', direct_checks=direct, rows=rows)


def run():
    deps = ('research_note_615.md', 'joint_subgroup_measure_source.py',
        'joint_subgroup_measure_source_results.json', 'research_note_677.md',
        'research_note_678.md', 'research_note_684.md', 'joint_local_source_lift.py',
        'round685_drafts/body_weight_limit_probe.py',
        'round685_drafts/body_weight_limit_probe_results.json',
        'round685_drafts/body_weight_limit_entry.md', 'round685_drafts/entry_checks.json')
    return dict(date='2026-10-02', round=685, tests_run=2, failures=0, errors=0,
        original_body_identity=body_identity(), original_nonzero_physical_support=physical_support(),
        dependency_hashes={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_conditional_all_source_weight_requires_matching=True,
            original_fully_S9_averaged_nonzero_support_used=True,
            finite_regulator_coefficient_nonzero_from_analytic_convergence=True,
            constant_or_curvature_only_subtraction_insufficient=True,
            all_counterterms_excluded=False, full_Haar_normalized_difference_proved=False,
            original_dynamic_Gauss_RP_decided=False,
            original_HF_continuum_or_quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=685, tests_run=2, all_checks_passed=True,
        physical_weights=[r['exact_full_S9_target_weight'] for r in result['original_nonzero_physical_support']['rows']])))
