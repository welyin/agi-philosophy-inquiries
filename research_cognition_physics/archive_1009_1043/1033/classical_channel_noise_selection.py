"""1033: conditional noise threshold for a two-qubit Markov channel.

This is a finite-dimensional calibration, not a gravitational experiment or an
autonomous implementation. --write creates results exclusively; default checks.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'classical_channel_noise_selection_results.json'
HISTORY = (
    'archive_935_955/research_note_945.md',
    'archive_956_989/research_note_989.md',
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_990_1008/1008/overall_operation_hypothesis_v2_2.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_1009_/1032/NEXT.md',
    'archive_1009_/1032/input_dependency_update_v0_21.md',
)
I = np.eye(2, dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
X = np.array([[0., 1.], [1., 0.]], complex)
Y = np.array([[0., -1j], [1j, 0.]], complex)
A, B = np.kron(Z, I), np.kron(I, Z)
AV, BV = np.diag(A).real, np.diag(B).real
PLUS = np.array([1., 1.])/math.sqrt(2)
MINUS = np.array([1., -1.])/math.sqrt(2)
PP = np.kron(PLUS, PLUS)
RHO = np.outer(PP, PP)
NULL = np.stack([np.kron(MINUS, PLUS), np.kron(PLUS, MINUS)], axis=1)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pt(rho):
    return rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def dissipator(c, rho):
    k = c.conj().T @ c
    return c @ rho @ c.conj().T - (k @ rho + rho @ k)/2


def generator(rho, g, ga, gb):
    h = g*A@B
    return -1j*(h@rho-rho@h)+ga*dissipator(A, rho)+gb*dissipator(B, rho)


def multiplier(g, ga, gb, t):
    energy = g*AV*BV
    return np.exp(t*(-1j*(energy[:, None]-energy[None, :])
                    -ga*(AV[:, None]-AV[None, :])**2/2
                    -gb*(BV[:, None]-BV[None, :])**2/2))


def feedback_generator(rho, c, f):
    h = (c.conj().T@f + f@c)/2
    return -1j*(h@rho-rho@h)+dissipator(c-1j*f, rho)


def one_way_multiplier(measured, target, kappa, force, dt):
    """Two actual local measurement outcomes and conditional remote unitaries."""
    assert 4*kappa*dt <= 1
    answer = np.zeros((4, 4), complex)
    for s in (-1, 1):
        diagonal = np.sqrt((1+2*s*math.sqrt(kappa*dt)*measured)/2)
        diagonal = diagonal*np.exp(-1j*s*force*math.sqrt(dt)*target)
        answer += np.outer(diagonal, diagonal.conj())
    return answer


def finite_feedback(g, ga, gb, t, steps):
    if g == 0:
        return multiplier(0, ga, gb, t)
    ka, kb = ga/2, g*g/(2*ga)
    excess_b = gb-g*g/ga
    assert excess_b >= -1e-14
    dt = t/steps
    fa = g/(2*math.sqrt(ka))
    fb = g/(2*math.sqrt(kb))
    step = one_way_multiplier(AV, BV, ka, fa, dt)
    step *= one_way_multiplier(BV, AV, kb, fb, dt)
    step *= multiplier(0, 0, max(0., excess_b), dt)
    return step**steps


def grid_cases():
    rows = []
    for g in (-1., -.25, 0., .25, 1.):
        for ga in (0., .25, 1., 4.):
            for gb in (0., .25, 1., 4.):
                block = NULL.conj().T @ pt(generator(RHO, g, ga, gb)) @ NULL
                expected = np.array([[ga, -1j*g], [1j*g, gb]])
                err = float(np.max(np.abs(block-expected)))
                assert err < 1e-13
                eig, vectors = np.linalg.eigh(expected)
                negative = ga*gb < g*g
                assert (eig[0] < -1e-13) == negative
                min_pt, min_cp = [], []
                for t in (.0001, .01, .2):
                    fm = multiplier(g, ga, gb, t)
                    min_pt.append(float(np.linalg.eigvalsh(pt(fm*RHO))[0]))
                    # Nonzero eigenvalues of the normalized Choi matrix are F/4.
                    min_cp.append(float(np.linalg.eigvalsh(fm/4)[0]))
                assert min(min_cp) > -1e-13
                if negative:
                    assert min_pt[0] < -1e-8
                else:
                    assert min(min_pt) > -1e-13
                rows.append(dict(g=g, gamma_A=ga, gamma_B=gb,
                    classical_threshold_satisfied=not negative,
                    tangent_smallest_eigenvalue=float(eig[0]),
                    tangent_residual=err, pt_min_eigenvalues=min_pt,
                    choi_min_eigenvalues=min_cp))
    return rows


def feedback_cases():
    rows = []
    for g, ga, gb in ((1., .5, 2.), (1., 1., 1.), (1., 2., .5),
                      (-1., 1., 2.), (.3, .4, .3), (0., 0., .5)):
        residual = 0.
        if g:
            ka, kb = ga/2, g*g/(2*ga)
            fa, fb = g/(2*math.sqrt(ka)), g/(2*math.sqrt(kb))
            for i in range(4):
                for j in range(4):
                    e = np.zeros((4, 4), complex); e[i, j] = 1
                    got = feedback_generator(e, math.sqrt(ka)*A, fa*B)
                    got += feedback_generator(e, math.sqrt(kb)*B, fb*A)
                    got += (gb-g*g/ga)*dissipator(B, e)
                    residual = max(residual, float(np.max(np.abs(got-generator(e, g, ga, gb)))))
        assert residual < 1e-13
        errors = []
        for n in (20, 80, 320, 1280):
            fm = finite_feedback(g, ga, gb, .2, n)
            err = float(np.max(np.abs(fm-multiplier(g, ga, gb, .2))))
            assert np.max(np.abs(np.diag(fm)-1)) < 3e-12
            assert np.linalg.eigvalsh(pt(fm*RHO))[0] > -1e-13
            errors.append(err)
        if g:
            assert all(errors[k+1] < .3*errors[k] for k in range(3))
        rows.append(dict(g=g, gamma_A=ga, gamma_B=gb,
                    generator_matrix_basis_residual=residual,
                    steps=[20, 80, 320, 1280],
                    max_multiplier_error=errors))
    return rows


def finite_witness():
    g, ga, gb, t = 1., .25, .25, .01
    _, v = np.linalg.eigh(np.array([[ga, -1j*g], [1j*g, gb]]))
    null_vector = NULL@v[:, 0]
    w = pt(np.outer(null_vector, null_vector.conj()))
    rho = multiplier(g, ga, gb, t)*RHO
    direct = float(np.trace(w@rho).real)
    coeffs, reconstructed = {}, np.zeros((4, 4), complex)
    for i, pa in enumerate((I, X, Y, Z)):
        for j, pb in enumerate((I, X, Y, Z)):
            op = np.kron(pa, pb)
            c = float(np.trace(w@op).real/4)
            if abs(c) > 1e-13:
                coeffs['IXYZ'[i]+'IXYZ'[j]] = c
                reconstructed += c*op
    assert np.max(np.abs(reconstructed-w)) < 1e-13
    # ||L||_{1->1} <= 3. Exp remainder <= (3t)^2/[2(1-3t)].
    analytic_upper = -F(3, 4)*F(1, 100) + F(3, 100)**2/(2*(1-F(3, 100)))
    assert direct <= float(analytic_upper) < -0.007
    epsilon = F(1, 1000)  # half trace distance; ||W||=1/2
    assert analytic_upper+epsilon < 0
    return dict(g=g, gamma_A=ga, gamma_B=gb, time=t,
        direct_witness=direct, exact_upper_bound=str(analytic_upper),
        upper_bound_float=float(analytic_upper),
        witness_operator_norm=float(np.max(np.abs(np.linalg.eigvalsh(w)))),
        local_pauli_coefficients=coeffs,
        illustrative_half_trace_error=str(epsilon),
        robust_upper_bound=str(analytic_upper+epsilon),
        actual_physical_error_certified=False)


def visibility_cases():
    rows = []
    for g, ga, gb in ((1., 1., 1.), (1., .5, 2.), (-1., 2., 1.),
                      (1., .25, .25), (0., .5, 0.)):
        t = .1
        p_a = np.kron(PLUS, np.array([1., 0.]))
        p_b = np.kron(np.array([1., 0.]), PLUS)
        ra = multiplier(g, ga, gb, t)*np.outer(p_a, p_a)
        rb = multiplier(g, ga, gb, t)*np.outer(p_b, p_b)
        va, vb = float(2*abs(ra[0, 2])), float(2*abs(rb[0, 1]))
        bound = math.exp(-4*abs(g)*t)
        assert abs(va-math.exp(-2*ga*t)) < 1e-14
        assert abs(vb-math.exp(-2*gb*t)) < 1e-14
        if ga*gb >= g*g:
            assert va*vb <= bound+1e-14
        rows.append(dict(g=g, gamma_A=ga, gamma_B=gb, time=t,
            phase_lift=4*g*t, visibility_A=va, visibility_B=vb,
            product=va*vb, classical_upper_bound=bound,
            threshold_satisfied=ga*gb >= g*g))
    # A finite, conservative error envelope; all numbers are illustrative.
    lower_a = lower_b = F(99, 100)-F(1, 1000)
    phi_lower = F(4, 100)-F(1, 1000)
    exp_upper = 1/(1+phi_lower)  # exp(-x) <= 1/(1+x), x >= 0
    gap = lower_a*lower_b-exp_upper
    assert gap > 0
    return dict(cases=rows, finite_error_example=dict(
        observed_visibility_each='99/100', visibility_error_each='1/1000',
        observed_phase_lift='1/25', phase_error='1/1000',
        exact_positive_gap_lower_bound=str(gap),
        physical_measurement_data=False))


def run():
    grid = grid_cases()
    return dict(round=1033, date='2026-10-08',
        new_calibration_groups=1, cumulative_research_groups=3810,
        new_adopted_cognitive_axioms=0, goal_completed=False,
        external_independent_agent_review_completed=False,
        effective_Markov_contract_assumed=True,
        universal_classical_gravity_no_go=False,
        autonomous_internal_implementation_certified=False,
        physical_gravity_entanglement_observed=False,
        code_sha256=sha(Path(__file__)),
        historical_source_sha256={name: sha(BASE/name) for name in HISTORY},
        grid=grid, feedback=feedback_cases(), finite_witness=finite_witness(),
        visibility=visibility_cases())


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, abs_tol=1e-11, rel_tol=1e-9), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true')
    args = parser.parse_args(); result = run()
    if args.write:
        with OUT.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False); f.write('\n')
    else:
        compare(result, json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1033, grid_cases=len(result['grid']),
        feedback_cases=len(result['feedback']), finite_witness=result['finite_witness'],
        finite_error_example=result['visibility']['finite_error_example']), ensure_ascii=False))
