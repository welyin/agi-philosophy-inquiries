"""737: original local complement of the fermionic response null direction.

Reuses583's frame geometry,623/735's mass coordinates and632's matched vacuum.
Checks local object maps and derivative order, not a full constrained solution.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_local_source_normalization as massmap
import joint_background_contact_matching as matched
import joint_continuum_source_spectrum as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_response_null_completion_results.json'
original = matched.matter.original


def maxabs(a):
    return float(np.max(np.abs(a)))


def local_complement_check():
    _, _, _, diagnostic, _ = old.original_data()
    vacuum = np.array([0., np.sqrt(matched.U0[0]), 0., 0., np.sqrt(matched.U0[1])])
    N = massmap.basis()
    rows = []
    for label, phi in (('old630_diagnostic_not_stationary', diagnostic),
                       ('old632_matched_stationary_vacuum', vacuum)):
        F = float(original.F(phi))
        x = phi / np.sqrt(F)
        J = massmap.jacobian(phi)
        B = np.linalg.inv(J)
        # Pull back the ORIGINAL target metric, independently of the closed form.
        G = B.T @ original.metric(phi) @ B
        G_formula = np.eye(5) - np.outer(x, x) / (6 + x @ x)
        A = float(6 - x @ G @ x)
        A_formula = 6 * F / original.M
        H_original = np.zeros((6, 6))
        H_original[:5, :5] = G
        H_original[-1, -1] = -6
        transform = np.eye(6)
        transform[:5, -1] = -x
        H_active = transform.T @ H_original @ transform
        schur = H_active[:5, :5] - np.outer(H_active[:5, -1], H_active[-1, :5]) / H_active[-1, -1]
        null = np.r_[-x, 1.]
        null_mass = -np.einsum('a,aij->ij', x, N) + massmap.mass(phi)
        error = max(maxabs(G-G_formula), abs(A-A_formula),
                    maxabs(schur-np.eye(5)), abs(null@H_original@null+A), maxabs(null_mass))
        assert A > 0 and error < 2e-14
        rows.append(dict(background=label, F=F, mass_coordinate_norm=float(np.linalg.norm(x)),
                         complementary_principal_magnitude=A,
                         original_metric_to_active_block_error=error,
                         active_schur_eigenvalues=np.linalg.eigvalsh(schur).tolist(),
                         compensated_mass_error=maxabs(null_mass)))

    # The actual matched stationary background has TWO coupled neutral radials.
    # A common scaling direction need not be an invariant dynamical subspace.
    F = float(original.F(vacuum))
    K = original.metric(vacuum)
    v = vacuum / np.sqrt(vacuum @ K @ vacuum)
    _, _, W_hessian = matched.W_jets(matched.U0)
    D = np.diag(2*np.sqrt(matched.U0))
    closure = []
    for epsilon in (0., 1.):
        H = np.zeros((5, 5))
        H[np.ix_([1, 4], [1, 4])] = D @ (matched.L/2 + epsilon*W_hessian) @ D / F**2
        acceleration = np.linalg.solve(K, H @ v)
        eigenvalue = float(v @ H @ v)
        transverse = acceleration - eigenvalue*v
        defect = float(np.sqrt(transverse @ K @ transverse))
        assert defect > .001
        closure.append(dict(loop_multiplier=epsilon,
                            common_scaling_K_unit_norm=float(v@K@v),
                            radial_rayleigh_value=eigenvalue,
                            transverse_acceleration_K_norm=defect,
                            transverse_components=transverse.tolist()))
    # Only the rank of the finite mass map is certified here, not its full
    # continuum mixed-mass response or its physical gauge-reduced rank.
    gram = np.einsum('aij,bij->ab', N.conj(), N).real
    assert np.min(np.linalg.eigvalsh(gram)) > 0
    q = old.QSTAR
    b = float(1/np.sqrt(6)/np.tanh(q/np.sqrt(6)))
    radial_A = 6 - 1/b**2
    assert abs(radial_A-rows[0]['complementary_principal_magnitude']) < 3e-14
    return dict(original_backgrounds=rows,matched_vacuum_radial_nonclosure=closure,
                finite_mass_map_Gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                original630_b=b,original630_null_coefficient=radial_A,
                old583_Schur_tool_reused=True,
                no_claim_that_the_seven_source_block_is_a_closed_model=True)


def curvature_order_check():
    # Independent nonlinear evaluation of sqrt(g) R^2 on g=e^(2 sigma) delta,
    # with sigma depending on one periodic coordinate. No boundary term is
    # discarded in R^2; the Ricci scalar itself is evaluated explicitly.
    t = np.arange(2048)*2*np.pi/2048
    n = 3
    rows = []
    for eta in (1e-3, 5e-4):
        actions = []
        for sign in (-1., 1.):
            sigma = sign*eta*np.cos(n*t)
            ds = -sign*eta*n*np.sin(n*t)
            dds = -sign*eta*n*n*np.cos(n*t)
            R = -6*np.exp(-2*sigma)*(dds+ds*ds)
            actions.append(float(2*np.pi*np.mean(np.exp(4*sigma)*R*R)))
        observed = sum(actions)/eta**2
        expected = 72*np.pi*n**4
        exact_finite_eta = expected + 54*np.pi*n**4*eta**2
        assert abs(observed-exact_finite_eta) < 2e-11
        rows.append(dict(amplitude=eta,nonlinear_second_variation=observed,
                         predicted_quadratic_hessian=expected,
                         finite_amplitude_remainder=observed-expected,
                         exact_nonlinear_identity_error=abs(observed-exact_finite_eta)))
    assert abs(rows[0]['finite_amplitude_remainder']/rows[1]['finite_amplitude_remainder']-4) < 2e-6
    # Analytic sequence, not a numerical proof of unboundedness:
    # sigma_n=cos(nt)/n^2 has bounded C^2 norm but fourth derivative norm n^2.
    witness = [dict(n=n,second_derivative_sup=1.,fourth_derivative_sup=float(n*n))
               for n in (2, 4, 8, 16)]
    # Original leading null coefficient and the exact two alternative causal
    # eliminations, calibrated on prescribed smooth functions, NOT a solver.
    b = float(1/np.sqrt(6)/np.tanh(old.QSTAR/np.sqrt(6)))
    c = 1/b**2
    A = 6-c
    u = np.sin(t)
    v2 = np.cos(2*t)
    g2 = -c*u-A*v2
    recovered2 = -(g2+c*u)/A
    alpha = .02  # arbitrary nonzero local coefficient fixture, not a prediction
    r = 72*alpha
    sigma2 = -4*np.cos(2*t)
    v4 = 16*np.cos(2*t)
    g4 = r*v4-c*u-A*sigma2
    recovered4 = (g4+c*u+A*sigma2)/r
    assert max(maxabs(recovered2-v2),maxabs(recovered4-v4)) < 2e-14
    return dict(nonlinear_R_squared_checks=rows,
                null_channel_derivative_order_witness=witness,
                R_squared_Hessian_coefficient_per_alpha=72,
                local_elimination_identity_errors=[maxabs(recovered2-v2),maxabs(recovered4-v4)],
                nonzero_alpha_fixture=alpha,
                zero_alpha_branch_requires_nonzero_effective_second_order_coefficient=True,
                no_uniform_alpha_to_zero_or_physical_branch_claim=True,
                lapse_and_shift_constraints_not_solved_by_these_checks=True)


def run():
    names = ('local_complement_check', 'curvature_order_check')
    results = {name: globals()[name]() for name in names}
    deps = ('research_note_570.md','research_note_583.md','research_note_601.md',
            'research_note_623.md','research_note_630.md','research_note_631.md',
            'research_note_632.md','research_note_735.md','research_note_736.md',
            'joint_local_source_normalization.py','joint_background_contact_matching.py',
            'joint_continuum_source_spectrum.py','joint_curved_quantum_source.py')
    return dict(round=737,tests_run=2,failures=0,errors=0,checks=list(names),results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original local scalar/conformal complement mapped to the fermionic null direction; matched-vacuum radial nonclosure and the R-squared derivative-order distinction. Conditional short-time local-complement criteria, not an unconstrained seven-source model or nonlinear semiclassical spacetime.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
