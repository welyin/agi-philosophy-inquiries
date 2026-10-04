"""763: same-algebra linear initial response, not a quantum-GR state construction.

The matrix-valued energy profiles below are DECLARED diagnostic sources.
They use original neutral CAR matrices on the original 753 background, but
are not a computed continuum stress tensor or a graph-to-continuum map.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as bg
import joint_loop_scale_transport as packet

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_operator_constraint_lift_results.json'
geo = bg.geo


def norm(a):
    return float(np.linalg.norm(a, 2))


def maximum(a):
    return float(np.max(np.abs(a)))


def constraint_lift():
    q, psi, tensor, base, _ = bg.completed(15, 1.)
    _, _, _, B, M, vp, vm, _ = packet.original_data()
    A = np.sum(tensor*tensor, axis=(-1, -2)) + q['pKp']
    potential = (5*q['C']*psi**4 - q['B'] + 7*A*psi**-8
                 + 6*q['Y']*psi**-4)
    positive_certificate = 4*q['C']*psi**5+8*A*psi**-7+8*q['Y']*psi**-3
    op = lambda v: -8*geo.laplace(v)+potential*v
    k2 = np.sum(geo.waves(15)**2, axis=-1)
    pre = lambda v: geo.ifft(geo.fft(v)/(8*k2+potential.mean()))
    # Two linearly independent, explicitly supplied INITIAL energy profiles.
    profiles = (np.ones_like(psi), q['f']['s'])
    solutions = []
    rows = []
    for profile in profiles:
        rhs = 2*psi**5*profile
        u, iterations = geo.cg(op, rhs, pre, tol=2e-13)
        solutions.append(u)
        rows.append(dict(iterations=iterations, linear_residual=maximum(op(u)-rhs)))
    a, b = solutions
    assert a.min() > 0
    ratio = b/a
    imin, imax = np.unravel_index(np.argmin(ratio), ratio.shape), np.unravel_index(np.argmax(ratio), ratio.shape)
    coefficients = np.array([[a[imin], b[imin]], [a[imax], b[imax]]])
    wedge = float(np.linalg.det(coefficients))
    U, V = [row[0]*B+row[1]*M for row in coefficients]
    comm = U@V-V@U
    assert maximum(comm-wedge*(B@M-M@B)) < 2e-15
    assert norm(comm) > 1e-7
    assert max(r['linear_residual'] for r in rows) < 2e-10
    assert maximum(op(psi)-positive_certificate) < 2e-8
    # Independent scalar nonlinear checks for two legal neutral input states.
    # These are state-dependent classical diagnostics, NOT a quantum channel.
    fd = []
    for state in (vp, vm):
        eB = float(np.vdot(state, B@state).real)
        eM = float(np.vdot(state, M@state).real)
        rho = profiles[0]*eB+profiles[1]*eM
        tangent = a*eB+b*eM
        errors = []
        for eta in (.02, .01, .005):
            updated = dict(q)
            updated['C'] = q['C']-2*eta*rho
            shifted, _, info = geo.solve_hamiltonian(updated, initial=psi)
            error = maximum((shifted-psi)/eta-tangent)
            errors.append(dict(amplitude=eta, derivative_error=error,
                               error_over_amplitude=error/eta,
                               Hamiltonian_residual=info['original_Hamiltonian_residual']))
        assert errors[-1]['derivative_error'] < errors[0]['derivative_error']*.3
        fd.append(dict(B_mean=eB, M_mean=eM, rows=errors))
    report = dict(N=15, color_amplitude=1., Fock_dimension=16,
                  background_residual=base['original_Hamiltonian_residual'],
                  positive_certificate_min=float(positive_certificate.min()),
                  positive_certificate_error=maximum(op(psi)-positive_certificate),
                  profile_2_range=[float(profiles[1].min()), float(profiles[1].max())],
                  scalar_solutions=rows, chosen_indices=[[int(i) for i in imin], [int(i) for i in imax]],
                  response_coefficients=coefficients.tolist(), wedge=wedge,
                  original_B_M_commutator_norm=norm(B@M-M@B),
                  response_commutator_norm=norm(comm),
                  commutator_identity_error=maximum(comm-wedge*(B@M-M@B)),
                  Gauss_and_momentum_unchanged_in_conformal_initial_sector=True,
                  nonlinear_mean_diagnostics=fd,
                  complete_continuum_quantum_source_computed=False)
    return report, (B, M, U, V, vp, vm)


def same_algebra(data):
    B, M, U, V, _, _ = data
    # A state that explicitly detects the antisymmetric response product.
    W = (U@V-V@U)/1j
    eigen, vec = np.linalg.eigh(W)
    v = vec[:, np.argmax(eigen)]
    rho = .9*np.outer(v, v.conj())+.1*np.eye(16)/16
    words = [np.eye(16), B, M, U, V, U@V, V@U]
    gram = np.array([[np.trace(rho@x.conj().T@y) for y in words] for x in words])
    herm_error = maximum(gram-gram.conj().T)
    min_eigen = float(np.linalg.eigvalsh((gram+gram.conj().T)/2).min())
    ordered_gap = np.trace(rho@(U@V-V@U))
    assert herm_error < 1e-13 and min_eigen > -1e-12
    assert abs(ordered_gap.imag) > 1e-7 and abs(ordered_gap.real) < 1e-13
    # The same linear identity survives all matrix element/state evaluations.
    same_moment_error = maximum(gram[5, 5]-np.trace(rho@V@U@U@V))
    assert same_moment_error < 1e-13
    return dict(state_min_eigenvalue=float(np.linalg.eigvalsh(rho).min()),
                ordered_moment_Gram_size=len(words), Gram_hermiticity_error=herm_error,
                Gram_min_eigenvalue=min_eigen, original_state_modified=False,
                antisymmetric_response_expectation_imag=float(ordered_gap.imag),
                product_identity_error=same_moment_error,
                products_are_evaluated_before_any_compression=True,
                independent_geometry_factor_constructed=False,
                scope='Positive same-algebra moments of declared initial response operators; '
                      'not intrinsic graviton covariance, physical geometry instruments or a full evolution.')


def independent_factor_witness(data):
    _, M, _, _, vp, vm = data
    basis = np.column_stack((vp, vm))
    logical_source = basis.conj().T@M@basis
    # Known information/disturbance calibration. This encoder is not the
    # original instrument and is not claimed to manufacture a gravitational state.
    V = np.zeros((4, 2), complex)
    V[0, 0] = 1
    V[3, 1] = 1
    plus = np.ones(2)/np.sqrt(2)
    initial = np.outer(plus, plus)
    joint = V@initial@V.conj().T
    reduced = np.trace(joint.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    distance = .5*float(np.abs(np.linalg.eigvalsh(reduced-initial)).sum())
    gap = float((logical_source[0, 0]-logical_source[1, 1]).real)
    assert maximum(basis.conj().T@basis-np.eye(2)) < 1e-14
    assert maximum(V.conj().T@V-np.eye(2)) < 1e-14
    assert maximum(logical_source-np.diag(np.diag(logical_source))) < 1e-14
    assert gap > .1 and abs(distance-.5) < 1e-14
    return dict(original_source_eigenvalue_gap=gap,
                source_basis_states_preserved=True,
                coherent_input_trace_distance=distance,
                matrix_embedding_isometry_error=maximum(V.conj().T@V-np.eye(2)),
                calibration_only_not_physical_geometry_or_original_instrument=True,
                scope='Illustrates the exact full-state assignment obstruction; restricted '
                      'commuting input states and relational encodings remain allowed.')


def run():
    lifted, data = constraint_lift()
    result = dict(round=763, tests_run=3, failures=0, errors=0,
                  operator_initial_lift=lifted, same_algebra=same_algebra(data),
                  independent_factor=independent_factor_witness(data))
    dependencies = ('research_note_731.md', 'research_note_753.md',
                    'research_note_754.md', 'research_note_756.md', 'research_note_762.md',
                    'joint_reference_constraint_strata.py', 'joint_loop_scale_transport.py',
                    'joint_gauss_einstein_initial_data.py', 'joint_irreducible_source_completion.py')
    result['dependency_hashes'] = {
        p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in dependencies}
    result['scope'] = ('Conditional operator-valued linear initial response in the SAME source algebra. '
                       'Numerics use declared energy profiles and original CAR coefficients on background753. '
                       'Not a graph-to-continuum map, full quantum constraints, canonical gravity quantization, '
                       'autonomous measurement or a replacement for intrinsic background fluctuations.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(payload)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('round', 'tests_run', 'failures', 'errors')}))
