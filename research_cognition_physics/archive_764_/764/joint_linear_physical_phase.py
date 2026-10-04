"""764: derivative-core rank and symplectic/CCR completion diagnostics.

These matrices calibrate analytic identities on the universal diagonal
derivative core and an explicitly chosen finite symplectic section. They are NOT a
discretization of all variable-coefficient constraints on background753.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as bg

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_linear_physical_phase_results.json'
NQ, NC = 47, 16


def maximum(a):
    return float(np.max(np.abs(a)))


def symplectic(n):
    return np.block([[np.zeros((n, n)), np.eye(n)],
                     [-np.eye(n), np.zeros((n, n))]])


J = symplectic(NQ)


def symmetric_basis():
    basis = []
    for i in range(3):
        a = np.zeros((3, 3))
        a[i, i] = 1
        basis.append(a)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        a = np.zeros((3, 3))
        a[i, j] = a[j, i] = 1/np.sqrt(2)
        basis.append(a)
    return np.array(basis)


T = symmetric_basis()


def principal(k):
    """Normalized real diagonal derivative core, row Fourier phases removed.

    q=(metric[6], twelve spatial connections[36], scalars[5]).
    Rows: Hamiltonian, three momenta, twelve Gauss constraints.
    Background-dependent derivative/mixing terms are NOT sampled here.
    This is not the entire principal symbol of the coupled PDE.
    """
    k = np.asarray(k, float)
    k /= np.linalg.norm(k)
    L = np.zeros((NC, 2*NQ))
    L[0, :6] = np.einsum('aij,ij->a', T, np.outer(k, k)-np.eye(3))
    for i in range(3):
        L[1+i, NQ:NQ+6] = -2*np.einsum('aj,j->a', T[:, i, :], k)
    for a in range(12):
        L[4+a, NQ+6+3*a:NQ+9+3*a] = k
    return L, k


def tt_vectors(k):
    axis = np.eye(3)[np.argmin(abs(k))]
    e = axis-k*np.dot(k, axis)
    e /= np.linalg.norm(e)
    f = np.cross(k, e)
    plus = (np.outer(e, e)-np.outer(f, f))/np.sqrt(2)
    cross = (np.outer(e, f)+np.outer(f, e))/np.sqrt(2)
    U = np.zeros((2*NQ, 4))
    for i, tensor in enumerate((plus, cross)):
        components = np.einsum('aij,ij->a', T, tensor)
        U[:6, i] = components
        U[NQ:NQ+6, 2+i] = components
    return U


def principal_rank_check():
    rows = []
    for wave in ((1, 0, 0), (0, 0, 2), (1, 2, -3), (.13, -.71, .43)):
        L, k = principal(wave)
        tt = tt_vectors(k)
        rank = int(np.linalg.matrix_rank(L))
        assert rank == 16 and maximum(L@J@L.T) < 1e-14
        assert maximum(L@tt) < 1e-14
        assert maximum(tt.T@J@tt-symplectic(2)) < 1e-14
        rows.append(dict(wave=list(wave), constraint_rank=rank,
                         first_class_principal_error=maximum(L@J@L.T),
                         TT_constraint_error=maximum(L@tt),
                         TT_canonical_pair_error=maximum(tt.T@J@tt-symplectic(2))))
    q, psi, _, info, _ = bg.completed(12, 1.)
    return dict(rows=rows, canonical_configuration_components=NQ,
                physical_phase_principal_dimension=2*NQ-2*NC,
                physical_configuration_principal_count=NQ-NC,
                matter_only_configuration_principal_count=3*12+5-12,
                gravity_configuration_principal_difference=2,
                original_background_positive_F_min=float(q['F'].min()),
                original_background_psi_range=[float(psi.min()), float(psi.max())],
                original_background_constraint_residual=info['original_Hamiltonian_residual'],
                actual_full_constraint_matrix_discretized=False,
                omitted_background_dependent_derivative_terms=True,
                scope='Universal diagonal derivative-core rank for nonzero covectors. '
                      'Not the full coupled principal symbol, a finite wavelength solution, '
                      'or a simulation of the variable background.')


def projection_check():
    L, _ = principal((1, 2, -3))
    R0 = L.T@np.linalg.inv(L@L.T)
    rng = np.random.default_rng(764)
    R = R0+(np.eye(2*NQ)-R0@L)@(.03*rng.normal(size=(2*NQ, NC)))
    K = J@L.T
    skew = R.T@J@R
    S = R+.5*K@skew
    P = np.eye(2*NQ)-S@L+K@S.T@J
    naive = np.eye(2*NQ)-R@L
    errors = dict(original_right_inverse=maximum(L@R-np.eye(NC)),
                  corrected_right_inverse=maximum(L@S-np.eye(NC)),
                  isotropic_complement=maximum(S.T@J@S),
                  constraints=maximum(L@P),
                  gauge_section=maximum(S.T@J@P),
                  gauge_removed=maximum(P@K),
                  complement_removed=maximum(P@S),
                  projector=maximum(P@P-P),
                  symplectic_adjoint=maximum(P.T@J-J@P))
    assert max(errors.values()) < 5e-13
    assert maximum(skew) > .01
    assert int(np.linalg.matrix_rank(P, tol=1e-10)) == 62
    assert int(np.linalg.matrix_rank(naive, tol=1e-10)) == 78
    assert maximum(naive@K-K) < 1e-13
    # Arbitrary in-constraint data differ from their section by a pure gauge.
    v = naive@rng.normal(size=2*NQ)
    gauge_parameter = S.T@J@v
    orbit_error = maximum(P@v-v-K@gauge_parameter)
    assert orbit_error < 1e-13
    report = dict(errors=errors,
                  original_complement_skew_max=maximum(skew),
                  reduced_phase_rank=int(np.linalg.matrix_rank(P, tol=1e-10)),
                  naive_constraint_projection_rank=int(np.linalg.matrix_rank(naive, tol=1e-10)),
                  naive_projection_keeps_all_gauge_directions=True,
                  gauge_orbit_representative_error=orbit_error,
                  right_inverse_is_a_declared_finite_calibration=True)
    return report, (L, P)


def gaussian_ccr_check(data):
    L, P = data
    # Same position/momentum variance relation as a Gaussian packet,
    # extended here to ALL kinematical field coordinates as an explicit choice.
    c = np.linspace(.4, 1.7, NQ)
    V0 = np.diag(np.r_[c, .25/c])
    V = P@V0@P.T
    physical_J = P@J@P.T
    physical_quantum_form = V+.5j*physical_J
    wrong_quantum_form = V+.5j*J
    eig = np.linalg.eigvalsh(physical_quantum_form)
    wrong = np.linalg.eigvalsh(wrong_quantum_form)
    assert eig.min() > -1e-12
    assert wrong.min() < -.1
    constraint_variance = maximum(L@V@L.T)
    assert constraint_variance < 1e-13
    # Finite Weyl positivity check, including the CCR phase.
    rng = np.random.default_rng(765)
    labels = .1*rng.normal(size=(9, 2*NQ))
    gram = np.empty((len(labels), len(labels)), complex)
    for i, x in enumerate(labels):
        for j, y in enumerate(labels):
            d = y-x
            gram[i, j] = np.exp(.5j*x@physical_J@y-.5*d@V@d)
    gram_min = float(np.linalg.eigvalsh(gram).min())
    assert maximum(gram-gram.conj().T) < 1e-13 and gram_min > -1e-12
    return dict(physical_uncertainty_form_min=float(eig.min()),
                wrong_unreduced_CCR_uncertainty_form_min=float(wrong.min()),
                constraint_variance_residual=constraint_variance,
                physical_commutator_rank=int(np.linalg.matrix_rank(physical_J, tol=1e-10)),
                Weyl_Gram_min=gram_min, Weyl_Gram_size=len(labels),
                auxiliary_covariance_is_new_preparation_input=True,
                continuum_Hadamard_state_claimed=False,
                full_nonlinear_constraints_claimed=False)


def run():
    rank = principal_rank_check()
    reduced, data = projection_check()
    result = dict(round=764, tests_run=3, failures=0, errors=0,
                  principal_rank=rank, symplectic_section=reduced,
                  positive_state_and_CCR=gaussian_ccr_check(data))
    dependencies = ('research_note_568.md', 'research_note_573.md',
                    'research_note_574.md', 'research_note_590.md',
                    'research_note_731.md', 'research_note_753.md',
                    'research_note_754.md', 'research_note_758.md',
                    'research_note_762.md', 'research_note_763.md',
                    'joint_reference_constraint_strata.py')
    result['dependency_hashes'] = {
        name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dependencies}
    result['scope'] = ('Full classical linear constraint reduction and a positive regular CCR-state '
                       'construction are analytic, conditional on the inherited smooth canonical model '
                       'and its right inverse. Numerical checks are derivative-core/finite-symplectic '
                       'calibrations, NOT all coupled constraint modes, Hadamard renormalization, original '
                       'instrument realization, a Q-to-continuum equivalence, or quantum GR completion.')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--write-results', action='store_true')
    args = p.parse_args()
    result = run()
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(payload)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('round', 'tests_run', 'failures', 'errors')}))
