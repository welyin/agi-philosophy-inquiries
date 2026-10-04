"""735: object checks for the local normalization bridge in note735.

These are complete-species finite-symbol checks, not computations of continuum
anomaly coefficients or of a renormalized stress tensor.  The existence and
scope of the continuum prescription require the analytic argument in the note.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_dynamic_continuum_reference as ref
import joint_relative_source_development as prior
import joint_retarded_reference_response as response

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_local_source_normalization_results.json'


def maxabs(x):
    return float(np.max(np.abs(x)))


def mass(phi):
    return ref.old.bdg(*ref.old.matter.mass_matrices(phi))


def basis():
    return np.array([np.sqrt(ref.old.matter.original.F(e)) * mass(e)
                     for e in np.eye(5)])


def jacobian(phi):
    F = ref.old.matter.original.F(phi)
    return np.eye(5) / np.sqrt(F) + np.outer(phi, phi) / (6 * F**1.5)


def jacobian_jet(phi, v):
    F = ref.old.matter.original.F(phi)
    return ((phi @ v) * np.eye(5) + np.outer(v, phi) + np.outer(phi, v)) / (6 * F**1.5) + (phi @ v) * np.outer(phi, phi) / (12 * F**2.5)


def inverse_mass_coordinates(x):
    return x * np.sqrt(2 / (1 + x @ x / 6))


def scalar_generator(phi, kind, index):
    result = np.zeros(5)
    X = phi[:2] + 1j * phi[2:4]
    if kind == 'weak':
        dX = 1j * ref.SIG[index] @ X / 2
        result[:4] = np.r_[dX.real, dX.imag]
    elif kind == 'circle':
        dX = 3j * X
        result[:4] = np.r_[dX.real, dX.imag]
    return result


def mass_covariant_jet_check():
    p = ref.point()
    phi, v = p['phi'], p['v']
    N = basis()
    F = ref.old.matter.original.F(phi)
    x = phi / np.sqrt(F)
    J = jacobian(phi)
    dx = J @ v
    ddx = jacobian_jet(phi, v) @ v
    M = np.einsum('a,aij->ij', x, N)
    dM = np.einsum('a,aij->ij', dx, N)
    ddM = np.einsum('a,aij->ij', ddx, N)
    h = 2e-4
    fd1 = (mass(phi + h*v) - mass(phi - h*v)) / (2*h)
    fd2 = (mass(phi + h*v) - 2*mass(phi) + mass(phi - h*v)) / h**2
    jet_error = max(maxabs(M - mass(phi)), maxabs(dM - fd1), maxabs(ddM - fd2))
    # The background connection and its actual time derivative do not commute
    # with the mass.  Compare the differentiated covariant jet independently.
    A = ref.gauge_h(p['a'], p['a0'])[0]
    dA = ref.gauge_h(p['da'], p['da0'])[0]
    A = ref.old.bdg(A, np.zeros_like(A))
    dA = ref.old.bdg(dA, np.zeros_like(dA))
    def comm(a, b):
        return a @ b - b @ a
    mixed = ddM - 1j * (comm(dA, M) + comm(A, dM))
    def covariant(s):
        MM, D = prior.mass_jet(phi+s*v, v)
        return D - 1j*comm(A+s*dA, MM)
    mixed_fd = (covariant(h) - covariant(-h)) / (2*h)
    mixed_error = maxabs(mixed - mixed_fd)
    equivariance_error = 0.
    for kind, index, H in prior.physical_generators():
        Q = ref.old.bdg(H, np.zeros_like(H))
        k = scalar_generator(phi, kind, index)
        dk = scalar_generator(v, kind, index)
        K = np.einsum('a,aij->ij', J @ k, N)
        dK = np.einsum('a,aij->ij', jacobian_jet(phi, v) @ k + J @ dk, N)
        equivariance_error = max(equivariance_error,
                                maxabs(K-1j*comm(Q, M)),
                                maxabs(dK-1j*comm(Q, dM)))
    assert jet_error < 3e-8 and mixed_error < 3e-8 and equivariance_error < 2e-12
    witness = maxabs(comm(A, M))
    assert witness > 1e-3 and np.linalg.norm(ddx) > 1e-8
    return dict(mass_dimension=64,original_generators=12,
                first_and_second_mass_jet_error=jet_error,
                connection_mass_mixed_jet_error=mixed_error,
                differentiated_spurion_equivariance_error=equivariance_error,
                noncommuting_connection_mass_witness=witness,
                second_mass_coordinate_jet_norm=float(np.linalg.norm(ddx)),
                coordinate_map_is_inherited_from_round623=True,
                no_commuting_mass_or_zero_Yukawa_reduction=True)


def actual_record_source_pullback_check():
    t = -.07
    P, dP, _ = response.flow(t)
    post = P + response.family.record_change(P)
    dpost = dP + response.family.record_change(dP)
    _, phi, _, _, _ = ref.collar(t, .08)
    w = response.entry.pulse(t)
    dphi = w * phi
    N = [prior.assemble([n, n]) for n in basis()]
    sx = np.array([prior.source(post, n) for n in N])
    dsx = np.array([prior.source(dpost, n) for n in N])
    J = jacobian(phi)
    dJ = jacobian_jet(phi, dphi)
    source = J.T @ sx
    state = J.T @ dsx
    contact = dJ.T @ sx
    independent = []
    for e in np.eye(5):
        _, D = prior.mass_jet(phi, e)
        independent.append(prior.source(post, prior.assemble([D, D])))
    source_error = maxabs(source - independent)
    h = 2e-5
    sources = []
    for gamma in (h, -h):
        Q, _, _ = response.flow(t, gamma)
        Q += response.family.record_change(Q)
        phig = np.exp(gamma*w) * phi
        sources.append(np.array([prior.source(Q, prior.assemble([prior.mass_jet(phig,e)[1]]*2)) for e in np.eye(5)]))
    fd = (sources[0] - sources[1]) / (2*h)
    derivative_error = maxabs(fd - state - contact)
    exchange_error = 0.
    for kind, index, _ in prior.physical_generators():
        kphi = scalar_generator(phi, kind, index)
        kx = scalar_generator(phi / np.sqrt(ref.old.matter.original.F(phi)), kind, index)
        exchange_error = max(exchange_error, abs(source @ kphi - sx @ kx))
    velocity = ref.point()['v']
    exchange_error = max(exchange_error, abs(source @ velocity - sx @ (J @ velocity)))
    # Exact eigenvalue identity exposes the non-uniform boundary; a numerical
    # check of that identity is not a scan or a claim about a physical singularity.
    F = ref.old.matter.original.F(phi)
    expected = np.array([F**-.5]*4 + [2*F**-1.5])
    eigen_error = maxabs(np.linalg.eigvalsh(J) - expected)
    assert source_error < 2e-12 and derivative_error < 3e-7
    assert exchange_error < 2e-12 and eigen_error < 2e-12
    assert maxabs(contact) > 1e-2 and maxabs(state) > 1e-6
    return dict(original_dimension=128,original_record_retained=True,
                scalar_sources_phi=source.tolist(),scalar_sources_x=sx.tolist(),
                source_pullback_error=source_error,
                history_state_response=state.tolist(),coordinate_contact=contact.tolist(),
                complete_derivative=(state+contact).tolist(),independent_difference=fd.tolist(),
                complete_derivative_error=derivative_error,
                omission_defect=maxabs(fd-state),
                gauge_and_time_exchange_pullback_error=float(exchange_error),
                jacobian_eigenvalue_identity_error=eigen_error,
                local_field_reparametrization_not_a_new_reference_preparation=True)


def run():
    names = ('mass_covariant_jet_check', 'actual_record_source_pullback_check')
    results = {name:globals()[name]() for name in names}
    dependencies = ('research_note_623.md','research_note_628.md','research_note_734.md',
                    'joint_dynamic_continuum_reference.py','joint_relative_source_development.py',
                    'joint_retarded_reference_response.py','joint_retarded_reference_response_results.json',
                    'round735_drafts/counterterm_scope_entry_results.json',
                    'round735_drafts/gauge_history_scope_probe_results.json')
    return dict(round=735,tests_run=2,failures=0,errors=0,checks=list(names),results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
                scope='Finite complete-species object checks for a conditional one-fermion-loop local normalization bridge: noncommuting mass and connection jets and the same actual record with full nonlinear scalar-source contacts. No numerical anomaly cancellation, continuum stress values, all-loop construction, or nonlinear semiclassical solution is claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as handle:
            handle.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
