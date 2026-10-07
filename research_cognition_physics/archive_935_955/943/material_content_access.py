"""943: content access for the actual 929 history code in the 942 material.

Reuse the old spectral lemma.  A symmetry-preserving bulk probe sees energy
blocks, not the logical multiplicity.  A new, explicit flavour matrix can
export content while still commuting with the free mass.  The finite probe
below is a witness, not a simulation of the full Einstein--matter action.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'material_content_access_results.json'
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def norm(a):
    return float(np.linalg.norm(a, 2))


def expm_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def trace_distance(a, b):
    return float(np.sum(abs(np.linalg.eigvalsh(a-b))) / 2)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run():
    spec = importlib.util.spec_from_file_location('protocol929', STAGE/'929/joint_protocol_selection.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    gates, hist = old.history(0., (2, 3))
    output = old.expected(0., (2, 3))
    inherited_error = norm(hist[-1]-output)
    # Actual event witnesses from both rounds, independent of the corrected fault labels.
    event00 = ((old.IND >> 1) & 1 == 0) & ((old.IND >> 7) & 1 == 0)
    effect00 = output[event00].conj().T @ output[event00]
    probabilities = np.real(np.diag(effect00))
    assert inherited_error < 1e-13
    assert norm(effect00-np.diag([.64, .04])) < 1e-13

    n = 30
    coupling = np.array([.5*np.sqrt((j+1)*(n-j)) for j in range(n)])
    hc = 15*np.eye(n+1)+np.diag(coupling, 1)+np.diag(coupling, -1)
    energies, v = np.linalg.eigh(hc)
    weights = abs(v[0])**2
    analytic = np.array([math.comb(n, j)/2**n for j in range(n+1)])
    spectral_error = float(np.max(abs(weights-analytic)))
    # Applying the full gate-weighted Hamiltonian to actual history columns.
    rng = np.random.default_rng(943)
    z = rng.normal(size=n+1)+1j*rng.normal(size=n+1)
    z /= np.linalg.norm(z)
    state = z[:, None, None]*hist
    applied = 15*state.copy()
    for j, gate in enumerate(gates):
        applied[j+1] += coupling[j]*old.act(state[j], gate)
        applied[j] += coupling[j]*old.act(state[j+1], gate, inverse=True)
    intertwiner_error = norm((applied-(hc@z)[:, None, None]*hist).reshape(-1, 2))
    assert spectral_error < 1e-13 and intertwiner_error < 1e-12

    # A genuinely interacting two-qubit neutral probe, with noncommuting probe terms.
    # The conditional source matrices are functions of the actual mass; this is a
    # finite symmetry-preserving witness, not an asserted gravitational Hamiltonian.
    hb = .31*np.kron(Z, I)+.23*np.kron(I, X)+.17*np.kron(Y, Z)
    b1, b2 = np.kron(X, I), np.kron(Z, Y)
    tau = np.array([1., 0., 0., 0.], dtype=complex)
    time, rest_mass = .73, 100.
    us, probe_h = [], []
    for e in energies:
        mass = rest_mass+e
        he = hb+.9*(mass/115)*b1+.7*((mass/115)**2)*b2
        probe_h.append(he)
        us.append(expm_h(he, time))
    us = np.asarray(us)
    vectors = us @ tau
    sigma = np.einsum('e,ea,eb->ab', weights, vectors, vectors.conj())
    # First finish the original process at pi, then let the invariant probe interact.
    # Dressed histories are used explicitly, not replaced by an identity data circuit.
    phases = np.exp(-1j*(np.pi+time)*energies)
    clock_probe = np.einsum('le,e,ea->la', v, v[0].conj()*phases, vectors)
    inputs = [np.array([1., 0.]), np.array([0., 1.]),
              np.array([1., 1.])/np.sqrt(2), np.array([1., 1j])/np.sqrt(2)]
    output_errors = []
    normalization_error = 0.
    for psi in inputs:
        physical = np.einsum('lwq,q->lw', hist, psi)
        amplitude = physical[:, :, None]*clock_probe[:, None, :]
        flat = amplitude.reshape(-1, 4)
        actual = flat.T @ flat.conj()
        output_errors.append(trace_distance(actual, sigma))
        normalization_error = max(normalization_error, float(abs(np.trace(actual)-1)))
    # Include an old passive reference entangled with the unknown logical input.
    bell = hist[:, :, :, None]*clock_probe[:, None, None, :]/np.sqrt(2)
    flat = bell.reshape(-1, 8)
    reference_probe = flat.T @ flat.conj()
    factor_error = trace_distance(reference_probe, np.kron(I/2, sigma))
    conditional_gap = trace_distance(np.outer(vectors[0], vectors[0].conj()),
                                     np.outer(vectors[-1], vectors[-1].conj()))
    unitary_error = max(norm(u.conj().T@u-np.eye(4)) for u in us)
    probe_internal_exchange = max(norm(h@hb-hb@h) for h in probe_h)
    assert max(output_errors) < 1e-12 and factor_error < 1e-12
    assert normalization_error < 1e-12 and unitary_error < 1e-12
    assert conditional_gap > .1 and probe_internal_exchange > .1

    # Compression of all central spectral effects, using several actual history
    # eigenvectors to check the analytic lemma in physical (gate-dressed) coordinates.
    code0 = hist[0]
    spectral_compression_error = 0.
    for k in (0, 5, 15, 24, 30):
        overlap = v[0, k].conj()*hist[0].conj().T@code0
        compressed = overlap.conj().T@overlap
        spectral_compression_error = max(spectral_compression_error, norm(compressed-weights[k]*I))
    assert spectral_compression_error < 1e-13

    # New portal matrix: B=Gamma(I_clock tensor Z_input tensor I_rest)Gamma^dagger.
    # On the actual two-column history sector, it is represented exactly by I tensor Z.
    h_restricted = np.kron(hc, I)
    portal = np.kron(np.eye(n+1), Z)
    allowed_flavour_rotation = np.kron(np.eye(n+1), expm_h(X, .37))
    commutes_mass = norm(h_restricted@portal-portal@h_restricted)
    breaks_commutant = norm(portal@allowed_flavour_rotation-allowed_flavour_rotation@portal)
    assert commutes_mass < 1e-13 and breaks_commutant > .5
    assert norm(portal@portal-np.eye(2*(n+1))) < 1e-13
    # In the fixed positive-energy momentum sector, the scalar bilinear projects
    # to (M/E) B because B commutes with M.  All masses use the original clock spectrum.
    momentum = 20.
    masses = rest_mass+energies
    kappa = float(np.sum(weights*masses/np.sqrt(masses**2+momentum**2)))
    scalar_source_effect = kappa*Z
    assert .9 < kappa < 1
    # Exact finite pointer witness for the extra access matrix; no field-delivery claim.
    gT = np.pi/8
    pointer_u = expm_h(np.kron(Z, Y), gT)
    plusx = (I+X)/2
    pointer_probabilities = []
    for psi in inputs[:2]:
        joint = (pointer_u@np.kron(psi, np.array([1., 0.]))).reshape(2, 2)
        ptr = joint.T@joint.conj()
        pointer_probabilities.append(float(np.real(np.trace(ptr@plusx))))
    expected_gap = float(np.sin(2*gT))
    pointer_gap = abs(pointer_probabilities[0]-pointer_probabilities[1])
    assert abs(pointer_gap-expected_gap) < 1e-13
    paths = [Path(__file__), STAGE/'929/joint_protocol_selection.py',
             STAGE/'929/joint_protocol_selection_results.json',
             STAGE/'research_note_935.md', STAGE/'research_note_939.md',
             STAGE/'942/dirac_material_embedding_results.json']
    return dict(round=943, date='2026-10-07', all_scientific_checks_passed=True,
        actual_929_history=dict(full_declared_flavour_dimension=4063232,
            computed_history_shape=list(hist.shape), bounded_fault_sector=[2, 3],
            original_output_isometry_error=inherited_error,
            actual_gate_weighted_intertwiner_error=intertwiner_error,
            event_00_probabilities=probabilities.tolist(), event_probability_gap=float(np.ptp(probabilities)),
            inherited_spectral_weights_error=spectral_error,
            sampled_spectral_compression_error=spectral_compression_error),
        invariant_bulk_probe=dict(probe_dimension=4, duration=time,
            input_output_trace_distances=output_errors, passive_reference_factorization_error=factor_error,
            normalization_error=normalization_error, unitary_error=unitary_error,
            extreme_mass_conditional_probe_trace_distance=conditional_gap,
            nonzero_probe_free_energy_commutator=probe_internal_exchange,
            actual_dressed_history_used=True, full_Einstein_evolution_simulated=False),
        new_flavour_portal=dict(matrix_norm=norm(portal), mass_commutator_error=commutes_mass,
            original_flavour_group_commutator_norm=breaks_commutant,
            momentum=momentum, scalar_source_compression_coefficient=kappa,
            scalar_source_eigenvalues=np.linalg.eigvalsh(scalar_source_effect).tolist(),
            finite_pointer_probabilities=pointer_probabilities, finite_pointer_probability_gap=pointer_gap,
            new_Yukawa_matrix_is_additional_physical_input=True,
            inherited_p_field_can_be_used=True, actual_p_to_R_field_readout_proved=False),
        analytic_scope=dict(old_935_spectral_lemma_reused=True, old_939_symmetric_access_lemma_reused=True,
            no_extra_chi_flavour_reference_and_fixed_common_motion_required=True,
            theorem_conditional_on_symmetry_preserving_dynamics_and_well_defined_process=True,
            candidate_access_failure_not_whole_program_failure=True,
            all_information_unobservable_or_all_GR_couplings_blind_claimed=False,
            all_cognitive_subjects_must_be_spatially_separate_claimed=False),
        candidate_refinement_stopped_after_decisive_access_test=True,
        full_goal_completed=False,
        source_hashes={str(q.relative_to(ROOT)):sha(q) for q in paths})


def compare(x, y):
    if isinstance(x, dict):
        assert x.keys() == y.keys()
        for k in x:
            compare(x[k], y[k])
    elif isinstance(x, list):
        assert len(x) == len(y)
        for a, b in zip(x, y):
            compare(a, b)
    elif isinstance(x, float):
        assert abs(x-y) < 1e-10, (x, y)
    else:
        assert x == y, (x, y)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf-8') as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write('\n')
    else:
        compare(result, json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}, ensure_ascii=False, indent=2))
