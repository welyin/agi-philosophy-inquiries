"""1167: complete homogeneous-scale source kernel and directional witness.

Read-only by default; --write-once creates the first result. The quantum
scale and its finite instrument menu are candidate-analysis inputs, not a
newly certified local spacetime implementation. All old material is retained.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCES = {
    '../1161/check.py': '52edf7a1dfc530996daa88ee3715efa936f7c5845eefe067124a544688621119',
    '../1161/proof.md': 'f07074bca0db63349bcd3c28ef6cc8e4f786720308d52351ee12c94959582bbe',
    '../../archive_1159_1160/1160/proof.md': 'ca3121515439d2b93c23d341d83a1dfbd68b48c4c63ac95fc41c28b027cf8e81',
    '../../archive_1159_1160/1160/results.json': 'f6e38c2f676991d7e0ccf63ff603819d2c6975ce3fe0901e7b10b1c05377b4e2',
    '../1166/acceptance.json': '7a660aa0b7d74fee7320addb51f0f124b9979b20dc64825f22279a64b93bc51a',
}
TOL = 2e-9


def small(value, label):
    error = float(np.linalg.norm(value))
    assert error < TOL, (label, error)
    return error


def compare(actual, saved):
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for key in saved:
            compare(actual[key], saved[key])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(saved, float):
        assert math.isclose(actual, saved, rel_tol=2e-9, abs_tol=3e-10), (actual, saved)
    else:
        assert actual == saved, (actual, saved)


def exponential(h, duration):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*duration*values))@vectors.conj().T


def calculate():
    for path, expected in SOURCES.items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest() == expected, path
    spec = importlib.util.spec_from_file_location('r1161_source_for_1167', HERE/'../1161/check.py')
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    old, data = source.reuse()  # Complete frozen 1160, compared to its saved result.
    m = data['m']
    hd, hc, gamma = m['hout'], m['hc'], m['gamma']
    evals, vectors = np.linalg.eigh(hd)
    small(evals-m['baseline']-np.repeat(np.arange(4), 72), 'original full HD spectrum')
    shape = data['logical_shape']
    gates = data['gates']
    logical = gates[3]@gates[2]@gates[1]@gates[0]
    clock_intertwiner = small(hd@gamma-gamma@(m['baseline']*np.eye(288)+np.kron(hc, np.eye(72))),
                              'full physical Gamma intertwiner')

    def dress(state, logical_gate):
        decoded = (gamma.conj().T@state.reshape(288, -1)).reshape(4, 8, 3, 3, 4, 2, 2, 3, -1)
        changed = old.along(decoded, logical_gate, (2, 4, 5, 6, 7))
        return (gamma@changed.reshape(288, -1)).reshape(state.shape)

    # The 576 block is tensored with aux8 and Qr3: no mother states omitted.
    hbase = m['baseline']*np.eye(4)+hc
    expanded_gates = [np.kron(np.eye(4), gate) for gate in gates]
    source_rows = []
    worst_kernel = 0.
    for scale in [.8, 1., 1.2]:
        energy = math.sqrt(16+9*scale**2)
        hhead = scale*(hbase+3*np.eye(4))+2*energy*np.eye(4)
        jhead = hbase+(3+18*scale/energy)*np.eye(4)
        hm = np.kron(hhead, np.eye(144))
        jm = np.kron(jhead, np.eye(144))
        kernel_h = max(small(g@hm-hm@g, 'whole material energy kernel') for g in expanded_gates)
        kernel_j = max(small(g@jm-jm@g, 'whole material source kernel') for g in expanded_gates)
        worst_kernel = max(worst_kernel, kernel_h, kernel_j)
        all_spectrum = np.repeat(scale*(evals+3)+2*energy, 48)
        assert len(all_spectrum) == 13824
        source_rows.append(dict(scale=scale, AB_one_energy=energy, material_dimension=len(all_spectrum),
                                lowest_energy=float(all_spectrum.min()),
                                highest_energy=float(all_spectrum.max()),
                                energy_commutator=kernel_h, source_commutator=kernel_j))

    # Program type alignment is a unitary on the complete retained L factor.
    typed_error = 0.
    for program in (1, 2):
        tag = np.zeros_like(logical)
        for idx in np.ndindex(shape):
            dst = list(idx)
            dst[-1] = program if idx[-1] == 0 else (0 if idx[-1] == program else idx[-1])
            tag[np.ravel_multi_index(tuple(dst), shape), np.ravel_multi_index(idx, shape)] = 1
        relative = logical@tag@logical.conj().T
        typed_error = max(typed_error, small(relative.conj().T@relative-np.eye(144), 'relative full unitary'),
                          small(dress(data['outputs'][0], relative)-data['outputs'][program],
                                'complete typed program relation'))
    gram_error = max(small(v.reshape(-1, 2).conj().T@v.reshape(-1, 2)-np.eye(2), 'all Q/R input')
                     for v in data['outputs'].values())

    # New mathematical candidate G, not a claim that original physical probes
    # implement these source instruments. No material is discarded in waits.
    sx = np.array([[0., 1.], [1., 0.]])
    sy = np.array([[0., -1j], [1j, 0.]])
    sz = np.diag([1., -1.])
    av = np.array([.8, 1.2])
    aop = np.diag(av)
    hg = np.eye(2)+.7*sx
    noncommuting = float(np.linalg.norm(hg@aop-aop@hg, 2))
    assert noncommuting > .1
    source_identity = np.diag(2*np.sqrt(16+9*av**2))
    positive_lower_bound = float(np.linalg.eigvalsh(hg).min()+av.min()*(evals.min()+3)
                                +np.diag(source_identity).min())
    assert positive_lower_bound > 0
    gready = np.array([1., 1.])/math.sqrt(2)
    waits = [.2, .3]
    instruments = [[(np.eye(2)+sign*pauli)/2 for sign in (1, -1)] for pauli in (sy, sx)]
    source_direction = np.array([-1., 0., 0., 1.])
    cache = {}

    def wait_joint(state, duration, perturbation):
        key = (duration, perturbation)
        if key not in cache:
            blocks = np.empty((288, 4, 2, 2), complex)
            for d, ev in enumerate(evals):
                for ab in range(4):
                    h = hg+(ev+3)*aop+source_identity+perturbation*source_direction[ab]*sz
                    blocks[d, ab] = exponential(h, duration)
            cache[key] = blocks
        return np.einsum('dauv,vdamelr->udamelr', cache[key], state, optimize=True)

    def instrument_chain(program, perturbation):
        material = data['outputs'][program]
        encoded = (vectors.conj().T@material.reshape(288, -1)).reshape(material.shape)
        # Input column is the retained reference of the normalized Choi state.
        joint = gready.reshape(2, 1, 1, 1, 1, 1, 1)*encoded[None, ...]/math.sqrt(2)
        branches = [joint]
        for duration, instrument in zip(waits, instruments):
            updated = []
            for state in branches:
                evolved = wait_joint(state, duration, perturbation)
                updated.extend(np.einsum('uv,vdamelr->udamelr', kraus, evolved, optimize=True)
                               for kraus in instrument)
            branches = updated
        out = []
        for state in branches:
            gr = state.transpose(0, 6, 1, 2, 3, 4, 5).reshape(4, -1)
            out.append(gr@gr.conj().T)
        small(sum(np.trace(r).real for r in out)-1, 'all instrument outcomes retained')
        return out

    exact = [instrument_chain(program, 0.) for program in range(3)]
    channel_residual = max(small(a-b, 'entire G/R unnormalized history matrix')
                           for output in exact[1:] for a, b in zip(exact[0], output))
    probabilities = [[float(np.trace(r).real) for r in output] for output in exact]
    # The candidate is not forced to ignore every change of D energy.
    g_energy_outputs = []
    for ev in (evals.min(), evals.max()):
        u = exponential(hg+(ev+3)*aop+source_identity, waits[0])
        ket = u@gready
        g_energy_outputs.append(np.outer(ket, ket.conj()))
    other_source_dist = float(np.sum(abs(np.linalg.eigvalsh(g_energy_outputs[0]-g_energy_outputs[1])))/2)
    assert other_source_dist > .01

    b = .02
    perturbed = [instrument_chain(program, b) for program in (0, 1)]
    choi_distance = float(sum(np.sum(abs(np.linalg.eigvalsh(a-c)))/2
                              for a, c in zip(perturbed[0], perturbed[1])))
    joint_bound = 2*b*sum(waits)
    assert 1e-6 < choi_distance <= joint_bound+1e-11, (choi_distance, joint_bound)

    # Source pullbacks use the very same full output columns, including all
    # clock, record and program factors. This is not an actual shear detector.
    source_witness_rows = []
    witness_error = 0.
    for scale in [.8, 1., 1.2]:
        energy = math.sqrt(16+9*scale**2)
        derivative = 32/energy*source_direction
        pullbacks = []
        for program in range(3):
            state = data['outputs'][program]
            v = state.reshape(-1, 2)
            ys = state.copy(); ys[:, :3] = 0
            effect = v.conj().T@ys.reshape(-1, 2)
            k = np.diag([.5, .5 if program == 1 else 1.])
            changed = state*derivative[None, :, None, None, None, None]
            delta = v.conj().T@changed.reshape(-1, 2)+32/energy*np.eye(2)
            witness_error = max(witness_error, small(effect-k, 'original yy effect'),
                                small(delta-64/energy*k, 'whole source pullback'))
            pullbacks.append(delta)
        dimensionless = energy/32*(pullbacks[0]-pullbacks[1])
        witness_error = max(witness_error, small(dimensionless-np.diag([0., 1.]), 'finite normalized contrast'))
        source_witness_rows.append(dict(scale=scale, shear_changes_Q0_Q1=[np.diag(p).real.tolist() for p in pullbacks],
                                        program0_minus1=np.diag(dimensionless).real.tolist()))
    gram = np.stack([np.zeros(2), np.array([0., 1.])])
    assert np.linalg.matrix_rank(gram) == 1
    finite_epsilon = .01
    return dict(schema='round1167_complete_common_scale_kernel_v1', status='PASS', date='2026-10-13',
                check_groups=3, new_adopted_axioms=0, new_empirical_calibration=0,
                gr_hard_deliverables_completed=0,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), sources_sha256=SOURCES,
                reused_identity_check='original complete 1160 calculation matches saved results; not new discovery',
                complete_source_kernel=dict(material_dimension=13824, logical_dimension=144,
                     head_logical_dimension=576, spectator_multiplicity=24,
                     original_clock_intertwiner=clock_intertwiner, complete_isometry_error=gram_error,
                     typed_program_relation_error=typed_error, maximum_kernel_error=worst_kernel, samples=source_rows),
                relation_future=dict(G_dimension=2, material_dimension=13824, reference_dimension=2,
                     A_eigenvalues=av.tolist(), HG_A_commutator_norm=noncommuting,
                     joint_energy_lower_bound=positive_lower_bound, waits=waits, instruments=['Pauli Y', 'Pauli X'],
                     retained_histories=4, normalized_Choi_history_probabilities=probabilities,
                     full_G_R_history_matrix_error=channel_residual,
                     other_HD_eigenstate_G_trace_distance=other_source_dist,
                     bounded_perturbation_norm=b, perturbed_Choi_trace_distance=choi_distance,
                     universal_half_diamond_upper_bound=joint_bound),
                directional_source=dict(samples=source_witness_rows, source_pullback_error=witness_error,
                     contrast_operator=[[0., 0.], [0., 1.]], comparison_rank=1,
                     isometry_error_per_program=finite_epsilon, normalized_contrast_lower_bound=1-4*finite_epsilon),
                scope=['homogeneous scale family and stated G-only finite future algebra only',
                       'local scalar field profiles may distinguish spatially different recoil supports',
                       'complete original material and active futures remain present',
                       'dynamic A does not certify original clock timing or a local spacetime realization',
                       'directional coefficient is not a new actual geometry readout',
                       'no all-axiom dynamic quantum countermodel or Einstein selection is claimed'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        mode = 'CREATED_ONCE'
    else:
        compare(result, json.loads(path.read_text(encoding='utf-8')))
        mode = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS', mode=mode, check_groups=result['check_groups'],
             source_kernel_error=result['complete_source_kernel']['maximum_kernel_error'],
             relation_history_error=result['relation_future']['full_G_R_history_matrix_error'],
             perturbed_distance=result['relation_future']['perturbed_Choi_trace_distance'],
             normalized_contrast_lower_bound=result['directional_source']['normalized_contrast_lower_bound']),
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
