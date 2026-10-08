"""Finite local record / interacting-energy audit. NumPy only; no plotting."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'local_record_energy_selection_results.json'
I = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1, -1]).astype(complex)]
HISTORY = [
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_370_428/research_note_391.md',
    'archive_742_763/research_note_751.md',
    'archive_819_853/research_note_830.md',
    'archive_1009_/research_note_1031.md',
    'archive_1009_/research_note_1033.md',
    'archive_1009_/1035/NEXT.md',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hamiltonian(js):
    return sum(j * np.kron(p, p) for j, p in zip(js, PAULI))


def root_psd(a):
    w, v = np.linalg.eigh(a)
    assert w.min() > -1e-12
    return (v * np.sqrt(np.maximum(w, 0))) @ v.conj().T


def adjoint(ks, a):
    return sum(k.conj().T @ a @ k for k in ks)


def norm(a):
    return float(np.linalg.norm(a, 2))


def instrument_stats(js, groups):
    h = hamiltonian(js)
    ks = [k for group in groups for k in group]
    full = [np.kron(k, I) for k in ks]
    tp_error = norm(adjoint(ks, I) - I)
    effect = adjoint(groups[0], I)
    eig = np.linalg.eigvalsh(effect)
    c = float(eig[-1] - eig[0])
    coeff = np.array([[np.trace(p @ k) / 2 for p in [I] + PAULI] for k in ks])
    s = float(np.sum(abs(coeff[:, 1:]) ** 2))
    fe = float(np.sum(abs(coeff[:, 0]) ** 2))
    qs = 4 * (sum(j * j for j in js) - np.array(js) ** 2)
    comms = [h @ k - k @ h for k in full]
    gram = sum(a.conj().T @ a for a in comms)
    b = float(np.trace(gram).real / 4)
    b_pauli = float(np.sum(qs * np.sum(abs(coeff[:, 1:]) ** 2, axis=0)))
    energies, vecs = np.linalg.eigh(h)
    transitions = sum(abs(vecs.conj().T @ k @ vecs) ** 2 for k in full)
    b_spectral = float(np.sum(transitions * (energies[:, None] - energies[None, :]) ** 2) / 4)
    center = float((energies[-1] + energies[0]) / 2)
    hc = h - center * np.eye(4)
    d1 = adjoint(full, hc) - hc
    d2 = adjoint(full, hc @ hc) - hc @ hc
    gram_identity_error = norm(gram - (d2 - hc @ d1 - d1 @ hc))
    eps1, eps2 = norm(d1), norm(d2)
    q_min = float(min(qs))
    sharp_bound = q_min * (1 - np.sqrt(max(0., 1 - c * c))) / 2
    p, q = eig
    fixed_effect_fe_max = float((1 + np.sqrt(max(0., p*q)) +
                                np.sqrt(max(0., (1-p)*(1-q)))) / 2)
    assert tp_error < 2e-12
    assert abs(fe+s-1) < 2e-12
    assert abs(b-b_pauli) < 2e-11 and abs(b-b_spectral) < 2e-11
    assert gram_identity_error < 2e-11
    assert b+2e-11 >= sharp_bound
    assert fe <= fixed_effect_fe_max+2e-12
    assert b <= eps2+2*norm(hc)*eps1+2e-11
    return dict(c=c, effect_eigenvalues=eig.tolist(), traceless_weight=s,
                entanglement_fidelity=fe, fixed_effect_fidelity_bound=fixed_effect_fe_max,
                B=b, pauli_B=b_pauli, spectral_B=b_spectral,
                q=list(map(float, qs)), sharp_lower_bound=float(sharp_bound),
                bound_slack=float(b-sharp_bound), energy_center=center,
                epsilon1=eps1, epsilon2_centered=eps2,
                moment_envelope=eps2+2*norm(hc)*eps1,
                gram_identity_error=gram_identity_error, tp_error=tp_error)


def luders(c, axis):
    hi, lo = np.sqrt((1+c)/2), np.sqrt((1-c)/2)
    a, b = (hi+lo)/2, (hi-lo)/2
    return [[a*I+b*PAULI[axis]], [a*I-b*PAULI[axis]]]


def random_instrument(rng):
    a = rng.normal(size=(8, 2)) + 1j*rng.normal(size=(8, 2))
    q, _ = np.linalg.qr(a)
    ks = [q[2*i:2*i+2] for i in range(4)]
    return [ks[:2], ks[2:]]


def generator_adjoint(ls, a, k=None):
    v = np.zeros_like(a)
    if k is not None:
        v += 1j * (k @ a - a @ k)
    for l in ls:
        ll = l.conj().T @ l
        v += l.conj().T @ a @ l - (ll @ a + a @ ll)/2
    return v


def boundary_checks():
    # Mean-preserving, variance-changing valid qutrit channel.
    h = np.diag([-1., 0., 1.]).astype(complex)
    p = .4
    ks = [np.diag([1, np.sqrt(1-p), 1]).astype(complex)]
    for target in (0, 2):
        k = np.zeros((3, 3), complex)
        k[target, 1] = np.sqrt(p/2)
        ks.append(k)
    d1 = adjoint(ks, h)-h
    d2 = adjoint(ks, h@h)-h@h
    assert norm(d1) < 1e-13
    assert norm(d2-np.diag([0, p, 0])) < 1e-13
    # Continuous-time version, no unitary counterterm can cancel the Gram.
    ls = [k/np.sqrt(p) for k in ks[1:]]
    g1 = generator_adjoint(ls, h, h)
    g2 = generator_adjoint(ls, h@h, h)
    gram = sum((h@l-l@h).conj().T @ (h@l-l@h) for l in ls)
    assert norm(g1) < 1e-13
    assert norm(g2-np.diag([0, 1, 0])) < 1e-13
    assert norm(g2-h@g1-g1@h-gram) < 1e-13
    # Genuine surviving local sharp record in a commuting interaction.
    ising = instrument_stats((0., 0., 2.), luders(1, 2))
    assert ising['c'] > .999999 and ising['B'] < 1e-13
    # Dropping locality: energy-projective recording exactly preserves moments.
    hf = hamiltonian((1., 2., 3.))
    ev, v = np.linalg.eigh(hf)
    projectors = [np.outer(v[:, j], v[:, j].conj()) for j in range(4)]
    nonlocal_mean_error = norm(adjoint(projectors, hf)-hf)
    nonlocal_second_error = norm(adjoint(projectors, hf@hf)-hf@hf)
    assert max(nonlocal_mean_error, nonlocal_second_error) < 1e-12
    # Merely testing the one maximally mixed state misses source disturbance.
    local = [np.kron(k, I) for group in luders(.6, 2) for k in group]
    mixed_mean_drift = abs(np.trace(adjoint(local, hf)-hf).real/4)
    mixed_second_drift = abs(np.trace(adjoint(local, hf@hf)-hf@hf).real/4)
    assert max(mixed_mean_drift, mixed_second_drift) < 1e-12
    return dict(qutrit_p=p, qutrit_mean_defect=norm(d1),
                qutrit_second_defect=norm(d2), qutrit_middle_output_probabilities=[p/2, 1-p, p/2],
                GKLS_mean_defect=norm(g1), GKLS_second_defect=norm(g2),
                GKLS_gram_residual=norm(g2-h@g1-g1@h-gram),
                Ising_surviving_record=ising,
                nonlocal_energy_record_mean_defect=nonlocal_mean_error,
                nonlocal_energy_record_second_defect=nonlocal_second_error,
                maximally_mixed_first_and_second_drift=[mixed_mean_drift, mixed_second_drift])


def run():
    rng = np.random.default_rng(1037)
    arbitrary = [instrument_stats((1., 2., 3.), random_instrument(rng)) for _ in range(24)]
    # Deliberately biased effects and output feedback are included, not excluded.
    p, q = .8, .1
    biased = [[np.diag(np.sqrt([p, q])).astype(complex)],
              [PAULI[0] @ np.diag(np.sqrt([1-p, 1-q]))]]
    biased_row = instrument_stats((1., 2., 3.), biased)
    saturation = []
    for js in ((1., 2., 3.), (-1., 2., 3.), (1., 1., 1.)):
        qs = 4*(sum(j*j for j in js)-np.array(js)**2)
        axis = int(np.argmin(qs))
        for c in (0., .2, .6, 1.):
            row = instrument_stats(js, luders(c, axis))
            assert abs(row['bound_slack']) < 1e-11
            saturation.append(dict(J=list(js), axis=axis, **row))
    main = instrument_stats((1., 2., 3.), luders(.6, 2))
    assert abs(main['B']-2.) < 1e-12 and abs(main['traceless_weight']-.1) < 1e-13
    shift_errors = []
    ks = [np.kron(k, I) for group in luders(.6, 2) for k in group]
    h = hamiltonian((1., 2., 3.))
    for e in (-20., -1., 0., 5.):
        hc = h-e*np.eye(4)
        d1, d2 = adjoint(ks, hc)-hc, adjoint(ks, hc@hc)-hc@hc
        value = float(np.trace(d2-hc@d1-d1@hc).real/4)
        shift_errors.append(abs(value-2.))
    assert max(shift_errors) < 1e-11
    return dict(round=1037, date='2026-10-08', new_calibration_groups=1,
                new_adopted_cognitive_axioms=0, physical_implementation_certified=False,
                full_parent_mapping_certified=False, goal_completed=False,
                arbitrary_instruments=arbitrary, biased_feedback_instrument=biased_row,
                saturation_cases=saturation, exact_rational_case=main,
                energy_origin_residual=max(shift_errors), boundary_cases=boundary_checks(),
                finite_budget_example=dict(B_budget=1., max_record_contrast=float(np.sqrt(.19)),
                                           requested_contrast=.6, required_B=2.),
                historical_source_sha256={p:sha(BASE/p) for p in HISTORY})


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert np.isclose(a, b, rtol=2e-11, atol=2e-11), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with OUT.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write('\n')
    else:
        compare(result, json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1037, arbitrary_instruments=len(result['arbitrary_instruments']),
                         saturation_cases=len(result['saturation_cases']),
                         rational_case_B=result['exact_rational_case']['B'],
                         passed=True), ensure_ascii=False))
