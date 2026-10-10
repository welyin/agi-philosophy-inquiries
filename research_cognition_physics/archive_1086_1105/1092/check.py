"""Finite checks of bounded directed-type stabilization; default read-only."""
from pathlib import Path
import itertools
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent


def run():
    rng = np.random.default_rng(109220261010)
    groups = []

    def record(name, residuals, **data):
        maximum = float(max(residuals, default=0.0))
        assert maximum < 1e-11, (name, maximum)
        groups.append(dict(name=name, passed=True, maximum_residual=maximum, **data))

    def unitary(n):
        q, r = np.linalg.qr(rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)))
        return q @ np.diag(np.diag(r)/np.abs(np.diag(r)))

    stages = [frozenset(c) for n in range(1, 5)
              for c in itertools.combinations(range(4), n)]
    labels = {s: sorted({i % 3 for i in s}) for s in stages}
    gauge = {s: unitary(len(labels[s])) for s in stages}
    base = frozenset((0, 1, 2))

    def embedding(s, t):
        e = np.zeros((len(labels[t]), len(labels[s])), dtype=complex)
        for k, x in enumerate(labels[s]):
            e[labels[t].index(x), k] = 1
        return gauge[t].conj().T @ e @ gauge[s]

    pairs = [(s, t) for s in stages for t in stages if s <= t]
    triples = [(s, t, u) for s, t in pairs for u in stages if t <= u]
    record('isometric_inclusions', [np.linalg.norm(embedding(s, t).conj().T @
           embedding(s, t)-np.eye(len(labels[s]))) for s, t in pairs], pairs=len(pairs))
    record('directed_coherence', [np.linalg.norm(embedding(t, u) @ embedding(s, t)-
           embedding(s, u)) for s, t, u in triples], triples=len(triples))
    collapse = {}
    residuals = []
    upper_count = 0
    for s in stages:
        uppers = [t for t in stages if s | base <= t]
        maps = [embedding(base, t).conj().T @ embedding(s, t) for t in uppers]
        collapse[s] = maps[0]
        residuals += [np.linalg.norm(m-maps[0]) for m in maps]
        upper_count += len(uppers)
    record('upper_choice_independence', residuals, comparisons=upper_count)
    record('common_carrier_intertwining', [np.linalg.norm(collapse[t] @ embedding(s,t)-
           collapse[s]) for s,t in pairs])

    v = unitary(3)
    effect = v @ np.diag([0.12, 0.49, 0.88]) @ v.conj().T
    effects = {s: gauge[s].conj().T @ effect[np.ix_(labels[s], labels[s])] @ gauge[s]
               for s in stages}
    record('effects_natural_and_descend',
           [np.linalg.norm(embedding(s,t).conj().T @ effects[t] @ embedding(s,t)-effects[s])
            for s,t in pairs] +
           [np.linalg.norm(collapse[s].conj().T @ effects[base] @ collapse[s]-effects[s])
            for s in stages])

    phase = np.exp(1j*np.array([0.2, 0.7, -0.4]))
    phase_gates = {s: gauge[s].conj().T @ np.diag(phase[labels[s]]) @ gauge[s] for s in stages}

    def channel(s, rho, ref=1):
        z = np.kron(phase_gates[s], np.eye(ref))
        return 0.4*rho + 0.6*z @ rho @ z.conj().T

    residuals = []
    for s in stages:
        n = len(labels[s])*2
        psi = rng.normal(size=n)+1j*rng.normal(size=n)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        e = np.kron(collapse[s], np.eye(2))
        residuals.append(np.linalg.norm(channel(base, e @ rho @ e.conj().T, 2)-
                                       e @ channel(s, rho, 2) @ e.conj().T))
    record('full_channel_with_unknown_reference', residuals, reference_dimension=2)

    # Two equal-dimensional descriptions can disagree unless the effects are natural.
    first = np.diag([1.0, 0.0])
    second = np.array([[0.5, 0.5], [0.5, 0.5]])
    record('missing_naturality_diagnostic', [abs(first[0,0]-1), abs(second[0,0]-0.5)],
           same_state_probabilities=[1.0, 0.5], nature_of_check='violates_naturality_not_theorem')

    w = np.zeros((6, 4), complex)
    for p in range(2):
        for x in range(2):
            w[p*3+x+p, p*2+x] = 1
    record('controlled_move_image_vs_complete_type', [np.linalg.norm(w.conj().T @ w-np.eye(4))],
           input_dimension=4, image_rank=int(np.linalg.matrix_rank(w)),
           complete_output_dimension=6, same_composition_full_type_without_spare_capacity=False)
    u3 = np.zeros((3,3))
    for x in range(3):
        u3[(x+1)%3,x] = 1
    record('padding_is_not_global_translation', [np.linalg.norm(u3.T @ u3-np.eye(3))],
           valid_input_positions=[0,1], padded_input=2, padded_output=0,
           full_real_line_translation=False)

    # P,A,B,R are wholly unknown; W is separately supplied in |0>.
    join_gate = np.zeros((8,8),complex)
    for p, permutation in enumerate(([0,2,1,3], [1,0,2,3])):
        for col, row in enumerate(permutation):
            join_gate[4*p+row, 4*p+col] = 1
    insertion = np.zeros((48,24),complex)
    expected = np.zeros((48,24),complex)
    for p,a,b,r in itertools.product(range(2),range(2),range(2),range(3)):
        col = ((p*2+a)*2+b)*3+r
        insertion[((p*4+2*a)*2+b)*3+r,col] = 1
        expected[((p*4+a+p)*2+b)*3+r,col] = 1
    global_gate = np.kron(join_gate, np.eye(6))
    joined = global_gate @ insertion
    psi = rng.normal(size=24)+1j*rng.normal(size=24)
    psi /= np.linalg.norm(psi)
    recovered = insertion.conj().T @ global_gate.conj().T @ joined @ psi
    record('explicit_workspace_and_private_unknown_member',
           [np.linalg.norm(global_gate.conj().T @ global_gate-np.eye(48)),
            np.linalg.norm(joined-expected), np.linalg.norm(recovered-psi)],
           unknown_private_dimension=2, reference_dimension=3, added_prepared_workspace_dimension=2,
           input_prepared_workspace_is_unknown=False,
           mathematical_recovery_not_physical_rollback=True,
           autonomous_motion_and_full_CO_join_certified=False)

    counts = [{'valid_catalog_size': k, 'union_with_translate_size': k+1,
               'extra_orthogonal_port_required': True} for k in [1,2,4,8,16]]
    return {'schema':'round1092_bounded_directed_types_v1', 'round':1092, 'status':'PASS',
            'seed':109220261010, 'numpy_version':np.__version__, 'types':len(stages),
            'maximum_type_dimension':3, 'groups':groups, 'translation_catalogs':counts,
            'scope':{'new_universal_axioms':0, 'scientific_count_increment':0,
                     'proof_of_entire_new_conjecture':False, 'Lorentz_derived':False,
                     'bounded_directed_family_lemma_analytic_proof_in_proof_md':True,
                     'numerics_do_not_prove_infinite_quantifier':True,
                     'round1088_compactness_reused':True}}


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for left, right in zip(actual, saved):
            compare(left, right)
    elif isinstance(actual, float):
        assert abs(actual-saved) < 1e-12, (actual, saved)
    else:
        assert actual == saved, (actual, saved)


if __name__ == '__main__':
    result = run()
    output = HERE/'results.json'
    if '--write' in sys.argv:
        if output.exists():
            raise SystemExit('Refusing to overwrite existing results.')
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    else:
        compare(result, json.loads(output.read_text(encoding='utf8')))
    print(json.dumps({'status':'PASS', 'round':1092, 'groups':len(result['groups']),
        'maximum_residual':max(g['maximum_residual'] for g in result['groups']),
        'new_scientific_count':0, 'whole_conjecture_decided':False},ensure_ascii=False))
