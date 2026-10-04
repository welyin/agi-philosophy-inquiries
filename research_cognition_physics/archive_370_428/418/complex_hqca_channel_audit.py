"""Round 418: a vacuum-preserving complex extension of a local HQCA.

The finite tests use actual two-cell rules, not a pre-assumed history clock.
General input/reference and waiting-time claims are proved in the note.
"""
import argparse
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import platform
import unittest

import numpy as np

I2 = np.eye(2, dtype=complex)
T = np.diag([1, np.exp(1j * np.pi / 4)])
PHASE = T @ T
HAD = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
RY = np.array([[1, -1], [1, 1]], complex) / np.sqrt(2)
W = np.eye(4, dtype=complex)
W[2:, 2:] = RY
SWAP = np.eye(4, dtype=complex)[[0, 2, 1, 3]]
CNOT = np.eye(4, dtype=complex)[[0, 1, 3, 2]]
GATES = {'I': np.eye(4), 'W': W, 'S': SWAP, 'T': np.kron(T, I2)}
SYMBOLS = ['.', '>', 'I', 'W', 'S', 'T']


def evolution(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * t * values)) @ vectors.conj().T


def pair_action(g, vectors, n, site):
    """Act on two adjacent physical data bits, keeping all input columns."""
    a = np.asarray(vectors).reshape(2**site, 4, 2**(n-site-2), -1)
    return np.einsum('ab,lbrk->lark', g, a).reshape(2**n, -1)


def local_hamiltonian():
    # Cell ordering is (program, data) at each of the two sites.
    r = np.zeros((144, 144), complex)
    for name, gate in GATES.items():
        a = SYMBOLS.index(name)
        for hole in (0, 1):
            g = np.eye(4) if hole == 0 else gate
            for x in range(4):
                for y in range(4):
                    src = (2*hole+x//2)*12 + 2*a+x%2
                    dst = (2*a+y//2)*12 + 2*hole+y%2
                    r[dst, src] = g[y, x]
    return -(r+r.conj().T), r


@lru_cache(None)
def history(word=('I', 'S', 'I', 'T')):
    # N=2, K=2. Work sites 4,5 in an actual eight-cell chain.
    n, m, work = 8, 4, 4
    configs = list(combinations(range(n), m))
    ids = {c: i for i, c in enumerate(configs)}
    holes = ('.', '>', '.', '>')
    programs = []
    for c in configs:
        gates, empty = iter(word), iter(holes)
        programs.append(tuple(next(gates) if j in c else next(empty) for j in range(n)))
    seed = ids[tuple(range(m, n))]
    embed = np.zeros((2**n, 4), complex)
    for x in range(4):
        embed[x << (n-work-2), x] = 1
    edges = []
    adjacency = [[] for _ in configs]
    for i, program in enumerate(programs):
        for site in range(n-1):
            left, right = program[site:site+2]
            if left in ('.', '>') and right in GATES:
                new = list(program)
                new[site:site+2] = [right, left]
                occupied = tuple(j for j, s in enumerate(new) if s in GATES)
                j = ids[occupied]
                assert tuple(new) == programs[j]
                physical = np.eye(4) if left == '.' else GATES[right]
                active = left == '>' and site == work
                if left == '>' and site in (work-1, work+1):
                    assert right == 'I'  # the original guard identities
                logical = physical if active else np.eye(4)
                edges.append((i, j, site, physical, logical))
                adjacency[i].append((j, logical))
                adjacency[j].append((i, logical.conj().T))
    prefixes = {seed: np.eye(4, dtype=complex)}
    queue = [seed]
    loop_error = 0.
    for i in queue:
        for j, gate in adjacency[i]:
            candidate = gate @ prefixes[i]
            if j in prefixes:
                loop_error = max(loop_error, np.linalg.norm(candidate-prefixes[j]))
            else:
                prefixes[j] = candidate
                queue.append(j)
    assert len(prefixes) == len(configs)
    clock = np.zeros((len(configs), len(configs)))
    lifted = np.zeros((4*len(configs), 4*len(configs)), complex)
    physical_error = 0.
    for i, j, site, physical, logical in edges:
        clock[j, i] = clock[i, j] = -1
        lifted[4*j:4*j+4, 4*i:4*i+4] = -logical
        lifted[4*i:4*i+4, 4*j:4*j+4] = -logical.conj().T
        physical_error = max(physical_error, np.linalg.norm(
            pair_action(physical, embed @ prefixes[i], n, site) - embed @ prefixes[j]))
    target = GATES[word[3]] @ GATES[word[1]]
    done = [i for i, c in enumerate(configs) if c[-1] < m]
    done_error = max(np.linalg.norm(prefixes[i]-target) for i in done)
    return dict(configs=configs, seed=seed, prefixes=prefixes, clock=clock,
                lifted=lifted, target=target, done=done, edges=edges,
                loop_error=float(loop_error), physical_error=float(physical_error),
                done_error=float(done_error))


def trace_distance(a, b):
    d = a-b
    return float(np.sum(np.abs(np.linalg.eigvalsh((d+d.conj().T)/2)))/2)


def path_spectrum(length):
    j = np.arange(1, length+1)
    v = np.sqrt(2/(length+1)) * np.sin(np.pi * np.outer(j, j)/(length+1))
    return -2*np.cos(np.pi*j/(length+1)), v


def time_average_factor(values, horizon):
    d = values[:, None]-values[None, :]
    return np.exp(-.5j*horizon*d)*np.sinc(horizon*d/(2*np.pi))


@lru_cache(None)
def report():
    h, r = local_hamiltonian()
    number = np.diag([int(p in GATES) for p in SYMBOLS for _ in range(2)])
    left, right = np.kron(number, np.eye(12)), np.kron(np.eye(12), number)
    comm = h@left-left@h
    nested = comm@right-right@comm
    vacuum = np.eye(4)[:, 0]
    cnot_compiled = np.kron(PHASE, PHASE.conj().T) @ W @ W @ np.kron(I2, PHASE)
    h_compiled = RY @ np.linalg.matrix_power(T, 4)
    seed_one = np.vstack([np.zeros((2, 2)), np.eye(2)])
    histories = [history(), history(('I', 'W', 'I', 'T'))]
    a = histories[0]
    choi = np.eye(4).reshape(16)/2
    initial = np.zeros((70*4, 4), complex)
    initial[4*a['seed']:4*a['seed']+4] = np.eye(4)/2
    target_choi = np.kron(a['target'], np.eye(4)) @ choi
    channel_rows = []
    fermion_error = 0.
    for t in (.0, .7, 2.3, 4.8):
        amp = evolution(a['clock'], t)[:, a['seed']]
        full = evolution(a['lifted'], t) @ initial
        expected = np.concatenate([amp[i]*a['prefixes'][i]/2 for i in range(70)], axis=0)
        blocks = full.reshape(70, 16)
        reduced = blocks.T @ blocks.conj()
        mixture = sum(abs(amp[i])**2*np.outer(
            (a['prefixes'][i]/2).reshape(16), (a['prefixes'][i]/2).reshape(16).conj()) for i in range(70))
        p_done = float(np.sum(np.abs(amp[a['done']])**2))
        channel_rows.append(dict(time=t, joint_reference_error=float(np.linalg.norm(full-expected)),
            mixture_error=float(np.linalg.norm(reduced-mixture)), done_probability=p_done,
            choi_trace_error=trace_distance(reduced, np.outer(target_choi, target_choi.conj())),
            general_half_diamond_upper=1-p_done))
        single = np.zeros((8, 8))
        for j in range(7):
            single[j, j+1] = single[j+1, j] = -1
        one = evolution(single, t)
        dets = np.array([np.linalg.det(one[np.ix_(c, range(4, 8))]) for c in a['configs']])
        fermion_error = max(fermion_error, float(np.max(np.abs(dets-amp))))
    waiting = []
    for m, f in ((2, 4), (3, 8), (4, 16)):
        length = (f+2)*m
        values, vectors = path_spectrum(length)
        gap = float(np.min(np.diff(values)))
        exact_gap = 4*np.sin(3*np.pi/(2*(length+1)))*np.sin(np.pi/(2*(length+1)))
        horizon = 80*(length-1)/gap
        fac = time_average_factor(values, horizon)
        left_block = vectors[:f*m].T @ vectors[:f*m]
        occupied = vectors[f*m:].T @ vectors[f*m:]
        finite_mean = float(np.sum(left_block*occupied*fac).real)
        infinite_mean = float(np.sum(np.diag(left_block)*np.diag(occupied)))
        delta = 2*(length-1)/(horizon*gap)
        max_l1 = 0.
        for c in (0, length//2, length-1):
            coeff = np.outer(vectors[c], vectors[c])
            mean = np.einsum('xa,xb,ab->x', vectors, vectors, coeff*fac).real
            limiting = (vectors*vectors) @ (vectors[c]*vectors[c])
            max_l1 = max(max_l1, float(np.sum(abs(mean-limiting))))
        waiting.append(dict(M=m, f=f, length=length, particles=2*m, local_dimension=12,
            gap=gap, gap_formula_error=float(abs(gap-exact_gap)),
            gap_lower_bound=12/(length+1)**2, horizon=horizon, delta=delta,
            finite_mean_left=finite_mean, infinite_mean_left=infinite_mean,
            infinite_mean_formula=m*(2*f*m+1)/(length+1), max_sample_l1=max_l1,
            mean_error=float(abs(finite_mean-infinite_mean)),
            rigorous_average_success_lower=max(0., (f-2)*m/(length+1)-2*delta)))
    # An actual coherent measurement + feedback circuit, retaining its environment.
    instrument = np.zeros((8, 2), complex)
    instrument[0, 0] = instrument[4, 1] = 1
    instrument = pair_action(CNOT, instrument, 3, 0)
    instrument = pair_action(CNOT, instrument, 3, 1)
    # Controlled X from record site 1 to data site 0.
    instrument = pair_action(SWAP @ CNOT @ SWAP, instrument, 3, 0)
    kraus = []
    for e in range(2):
        kraus.append(instrument.reshape(2, 2, 2, 2)[:, :, e, :].reshape(4, 2))
    desired = [np.zeros((4, 2), complex) for _ in range(2)]
    desired[0][0, 0] = 1
    desired[1][1, 1] = 1
    joint = instrument/np.sqrt(2)  # maximally entangled input, last column = reference
    rho_dm_r = sum(np.outer((k/np.sqrt(2)).reshape(8), (k/np.sqrt(2)).reshape(8).conj()) for k in kraus)
    record_offdiagonal = rho_dm_r.reshape(2, 2, 2, 2, 2, 2)[:, 0, :, :, 1, :]
    return dict(round=418, scientific_baseline=417,
        scope='Fixed nearest-neighbor 12-state rule: all-input complex circuit channel and finite random-time completion; not a complete cognitive countermodel',
        local_rule=dict(dimension=12, matrix_dimension=144, norm=float(np.linalg.norm(h, 2)),
            hermitian_error=float(np.linalg.norm(h-h.conj().T)), r_square_error=float(np.linalg.norm(r@r)),
            adjacent_response_witness=float(np.linalg.norm(nested, 2)),
            vacuum_errors={k: float(np.linalg.norm(g@vacuum-vacuum)) for k,g in GATES.items()},
            invalid_hadamard_vacuum_leakage=float(np.sum(abs((np.kron(HAD,I2)@vacuum)[1:])**2))),
        compilation=dict(hadamard_error=float(np.linalg.norm(h_compiled-HAD)),
            cnot_error=float(np.linalg.norm(cnot_compiled-CNOT)),
            one_seed_return_error=float(np.linalg.norm(W@seed_one-seed_one@RY)),
            one_seed_is_extra_resource=True),
        history=dict(configurations=70, directed_forward_edges=len(a['edges']),
            maximum_loop_error=max(x['loop_error'] for x in histories),
            maximum_physical_rule_error=max(x['physical_error'] for x in histories),
            maximum_done_prefix_error=max(x['done_error'] for x in histories),
            common_program_graph_error=float(np.linalg.norm(histories[0]['clock']-histories[1]['clock'])),
            distinct_catalog_target_distance=float(np.linalg.norm(histories[0]['target']-histories[1]['target'])),
            fermion_amplitude_error=fermion_error),
        channels=channel_rows, waiting=waiting,
        instrument=dict(isometry_error=float(np.linalg.norm(instrument.conj().T@instrument-np.eye(2))),
            kraus_error=max(float(np.linalg.norm(k-d)) for k,d in zip(kraus,desired)),
            classical_record_offdiagonal=float(np.linalg.norm(record_offdiagonal)),
            joint_reference_norm=float(np.linalg.norm(joint)), all_environment_retained_before_reduction=True),
        same_fixed_finite_total_H_for_all_tasks=False, fixed_local_rule_for_all_chain_sizes=True,
        timing_and_final_access_internally_implemented=False, all_late_time_guarantee=False,
        no_postselection_in_output_error_bound=True, physical_dimension_derived=False,
        full_cognitive_countermodel_completed=False, phase_closure_triggered=False)


class AuditTests(unittest.TestCase):
    def test_01_actual_local_rule_and_vacuum_guard(self):
        r = report()['local_rule']
        self.assertAlmostEqual(r['norm'], 1)
        self.assertAlmostEqual(r['adjacent_response_witness'], 1)
        self.assertLess(max(r['hermitian_error'],r['r_square_error'],*r['vacuum_errors'].values()),1e-14)
        self.assertAlmostEqual(r['invalid_hadamard_vacuum_leakage'], .5)

    def test_02_complex_compilation_counts_seed(self):
        a=report()['compilation']
        self.assertLess(max(a['hadamard_error'],a['cnot_error'],a['one_seed_return_error']),2e-14)
        self.assertTrue(a['one_seed_is_extra_resource'])

    def test_03_every_local_edge_preserves_unknown_input(self):
        a=report()['history']
        self.assertEqual(a['configurations'],70)
        self.assertLess(max(a['maximum_loop_error'],a['maximum_physical_rule_error'],a['maximum_done_prefix_error']),2e-14)

    def test_04_same_local_hardware_different_programs(self):
        a=report()['history']
        self.assertEqual(a['common_program_graph_error'],0)
        self.assertGreater(a['distinct_catalog_target_distance'],1)

    def test_05_joint_reference_evolution_and_channel_bound(self):
        for a in report()['channels']:
            self.assertLess(max(a['joint_reference_error'],a['mixture_error']),2e-13)
            self.assertLessEqual(a['choi_trace_error'],a['general_half_diamond_upper']+2e-13)

    def test_06_independent_fermionic_amplitudes(self):
        self.assertLess(report()['history']['fermion_amplitude_error'],2e-13)

    def test_07_explicit_finite_spectral_waiting_bound(self):
        for a in report()['waiting']:
            self.assertLess(a['gap_formula_error'],2e-14)
            self.assertGreaterEqual(a['gap'],a['gap_lower_bound']-2e-14)
            self.assertAlmostEqual(a['infinite_mean_left'],a['infinite_mean_formula'])
            self.assertLessEqual(a['max_sample_l1'],a['delta']+1e-13)
            self.assertLessEqual(a['mean_error'],2*a['M']*a['delta']+1e-13)
            self.assertGreater(a['rigorous_average_success_lower'],0)

    def test_08_internal_instrument_and_reference_record(self):
        a=report()['instrument']
        self.assertLess(max(a['isometry_error'],a['kraus_error'],a['classical_record_offdiagonal']),2e-14)
        self.assertAlmostEqual(a['joint_reference_norm'],1)
        self.assertFalse(report()['timing_and_final_access_internally_implemented'])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(AuditTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    out=dict(report(),checks=dict(run=result.testsRun,failures=len(result.failures),errors=len(result.errors)),
             runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with Path(__file__).with_name('complex_hqca_channel_audit_results.json').open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))
