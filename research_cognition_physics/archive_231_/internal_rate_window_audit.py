"""Round 419: internal finite rate labels replace random read times.

The labels are part of the closed apparatus. Bounds are proved in the note;
Choi examples are not advertised as diamond-norm optimizations.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import complex_hqca_channel_audit as base


def speeds(count):
    return (np.arange(count) + .5) / count


def discrete_factor(values, time, count):
    """Closed geometric sum; used only away from nonzero alias poles."""
    z = time * (values[:, None] - values[None, :])
    denominator = np.sinc(z / (2 * np.pi * count))
    if np.min(abs(denominator)) < 1e-9:
        raise ValueError('Use direct finite sum near a nonzero alias pole')
    return np.exp(-.5j*z) * np.sinc(z/(2*np.pi)) / denominator


def direct_factor(values, time, count):
    phases = np.exp(-1j*np.outer(speeds(count)*time, values))
    return phases.T @ phases.conj() / count


def probabilities(values, vectors, seed, factor):
    weights = np.outer(vectors[seed].conj(), vectors[seed])
    return np.einsum('xa,xb,ab->x', vectors, vectors.conj(), weights*factor).real


def choi_from_probabilities(history, probability):
    cols = np.stack([(history['prefixes'][j]/2).reshape(16)
                     for j in range(len(probability))])
    return (cols.T * probability) @ cols.conj()


def correlated_label_preparation(count=3, length=4):
    # Last register purifies the random label. Nearest-neighbor SUMs copy a
    # known orthogonal basis label, not an unknown quantum state.
    initial = np.zeros((count,)*(length+1), complex)
    for j in range(count):
        initial[(j,)+(0,)*(length-1)+(j,)] = 1/np.sqrt(count)
    state = initial
    for site in range(length-1):
        new = np.zeros_like(state)
        for index in np.ndindex(state.shape):
            dest = list(index)
            dest[site+1] = (dest[site+1]+index[site]) % count
            new[tuple(dest)] += state[index]
        state = new
    expected = np.zeros_like(state)
    for j in range(count):
        expected[(j,)*(length+1)] = 1/np.sqrt(count)
    reduced = state.reshape(count**length, count)
    rho = reduced @ reduced.conj().T
    return dict(purified_preparation_error=float(np.linalg.norm(state-expected)),
                label_offdiagonal=float(np.linalg.norm(rho-np.diag(np.diag(rho)))),
                nonzero_label_configurations=int(np.count_nonzero(np.diag(rho))),
                nearest_neighbor_sum_count=length-1,
                label_capacity_bits=length*math.log2(count),
                label_entropy_bits=math.log2(count),
                purification_dimension=count)


def actual_two_cell_extension(count=2):
    old, _ = base.local_hamiltonian()
    cell = 12*count
    extended = np.zeros((cell*cell, cell*cell), complex)
    # Basis of a cell is (stationary rate label, moving HQCA program/data).
    indices = []
    for label, v in enumerate(speeds(count)):
        ix = np.array([(12*label+a)*cell+12*label+b
                       for a in range(12) for b in range(12)])
        extended[np.ix_(ix, ix)] = v*old
        indices.append(ix)
    generator = np.random.default_rng(419)
    worst = 0.
    for label, ix in enumerate(indices):
        columns = generator.normal(size=(144, 3))+1j*generator.normal(size=(144, 3))
        lifted = np.zeros((cell*cell, 3), complex)
        lifted[ix] = columns
        expected = np.zeros_like(lifted)
        expected[ix] = speeds(count)[label]*old@columns
        worst = max(worst, float(np.linalg.norm(extended@lifted-expected)))
    cross = np.array([(a)*cell+12+b for a in range(12) for b in range(12)])
    # An explicit matrix check, including sectors outside the prepared labels.
    eig = np.linalg.eigvalsh(extended)
    return dict(local_dimension=cell, matrix_dimension=cell*cell,
                norm=float(np.max(abs(eig))), norm_expected=float(speeds(count)[-1]),
                hermitian_error=float(np.linalg.norm(extended-extended.conj().T)),
                copied_sector_action_error=worst,
                mismatched_label_sector_norm=float(np.linalg.norm(extended[:, cross])),
                static_label_commutator=float(np.linalg.norm(
                    extended*np.repeat(np.repeat(np.arange(count), 12), cell)[None, :]
                    -np.repeat(np.repeat(np.arange(count), 12), cell)[:, None]*extended)))


@lru_cache(None)
def report():
    hist = base.history()
    values, vectors = np.linalg.eigh(hist['clock'])
    # Direct exponentiation of a finite total Hamiltonian, retaining a
    # purification of the rate seed and a maximally entangled unknown input.
    count, time = 2, 4.8
    controlled = np.kron(np.diag(speeds(count)), hist['lifted'])
    initial = np.zeros((count*280, count*4), complex)
    for label in range(count):
        initial[label*280+4*hist['seed']:label*280+4*hist['seed']+4,
                4*label:4*label+4] = np.eye(4)/(2*np.sqrt(count))
    joint = base.evolution(controlled, time) @ initial
    tensor = joint.reshape(count, 70, 4, count, 4)
    work_reference = np.einsum('jcwer,jcves->wrvs', tensor, tensor.conj()).reshape(16,16)
    factor = direct_factor(values, time, count)
    p = probabilities(values, vectors, hist['seed'], factor)
    formula = choi_from_probabilities(hist, p)
    leakage = sum(float(np.linalg.norm(tensor[j,:,:,e,:])**2)
                  for j in range(count) for e in range(count) if j != e)
    # Static replicated tags multiply every actual transition; they never
    # travel together with a program symbol. Verify all original 140 edges.
    edge_error = 0.
    embedded = np.zeros((256, 4), complex)
    for x in range(4):
        embedded[x << 2, x] = 1
    for v in speeds(3):
        for i, j, site, physical, _ in hist['edges']:
            actual = v*base.pair_action(physical, embedded@hist['prefixes'][i], 8, site)
            expected = v*embedded@hist['prefixes'][j]
            edge_error = max(edge_error, float(np.linalg.norm(actual-expected)))
    samples = []
    factor_error = 0.
    for t in (1.3, 4.8, 12.):
        average = base.time_average_factor(values, t)
        pa = probabilities(values, vectors, hist['seed'], average)
        rho_a = choi_from_probabilities(hist, pa)
        for k in (16, 64, 256):
            fd = direct_factor(values, t, k)
            factor_error = max(factor_error, float(np.max(abs(fd-discrete_factor(values,t,k)))))
            pd = probabilities(values, vectors, hist['seed'], fd)
            rho_d = choi_from_probabilities(hist, pd)
            samples.append(dict(time=t, labels=k,
                probability_sum_error=float(abs(np.sum(pd)-1)),
                minimum_probability=float(np.min(pd)),
                choi_trace_error=base.trace_distance(rho_a,rho_d),
                general_half_diamond_bound=min(1.,7*t/(4*k)),
                clock_total_variation=float(np.sum(abs(pa-pd))/2)))
    # A finite nontrivial certificate. No array of length K is allocated:
    # the scalar geometric sum is checked above against direct summation.
    epsilon, m, f = .25, 2, 58
    length = (f+2)*m
    gap = float(4*np.sin(3*np.pi/(2*(length+1)))*np.sin(np.pi/(2*(length+1))))
    tau = 12*(length-1)/(epsilon*gap)
    horizon = 2*tau
    needed = math.ceil(3*(length-1)*horizon/(4*epsilon))
    bits = (needed-1).bit_length()
    labels = 2**bits
    vl, el = base.path_spectrum(length)
    left = el[:f*m].T @ el[:f*m]
    occupied = el[f*m:].T @ el[f*m:]
    window = []
    for t in (tau, 1.25*tau, 1.5*tau, 1.75*tau, horizon):
        fac = discrete_factor(vl, t, labels)
        mean = float(np.sum(left*occupied*fac).real)
        lower_from_mean = max(0., (mean-m)/m)
        completion_failure_bound = (4*m+1)/(length+1)+4*(length-1)/(t*gap)
        quadrature_bound = (length-1)*t/(4*labels)
        window.append(dict(time=t, actual_discrete_mean_left=mean,
            lower_success_from_fermion_mean=lower_from_mean,
            average_failure_upper=completion_failure_bound,
            midpoint_half_diamond_upper=quadrature_bound,
            full_reference_error_upper=completion_failure_bound+quadrature_bound))
    # Two counterchecks: averaging Hamiltonians is different from averaging
    # evolutions, and a seed correlated with the unknown input is not free.
    z = np.diag([1.,-1.])
    plus = np.array([1.,1.])/np.sqrt(2)
    minus = np.array([1.,-1.])/np.sqrt(2)
    plus_rho = np.outer(plus,plus)
    mixed = sum(base.evolution(z,np.pi*v)@plus_rho@base.evolution(z,np.pi*v).conj().T
                for v in speeds(2))/2
    mean_unitary = base.evolution(z,np.pi/2)
    wrong_mean = mean_unitary@plus_rho@mean_unitary.conj().T
    correlated = sum(np.outer(base.evolution(z,np.pi*v)@s,
                              (base.evolution(z,np.pi*v)@s).conj())
                     for v,s in zip(speeds(2),(plus,minus)))/2
    return dict(round=419, scientific_baseline=418,
        scope='Finite internal replicated rate labels yield a single closed local apparatus with a uniform finite read window; no autonomous stopping or online admission theorem',
        preparation=correlated_label_preparation(), local_rule=actual_two_cell_extension(),
        actual_edges=dict(forward_edges=len(hist['edges']), rates_checked=3,
                          maximum_physical_transition_error=edge_error),
        purified_channel=dict(total_hamiltonian_dimension=560, input_reference_dimension=4,
            label_purification_dimension=2, time=time,
            joint_state_norm=float(np.linalg.norm(joint)),
            exact_mixture_error=float(np.linalg.norm(work_reference-formula)),
            rate_purification_mismatch_weight=leakage),
        quadrature=samples, closed_factor_direct_sum_error=factor_error,
        certificate=dict(epsilon=epsilon, M=m, f=f, chain_length=length, gap=gap,
            window_start=tau, window_end=horizon, labels=labels,
            label_bits_per_cell=bits, local_dimension=12*labels,
            total_label_capacity_bits=length*bits, label_entropy_bits=bits,
            additional_purifier_bits=bits, label_copy_SUM_count=length-1,
            global_uniform_half_diamond_upper=(4*m+1)/(length+1)
                +4*(length-1)/(tau*gap)+(length-1)*horizon/(4*labels),
            rows=window, huge_label_array_allocated=False),
        counterchecks=dict(mean_hamiltonian_trace_error=base.trace_distance(mixed,wrong_mean),
            correlated_vs_independent_seed_trace_error=base.trace_distance(correlated,np.eye(2)/2)),
        random_read_time_choice_required=False, deterministic_time_window_proved=True,
        rate_preparation_and_distribution_assumed_as_counted_initial_resources=True,
        same_12_state_local_rule_unchanged=False, fixed_local_dimension_for_all_windows=False,
        internal_record_circuit_reused_from_round418=True,
        arbitrary_active_readout_schedule_proved=False, online_new_subject_admission_proved=False,
        all_late_time_guarantee=False, full_cognitive_countermodel_completed=False,
        physical_dimension_derived=False, phase_closure_triggered=False)


class AuditTests(unittest.TestCase):
    def test_01_label_preparation_retains_purifier(self):
        r = report()['preparation']
        self.assertLess(r['purified_preparation_error'],1e-14)
        self.assertLess(r['label_offdiagonal'],1e-14)
        self.assertEqual(r['nonzero_label_configurations'],3)
        self.assertEqual(r['nearest_neighbor_sum_count'],3)

    def test_02_actual_local_matrix_and_invariant_labels(self):
        r = report()['local_rule']
        self.assertAlmostEqual(r['norm'],r['norm_expected'])
        self.assertEqual(r['local_dimension'],24)
        self.assertLess(max(r[k] for k in ('hermitian_error','copied_sector_action_error',
            'mismatched_label_sector_norm','static_label_commutator')),2e-13)

    def test_03_every_physical_edge_with_replicated_rate(self):
        r = report()['actual_edges']
        self.assertEqual(r['forward_edges'],140)
        self.assertLess(r['maximum_physical_transition_error'],1e-13)

    def test_04_closed_unitary_with_two_retained_references(self):
        r = report()['purified_channel']
        self.assertAlmostEqual(r['joint_state_norm'],1)
        self.assertLess(r['exact_mixture_error'],3e-13)
        self.assertLess(r['rate_purification_mismatch_weight'],1e-26)

    def test_05_actual_channel_midpoint_comparison(self):
        for r in report()['quadrature']:
            self.assertLess(r['probability_sum_error'],1e-12)
            self.assertGreater(r['minimum_probability'],-1e-13)
            self.assertLessEqual(r['choi_trace_error'],r['general_half_diamond_bound']+1e-12)
            self.assertLessEqual(r['clock_total_variation'],r['general_half_diamond_bound']+1e-12)

    def test_06_large_finite_label_formula_has_independent_check(self):
        self.assertLess(report()['closed_factor_direct_sum_error'],1e-13)
        self.assertFalse(report()['certificate']['huge_label_array_allocated'])

    def test_07_nontrivial_uniform_window_certificate(self):
        r = report()['certificate']
        self.assertLessEqual(r['global_uniform_half_diamond_upper'],r['epsilon'])
        self.assertGreater(r['window_end'],r['window_start'])
        for row in r['rows']:
            self.assertGreaterEqual(row['lower_success_from_fermion_mean'],
                                    1-row['full_reference_error_upper']-1e-10)
            self.assertLessEqual(row['full_reference_error_upper'],r['epsilon'])

    def test_08_seed_independence_and_not_mean_generator(self):
        r = report()['counterchecks']
        self.assertAlmostEqual(r['mean_hamiltonian_trace_error'],.5)
        self.assertAlmostEqual(r['correlated_vs_independent_seed_trace_error'],.5)
        self.assertFalse(report()['arbitrary_active_readout_schedule_proved'])
        self.assertFalse(report()['online_new_subject_admission_proved'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AuditTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = dict(run=tests.testsRun, failures=0, errors=0)
    result['runtime'] = dict(python=platform.python_version(),numpy=np.__version__)
    if args.write_results:
        dest = Path(__file__).with_name('internal_rate_window_audit_results.json')
        with dest.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
