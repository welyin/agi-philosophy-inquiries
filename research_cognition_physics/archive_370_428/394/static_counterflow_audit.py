"""Round 394: a static finite-range obstruction for round 393's transport W.

The infinite-chain theorem is analytic. These checks verify its exact wire
action, support-orbit mechanism, and finite-volume limitations. They do not
test the full computing rule (C tensor I) W or a prepared-clock code space.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('static_counterflow_audit_results.json')
TRACKS = ('t', 'a', 'b', 's', 'e')
GROUPS = {'t': 'right', 'e': 'right', 's': 'fast_left',
          'a': 'slow_left', 'b': 'slow_left'}
VELOCITIES = {'right': 2, 'fast_left': -2, 'slow_left': -1}


def content_destination(track, site):
    """Forward transport W, on the infinite line, with no wraparound."""
    return {'t': ('t', site+1), 'a': ('b', site-1),
            'b': ('a', site), 's': ('s', site-1),
            'e': ('e', site+1)}[track]


def content_power(track, site, steps):
    for _ in range(steps):
        track, site = content_destination(track, site)
    return track, site


def orbit_support(support, velocities, steps):
    return tuple((group, site+velocities[group]*steps)
                 for group, site in support)


def diameter(support):
    positions = [site for _, site in support]
    return max(positions)-min(positions)


def support_audit(include_stationary=False):
    """Exhaust supports in a three-cell window; operator labels are irrelevant.

    A nonidentity string may use several factors in one velocity group. Its
    escape depends only on occupied group/site pairs, so these supports suffice
    for this finite calibration. The general theorem uses arbitrary diameter.
    """
    velocities = dict(VELOCITIES)
    if include_stationary:
        velocities['clock'] = 0
    slots = tuple(product(velocities, range(-1, 2)))
    pure = mixed = 0
    latest_exit = 0
    initial_shrink_example = None
    for bits in range(1, 1 << len(slots)):
        support = tuple(slot for j, slot in enumerate(slots) if bits >> j & 1)
        if len({g for g, _ in support}) == 1:
            pure += 1
            assert all(diameter(orbit_support(support, velocities, n)) == diameter(support)
                       for n in (1, 5, 17))
        else:
            mixed += 1
            # Initial diameter is <= R=2. min velocity gap is 1. n>2R
            # ensures that a selected pair from distinct groups is >R apart.
            first = next(n for n in range(1, 6)
                         if diameter(orbit_support(support, velocities, n)) > 2)
            latest_exit = max(latest_exit, first)
            if diameter(orbit_support(support, velocities, 1)) < diameter(support):
                initial_shrink_example = support
    return dict(stationary_clock_included=include_stationary, total=pure+mixed,
                pure_group=pure, mixed_group=mixed, latest_exit=latest_exit,
                initial_shrink_example=initial_shrink_example)


def pauli_shift(word, cells):
    """Conjugation by a two-lane countershift, using the forward convention."""
    out = [0]*(2*cells)
    for i in range(cells):
        out[2*((i+1) % cells)] = word[2*i]
        out[2*((i-1) % cells)+1] = word[2*i+1]
    return tuple(out)


def local_words(cells):
    """All nonidentity Pauli strings supported in at most two adjacent cells."""
    words = set()
    for i in range(cells):
        axes = (2*i, 2*i+1, 2*((i+1) % cells), 2*((i+1) % cells)+1)
        for values in product(range(4), repeat=4):
            if not any(values):
                continue
            word = [0]*(2*cells)
            for axis, value in zip(axes, values):
                word[axis] = value
            words.add(tuple(word))
    return words


def local_commutant_audit(cells):
    words = local_words(cells)
    visited = set()
    survivor_orbits = []
    rejected = 0
    for initial in sorted(words):
        if initial in visited:
            continue
        orbit = []
        current = initial
        while current not in orbit:
            orbit.append(current)
            current = pauli_shift(current, cells)
        assert current == initial
        visited.update(orbit)
        if all(word in words for word in orbit):
            survivor_orbits.append(orbit)
        else:
            rejected += sum(word in words for word in orbit)
    survivors = [word for orbit in survivor_orbits for word in orbit]
    mixed_survivors = sum(any(word[::2]) and any(word[1::2]) for word in survivors)
    return dict(cells=cells, nonidentity_local_words=len(words),
                invariant_local_dimension=len(survivor_orbits),
                surviving_words=len(survivors), rejected_local_words=rejected,
                mixed_surviving_words=int(mixed_survivors))


def countershift_matrix(cells):
    dim = 2**(2*cells)
    out = np.zeros((dim, dim), complex)
    for source in range(dim):
        bits = [(source >> (2*cells-1-j)) & 1 for j in range(2*cells)]
        moved = pauli_shift(bits, cells)
        destination = sum(bit << (2*cells-1-j) for j, bit in enumerate(moved))
        out[destination, source] = 1
    return out


def exp_h(h):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*values))@vectors.conj().T


def finite_log_audit():
    """An exact finite-ring logarithm exists: this is NOT an infinite no-go test."""
    cells = 3
    u = countershift_matrix(cells)
    powers = [np.eye(len(u), dtype=complex)]
    for _ in range(1, cells):
        powers.append(powers[-1]@u)
    h = np.zeros_like(u)
    for k in range(cells):
        eigenvalue = np.exp(2j*np.pi*k/cells)
        projection = sum(eigenvalue**(-r)*powers[r] for r in range(cells))/cells
        h += -np.angle(eigenvalue)*projection
    paulis = (np.eye(2), np.array([[0, 1], [1, 0]]),
              np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    max_full_support = 0.
    witness = None
    # Enumerate all strings occupying all three cells. This is a finite
    # locality diagnostic of this chosen logarithm, not of every logarithm.
    for word in product(range(4), repeat=2*cells):
        if not all(word[2*i] or word[2*i+1] for i in range(cells)):
            continue
        operator = np.array([[1.]], complex)
        for label in word:
            operator = np.kron(operator, paulis[label])
        coefficient = abs(np.vdot(operator, h)/len(u))
        if coefficient > max_full_support+1e-14:
            max_full_support = float(coefficient)
            witness = word
    return dict(cells=cells, hilbert_dimension=len(u),
                hermiticity_error=float(np.linalg.norm(h-h.conj().T)),
                exponential_error=float(np.linalg.norm(exp_h(h)-u)),
                commutator_error=float(np.linalg.norm(h@u-u@h)),
                full_three_cell_coefficient=max_full_support,
                full_three_cell_witness=witness,
                finite_log_does_not_prove_infinite_local_generator=True)


@lru_cache(None)
def report():
    indices = {'right': Fraction(6**2), 'fast_left': Fraction(1, 3**2),
               'slow_left': Fraction(1, 4)}
    return dict(round=394,
                scope='Exact full-algebra transport W from round 393 on an infinite chain; uniformly bounded finite-range time-independent generators only. Not the full computing update, a finite-window approximation, or a prepared-clock code space.',
                squared_content_map={track: content_power(track, 0, 2) for track in TRACKS},
                velocity_group_dimensions={'right': 6, 'fast_left': 3, 'slow_left': 4},
                squared_indices={g: str(value) for g, value in indices.items()},
                total_squared_index=str(indices['right']*indices['fast_left']*indices['slow_left']),
                support_orbits=[support_audit(), support_audit(True)],
                finite_ring_local_commutants=[local_commutant_audit(n) for n in (5, 7)],
                finite_log=finite_log_audit(),
                static_finite_range_transport_generator_excluded_analytically=True,
                prepared_returning_clock_excluded=False,
                full_computing_rule_static_generator_excluded=False,
                spatial_dimension_selected=False)


class Audit(unittest.TestCase):
    def test_squared_transport_and_group_flux(self):
        for site in range(-4, 5):
            for track in TRACKS:
                self.assertEqual(content_power(track, site, 2),
                                 (track, site+VELOCITIES[GROUPS[track]]))
        data = report()
        self.assertEqual(data['total_squared_index'], '1')
        self.assertEqual(set(data['squared_indices'].values()), {'36', '1/9', '1/4'})

    def test_unbounded_mixed_support_mechanism(self):
        data = report()['support_orbits'][0]
        self.assertEqual((data['total'], data['pure_group'], data['mixed_group']), (511, 21, 490))
        self.assertLessEqual(data['latest_exit'], 5)
        self.assertIsNotNone(data['initial_shrink_example'])
        # A simple two-factor string has diameter 4n, with no periodic boundary.
        self.assertEqual(diameter(orbit_support((('right', 0), ('fast_left', 0)),
                                               VELOCITIES, 1000)), 4000)

    def test_full_algebra_identity_clock_does_not_allow_mixed_terms(self):
        data = report()['support_orbits'][1]
        self.assertEqual((data['total'], data['pure_group'], data['mixed_group']), (4095, 28, 4067))
        self.assertLessEqual(data['latest_exit'], 5)
        self.assertFalse(report()['prepared_returning_clock_excluded'])

    def test_finite_ring_commutant_has_only_same_velocity_terms(self):
        for data in report()['finite_ring_local_commutants']:
            self.assertEqual(data['invariant_local_dimension'], 24)
            self.assertEqual(data['mixed_surviving_words'], 0)
            self.assertEqual(data['surviving_words'], 24*data['cells'])
            self.assertEqual(data['surviving_words']+data['rejected_local_words'],
                             data['nonidentity_local_words'])

    def test_finite_log_exists_and_has_long_terms(self):
        data = report()['finite_log']
        for name in ('hermiticity_error', 'exponential_error', 'commutator_error'):
            self.assertLess(data[name], 2e-12)
        self.assertGreater(data['full_three_cell_coefficient'], .01)
        self.assertFalse(report()['full_computing_rule_static_generator_excluded'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
