"""Round 393: information-flow index and an explicit compensating transport.

Infinite-chain conclusions use the cited GNVW/RWW theorems. Finite rings below
check only wire identities, reference preservation, and local pulse endpoints.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('program_transport_index_audit_results.json')
TOL = 2e-12


def support_dimensions(width=1, origin=0):
    """Heisenberg preimages of two adjacent grouped cells, no periodic wrap."""
    atoms = [('t', i-1) for i in range(origin, origin+2*width)]
    atoms += [('q', j+1) for j in range(2*origin, 2*(origin+2*width))]
    left = right = 1
    for track, index in atoms:
        cell = index if track == 't' else index//2
        dim = 3 if track == 't' else 2
        if origin-width <= cell < origin+width:
            left *= dim
        elif origin+width <= cell < origin+3*width:
            right *= dim
        else:
            raise AssertionError('support escaped the nearest-neighbor blocks')
    d = 12**width
    assert left*right == d*d
    assert Fraction(left, d) == Fraction(d, right)
    return dict(width=width, cell_dimension=d, left_matrix_size=left,
                right_matrix_size=right, index=str(Fraction(left, d)))


def swap_matrix(d):
    out = np.zeros((d*d, d*d), complex)
    for a in range(d):
        for b in range(d):
            out[b*d+a, a*d+b] = 1
    return out


def exp_h(h):
    eig, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*eig))@vectors.conj().T


def source_gate():
    """The 12-dimensional controlled I/SWAP/G cell gate, not a compiler."""
    g = np.eye(4, dtype=complex)
    g[2:, 2:] = np.array([[1, -1], [1, 1]])/np.sqrt(2)
    u = np.zeros((12, 12), complex)
    u[:4, :4], u[4:8, 4:8], u[8:, 8:] = np.eye(4), swap_matrix(2), g
    return u


def source_hamiltonian():
    h = np.zeros((12, 12), complex)
    h[4:8, 4:8] = np.pi/2*(np.eye(4)-swap_matrix(2))
    h[10:12, 10:12] = np.pi/4*np.array([[0, -1j], [1j, 0]])
    return h


def entropy(rho):
    p = np.linalg.eigvalsh((rho+rho.conj().T)/2)
    p = p[p > 1e-13]
    return float(-np.sum(p*np.log(p)))


def reduced_pure(psi, dims, keep):
    order = list(keep)+[i for i in range(len(dims)) if i not in keep]
    a = psi.reshape(dims).transpose(order).reshape(math.prod(dims[i] for i in keep), -1)
    return a@a.conj().T


def crossing_state():
    # Order: left reference (3), right reference (2), left out (2), right out (3).
    psi = np.zeros((3, 2, 2, 3), complex)
    for t in range(3):
        for q in range(2):
            psi[t, q, q, t] = 1/np.sqrt(6)
    return psi.reshape(-1)


def entropy_flow(psi):
    dims = (3, 2, 2, 3)
    s = lambda k: entropy(reduced_pure(psi, dims, k))
    right = s((0,))+s((3,))-s((0, 3))
    left = s((1,))+s((2,))-s((1, 2))
    return dict(right_mutual_information=right, left_mutual_information=left,
                log_index=(right-left)/2)


def apply_gate(psi, dims, axes, gate):
    order = list(axes)+[i for i in range(len(dims)) if i not in axes]
    moved = psi.reshape(dims).transpose(order)
    out = (gate@moved.reshape(gate.shape[0], -1)).reshape(tuple(dims[i] for i in order))
    return out.transpose(np.argsort(order)).reshape(-1)


def layout(n, reference=1):
    # Each cell: original (t,a,b), additional (s,e), reference is idle.
    names = [(track, i) for i in range(n) for track in ('t', 'a', 'b', 's', 'e')]
    dims = tuple(3 if track in ('t', 's') else 2 for track, _ in names)
    return names, dims+(reference,)


def desired_sources(n):
    """For each output factor, identify its old input factor."""
    names, _ = layout(n)
    where = {name: i for i, name in enumerate(names)}
    source = []
    for track, i in names:
        previous = {'t': ('t', (i-1) % n), 'a': ('b', i),
                    'b': ('a', (i+1) % n), 's': ('s', (i+1) % n),
                    'e': ('e', (i-1) % n)}[track]
        source.append(where[previous])
    return source


def swap_layers(n):
    names, _ = layout(n)
    at = {name: i for i, name in enumerate(names)}
    # L1: t<->s and a<->b. L2: b<->e. K: s_i<->t_(i+1), b_i<->e_(i+1).
    first = [(at['t', i], at['s', i]) for i in range(n)]
    first += [(at['a', i], at['b', i]) for i in range(n)]
    second = [(at['b', i], at['e', i]) for i in range(n)]
    third = [(at['s', i], at['t', (i+1) % n]) for i in range(n)]
    third += [(at['b', i], at['e', (i+1) % n]) for i in range(n)]
    return (first, second, third)


def layer_sources(n):
    source = list(range(5*n))
    for layer in swap_layers(n):
        for a, b in layer:
            source[a], source[b] = source[b], source[a]
    return source


def route(psi, dims, source):
    return psi.reshape(dims).transpose(source+[len(source)]).reshape(-1)


def circuit(psi, n, reference=1, pulses=False):
    names, dims = layout(n, reference)
    out = psi.copy()
    for layer in swap_layers(n):
        for a, b in layer:
            if pulses:
                d = dims[a]
                out = apply_gate(out, dims, (a, b), exp_h(np.pi/2*(np.eye(d*d)-swap_matrix(d))))
            else:
                out = out.reshape(dims).swapaxes(a, b).reshape(-1)
    c = exp_h(source_hamiltonian()) if pulses else source_gate()
    for i in range(n):
        out = apply_gate(out, dims, (5*i, 5*i+1, 5*i+2), c)
    return out


def target(psi, n, reference=1):
    _, dims = layout(n, reference)
    out = route(psi, dims, desired_sources(n))
    for i in range(n):
        out = apply_gate(out, dims, (5*i, 5*i+1, 5*i+2), source_gate())
    return out


def random_state(d, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=d)+1j*rng.normal(size=d)
    return v/np.linalg.norm(v)


def vacuum_reference_input(n=2):
    _, dims = layout(n, 12**n)
    psi = np.zeros(dims, complex)
    for flat, values in enumerate(np.ndindex(*((3, 2, 2)*n))):
        index = []
        for i in range(n):
            index += list(values[3*i:3*i+3])+[0, 0]
        psi[tuple(index+[flat])] = 12**(-n/2)
    return psi.reshape(-1), dims


def report():
    _, dims = layout(2, 3)
    psi = random_state(math.prod(dims), 393)
    pulse_error = float(np.linalg.norm(circuit(psi, 2, 3, True)-target(psi, 2, 3)))
    reference, rdims = vacuum_reference_input()
    actual = circuit(reference, 2, rdims[-1])
    expected = target(reference, 2, rdims[-1])
    aux = (3, 4, 8, 9)
    aux_rho = reduced_pure(actual, rdims, aux)
    return dict(round=393,
                scope='Conditional audit of a specified 1D infinite-chain QCA; no dimension selection, no autonomous time-independent Hamiltonian construction, no full physically-universal compiler claim.',
                original_supports=[support_dimensions(w) for w in (1, 2, 3, 5)],
                boundary_choi=entropy_flow(crossing_state()),
                counterflow=dict(index='2/3', original_index='3/2', combined_index='1',
                                 minimum_uniform_auxiliary_cell_dimension_for_product_completion=6,
                                 original_cell_dimension=12, enlarged_cell_dimension=72),
                wire_identity_rings=[dict(cells=n, exact=layer_sources(n) == desired_sources(n)) for n in (2, 3, 4, 7)],
                pulse_reference_vector_error=pulse_error,
                maximally_entangled_reference_vector_error=float(np.linalg.norm(actual-expected)),
                blank_auxiliary_return_probability=float(aux_rho[0, 0].real),
                pulse_stages=4, largest_single_swap_generator_norm=float(np.pi),
                constant_autonomous_hamiltonian_derived=False,
                finite_ring_used_to_prove_infinite_chain_obstruction=False)


class Checks(unittest.TestCase):
    def test_01_source_gate_and_independent_generator(self):
        u = source_gate()
        np.testing.assert_allclose(u.conj().T@u, np.eye(12), atol=TOL)
        np.testing.assert_allclose(exp_h(source_hamiltonian()), u, atol=TOL)

    def test_02_support_algebra_sizes_and_regrouping(self):
        for w in range(1, 9):
            for origin in (-17, -1, 0, 5):
                s = support_dimensions(w, origin)
                self.assertEqual(s['index'], '3/2')
                self.assertEqual(s['left_matrix_size'], 3**(w+1)*2**(2*w-1))
        self.assertEqual(support_dimensions()['right_matrix_size'], 8)

    def test_03_independent_choi_entropy_direction(self):
        f = entropy_flow(crossing_state())
        self.assertAlmostEqual(f['right_mutual_information'], 2*np.log(3))
        self.assertAlmostEqual(f['left_mutual_information'], 2*np.log(2))
        self.assertAlmostEqual(np.exp(f['log_index']), 1.5)

    def test_04_boundary_unitary_decoration_does_not_change_flow(self):
        for seed in (11, 21, 31):
            rng = np.random.default_rng(seed)
            x = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
            u = exp_h(x+x.conj().T)
            psi = apply_gate(crossing_state(), (3, 2, 2, 3), (2, 3), u)
            self.assertAlmostEqual(entropy_flow(psi)['log_index'], np.log(1.5))

    def test_05_counterflow_dimension_and_integer_divisibility(self):
        self.assertEqual(Fraction(3, 2)*Fraction(2, 3), 1)
        for d in range(1, 25):
            for group in (1, 2, 5):
                self.assertEqual(d**group % 2 == 0 and d**group % 3 == 0, d % 6 == 0)
        # Static auxiliary factors multiply the index by 1, not by 2/3.
        self.assertNotEqual(Fraction(3, 2)*1, 1)

    def test_06_exact_wires_and_disjoint_factor_pulses(self):
        for n in (2, 3, 4, 7):
            _, dims = layout(n)
            self.assertEqual(layer_sources(n), desired_sources(n))
            for layer in swap_layers(n):
                used = [a for pair in layer for a in pair]
                self.assertEqual(len(used), len(set(used)))
                self.assertTrue(all(dims[a] == dims[b] for a, b in layer))

    def test_07_arbitrary_complex_correlated_input(self):
        for n, r in ((2, 3), (3, 2)):
            _, dims = layout(n, r)
            psi = random_state(math.prod(dims), n+r)
            np.testing.assert_allclose(circuit(psi, n, r), target(psi, n, r), atol=TOL)

    def test_08_full_input_choi_and_blank_return(self):
        psi, dims = vacuum_reference_input()
        actual = circuit(psi, 2, dims[-1])
        np.testing.assert_allclose(actual, target(psi, 2, dims[-1]), atol=TOL)
        aux_rho = reduced_pure(actual, dims, (3, 4, 8, 9))
        self.assertAlmostEqual(aux_rho[0, 0].real, 1.)
        self.assertAlmostEqual(np.trace(aux_rho@aux_rho).real, 1.)

    def test_09_counterflow_is_not_identity_on_unknown_auxiliary(self):
        _, dims = layout(3, 2)
        psi = np.zeros(dims, complex)
        base = [0]*len(dims)
        psi[tuple(base)] = 1/np.sqrt(2)
        base[4] = base[-1] = 1  # e_0 entangled with idle reference.
        psi[tuple(base)] = 1/np.sqrt(2)
        moved = route(psi.reshape(-1), dims, desired_sources(3))
        # The Bell pair moves from auxiliary e_0 to e_1.
        local = reduced_pure(moved, dims, (9, len(dims)-1))
        phi = np.array([1, 0, 0, 1])/np.sqrt(2)
        self.assertAlmostEqual(float((phi@local@phi).real), 1.)
        self.assertAlmostEqual(reduced_pure(moved, dims, (4,))[0, 0].real, 1.)

    def test_10_bounded_local_pulse_endpoint_with_reference(self):
        for d in (2, 3):
            swap = swap_matrix(d)
            h = np.pi/2*(np.eye(d*d)-swap)
            np.testing.assert_allclose(exp_h(h), swap, atol=TOL)
            self.assertAlmostEqual(np.linalg.norm(h, 2), np.pi)
        _, dims = layout(2, 3)
        psi = random_state(math.prod(dims), 393)
        self.assertLess(np.linalg.norm(circuit(psi, 2, 3, True)-target(psi, 2, 3)), TOL)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    check = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not check.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = dict(run=check.testsRun, failures=len(check.failures), errors=len(check.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
