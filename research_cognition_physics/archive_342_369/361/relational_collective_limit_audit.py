"""Round 361: site-uniform collective limit on a specified weighted graph.

N counts qubits per site; M counts graph sites. Bounds concern the specified
product-state family and collective means, not arbitrary-input channels.
"""
from math import ceil, comb, factorial
import unittest
import numpy as np
from growing_stream_audit import main
from collective_poisson_limit_audit import coherent, collective_full, tensor_qubits


EZ = np.array([0., 0., 1.])


def graph(kind, sites, budget=.5):
    if sites < 1 or budget < 0 or (kind == "cycle" and sites < 3):
        raise ValueError("Need sites>=1, budget>=0 and at least three sites for a cycle.")
    j = np.zeros((sites, sites))
    if kind == "path":
        pairs = [(i, i+1) for i in range(sites-1)]
        weight = budget/2 if sites > 2 else budget
    elif kind == "cycle":
        pairs = [(i, (i+1) % sites) for i in range(sites)]
        weight = budget/2
    elif kind == "star":
        pairs = [(0, i) for i in range(1, sites)]
        weight = budget/max(1, sites-1)
    elif kind == "complete":
        pairs = [(i, k) for i in range(sites) for k in range(i+1, sites)]
        weight = budget/max(1, sites-1)
    else:
        raise ValueError("Unknown graph.")
    for i, k in pairs:
        j[i, k] = j[k, i] = weight
    return j


def initial_vectors(sites):
    v = np.arange(sites)
    z = .3*np.cos(.7*v)
    phase = .2*np.sin(.5*v)-.1
    radius = np.sqrt(1-z*z)
    return np.column_stack((radius*np.cos(phase), radius*np.sin(phase), z))


def product_state(n, vectors):
    value = np.array(1.+0j)
    for vector in vectors:
        site = coherent(n, vector[2], np.arctan2(vector[1], vector[0]))
        value = np.multiply.outer(value, site)
    return value


def spin_action(state, site, component):
    n = state.shape[0]-1
    moved = np.moveaxis(state, site, 0)
    extension = (1,)*(state.ndim-1)
    if component == 2:
        out = ((2*np.arange(n+1)-n)/n).reshape((n+1,)+extension)*moved
    else:
        hop = (np.sqrt(np.arange(1, n+1)*np.arange(n, 0, -1))/n)
        hop = hop.reshape((n,)+extension)
        out = np.zeros_like(moved)
        lower, upper = (1., 1.) if component == 0 else (-1j, 1j)
        out[1:] += lower*hop*moved[:-1]
        out[:-1] += upper*hop*moved[1:]
    return np.moveaxis(out, 0, site)


def h_action(state, omega, j):
    n = state.shape[0]-1
    out = np.zeros_like(state)
    for v in range(len(omega)):
        out += n*omega[v]*spin_action(state, v, 0)
        for w in range(v+1, len(omega)):
            if j[v, w] != 0:
                out += n*j[v, w]*spin_action(spin_action(state, w, 2), v, 2)
    return out


def propagate(state, time, omega, j, scaled_max=.5, order=18):
    n = state.shape[0]-1
    norm_upper = n*(sum(abs(omega))+np.triu(abs(j), 1).sum())
    steps = max(1, int(ceil(abs(time)*norm_upper/scaled_max)))
    dt = time/steps
    x = abs(dt)*norm_upper
    local_error = np.exp(x)*x**(order+1)/factorial(order+1)
    total_error = steps*local_error*np.exp(steps*local_error)
    value = state.copy()
    for _ in range(steps):
        term = value.copy()
        update = value.copy()
        for degree in range(1, order+1):
            term = (-1j*dt/degree)*h_action(term, omega, j)
            update += term
        value = update
    return value, {"substeps": steps, "norm_upper": float(norm_upper),
                   "exact_arithmetic_truncation_upper": float(total_error),
                   "norm_error": float(abs(np.vdot(value, value).real-1))}


def classical_rhs(value, omega, j):
    field = np.zeros_like(value)
    field[:, 0] = omega
    field[:, 2] = j @ value[:, 2]
    return 2*np.cross(field, value)


def classical(initial, time, omega, j, max_step=.002):
    value = initial.copy()
    steps = max(1, int(ceil(abs(time)/max_step)))
    dt = time/steps
    for _ in range(steps):
        a = classical_rhs(value, omega, j)
        b = classical_rhs(value+dt*a/2, omega, j)
        c = classical_rhs(value+dt*b/2, omega, j)
        d = classical_rhs(value+dt*c, omega, j)
        value += dt*(a+2*b+2*c+d)/6
    return value


def classical_energy(value, omega, j):
    return float(omega @ value[:, 0]+.5*value[:, 2] @ j @ value[:, 2])


def quantum_moments(state):
    means, variances = [], []
    for v in range(state.ndim):
        mean, var = [], []
        for component in range(3):
            applied = spin_action(state, v, component)
            m = np.vdot(state, applied).real
            mean.append(m)
            var.append(np.vdot(applied, applied).real-m*m)
        means.append(mean)
        variances.append(var)
    return np.array(means), np.array(variances)


def exp_symmetric(matrix):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors*np.exp(values)) @ vectors.T


def budgets(n, time, j):
    a = abs(j)
    incident = a.sum(axis=1)
    kappa = max(incident, default=0.)
    b = 2*(np.diag(incident)+a)
    local = exp_symmetric(time*b) @ np.full(len(j), 2/n)
    return {"local_mean_square": local,
            "uniform_mean_square": float(2*np.exp(4*kappa*time)/n),
            "uniform_mean_bias": float(2*np.exp(2*kappa*time)*np.expm1(2*kappa*time)/n),
            "kappa": float(kappa), "comparison_matrix": b}


def fluctuation_derivatives(state, target, omega, j):
    hstate = h_action(state, omega, j)
    target_dot = classical_rhs(target, omega, j)
    direct = np.zeros(state.ndim)
    cancelled = np.zeros(state.ndim)
    d = np.zeros(state.ndim)
    for v in range(state.ndim):
        deltas = [spin_action(state, v, i)-target[v, i]*state for i in range(3)]
        d[v] = sum(np.vdot(x, x).real for x in deltas)
        for i in range(3):
            delta_h = spin_action(hstate, v, i)-target[v, i]*hstate
            direct[v] += 2*np.vdot(deltas[i], -1j*delta_h).real
            direct[v] -= 2*target_dot[v, i]*np.vdot(state, deltas[i]).real
        direction = np.cross(EZ, target[v])
        transverse = sum(direction[i]*deltas[i] for i in range(3))
        for w in range(state.ndim):
            dz = spin_action(state, w, 2)-target[w, 2]*state
            cancelled[v] += 4*j[v, w]*np.vdot(dz, transverse).real
    return direct, cancelled, d


def full_tensor_means(n, initial, time, omega, j):
    size = 2**n
    single = [2*op/n for op in collective_full(n)]
    sites = len(initial)

    def embedded(op, position):
        out = np.array([[1.+0j]])
        for v in range(sites):
            out = np.kron(out, op if v == position else np.eye(size))
        return out

    operators = [[embedded(single[i], v) for i in range(3)] for v in range(sites)]
    h = np.zeros((size**sites, size**sites), complex)
    psi = np.array([1.+0j])
    for v in range(sites):
        psi = np.kron(psi, tensor_qubits(
            n, initial[v, 2], np.arctan2(initial[v, 1], initial[v, 0])))
        h += n*omega[v]*operators[v][0]
        for w in range(v+1, sites):
            h += n*j[v, w]*operators[v][2] @ operators[w][2]
    energies, basis = np.linalg.eigh(h)
    psi = basis @ (np.exp(-1j*time*energies)*(basis.conjugate().T @ psi))
    return np.array([[np.vdot(psi, op @ psi).real for op in row] for row in operators])


def tail_series(value, distance):
    term = value**distance/factorial(distance)
    result = term
    for degree in range(distance+1, distance+150):
        term *= value/degree
        result += term
        if term < 1e-17*max(result, 1e-300):
            break
    return float(result)


def initial_perturbation(initial, site=0, angle=.8):
    changed = initial.copy()
    c, s = np.cos(angle), np.sin(angle)
    changed[site] = np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]]) @ changed[site]
    return changed


def propagation_case(sites=7, time=.8):
    j = graph("path", sites)
    omega = np.full(sites, .6)
    start = initial_vectors(sites)
    changed = initial_perturbation(start)
    base = classical(start, time, omega, j)
    other = classical(changed, time, omega, j)
    base_fine = classical(start, time, omega, j, max_step=.001)
    other_fine = classical(changed, time, omega, j, max_step=.001)
    delta0 = np.linalg.norm(start-changed, axis=1)
    comparison = exp_symmetric(2*time*abs(j)) @ delta0
    kappa = max(abs(j).sum(axis=1))
    return {"M_sites": sites, "time": time,
            "initial_difference_at_source": float(delta0[0]),
            "actual_classical_response": np.linalg.norm(other-base, axis=1).tolist(),
            "response_vector_refinement_difference": np.linalg.norm(
                (other-base)-(other_fine-base_fine), axis=1).tolist(),
            "weighted_path_bound": comparison.tolist(),
            "row_sum_factorial_bound": [delta0[0]*tail_series(2*kappa*time, r)
                                         for r in range(sites)]}


def quantum_case(n, kind="path", sites=3, time=.8):
    j = graph(kind, sites)
    omega = .55*np.cos(.6*np.arange(sites))
    initial = initial_vectors(sites)
    state = product_state(n, initial)
    final, info = propagate(state, time, omega, j)
    means, variance = quantum_moments(final)
    target = classical(initial, time, omega, j)
    bias = np.linalg.norm(means-target, axis=1)
    d = variance.sum(axis=1)+bias*bias
    bound = budgets(n, time, j)
    energy_error = abs(np.vdot(final, h_action(final, omega, j)).real-
                       np.vdot(state, h_action(state, omega, j)).real)/n
    return {"N_per_site": n, "M_sites": sites, "graph": kind, "time": time,
            "incident_absolute_weight_max": bound["kappa"],
            "quantum_means": means.tolist(), "classical_means": target.tolist(),
            "mean_bias_each_site": bias.tolist(), "mean_square_each_site": d.tolist(),
            "local_mean_square_bounds": bound["local_mean_square"].tolist(),
            "uniform_mean_square_bound": bound["uniform_mean_square"],
            "uniform_mean_bias_bound": bound["uniform_mean_bias"],
            "energy_per_N_error": float(energy_error), "numerical_action": info}


def binomial_single_site_failure(n, epsilon=.5):
    failures = sum(comb(n, k) for k in range(n+1) if abs(2*k/n-1) > epsilon)
    return failures/(2**n)


def any_site_failure(n, sites, epsilon=.5):
    p = binomial_single_site_failure(n, epsilon)
    if p == 1:
        return 1.
    return float(-np.expm1(sites*np.log1p(-p)))


def report():
    growth = []
    for n in (16, 32, 64):
        p = binomial_single_site_failure(n)
        sites = int(ceil(1/p))
        growth.append({"N_per_site": n, "single_site_failure": p,
                       "M_sites_ceil_inverse_probability": sites,
                       "probability_any_site_fails": any_site_failure(n, sites)})
    return {
        "round": 361,
        "scope": {
            "proved": ["Specified graph collective-spin dynamics has sitewise mean-square comparison D<=exp(Bt)D0.",
                       "Bounded incident absolute weights give site-uniform O(1/N) concentration and mean bias independent of M at fixed times.",
                       "Classical initial-data dependence has a weighted graph-path factorial tail; this is not a strict causal cone.",
                       "Single-site convergence does not imply vanishing probability of any failure as M also grows."],
            "inputs": ["Undirected weighted relation graph J with zero diagonal",
                       "N qubits per graph site and uniform incident absolute-weight bound",
                       "Fixed transverse on-site fields, scaled all-to-all couplings between adjacent cells",
                       "Specified identical pure qubits inside each site and product across sites"],
            "not_proved": ["Arbitrary-input diamond convergence", "Full many-body product-state convergence",
                           "Selection of graph or physical dimension", "Lorentz symmetry, metric, HDA or Einstein dynamics",
                           "General existence or uniqueness of an infinite-volume quantum limit"]},
        "quantum_scaling": [quantum_case(n) for n in (2, 4, 8, 16)],
        "same_incident_budget_different_graphs": [
            quantum_case(4, kind, 4, .5) for kind in ("path", "cycle", "star", "complete")],
        "classical_graph_propagation": propagation_case(),
        "two_limits_exact_initial_measurement": growth,
        "fixed_N16_increasing_M": [{"M_sites": m, "any_site_failure": any_site_failure(16, m)}
                                  for m in (1, 10, 100, 1000)],
        "formulae": {"A": "|J|", "B": "2*(diag(A*1)+A)",
                     "mean_square_uniform": "2*exp(4*kappa*t)/N",
                     "mean_bias_uniform": "2*(exp(4*kappa*t)-exp(2*kappa*t))/N",
                     "classical_dependence": "r(t)<=exp(2*A*t)*r(0)",
                     "all_sites_Chebyshev_sufficient": "min(1,2*M*exp(4*kappa*t)/(N*epsilon^2))"},
        "resource_account": {"physical_qubits": "M*N",
                             "full_dimension": "2^(M*N)", "symmetric_dimension": "(N+1)^M",
                             "microscopic_cross_pairs": "|E|*N^2",
                             "pair_coupling": "J_vw/N",
                             "norm_upper_per_qubit": "Omega_max+kappa/2"},
        "source": ["https://arxiv.org/html/2403.17163", "https://arxiv.org/html/0902.0025"],
        "numerical_scope": ("Symmetric-sector evolution is exact within the invariant sector except "
                            "Taylor truncation and floating-point error. Taylor bounds certify exact "
                            "arithmetic only; full small tensors, RK4 refinement, energy and norm "
                            "checks provide independent numerical cross-checks."),
    }


class Checks(unittest.TestCase):
    def test_01_full_quantum_tensor_against_symmetric_graph_dynamics(self):
        for n, sites in ((1, 3), (2, 3), (1, 4)):
            initial = initial_vectors(sites)
            j = graph("path", sites)
            omega = .55*np.cos(.6*np.arange(sites))
            final, _ = propagate(product_state(n, initial), .6, omega, j)
            actual = quantum_moments(final)[0]
            expected = full_tensor_means(n, initial, .6, omega, j)
            np.testing.assert_allclose(actual, expected, atol=4e-14)

    def test_02_initial_sitewise_budget_does_not_grow_with_M(self):
        for n, sites in ((2, 2), (3, 4), (2, 6)):
            initial = initial_vectors(sites)
            mean, var = quantum_moments(product_state(n, initial))
            np.testing.assert_allclose(mean, initial, atol=2e-14)
            np.testing.assert_allclose(var.sum(axis=1), 2/n, atol=3e-14)

    def test_03_exact_each_site_fluctuation_cancellation_on_correlated_states(self):
        rng = np.random.default_rng(361)
        psi = rng.normal(size=(3, 3, 3))+1j*rng.normal(size=(3, 3, 3))
        psi /= np.linalg.norm(psi)
        j = graph("path", 3)
        j[1, 2] = j[2, 1] = -.2
        omega = np.array([.7, -.4, .2])
        direct, cancelled, d = fluctuation_derivatives(psi, initial_vectors(3), omega, j)
        np.testing.assert_allclose(direct, cancelled, atol=3e-14)
        incident = abs(j).sum(axis=1)
        comparison = 2*(np.diag(incident)+abs(j)) @ d
        self.assertTrue(np.all(abs(direct) <= comparison+3e-14))

    def test_04_local_comparison_and_uniform_bias_bounds(self):
        for n in (2, 5, 9):
            row = quantum_case(n)
            self.assertTrue(np.all(np.array(row["mean_square_each_site"]) <=
                                   np.array(row["local_mean_square_bounds"])+1e-12))
            self.assertLessEqual(max(row["mean_bias_each_site"]),
                                 row["uniform_mean_bias_bound"]+1e-12)

    def test_05_isolated_cell_has_unchanged_variance_and_zero_classical_bias(self):
        j = graph("path", 3)
        j[1, 2] = j[2, 1] = 0.
        omega = np.array([.4, -.3, .7])
        initial = initial_vectors(3)
        final, _ = propagate(product_state(5, initial), .8, omega, j)
        means, var = quantum_moments(final)
        target = classical(initial, .8, omega, j)
        np.testing.assert_allclose(means[2], target[2], atol=2e-12)
        self.assertAlmostEqual(var[2].sum(), 2/5, places=12)
        self.assertAlmostEqual(budgets(5, .8, j)["local_mean_square"][2], 2/5, places=13)

    def test_06_classical_length_energy_and_time_refinement(self):
        j = graph("cycle", 5)
        omega = .5*np.cos(np.arange(5))
        initial = initial_vectors(5)
        final = classical(initial, .8, omega, j)
        refined = classical(initial, .8, omega, j, max_step=.001)
        np.testing.assert_allclose(final, refined, atol=3e-12)
        np.testing.assert_allclose(np.linalg.norm(final, axis=1), 1., atol=3e-12)
        self.assertAlmostEqual(classical_energy(final, omega, j),
                               classical_energy(initial, omega, j), places=11)

    def test_07_weighted_classical_path_bound(self):
        row = propagation_case()
        actual = np.array(row["actual_classical_response"])
        matrix_bound = np.array(row["weighted_path_bound"])
        tail_bound = np.array(row["row_sum_factorial_bound"])
        self.assertTrue(np.all(actual <= matrix_bound+2e-12))
        self.assertTrue(np.all(matrix_bound <= tail_bound+2e-12))
        self.assertGreater(actual[2], 1e-7)

    def test_08_weighted_powers_vanish_before_graph_distance(self):
        a = graph("path", 7)
        for distance in range(1, 7):
            for power in range(distance):
                self.assertEqual(np.linalg.matrix_power(a, power)[distance, 0], 0.)
            self.assertGreater(np.linalg.matrix_power(a, distance)[distance, 0], 0.)

    def test_09_quantum_response_is_bounded_by_classical_response_plus_bias(self):
        n, sites, time = 6, 3, .6
        j = graph("path", sites)
        omega = np.full(sites, .6)
        initial = initial_vectors(sites)
        changed = initial_perturbation(initial)
        responses = []
        for start in (initial, changed):
            final, _ = propagate(product_state(n, start), time, omega, j)
            responses.append(quantum_moments(final)[0])
        actual = np.linalg.norm(responses[0]-responses[1], axis=1)
        delta0 = np.linalg.norm(initial-changed, axis=1)
        bound = exp_symmetric(2*time*abs(j)) @ delta0
        bound += 2*budgets(n, time, j)["uniform_mean_bias"]
        self.assertTrue(np.all(actual <= bound+1e-12))

    def test_10_different_graphs_retain_the_same_incident_budget(self):
        means = []
        for kind in ("path", "cycle", "star", "complete"):
            row = quantum_case(3, kind, 4, .5)
            self.assertAlmostEqual(row["incident_absolute_weight_max"], .5)
            means.append(np.array(row["quantum_means"]))
        self.assertGreater(np.linalg.norm(means[0]-means[1]), .01)
        self.assertGreater(np.linalg.norm(means[2]-means[3]), .01)

    def test_11_negative_weights_enter_the_bound_by_absolute_value(self):
        initial = initial_vectors(3)
        j = graph("path", 3)
        j[1, 2] = j[2, 1] = -.25
        omega = np.array([.7, -.4, .2])
        final, _ = propagate(product_state(4, initial), .7, omega, j)
        mean, var = quantum_moments(final)
        target = classical(initial, .7, omega, j)
        actual = var.sum(axis=1)+np.sum((mean-target)**2, axis=1)
        self.assertTrue(np.all(actual <= budgets(4, .7, j)["local_mean_square"]+1e-12))

    def test_12_quantum_norm_energy_reverse_and_Taylor_refinement(self):
        initial = product_state(5, initial_vectors(3))
        j = graph("path", 3)
        omega = np.array([.5, -.3, .7])
        final, info = propagate(initial, .8, omega, j)
        finer, _ = propagate(initial, .8, omega, j, scaled_max=.25, order=20)
        returned, _ = propagate(final, -.8, omega, j)
        np.testing.assert_allclose(final, finer, atol=3e-14)
        np.testing.assert_allclose(returned, initial, atol=3e-14)
        self.assertLess(info["exact_arithmetic_truncation_upper"], 1e-18)
        self.assertLess(info["norm_error"], 1e-13)
        self.assertAlmostEqual(np.vdot(initial, h_action(initial, omega, j)).real,
                               np.vdot(final, h_action(final, omega, j)).real, places=12)

    def test_13_exact_binomial_initial_counterexample(self):
        self.assertEqual(binomial_single_site_failure(16), 1394/65536)
        probabilities = [any_site_failure(16, sites) for sites in (1, 10, 100, 1000)]
        self.assertTrue(np.all(np.diff(probabilities) > 0))
        self.assertGreater(probabilities[-1], .999)
        for n in (16, 32, 64):
            p = binomial_single_site_failure(n)
            sites = int(ceil(1/p))
            self.assertGreater(any_site_failure(n, sites), 1-np.exp(-1))

    def test_14_uniform_comparison_constant_is_independent_of_sites(self):
        for sites in (4, 16, 64):
            for kind in ("path", "star", "complete"):
                bound = budgets(32, .8, graph(kind, sites))
                self.assertAlmostEqual(bound["kappa"], .5, places=14)
                self.assertAlmostEqual(bound["uniform_mean_square"], 2*np.exp(1.6)/32)
                self.assertTrue(np.all(bound["local_mean_square"] <=
                                       bound["uniform_mean_square"]+1e-12))

    def test_15_microscopic_resource_norm_budget(self):
        n, sites = 2, 3
        j = graph("path", sites)
        omega = np.array([.5, -.3, .7])
        dimension = (n+1)**sites
        basis = np.eye(dimension, dtype=complex)
        h = np.column_stack([h_action(basis[:, i].reshape((n+1,)*sites), omega, j).ravel()
                             for i in range(dimension)])
        np.testing.assert_allclose(h, h.conjugate().T, atol=1e-15)
        upper = n*(sum(abs(omega))+np.triu(abs(j), 1).sum())
        self.assertLessEqual(np.linalg.norm(h, 2), upper+1e-13)
        self.assertLessEqual(upper/(n*sites), max(abs(omega))+.5*max(abs(j).sum(axis=1)))


if __name__ == "__main__":
    main(__name__, "relational_collective_limit_audit", report)
