"""Round 374: transverse two-band nodes, degree, and controlled Weyl limits."""
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)


def pauli(v):
    return np.einsum("...i,ijk->...jk", v, SIGMA)


def unitary(h, t):
    eig, vec = np.linalg.eigh(h)
    return (vec*np.exp(-1j*t*eig)) @ vec.conjugate().T


def trace_distance(a, b):
    d = a-b
    return float(abs(np.linalg.eigvalsh((d+d.conjugate().T)/2)).sum()/2)


def lattice_h(k, mass=0.):
    k = np.asarray(k)
    x, y, z = np.moveaxis(k, -1, 0)
    return np.stack([np.sin(x), np.sin(y),
                     mass+np.cos(z)+2-np.cos(x)-np.cos(y)], axis=-1)


def lattice_jac(k):
    k = np.asarray(k)
    x, y, z = np.moveaxis(k, -1, 0)
    j = np.zeros(k.shape[:-1]+(3, 3))
    j[..., 0, 0], j[..., 1, 1] = np.cos(x), np.cos(y)
    j[..., 2, 0], j[..., 2, 1], j[..., 2, 2] = np.sin(x), np.sin(y), -np.sin(z)
    return j


def lattice_H(k, tilt=.2, mass=0.):
    return tilt*np.sin(np.asarray(k)[..., 0])[..., None, None]*I+pauli(lattice_h(k, mass))


def local_linear_H(q, node, tilt=.2):
    q = np.asarray(q)
    return tilt*q[..., 0, None, None]*I+pauli(q @ lattice_jac(node).T)


def exact_nodes(mass=0.):
    # Exhausts sin(kx)=sin(ky)=0 on the 3-torus, including other branches.
    nodes = []
    for x in (0., np.pi):
        for y in (0., np.pi):
            cosine = -mass-2+np.cos(x)+np.cos(y)
            if -1 <= cosine <= 1:
                angle = float(np.arccos(cosine))
                zs = [angle] if np.isclose(angle, 0) or np.isclose(angle, np.pi) else [angle, -angle]
                for z in zs:
                    k = np.array([x, y, z])
                    nodes.append((k, float(np.linalg.det(lattice_jac(k)))))
    return nodes


def sphere_grid(order=48):
    u, weight = np.polynomial.legendre.leggauss(order)
    phi = 2*np.pi*np.arange(2*order)/(2*order)
    uu, pp = np.meshgrid(u, phi, indexing="ij")
    root = np.sqrt(1-uu*uu)
    n = np.stack([root*np.cos(pp), root*np.sin(pp), uu], axis=-1)
    du = np.stack([-uu/root*np.cos(pp), -uu/root*np.sin(pp), np.ones_like(uu)], axis=-1)
    dp = np.stack([-root*np.sin(pp), root*np.cos(pp), np.zeros_like(uu)], axis=-1)
    weights = np.broadcast_to(weight[:, None]*2*np.pi/(2*order), uu.shape)
    # Outward orientation is (phi,u), since d_phi n cross d_u n = n.
    return n, du, dp, weights


def degree_flux(hfun, jfun, center, radius=.3, order=48):
    n, du, dp, weights = sphere_grid(order)
    points = np.asarray(center)+radius*n
    h, jac = hfun(points), jfun(points)
    hu = np.einsum("...ij,...j->...i", jac, radius*du)
    hp = np.einsum("...ij,...j->...i", jac, radius*dp)
    norms = np.linalg.norm(h, axis=-1)
    if np.min(norms) < 1e-10:
        raise ValueError("Enclosing surface must be gapped")
    density = np.einsum("...i,...i->...", h, np.cross(hp, hu))/norms**3
    return float(np.sum(weights*density)/(4*np.pi))


def double_h(k):
    x, y, z = np.moveaxis(np.asarray(k), -1, 0)
    return np.stack([x*x-y*y, 2*x*y, z], axis=-1)


def double_jac(k):
    k = np.asarray(k)
    x, y, _ = np.moveaxis(k, -1, 0)
    j = np.zeros(k.shape[:-1]+(3, 3))
    j[..., 0, 0], j[..., 0, 1] = 2*x, -2*y
    j[..., 1, 0], j[..., 1, 1], j[..., 2, 2] = 2*y, 2*x, 1
    return j


def block_diagonal(matrices):
    size = sum(len(m) for m in matrices)
    out = np.zeros((size, size), complex)
    cursor = 0
    for m in matrices:
        out[cursor:cursor+len(m), cursor:cursor+len(m)] = m
        cursor += len(m)
    return out


def linearization_case():
    node = np.array([0., 0., np.pi/2])
    rng = np.random.default_rng(374)
    momenta = rng.normal(size=(300, 3))
    momenta /= np.maximum(1., np.linalg.norm(momenta, axis=1))[:, None]
    momenta = np.vstack([momenta, np.eye(3), -np.eye(3)])
    tilt, time, cap = .2, 1.5, 1.
    constant = (abs(tilt)+np.sqrt(2))/2
    rows = []
    for spacing in (.4, .2, .1, .05):
        exact = lattice_H(node+spacing*momenta, tilt)/spacing
        linear = local_linear_H(momenta, node, tilt)
        herrors = [np.linalg.norm(a-b, 2) for a, b in zip(exact, linear)]
        uerrors = [np.linalg.norm(unitary(a, time)-unitary(b, time), 2)
                   for a, b in zip(exact, linear)]
        rows.append({"a": spacing, "K": cap, "time": time,
                     "sample_hamiltonian_error": float(max(herrors)),
                     "uniform_hamiltonian_bound": constant*spacing*cap**2,
                     "sample_unitary_error": float(max(uerrors)),
                     "uniform_half_diamond_bound": min(1., time*constant*spacing*cap**2)})
    # A coherent momentum superposition with a correlated, inactive R.
    selected = momenta[:4]
    spacing = .2
    u = block_diagonal([unitary(h, time) for h in lattice_H(node+spacing*selected, tilt)/spacing])
    v = block_diagonal([unitary(h, time) for h in local_linear_H(selected, node, tilt)])
    raw = rng.normal(size=24)+1j*rng.normal(size=24)
    psi = raw/np.linalg.norm(raw)
    exact_out, linear_out = np.kron(u, np.eye(3)) @ psi, np.kron(v, np.eye(3)) @ psi
    distance = trace_distance(np.outer(exact_out, exact_out.conj()),
                              np.outer(linear_out, linear_out.conj()))
    return {"taylor_constant": constant, "rows": rows, "reference_dimension": 3,
            "coherent_momentum_reference_distance": distance,
            "reference_case_half_diamond_bound": time*constant*spacing,
            "scope": "Every input and inactive reference on the fixed rescaled momentum band; "
                     "sample maxima do not replace the analytic supremum bound."}


def topological_case():
    rows = []
    for node, determinant in exact_nodes():
        flux = degree_flux(lattice_h, lattice_jac, node)
        rows.append({"node": node.tolist(), "determinant_Dh": determinant,
                     "chirality": int(np.sign(determinant)), "degree_quadrature": flux,
                     "upper_band_chern_for_A_i_u_du": -flux})
    perturb = np.array([.03, -.02, .04])
    radius = .3
    lower = radius-radius**2/np.sqrt(2)
    perturbed_degree = degree_flux(lambda k: lattice_h(k)+perturb,
                                  lattice_jac, np.array([0, 0, np.pi/2]), radius)
    multi_degree = degree_flux(double_h, double_jac, np.zeros(3), .5, order=96)
    return {"lattice_nodes": rows, "net_charge": sum(x["chirality"] for x in rows),
            "sphere_radius": radius, "analytic_unperturbed_boundary_gap_in_h": lower,
            "constant_vector_perturbation_norm": float(np.linalg.norm(perturb)),
            "perturbed_degree": perturbed_degree,
            "nontransverse_double_node_degree": multi_degree,
            "double_node_Dh_rank_at_origin": int(np.linalg.matrix_rank(double_jac(np.zeros(3))))}


def report():
    return {"round": 374,
            "scope": "Conditional smooth two-band crossing interface: codimension three, "
                     "degree protection, C2 linearization with arbitrary-reference band error. "
                     "Does not select physical space, node existence, filling, common metric, or GR.",
            "codimension": {"regular_zero_set_dimension": "d-3 when rank Dh=3",
                            "d2_mass_gap_at_m_0_3": .6,
                            "d4_model": "h=(q1,q2,q3), arbitrary q4: a nodal line",
                            "chemical_potential": "Crossing energy equals Fermi level only with an extra filling condition."},
            "topology": topological_case(), "controlled_linearization": linearization_case(),
            "periodic_mass_deformation": {"nodes_at_m_0": len(exact_nodes(0)),
                                         "nodes_at_m_0_3": len(exact_nodes(.3)),
                                         "nodes_at_m_1_2": len(exact_nodes(1.2)),
                                         "uniform_band_gap_lower_bound_m_1_2": .4},
            "extra_inputs": ["translation/momentum parameters and their dimension",
                             "C2 autonomous isolated two-band block",
                             "existence of a crossing and transverse rank-three Jacobian",
                             "near-node preparation and optional continuum clock/momentum scaling",
                             "chemical potential, tilt and anisotropic velocity matrix"],
            "sources": ["https://arxiv.org/html/1007.0016",
                        "https://arxiv.org/html/1105.5138",
                        "https://arxiv.org/pdf/1111.7309"]}


class Checks(unittest.TestCase):
    def test_01_general_hermitian_two_band_spectrum(self):
        rng = np.random.default_rng(1374)
        for _ in range(30):
            h = rng.normal(size=3)
            h0 = rng.normal()
            np.testing.assert_allclose(np.linalg.eigvalsh(h0*I+pauli(h)),
                                       [h0-np.linalg.norm(h), h0+np.linalg.norm(h)], atol=2e-14)

    def test_02_two_dimensional_mass_opens_gap(self):
        for x, y in [(0, 0), (.2, -.1), (1, 1)]:
            h = pauli([x, y, .3])
            np.testing.assert_allclose(np.linalg.eigvalsh(h),
                                       np.array([-1, 1])*np.sqrt(x*x+y*y+.09))

    def test_03_four_dimensional_transverse_zero_is_line(self):
        jac = np.column_stack([np.eye(3), np.zeros(3)])
        self.assertEqual(np.linalg.matrix_rank(jac), 3)
        for fourth in [-3., 0., 2.]:
            q = np.array([0, 0, 0, fourth])
            np.testing.assert_array_equal(jac @ q, np.zeros(3))
            np.testing.assert_allclose(np.linalg.eigvalsh(fourth*.2*I+pauli(jac @ q)),
                                       [fourth*.2, fourth*.2])

    def test_04_periodic_node_classification_and_opposite_charges(self):
        nodes = exact_nodes()
        self.assertEqual(len(nodes), 2)
        for k, determinant in nodes:
            np.testing.assert_allclose(lattice_h(k), np.zeros(3), atol=3e-16)
            self.assertAlmostEqual(abs(determinant), 1.)
        self.assertAlmostEqual(sum(d for _, d in nodes), 0.)

    def test_05_jacobian_against_finite_differences(self):
        k = np.array([.17, -.12, 1.48])
        eps = 1e-5
        numerical = np.column_stack([(lattice_h(k+eps*v)-lattice_h(k-eps*v))/(2*eps)
                                      for v in np.eye(3)])
        np.testing.assert_allclose(numerical, lattice_jac(k), atol=5e-10)

    def test_06_degree_matches_det_for_general_invertible_E(self):
        e = np.array([[1.2, .2, 0], [0, .9, .1], [.1, 0, 1.4]])
        for reflection in [np.eye(3), np.diag([1, 1, -1])]:
            a = reflection @ e
            flux = degree_flux(lambda q: q @ a.T,
                               lambda q: np.broadcast_to(a, np.asarray(q).shape[:-1]+(3, 3)),
                               np.zeros(3))
            self.assertAlmostEqual(flux, np.sign(np.linalg.det(a)), places=11)

    def test_07_periodic_node_degree_and_small_perturbation(self):
        row = topological_case()
        for node in row["lattice_nodes"]:
            self.assertAlmostEqual(node["degree_quadrature"], node["chirality"], places=11)
        self.assertLess(row["constant_vector_perturbation_norm"],
                        row["analytic_unperturbed_boundary_gap_in_h"])
        self.assertAlmostEqual(row["perturbed_degree"], -1., places=11)

    def test_08_upper_band_curvature_sign_from_projector(self):
        n, du, dp, _ = sphere_grid(6)
        for index in [(0, 0), (2, 4), (4, 8)]:
            p = (I+pauli(n[index]))/2
            pp, pu = pauli(dp[index])/2, pauli(du[index])/2
            curvature = (1j*np.trace(p @ (pp @ pu-pu @ pp))).real
            density = -.5*np.dot(n[index], np.cross(dp[index], du[index]))
            self.assertAlmostEqual(curvature, density, places=14)
            self.assertAlmostEqual(curvature, -.5, places=14)

    def test_09_degree_protection_does_not_imply_full_rank(self):
        self.assertEqual(np.linalg.matrix_rank(double_jac(np.zeros(3))), 1)
        self.assertAlmostEqual(degree_flux(double_h, double_jac, np.zeros(3), .5, 96), 2., places=9)
        eps = .04
        for sign in [-1, 1]:
            point = np.array([sign*np.sqrt(eps), 0, 0])
            np.testing.assert_allclose(double_h(point)-[eps, 0, 0], np.zeros(3), atol=1e-15)
            self.assertGreater(np.linalg.det(double_jac(point)), 0)

    def test_10_exact_C2_remainder_bound_and_scaled_unitaries(self):
        row = linearization_case()
        for sample in row["rows"]:
            self.assertLessEqual(sample["sample_hamiltonian_error"],
                                 sample["uniform_hamiltonian_bound"]+1e-14)
            self.assertLessEqual(sample["sample_unitary_error"],
                                 sample["uniform_half_diamond_bound"]+1e-14)

    def test_11_coherent_momentum_and_reference_channel(self):
        row = linearization_case()
        self.assertGreater(row["coherent_momentum_reference_distance"], .001)
        self.assertLessEqual(row["coherent_momentum_reference_distance"],
                             row["reference_case_half_diamond_bound"])

    def test_12_gap_after_periodic_node_annihilation(self):
        self.assertEqual(len(exact_nodes(.3)), 2)
        self.assertEqual(len(exact_nodes(1.2)), 0)
        for k in [np.zeros(3), np.array([0, 0, np.pi]), np.array([.3, .8, 2.7])]:
            self.assertGreaterEqual(np.linalg.norm(lattice_h(k, 1.2)), .2-1e-14)

    def test_13_scalar_offset_and_tilt_do_not_change_crossing(self):
        k0 = np.array([0., 0., np.pi/2])
        np.testing.assert_allclose(np.linalg.eigvalsh(.7*I+lattice_H(k0, tilt=2.)),
                                   [.7, .7], atol=3e-16)
        e = np.diag([1., 2., 3.])
        q = np.array([1., 0, 0])
        self.assertAlmostEqual(np.linalg.norm(e @ q), 1.)
        self.assertAlmostEqual(np.linalg.norm(e @ [0, 1, 0]), 2.)
        # Overtilted spectrum: both branches have the same sign along +x.
        self.assertTrue(np.all(np.linalg.eigvalsh(2*I+pauli(e @ q)) > 0))


if __name__ == "__main__":
    main(__name__, "stable_two_band_node_audit", report)
