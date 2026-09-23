"""Round 375: fixed-metric Clifford compatibility and operator identifiability.

The manifold, volume, irreducible rank-two symbol and (in the main theorem)
Levi-Civita connection are inputs. This is a local geometric classification,
not a derivation of spacetime, Maxwell theory or Einstein dynamics.
"""
import unittest
import numpy as np
from growing_stream_audit import main

I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)
ANTI = 1j*np.concatenate((I[None], PAULI))
E = np.array([[1., .25, 0.], [0., 1.6, .1], [.15, 0., .8]])
POINT = np.array([.2, -.3, .4])
HESS = np.array([[.2, .07, 0.], [.07, .4, 0.], [0., 0., .6]])


def sig(v):
    return np.einsum("...a,aij->...ij", np.asarray(v), PAULI)


def comm(a, b):
    return a@b-b@a


def real_vector(a):
    flat = np.asarray(a).ravel()
    return np.concatenate((flat.real, flat.imag))


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conjugate().T


def symbol(e):
    # A[j] = e[a,j] sigma[a].
    return sig(np.asarray(e).T)


def compatibility_map(a):
    return np.column_stack([real_vector(np.array([comm(b, aj) for aj in a]))
                            for b in ANTI])


def compatibility_residual(a, da, gamma, omega):
    # gamma[i,j,k] = Gamma^j_{i k}; da[i,j] = partial_i A^j.
    return np.array([[da[i, j]+comm(omega[i], a[j])+
                      sum(gamma[i, j, k]*a[k] for k in range(3))
                      for j in range(3)] for i in range(3)])


def solve_compatibility(a, da, gamma):
    matrix = compatibility_map(a)
    rows = []
    for i in range(3):
        target = -da[i]-np.einsum("jk,kab->jab", gamma[i], a)
        coeff = np.linalg.lstsq(matrix, real_vector(target), rcond=None)[0]
        rows.append(np.einsum("a,aij->ij", coeff, ANTI))
    omega = np.array(rows)
    return omega, int(np.linalg.matrix_rank(matrix)), float(np.max(abs(
        compatibility_residual(a, da, gamma, omega))))


def conformal_geometry(x):
    # g_ij=exp(2f) delta_ij, f=x.HESS.x/2.
    x = np.asarray(x)
    f, grad = float(x@HESS@x/2), HESS@x
    e = np.exp(-f)*np.eye(3)
    a = symbol(e)
    da = -grad[:, None, None, None]*a[None]
    gamma = np.zeros((3, 3, 3))
    dgamma = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                gamma[i, j, k] = ((j == i)*grad[k]+(j == k)*grad[i]-
                                  (i == k)*grad[j])
                for d in range(3):
                    dgamma[d, i, j, k] = (
                        (j == i)*HESS[k, d]+(j == k)*HESS[i, d]-
                        (i == k)*HESS[j, d])
    omega = np.array([.5j*sig(np.cross(axis, grad)) for axis in np.eye(3)])
    domega = np.array([[.5j*sig(np.cross(axis, HESS[:, d]))
                       for axis in np.eye(3)] for d in range(3)])
    return f, grad, e, a, da, gamma, dgamma, omega, domega


def curvature(connection, derivative=None):
    if derivative is None:
        derivative = np.zeros((3,)+connection.shape, dtype=complex)
    return np.array([[derivative[i, j]-derivative[j, i]+
                      comm(connection[i], connection[j]) for j in range(3)]
                     for i in range(3)])


def curvature_compatibility(a, spin_f, tangent_r):
    return np.array([[[comm(spin_f[i, j], a[k])+
                       sum(tangent_r[i, j, k, ell]*a[ell] for ell in range(3))
                       for k in range(3)] for j in range(3)] for i in range(3)])


def scalar_curvature(g_inverse, tangent_r):
    # R[i,j,k,l] = R^k_{l i j}; Ric[l,j]=sum_i R^i_{l i j}.
    ricci = np.array([[sum(tangent_r[i, j, i, ell] for i in range(3))
                       for j in range(3)] for ell in range(3)])
    return float(np.einsum("ij,ij", g_inverse, ricci).real)


def contraction(a, omega):
    return -1j*sum(a[j]@omega[j] for j in range(3))


def potential_coefficients(v):
    phi = float(np.trace(v).real/2)
    b = np.array([np.trace(p@v).real/2 for p in PAULI])
    return phi, b


def absorb_potential(e, central, v):
    phi, b = potential_coefficients(v)
    return np.asarray(central)+np.linalg.solve(e, b), phi


def su2_connection(b):
    return 1j*sig(np.asarray(b))


def induced_tangent_connection(b):
    # -Gamma_i is the generator of the Bloch-vector parallel transport.
    gamma = np.zeros((3, 3, 3))
    for i in range(3):
        for k in range(3):
            gamma[i, :, k] = -2*np.cross(b[i], np.eye(3)[k])
    return gamma


def torsion(gamma):
    # output[k,i,j] = Gamma^k_{i j}-Gamma^k_{j i}.
    return np.array([[[gamma[i, k, j]-gamma[j, k, i] for j in range(3)]
                      for i in range(3)] for k in range(3)])


def traceless_symmetric_basis():
    result = [np.diag([1., -1., 0.])/np.sqrt(2),
              np.diag([1., 1., -2.])/np.sqrt(6)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        b = np.zeros((3, 3))
        b[i, j] = b[j, i] = 1/np.sqrt(2)
        result.append(b)
    return np.array(result)


def contraction_map():
    return np.column_stack([real_vector(contraction(PAULI, su2_connection(
        np.eye(9)[k].reshape(3, 3)))) for k in range(9)])


def finite_mode_generator(b):
    # Exact reducing Fourier-mode subspace of the constant-coefficient operator
    # on a flat torus with an explicitly chosen periodic spin structure.
    momenta = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 1], [-1, 2, 0]])
    result = np.zeros((2*len(momenta), 2*len(momenta)), complex)
    w = contraction(PAULI, su2_connection(b))
    for n, k in enumerate(momenta):
        result[2*n:2*n+2, 2*n:2*n+2] = sig(k)+w
    return result


def choi_state(u):
    # Acting on one half of a normalized maximally entangled reference.
    n = len(u)
    vector = u.reshape(n*n)/np.sqrt(n)
    return np.outer(vector, vector.conjugate())


def blind_connection_case():
    b = np.diag([1., 2., -3.])
    omega = su2_connection(b)
    gamma = induced_tangent_connection(b)
    f, r = curvature(omega), curvature(gamma)
    h0, h1 = finite_mode_generator(np.zeros((3, 3))), finite_mode_generator(b)
    u0, u1 = unitary(h0, .37), unitary(h1, .37)
    return {
        "b": b.tolist(),
        "contraction_norm": opnorm(contraction(PAULI, omega)),
        "compatibility_residual": float(np.max(abs(compatibility_residual(
            PAULI, np.zeros((3, 3, 2, 2)), gamma, omega)))),
        "torsion_max_component": float(np.max(abs(torsion(gamma)))),
        "spin_curvature_xy_yz_zx_operator_norms": [opnorm(f[i, j])
                                                 for i, j in ((0, 1), (1, 2), (2, 0))],
        "tangent_curvature_max_component": float(np.max(abs(r))),
        "spin_vector_curvature_compatibility_residual": float(np.max(abs(
            curvature_compatibility(PAULI, f, r)))),
        "Levi_Civita_Riemann_curvature": 0.,
        "full_mode_H_difference": opnorm(h1-h0),
        "propagator_difference": opnorm(u1-u0),
        "maximally_entangled_reference_output_difference": opnorm(choi_state(u1)-choi_state(u0)),
        "interpretation": "These are different torsionful connections, but the same Weyl operator and every propagation channel; without an independently supplied transport readout they are not different observable propagation geometries."
    }


def oscillator_annihilator(size):
    result = np.zeros((size, size), complex)
    for n in range(1, size):
        result[n-1, n] = np.sqrt(n)
    return result


def landau_hamiltonian(size, magnetic, kz):
    # Pi=-i partial+a, a_y=-B x; [Pi_x,Pi_y]=iB.
    ann = oscillator_annihilator(size)
    px = np.sqrt(magnetic/2)*(ann+ann.conjugate().T)
    py = -1j*np.sqrt(magnetic/2)*(ann-ann.conjugate().T)
    return (np.kron(PAULI[0], px)+np.kron(PAULI[1], py)+
            np.kron(PAULI[2], kz*np.eye(size)))


def landau_block_embedding(size, n):
    if not 1 <= n < size:
        raise ValueError("Use a complete physical two-state block below the cutoff.")
    embedding = np.zeros((2*size, 2), complex)
    embedding[n, 0] = 1
    embedding[size+n-1, 1] = 1
    return embedding


def landau_case(magnetic=.7, kz=.4, n=2, size=9, time=.6):
    full = landau_hamiltonian(size, magnetic, kz)
    q = landau_block_embedding(size, n)
    block = q.conjugate().T@full@q
    actual = unitary(full, time)@q[:, 0]
    probability = float(abs(np.vdot(q[:, 1], actual))**2)
    energy = np.sqrt(kz*kz+2*magnetic*n)
    return {
        "B": magnetic, "kz": kz, "landau_index": n, "oscillator_cutoff": size,
        "metric": "flat Euclidean spatial metric; fixed Weyl principal symbol sigma.k",
        "central_connection_a": "(0,-B*x,0)",
        "central_curvature_f_xy": -magnetic,
        "block_eigenvalues": np.linalg.eigvalsh(block).tolist(),
        "analytic_eigenvalues": [-float(energy), float(energy)],
        "complete_block_leakage": float(np.linalg.norm(full@q-q@block)),
        "time": time,
        "fixed_readout_transition_probability": probability,
        "analytic_transition_probability": float(2*magnetic*n/energy**2*np.sin(time*energy)**2),
        "upper_zero_mode_energy": kz,
        "cutoff_boundary_warning": "The truncated CCR fails on the highest oscillator state. The artificial unpaired lower top state is excluded; only the exact reducing blocks n<size and physical upper n=0 mode are certified.",
        "physical_scope": "A supplied nonzero central curvature can change the spectrum without changing the principal metric. No field equation, charge value, laboratory field source, global degeneracy or chiral anomaly is derived."
    }


def report():
    f, grad, e, a, da, gamma, dgamma, omega, domega = conformal_geometry(POINT)
    recovered, rank, residual = solve_compatibility(a, da, gamma)
    tangent_r, spin_f = curvature(gamma, dgamma), curvature(omega, domega)
    r_actual = scalar_curvature(e.T@e, tangent_r)
    r_expected = np.exp(-2*f)*(-4*np.trace(HESS)-2*grad@grad)
    v = .3*I+sig([.2, -.4, .1])
    central = np.array([.1, -.2, .05])
    effective, phi = absorb_potential(E, central, v)
    before = sig(E@central)+v
    after = sig(E@effective)+phi*I
    return {
        "round": 375,
        "scope": {
            "baseline": "Frozen through 373; independent of 374.",
            "claim": "A fixed metric connection and irreducible rank-two Clifford symbol determine a Hermitian spin connection up to an imaginary scalar one-form. Connection coefficients and full propagation operators have different identifiability once additional potentials or torsion are permitted.",
            "inputs": ["oriented smooth three-dimensional spatial patch", "positive spatial metric, its volume and nondegenerate rank-two Clifford symbol", "Levi-Civita connection fixed for the main classification", "compatibility of spinor and vector parallel transport", "extra choices of domain, spin structure, source fields and measurements for witnesses"],
            "not_proved": ["dimension three from FUCP", "Levi-Civita or torsionlessness from probability conservation", "an operationally unique geometry from an unspecified notion of transport", "Maxwell equations, standard-model group or field sources", "spacetime or Einstein dynamics from cognition"],
            "novelty": "Goes beyond round349's arbitrary Hermitian potential: classifies the fixed-connection center, removes duplicate a/V parameterizations, and exhibits a five-dimensional torsionful connection kernel invisible to the entire Dirac propagation operator.",
            "no_double_counting": "D_a+V=D_(a+E^-1 b)+phi I. An overall spatially constant scalar energy changes a propagator phase but not its quantum channel."
        },
        "fixed_metric_connection": {
            "point": POINT.tolist(), "conformal_f": f, "gradient_f": grad.tolist(),
            "rank_per_direction_anti_Hermitian_4dim": rank,
            "nullity_per_direction": 4-rank,
            "recovered_spin_connection_max_error": float(np.max(abs(recovered-omega))),
            "compatibility_max_error": residual,
            "spin_vector_curvature_compatibility_max_error": float(np.max(abs(
                curvature_compatibility(a, spin_f, tangent_r)))),
            "scalar_curvature_direct": r_actual,
            "scalar_curvature_conformal_formula": float(r_expected)
        },
        "central_contraction": {
            "singular_values": np.linalg.svd(E, compute_uv=False).tolist(),
            "a_original": central.tolist(), "a_effective": effective.tolist(),
            "scalar_phi": phi, "potential_reallocation_operator_error": opnorm(before-after)
        },
        "unfixed_metric_connection": {
            "real_domain_dimension": 9,
            "Dirac_contraction_rank": int(np.linalg.matrix_rank(contraction_map())),
            "kernel_dimension": 9-int(np.linalg.matrix_rank(contraction_map())),
            "kernel": "symmetric trace-free b; all five basis elements checked",
            "explicit_case": blind_connection_case()
        },
        "physical_central_curvature": [landau_case(magnetic=b) for b in (.4, .7)],
        "sources": [
            "https://arxiv.org/pdf/math/0101155",
            "https://arxiv.org/html/1401.3160",
            "https://arxiv.org/pdf/1608.06632"
        ]
    }


class Checks(unittest.TestCase):
    def test_01_nondegenerate_symbol_Clifford_identity(self):
        a, ginv = symbol(E), E.T@E
        for i in range(3):
            for j in range(3):
                np.testing.assert_allclose(a[i]@a[j]+a[j]@a[i], 2*ginv[i,j]*I, atol=2e-15)
        self.assertGreater(np.linalg.det(E), 0)

    def test_02_compatibility_kernel_is_exactly_center(self):
        for e in (E, np.eye(3), np.diag([.2, 2., 4.])):
            matrix = compatibility_map(symbol(e))
            self.assertEqual(np.linalg.matrix_rank(matrix), 3)
            np.testing.assert_array_equal(matrix[:, 0], 0)
            self.assertEqual(np.linalg.matrix_rank(matrix[:, 1:]), 3)

    def test_03_curved_LC_compatible_spin_connection_independently_recovered(self):
        for x in (POINT, np.array([-.1, .4, .2])):
            _, _, _, a, da, gamma, _, omega, _ = conformal_geometry(x)
            recovered, rank, residual = solve_compatibility(a, da, gamma)
            self.assertEqual(rank, 3)
            self.assertLess(residual, 2e-15)
            np.testing.assert_allclose(recovered, omega, atol=2e-15)
            central = 1j*np.array([.3, -.1, .7])[:, None, None]*I
            self.assertLess(np.max(abs(compatibility_residual(a, da, gamma, omega+central))), 2e-15)

    def test_04_metric_torsionless_and_actual_curvature_checks(self):
        f, grad, e, a, da, gamma, dgamma, omega, domega = conformal_geometry(POINT)
        metric = np.exp(2*f)*np.eye(3)
        for i in range(3):
            np.testing.assert_allclose(2*grad[i]*metric,
                gamma[i].T@metric+metric@gamma[i], atol=2e-15)
        np.testing.assert_array_equal(torsion(gamma), 0)
        r, spin_f = curvature(gamma, dgamma), curvature(omega, domega)
        self.assertLess(np.max(abs(curvature_compatibility(a, spin_f, r))), 2e-15)
        scalar = scalar_curvature(e.T@e, r)
        self.assertAlmostEqual(scalar, np.exp(-2*f)*(-4*np.trace(HESS)-2*grad@grad), places=13)
        self.assertGreater(abs(scalar), 1)

    def test_05_weighted_formal_adjoint_matches_contracted_LC(self):
        f, grad, _, a, da, _, _, omega, _ = conformal_geometry(POINT)
        lower = contraction(a, omega)
        weighted_divergence = sum(da[i, i]+3*grad[i]*a[i] for i in range(3))
        np.testing.assert_allclose(lower, -.5j*weighted_divergence, atol=2e-15)
        np.testing.assert_allclose(lower, -1j*np.exp(-f)*sig(grad), atol=2e-15)

    def test_06_general_local_SU2_frame_covariance_of_compatibility(self):
        _, _, _, a, da, gamma, _, omega, _ = conformal_geometry(POINT)
        u = unitary(sig([.3, -.2, .4]), .47)
        tangent = 1j*sig(np.array([[.1,.2,.3],[-.2,.1,0],[.4,0,-.1]]))
        aprime = np.array([u@aj@u.conjugate().T for aj in a])
        daprime = np.array([[u@da[i,j]@u.conjugate().T+comm(tangent[i], aprime[j])
                            for j in range(3)] for i in range(3)])
        oprime = np.array([u@omega[i]@u.conjugate().T-tangent[i] for i in range(3)])
        self.assertLess(np.max(abs(compatibility_residual(aprime, daprime, gamma, oprime))), 2e-15)

    def test_07_central_contraction_is_injective_and_gives_metric_norm(self):
        for alpha in ([.1, -.3, .5], [1., 0., 0.], [0., 0., 0.]):
            alpha = np.asarray(alpha)
            contracted = contraction(symbol(E), 1j*alpha[:, None, None]*I)
            np.testing.assert_allclose(contracted, sig(E@alpha), atol=2e-15)
            self.assertAlmostEqual(opnorm(contracted)**2, float(alpha@(E.T@E)@alpha), places=13)

    def test_08_potential_connection_reallocation_preserves_whole_operator(self):
        central = np.array([.1, -.2, .05])
        v = .3*I+sig([.2, -.4, .1])
        effective, phi = absorb_potential(E, central, v)
        for k in ([0.,0.,0.], [.3,.4,-.1], [-1.,.7,.2]):
            first = sig(E@(np.asarray(k)+central))+v
            second = sig(E@(np.asarray(k)+effective))+phi*I
            np.testing.assert_allclose(first, second, atol=2e-15)
            np.testing.assert_allclose(unitary(first,.37),unitary(second,.37),atol=2e-15)

    def test_09_global_scalar_energy_does_not_change_reference_channel(self):
        h = finite_mode_generator(np.zeros((3,3)))
        u, uphi = unitary(h,.41), unitary(h+.7*np.eye(len(h)),.41)
        self.assertGreater(opnorm(u-uphi), .2)
        np.testing.assert_allclose(choi_state(u), choi_state(uphi), atol=2e-16)

    def test_10_induced_torsionful_connection_is_metric_and_Clifford_compatible(self):
        b = np.array([[1.,.2,-.1],[.2,2.,.3],[-.1,.3,-3.]])
        gamma, omega = induced_tangent_connection(b), su2_connection(b)
        for gi in gamma:
            np.testing.assert_array_equal(gi+gi.T, 0)
        self.assertLess(np.max(abs(compatibility_residual(
            PAULI,np.zeros((3,3,2,2)),gamma,omega))), 2e-15)
        self.assertGreater(np.max(abs(torsion(gamma))), 1)

    def test_11_exact_five_dimensional_invisible_torsion_kernel(self):
        matrix = contraction_map()
        self.assertEqual(np.linalg.matrix_rank(matrix), 4)
        basis = traceless_symmetric_basis()
        self.assertEqual(np.linalg.matrix_rank(basis.reshape(5,9)), 5)
        for b in basis:
            np.testing.assert_allclose(contraction(PAULI,su2_connection(b)),0,atol=2e-15)
        b = np.array([[.4,.3,-.2],[-.1,.6,.7],[.2,.1,.5]])
        w = contraction(PAULI,su2_connection(b))
        axial = sum((b[i,j]*np.cross(np.eye(3)[i],np.eye(3)[j])
                     for i in range(3) for j in range(3)), start=np.zeros(3))
        np.testing.assert_allclose(w,np.trace(b)*I+1j*sig(axial),atol=2e-15)

    def test_12_nonzero_torsion_curvature_but_zero_Levi_Civita_curvature(self):
        row = blind_connection_case()
        self.assertEqual(row["contraction_norm"],0)
        np.testing.assert_allclose(row["spin_curvature_xy_yz_zx_operator_norms"],[4.,12.,6.],atol=1e-14)
        self.assertGreater(row["tangent_curvature_max_component"],0)
        self.assertEqual(row["Levi_Civita_Riemann_curvature"],0)
        self.assertLess(row["spin_vector_curvature_compatibility_residual"],2e-15)

    def test_13_same_full_generator_and_propagation_with_entangled_reference(self):
        h0 = finite_mode_generator(np.zeros((3,3)))
        for b in (np.diag([1.,2.,-3.]), traceless_symmetric_basis()[2]):
            h1 = finite_mode_generator(b)
            np.testing.assert_allclose(h1,h0,atol=2e-15)
            np.testing.assert_allclose(h1,h1.conjugate().T,atol=2e-15)
            u0,u1 = unitary(h0,.37),unitary(h1,.37)
            np.testing.assert_allclose(u0,u1,atol=2e-15)
            np.testing.assert_allclose(choi_state(u0),choi_state(u1),atol=2e-16)

    def test_14_Landau_blocks_and_transition_are_independently_diagonalized(self):
        for magnetic in (.4,.7):
            for n in (1,2,4):
                size,kz = 9,.4
                h = landau_hamiltonian(size,magnetic,kz)
                q = landau_block_embedding(size,n)
                block = q.conjugate().T@h@q
                expected = np.array([[kz,np.sqrt(2*magnetic*n)],[np.sqrt(2*magnetic*n),-kz]])
                np.testing.assert_allclose(block,expected,atol=2e-15)
                np.testing.assert_allclose(h@q,q@block,atol=2e-15)
                row = landau_case(magnetic,kz,n,size)
                np.testing.assert_allclose(row["block_eigenvalues"],row["analytic_eigenvalues"],atol=2e-15)
                self.assertAlmostEqual(row["fixed_readout_transition_probability"],
                                       row["analytic_transition_probability"],places=14)

    def test_15_Landau_chiral_mode_and_finite_cutoff_boundary_are_distinct(self):
        size,magnetic,kz = 8,.7,.4
        ann = oscillator_annihilator(size)
        ccr = comm(ann,ann.conjugate().T)
        np.testing.assert_allclose(ccr,np.eye(size)-size*np.diag([0.]*(size-1)+[1.]),atol=2e-15)
        h = landau_hamiltonian(size,magnetic,kz)
        physical = np.eye(2*size)[:,0]
        artificial = np.eye(2*size)[:,-1]
        np.testing.assert_allclose(h@physical,kz*physical,atol=2e-15)
        np.testing.assert_allclose(h@artificial,-kz*artificial,atol=2e-15)
        for n in (1,2,4):
            r1=landau_case(magnetic,kz,n,size)
            r2=landau_case(magnetic,kz,n,size+3)
            np.testing.assert_allclose(r1["block_eigenvalues"],r2["block_eigenvalues"],atol=2e-15)
        with self.assertRaises(ValueError):
            landau_block_embedding(size,size)

    def test_16_central_curvature_changes_spectral_gap_not_principal_symbol(self):
        first,second = landau_case(.4),landau_case(.7)
        self.assertEqual(first["metric"],second["metric"])
        self.assertNotEqual(first["central_curvature_f_xy"],second["central_curvature_f_xy"])
        self.assertGreater(second["block_eigenvalues"][1]-first["block_eigenvalues"][1],.3)
        self.assertGreater(abs(second["fixed_readout_transition_probability"]-
                               first["fixed_readout_transition_probability"]),.1)


if __name__ == "__main__":
    main(__name__,"clifford_connection_compatibility_audit",report)

