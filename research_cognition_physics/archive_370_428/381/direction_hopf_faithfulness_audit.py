"""Round 381: faithful direction labels, Hopf fibres, and source resources.

The S^3 label space is an explicit candidate boundary, not emergent space.
Analytic proofs in research_note_381.md carry universal quantifiers.  This
file cross-checks the representations, instruments, and sharp recovery bound.
Only NumPy and the existing standard-library test/result runner are used.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

from growing_stream_audit import main


I2 = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
E0 = np.array([1., 0.], complex)
P0 = np.outer(E0, E0.conj())


def normalize(z):
    z = np.asarray(z, complex)
    return z / np.linalg.norm(z)


def projector(z):
    z = np.asarray(z, complex)
    return np.outer(z, z.conj())


def bloch(z):
    z = np.asarray(z, complex)
    return np.real(np.einsum("i,aij,j->a", z.conj(), PAULI, z))


def su2(z):
    """Smooth SU(2) representative determined by its full first column z."""
    a, b = np.asarray(z, complex)
    return np.array([[a, -b.conjugate()], [b, a.conjugate()]])


def controlled(u):
    out = np.eye(4, dtype=complex)
    out[2:, 2:] = u
    return out


def source(z, visibility=1.):
    """CQ order.  This is a newly specified coherent source, not a map of P_z."""
    z = np.asarray(z, complex)
    rho = np.zeros((4, 4), complex)
    rho[:2, :2] = P0 / 2
    rho[2:, 2:] = projector(z) / 2
    rho[:2, 2:] = visibility * np.outer(E0, z.conj()) / 2
    rho[2:, :2] = rho[:2, 2:].conj().T
    return rho


def trace_distance(a, b):
    delta = np.asarray(a) - np.asarray(b)
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2))) / 2)


def trace_last(rho, dim_a, dim_b):
    return np.einsum("abcb->ac", rho.reshape(dim_a, dim_b, dim_a, dim_b))


def trace_first(rho, dim_a, dim_b):
    return np.einsum("abad->bd", rho.reshape(dim_a, dim_b, dim_a, dim_b))


def sqrt_psd(a):
    eig, basis = np.linalg.eigh((a+a.conj().T)/2)
    return (basis * np.sqrt(np.maximum(eig, 0.))) @ basis.conj().T


def helstrom_instrument(a, b):
    eig, basis = np.linalg.eigh(a-b)
    positive = basis[:, eig > 1e-12]
    effect = positive @ positive.conj().T
    return effect, [sqrt_psd(effect), sqrt_psd(np.eye(len(a))-effect)]


def apply_kraus(kraus, rho):
    return sum(k @ rho @ k.conj().T for k in kraus)


def recovery_kraus():
    """R(rho) = (|0><0| x P0 Tr rho + |1><1| x rho)/2."""
    out = []
    for j in range(2):
        k = np.zeros((4, 2), complex)
        k[0, j] = 1/np.sqrt(2)
        out.append(k)
    k = np.zeros((4, 2), complex)
    k[2:, :] = I2/np.sqrt(2)
    return out + [k]


def sections(theta, phi):
    c, s = np.cos(theta/2), np.sin(theta/2)
    north = np.array([c, np.exp(1j*phi)*s])
    south = np.array([np.exp(-1j*phi)*c, s])
    return north, south


def reference_observables():
    """Four norm-one observables, ordered Re z0, Im z0, Re z1, Im z1."""
    values = []
    for j in range(2):
        a = np.zeros((4, 4), complex)
        a[0, 2+j] = 1
        values.extend([a+a.conj().T, -1j*a+1j*a.conj().T])
    return np.asarray(values)


def reference_probabilities(z, visibility):
    rho = source(z, visibility)
    return np.array([np.trace((np.eye(4)+b) @ rho).real/2 for b in reference_observables()])


def decode_reference(probabilities, visibility):
    if visibility <= 0:
        raise ValueError("A calibrated nonzero visibility is necessary.")
    real = (2*np.asarray(probabilities)-1)/visibility
    return real[0::2] + 1j*real[1::2]


def visible_dilation(z, visibility):
    """C,Q,E; E starts blank and is retained until the declared partial trace."""
    rotation = np.array([[visibility, -np.sqrt(1-visibility**2)],
                         [np.sqrt(1-visibility**2), visibility]], complex)
    joint_u = np.eye(8, dtype=complex)
    joint_u[4:, 4:] = np.kron(I2, rotation)
    before = np.kron(source(z), P0)
    return joint_u @ before @ joint_u.conj().T, joint_u


def tensor_power(a, n):
    out = np.array([[1.]], complex)
    for _ in range(n):
        out = np.kron(out, a)
    return out


def random_channel(input_dim, rng):
    """A generic CPTP Q^n -> CQ by independent Stinespring QR construction."""
    env_dim = max(2, math.ceil(input_dim/4))
    raw = rng.normal(size=(4*env_dim, input_dim)) + 1j*rng.normal(size=(4*env_dim, input_dim))
    isometry, _ = np.linalg.qr(raw)
    return [isometry.reshape(4, env_dim, input_dim)[:, e, :] for e in range(env_dim)]


def density(dim, rng):
    a = rng.normal(size=(dim, dim)) + 1j*rng.normal(size=(dim, dim))
    rho = a @ a.conj().T
    return rho / np.trace(rho)


def real_coordinates(z):
    return np.asarray(z, complex).view(float)


def tangent_basis(z):
    _, _, vh = np.linalg.svd(real_coordinates(z).reshape(1, 4))
    return vh[1:].T


def hopf_differential(z):
    columns = []
    for real in np.eye(4):
        dz = real[0::2] + 1j*real[1::2]
        columns.append(2*np.real(np.einsum("i,aij,j->a", np.conj(z), PAULI, dz)))
    return np.asarray(columns).T


def tetrahedral_certificate():
    """Exactly the same four qubit effects can be lifted to S^3 labels."""
    ns = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1.]])/np.sqrt(3)
    zs = [np.array([np.sqrt((1+n[2])/2), (n[0]+1j*n[1])/np.sqrt(2*(1+n[2]))]) for n in ns]
    zs = [np.exp(1j*f)*z for f, z in zip([.2, 1.3, -.7, 2.1], zs)]
    probes = [(I2+sum(r[j]*PAULI[j] for j in range(3)))/2 for r in np.r_[np.zeros((1, 3)), .75*np.eye(3)]]
    probabilities = np.array([[np.trace(projector(z) @ rho).real for rho in probes] for z in zs])
    differences = probabilities[1:, 1:]-probabilities[0, 1:]-probabilities[1:, :1]+probabilities[0, 0]
    estimated = differences/.75
    return {"labels": zs, "probabilities": probabilities, "contrast": estimated,
            "singular_values": np.linalg.svd(estimated, compute_uv=False)}


def report():
    z = normalize([.3+.4j, .7-.2j])
    phi, v = 1.2, .65
    cert = tetrahedral_certificate()
    n, s = sections(np.pi/2, np.pi)
    phases = np.linspace(0, 2*np.pi, 65)
    transitions = np.array([np.vdot(*sections(np.pi/2, p)) for p in phases])
    winding = np.sum(np.angle(transitions[1:]*transitions[:-1].conj()))/(2*np.pi)
    errors = []
    for phase in phases:
        errors.append(trace_distance(apply_kraus(recovery_kraus(), projector(z)),
                                     source(np.exp(1j*phase)*z, v)))
    noise = np.array([1., -1., .7, -.4])*1e-4
    fitted = decode_reference(reference_probabilities(z, v)+noise, v)
    tangent = tangent_basis(z)
    return {
        "round": 381,
        "scope": {
            "baseline": "frozen through round 380; no concurrent-round dependence",
            "result": "Removing label injectivity leaves an S3 boundary with all other round-380 direction hypotheses; coherent-source recovery from fibre-blind qubit data has sharp minimax error v/2.",
            "extra_inputs": ["S3 declared as a candidate complete boundary of a four-dimensional neighbourhood",
                             "full qubit direction effects and SU2 redirection action",
                             "a separately supplied, phase-calibrated controlled preparation source",
                             "blank registers, relative-phase reference, known visibility and repeatable preparations"],
            "not_claimed": ["a four-dimensional universe", "a measurement of an isolated state's global phase",
                            "coherent control can be obtained from an unknown gate channel",
                            "finite samples certify global label injectivity or topology",
                            "space is generated by the cognitive axioms"]},
        "candidate_boundary_dimension": 3,
        "candidate_neighbourhood_dimension": 4,
        "qubit_effect_orbit_dimension": 2,
        "tetrahedral_probability_table": cert["probabilities"].tolist(),
        "contrast_singular_values": cert["singular_values"].tolist(),
        "contrast_minimum_margin_exact_data": float(cert["singular_values"][-1]),
        "local_hopf_singular_values": np.linalg.svd(hopf_differential(z) @ tangent, compute_uv=False).tolist(),
        "relative_source_singular_values": np.linalg.svd(v*tangent, compute_uv=False).tolist(),
        "north_south_transition_winding": float(winding),
        "patch_switch_at_phi_pi": {"qubit_distance": trace_distance(projector(n), projector(s)),
                                  "coherent_source_distance": trace_distance(source(n), source(s))},
        "phase_discrimination": {"visibility": v, "phase_difference": phi,
                                 "qubit_distance": trace_distance(projector(z), projector(np.exp(1j*phi)*z)),
                                 "source_distance_matrix": trace_distance(source(z, v), source(np.exp(1j*phi)*z, v)),
                                 "source_distance_formula": v*abs(np.sin(phi/2))},
        "three_interfaces": {"density_source": "S3/U1 = S2",
                             "ordinary_su2_channel": "S3/{+1,-1} = RP3; generic phase is detectable",
                             "phase_calibrated_controlled_source": "injective S3 for v > 0",
                             "ordinary_gate_choi_distance_at_phase_phi": abs(np.sin(phi)),
                             "ordinary_gate_antipodal_distance": 0.},
        "optimal_recovery": {"analytic_minimax_error": v/2,
                             "explicit_channel_min_error_on_65_checks": min(errors),
                             "explicit_channel_max_error_on_65_checks": max(errors),
                             "antipodal_source_separation": trace_distance(source(z, v), source(-z, v)),
                             "input_copies_cannot_remove_same_fibre_identity": True},
        "four_reference_readouts": {"visibility": v, "probability_error_cap": 1e-4,
                                    "spinor_fit_error": float(np.linalg.norm(fitted-z)),
                                    "deterministic_bound": 4e-4/v,
                                    "normalization_of_raw_estimate_not_assumed": True},
        "source_resource_boundary": "The conditional source is supplied as a physical I2 direct-sum V_z, not synthesized from P_z or Ad(V_z). Any unknown input reference is included in unitary checks; source tomography assumes newly prepared blank Q.",
        "primary_sources": ["https://arxiv.org/pdf/quant-ph/0108137", "https://arxiv.org/pdf/1309.7976"],
    }


class Checks(unittest.TestCase):
    def test_01_effects_and_exact_fibres(self):
        z = normalize([.3+.4j, .7-.2j])
        e = projector(z)
        np.testing.assert_allclose(e @ e, e, atol=2e-16)
        np.testing.assert_allclose(e, (I2+np.einsum("a,aij->ij", bloch(z), PAULI))/2, atol=2e-16)
        self.assertAlmostEqual(np.linalg.norm(bloch(z)), 1.)
        for phi in (0., .4, np.pi, 4.1):
            np.testing.assert_allclose(projector(np.exp(1j*phi)*z), e, atol=3e-16)
        self.assertGreater(np.linalg.norm(z-(-z)), 1.9)

    def test_02_transitive_whole_family_su2_covariance(self):
        z, w = normalize([1+2j, 3-1j]), normalize([2-.2j, -1+1j])
        u = su2(w) @ su2(z).conj().T
        np.testing.assert_allclose(u @ z, w, atol=2e-16)
        np.testing.assert_allclose(u.conj().T @ u, I2, atol=4e-16)
        self.assertAlmostEqual(np.linalg.det(u), 1.)
        for a in (z, w, E0, normalize([1, 1j])):
            np.testing.assert_allclose(projector(u @ a), u @ projector(a) @ u.conj().T, atol=3e-16)

    def test_03_finite_contrast_certificate_does_not_detect_fibre(self):
        row = tetrahedral_certificate()
        np.testing.assert_allclose(row["singular_values"], [2/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)], atol=3e-16)
        replacement = [np.exp(1j*f)*z for f, z in zip([.4, -.3, 1.7, .2], row["labels"])]
        for z, w in zip(row["labels"], replacement):
            np.testing.assert_allclose(projector(z), projector(w), atol=3e-16)
        self.assertEqual(np.linalg.matrix_rank(row["contrast"]), 3)

    def test_04_local_sections_transition_and_nonzero_winding(self):
        phases = np.linspace(0, 2*np.pi, 65)
        gs = []
        for phi in phases:
            north, south = sections(np.pi/2, phi)
            np.testing.assert_allclose(south, np.exp(-1j*phi)*north, atol=2e-16)
            np.testing.assert_allclose(projector(north), projector(south), atol=2e-16)
            gs.append(np.vdot(north, south))
        winding = np.sum(np.angle(np.array(gs[1:])*np.conj(gs[:-1])))/(2*np.pi)
        self.assertAlmostEqual(winding, -1.)
        # Chart regularity at its own pole, not a false global smooth section.
        for phi in (0., .7, 2.3):
            np.testing.assert_allclose(sections(0., phi)[0], E0, atol=1e-16)
            np.testing.assert_allclose(sections(np.pi, phi)[1], [0, 1], atol=1e-16)

    def test_05_controlled_preparation_and_extended_family_covariance(self):
        z, w = normalize([1+2j, 3-1j]), normalize([2-.2j, -1+1j])
        c = controlled(su2(z))
        blank = np.kron(np.array([1, 1])/np.sqrt(2), E0)
        np.testing.assert_allclose(c.conj().T @ c, np.eye(4), atol=3e-16)
        np.testing.assert_allclose(projector(c @ blank), source(z), atol=3e-16)
        u = su2(w)
        for v in (0., .4, 1.):
            np.testing.assert_allclose(controlled(u) @ source(z, v) @ controlled(u).conj().T,
                                       source(u @ z, v), atol=3e-16)

    def test_06_patch_change_requires_corresponding_relative_phase_control(self):
        for phi in (.2, 1., np.pi):
            n, s = sections(1.1, phi)
            mismatch = trace_distance(source(n), source(s))
            self.assertAlmostEqual(mismatch, abs(np.sin(phi/2)))
            correction = np.diag([1, 1, np.exp(1j*phi), np.exp(1j*phi)])
            np.testing.assert_allclose(correction @ source(s) @ correction.conj().T, source(n), atol=3e-16)

    def test_07_visibility_has_an_explicit_unitary_environment(self):
        z = normalize([.3+.4j, .7-.2j])
        for v in (0., .3, .8, 1.):
            joint, unitary = visible_dilation(z, v)
            np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(8), atol=2e-16)
            np.testing.assert_allclose(trace_last(joint, 4, 2), source(z, v), atol=2e-16)
            self.assertGreaterEqual(np.linalg.eigvalsh(source(z, v)).min(), -1e-15)
            self.assertAlmostEqual(np.trace(source(z, v)).real, 1.)

    def test_08_phase_discrimination_and_complete_helstrom_instrument(self):
        z = normalize([.3+.4j, .7-.2j])
        for v, phi in ((1., np.pi), (.65, 1.2), (.3, .8), (0., 1.)):
            a, b = source(z, v), source(np.exp(1j*phi)*z, v)
            expected = v*abs(np.sin(phi/2))
            self.assertAlmostEqual(trace_distance(a, b), expected)
            effect, kraus = helstrom_instrument(a, b)
            np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(4), atol=2e-15)
            self.assertAlmostEqual(np.trace(effect @ (a-b)).real, expected)

    def test_09_four_relative_phase_readouts_recover_full_label_with_error_bound(self):
        z, v, eps = normalize([.3+.4j, .7-.2j]), .65, 1e-4
        for observable in reference_observables():
            eig = np.linalg.eigvalsh((np.eye(4)+observable)/2)
            self.assertGreaterEqual(eig.min(), 0.)
            self.assertLessEqual(eig.max(), 1.)
        probabilities = reference_probabilities(z, v)
        np.testing.assert_allclose(decode_reference(probabilities, v), z, atol=3e-16)
        fitted = decode_reference(probabilities+eps*np.array([1, -1, 1, -1]), v)
        self.assertAlmostEqual(np.linalg.norm(fitted-z), 4*eps/v)
        with self.assertRaises(ValueError):
            decode_reference(probabilities, 0.)

    def test_10_unknown_reference_does_not_make_ordinary_global_phase_visible(self):
        rng = np.random.default_rng(381)
        rho_qr = density(4, rng)
        z = normalize([.3+.4j, .7-.2j])
        u = su2(z)
        np.testing.assert_allclose(su2(-z), -u, atol=1e-16)
        phi = 1.2
        shifted = su2(np.exp(1j*phi)*z)
        np.testing.assert_allclose(shifted, u @ np.diag([np.exp(1j*phi), np.exp(-1j*phi)]), atol=2e-16)
        bell = np.array([1., 0, 0, 1])/np.sqrt(2)
        choi = projector(np.kron(u, I2) @ bell)
        shifted_choi = projector(np.kron(shifted, I2) @ bell)
        self.assertAlmostEqual(trace_distance(choi, shifted_choi), abs(np.sin(phi)))
        plus = projector(np.array([1, 1])/np.sqrt(2))
        original = np.kron(plus, rho_qr)
        ordinary = np.kron(u, I2)
        np.testing.assert_allclose(ordinary @ rho_qr @ ordinary.conj().T,
                                   (-ordinary) @ rho_qr @ (-ordinary).conj().T, atol=1e-16)
        outputs = []
        for op in (u, -u):
            gate = np.kron(controlled(op), I2)
            output = gate @ original @ gate.conj().T
            outputs.append(output)
            np.testing.assert_allclose(gate.conj().T @ output @ gate, original, atol=4e-16)
            np.testing.assert_allclose(trace_first(output, 4, 2), trace_first(rho_qr, 2, 2), atol=3e-16)
        self.assertAlmostEqual(trace_distance(*outputs), 1.)

    def test_11_minimax_upper_bound_is_a_single_cptp_recovery_for_all_labels(self):
        kraus = recovery_kraus()
        np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), I2, atol=3e-16)
        for z in (E0, normalize([1, 1]), normalize([.3+.4j, .7-.2j])):
            recovered = apply_kraus(kraus, projector(z))
            np.testing.assert_allclose(recovered, source(z, 0.), atol=2e-16)
            for v, phi in ((.3, .1), (.65, 1.2), (1., np.pi)):
                self.assertAlmostEqual(trace_distance(recovered, source(np.exp(1j*phi)*z, v)), v/2)

    def test_12_antipodal_lower_bound_for_multiple_copies_and_generic_channels(self):
        rng = np.random.default_rng(1381)
        z, v = normalize([.3+.4j, .7-.2j]), .65
        for n in (1, 2, 3):
            inp = tensor_power(projector(z), n)
            np.testing.assert_allclose(inp, tensor_power(projector(-z), n), atol=1e-16)
            kraus = random_channel(2**n, rng)
            np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(2**n), atol=6e-16)
            out = apply_kraus(kraus, inp)
            errors = [trace_distance(out, source(sign*z, v)) for sign in (1, -1)]
            self.assertGreaterEqual(sum(errors)+1e-15, v)
            self.assertGreaterEqual(max(errors)+1e-15, v/2)

    def test_13_label_leakage_changes_the_lower_bound_only_by_contractivity(self):
        z, v = E0, .65
        # Unlike the ideal Hopf source, this explicitly leaks a small label bit.
        a = projector(z)
        b = .96*a+.04*projector([0, 1])
        delta = trace_distance(a, b)
        self.assertAlmostEqual(delta, .04)
        for n in (1, 2, 3):
            qa, qb = tensor_power(a, n), tensor_power(b, n)
            self.assertLessEqual(trace_distance(qa, qb), n*delta+1e-15)
        rng = np.random.default_rng(2381)
        kraus = random_channel(2, rng)
        oa, ob = apply_kraus(kraus, a), apply_kraus(kraus, b)
        self.assertLessEqual(trace_distance(oa, ob), delta+1e-15)
        max_error = max(trace_distance(oa, source(z, v)), trace_distance(ob, source(-z, v)))
        self.assertGreaterEqual(max_error+1e-15, (v-delta)/2)

    def test_14_relative_interface_resolves_the_missing_tangent_direction(self):
        z, v = normalize([.3+.4j, .7-.2j]), .65
        tangent = tangent_basis(z)
        differential = hopf_differential(z)
        np.testing.assert_allclose(differential @ real_coordinates(1j*z), np.zeros(3), atol=3e-16)
        np.testing.assert_allclose(np.linalg.svd(differential @ tangent, compute_uv=False), [2, 2, 0], atol=8e-16)
        np.testing.assert_allclose(np.linalg.svd(v*tangent, compute_uv=False), [v]*3, atol=3e-16)
        direction = tangent[:, 0]
        h = 1e-6
        dz = direction[0::2]+1j*direction[1::2]
        numeric = (bloch(normalize(z+h*dz))-bloch(normalize(z-h*dz)))/(2*h)
        np.testing.assert_allclose(numeric, differential @ direction, atol=1e-9)


if __name__ == "__main__":
    main(__name__, "direction_hopf_faithfulness_audit", report)
