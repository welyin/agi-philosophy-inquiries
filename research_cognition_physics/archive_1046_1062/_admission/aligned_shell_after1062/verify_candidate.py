"""Aligned preparation shell: finite diagnostics, not a continuum proof.

Default compares the saved result read-only. --save-exclusive creates it once.
No author functions from earlier rounds are imported. Exact continuous bounds
are inherited with file hashes and checked here using Fraction. A rotating
finite quadrature only diagnoses the same algebraic covariance and poststate
identities; it is not the finite-mass continuous evolution.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results.json"
HISTORY = (
    "archive_370_428/research_note_382.md",
    "archive_370_428/research_note_383.md",
    "archive_370_428/research_note_384.md",
    "archive_370_428/research_note_386.md",
    "archive_370_428/research_note_425.md",
    "archive_467_530/research_note_523.md",
    "archive_1046_/1046/proof.md",
    "archive_1046_/1054/proof.md",
    "archive_1046_/1055/proof.md",
    "archive_1046_/1055/results.json",
)
PAULI = np.array([[[0, 1], [1, 0]],
                  [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)
I2 = np.eye(2, dtype=complex)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sdot(v):
    return np.einsum("j,jab->ab", v, PAULI)


def norm(a):
    return float(np.linalg.norm(a, 2))


def rational():
    radius = F(1, 200000)
    old_radius = F(1, 100000)
    m, mu, sigma = F(10**20), F(10**20, 2), F(1, 10**6)
    delta = F(3, 400000)
    old_p_lower = F(6889, 4000000)
    b_lower = old_p_lower/2-delta
    gap = 2*b_lower
    phase = F(3, 5)*mu*radius**2
    assert radius < old_radius
    assert b_lower == F(6829, 8000000)
    assert gap == F(6829, 4000000) > F(17, 10000)
    assert phase == 750000000
    # Uses the unchanged old full-domain RMS and K bounds, not new optimization.
    state_error = F(2*10**12)/m + F(2, 7)*F(13, 10**6)
    assert 2*state_error == F(1307, 175000000) < delta
    items = dict(radius=radius, old_task_radius=old_radius, mass=m,
                 sigma=sigma, transport_time=F(1), direction_read_time=F(1),
                 inherited_effect_error_upper=delta,
                 inherited_p_lower=old_p_lower,
                 aligned_abs_b_lower=b_lower, aligned_antipodal_gap_lower=gap,
                 continuous_endpoint_summary_error_lower=b_lower,
                 fixed_shell_force_phase=phase,
                 force_area=3*mu*radius,
                 total_comparison_state_error_upper=state_error)
    return {"exact": {k: str(v) for k, v in items.items()},
            "decimal": {k: float(v) for k, v in items.items()},
            "uniform_domain": "all 0 < radius < 1e-5; one displayed shell at 1/200000"}


def rotation(axis, theta):
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0.]])
    rot = (np.cos(theta)*np.eye(3) +
           (1-np.cos(theta))*np.outer(axis, axis) + np.sin(theta)*cross)
    spin = np.cos(theta/2)*I2-1j*np.sin(theta/2)*sdot(axis)
    err = max(norm(spin@sdot(e)@spin.conj().T-sdot(rot@e))
              for e in np.eye(3))
    assert err < 2e-15
    return rot, spin


def mesh():
    # Positive paired angular rule; radial entries diagnose covariance only.
    c, wc = np.polynomial.legendre.leggauss(6)
    k, pol, weight, omega = [], [], [], []
    for frequency, radial_weight in ((1.2, .5), (1.8, .5)):
        for ct, wct in zip(c, wc):
            st = np.sqrt(1-ct*ct)
            for j in range(8):
                phi = 2*np.pi*j/8
                direction = np.array([st*np.cos(phi), st*np.sin(phi), ct])
                e_theta = np.array([ct*np.cos(phi), ct*np.sin(phi), -st])
                e_phi = np.array([-np.sin(phi), np.cos(phi), 0.])
                for e in (e_theta, e_phi):
                    k.append(frequency*direction)
                    pol.append(e)
                    weight.append(radial_weight*wct/16)
                    omega.append(frequency)
    return (np.array(k), np.array(pol), np.array(weight), np.array(omega))


def matrices(k, pol, weight, omega, d):
    # Photon coordinates are mode first, atom spin second.
    blocks = np.array([np.sqrt(w)*np.exp(1j*kv@d)*sdot(e)
                       for kv, e, w in zip(k, pol, weight)])
    a = np.concatenate(list(blocks), axis=1)
    assert norm(a@a.conj().T-2*I2) < 2e-14

    def apply(z):
        ground, excited = z[:-2], z[-2:]
        return np.vstack((np.repeat(omega, 2)[:, None]*ground
                          + .1*a.conj().T@excited,
                          .1*a@ground+1.5*excited))
    return a, apply


def inject(f):
    return np.vstack((np.kron(f[:, None], I2), np.zeros((2, 2), complex)))


def evolve(apply, initial, t=1.):
    result = initial.copy()
    term = initial.copy()
    for order in range(1, 41):
        term = (-1j*t/order)*apply(term)
        result += term
    # Uniform norm <= 2 + sqrt(2)/10 < 2.15 for these finite diagnostics.
    tail = np.exp(2.15)*2.15**41/math.factorial(41)
    assert tail < 1e-34
    return result, float(tail)


def polarized_source(k, pol, weight, n, eps):
    cosine = k@n/np.linalg.norm(k, axis=1)
    win = np.where(cosine > .5, (2*cosine-1)**2, 0.)
    f = np.sqrt(weight)*win*(pol@eps)
    f /= np.linalg.norm(f)
    return f


def spin_output_action(z, u):
    return np.einsum("ab,mbj->maj", u, z.reshape(-1, 2, z.shape[1])).reshape(z.shape)


def finite_covariance():
    k, pol, weight, omega = mesh()
    n0 = np.array([0., 0., 1.])
    eps0 = np.array([1., 1j, 0.])/np.sqrt(2)
    r = 1/200000
    f = polarized_source(k, pol, weight, n0, eps0)
    a, h = matrices(k, pol, weight, omega, r*n0)
    j = inject(f)
    z, tail = evolve(h, j)
    effect = z[-2:].conj().T@z[-2:]
    av = float(np.trace(effect).real/2)
    bv = float(np.trace(effect@sdot(n0)).real/2)
    assert norm(effect-av*I2-bv*sdot(n0)) < 1e-14
    assert bv < -1e-4
    errors = []
    tests = (([1, 0, 0], .73), ([0, 1, 0], np.pi),
             ([1, 2, 3], 1.11), ([-2, 1, 4], -.8))
    for axis, theta in tests:
        rot, u = rotation(axis, theta)
        kr, pr = k@rot.T, pol@rot.T
        nr, er = rot@n0, rot@eps0
        fr = polarized_source(kr, pr, weight, nr, er)
        ar, hr = matrices(kr, pr, weight, omega, r*nr)
        zr, _ = evolve(hr, inject(fr))
        predicted = spin_output_action(z, u)@u.conj().T
        effr = zr[-2:].conj().T@zr[-2:]
        errors.append({
            "full_isometry_covariance": norm(zr-predicted),
            "effect_covariance": norm(effr-u@effect@u.conj().T),
            "effect_axis_formula": norm(effr-av*I2-bv*sdot(nr)),
            "isometry": norm(zr.conj().T@zr-I2)})
    maxima = {key: max(x[key] for x in errors) for key in errors[0]}
    assert max(maxima.values()) < 2e-13
    dark = float(np.linalg.norm(z[-2:]@np.array([1., 0.])))
    assert dark < 1e-14
    # Full output isometry, hence its tensor extension covers arbitrary R.
    # Bell is just one implementation diagnostic, not proof of the quantifier.
    bell_output = z.reshape(-1)/np.sqrt(2)
    assert abs(float(np.vdot(bell_output, bell_output).real)-1) < 2e-13
    return {"scope": "pinned two-radial-node comparison, co-rotated positive quadrature; no finite-mass simulation",
            "photon_modes": len(weight), "rotations": len(tests),
            "a_diagnostic": av, "b_diagnostic": bv,
            "dark_spin_amplitude_diagnostic": dark,
            "analytic_exponential_tail_bound_for_finite_matrix": tail,
            "max_errors": maxima}


def countercontrols():
    k, pol, weight, omega = mesh()
    ez = np.array([0., 0., 1.])
    eps = np.array([1., 1j, 0.])/np.sqrt(2)
    r = 1/200000
    # Same isotropic mixed preparation, three coherent component injections.
    source_components = [np.sqrt(weight)*pol[:, j] for j in range(3)]
    jin = np.concatenate([inject(f)/np.sqrt(2) for f in source_components], axis=1)
    effects, outputs = [], []
    for d in (r*ez, -r*ez, np.zeros(3)):
        _, h = matrices(k, pol, weight, omega, d)
        z, _ = evolve(h, jin)
        ke = [z[-2:, 2*j:2*j+2] for j in range(3)]
        effect = sum(v.conj().T@v for v in ke)
        effects.append(effect)
        outputs.append(ke)
    scalar_residual = max(norm(e-np.trace(e)*I2/2) for e in effects)
    antipodal = norm(effects[0]-effects[1])
    assert max(scalar_residual, antipodal) < 2e-14
    # Scalar effect is not a scalar-amplitude identity instrument.
    rho = np.array([[1, 0], [0, 0]], complex)
    post = sum(v@rho@v.conj().T for v in outputs[-1])
    post /= np.trace(post)
    expected = sum(s@rho@s for s in PAULI)/3
    channel_residual = norm(post-expected)
    assert channel_residual < 2e-13
    assert norm(post-rho) > .6
    # A fixed circular source on the small d ball stays near the same spin axis.
    f = polarized_source(k, pol, weight, ez, eps)
    fixed = []
    for d in (r*ez, -r*ez, np.array([r, 0, 0]), np.array([0, r, 0])):
        _, h = matrices(k, pol, weight, omega, d)
        z, _ = evolve(h, inject(f))
        e = z[-2:].conj().T@z[-2:]
        fixed.append(float(np.trace(e@sdot(ez)).real/2))
    assert max(fixed) < -1e-4
    # At the same zero displacement, different prepared source axes remain
    # distinct. This is a finite diagnostic of the analytic central obstruction.
    _, hz = matrices(k, pol, weight, omega, np.zeros(3))
    z0, _ = evolve(hz, inject(f))
    e0 = z0[-2:].conj().T@z0[-2:]
    rot, u = rotation([0, 1, 0], np.pi)
    kr, pr = k@rot.T, pol@rot.T
    fm = polarized_source(kr, pr, weight, rot@ez, rot@eps)
    _, hm = matrices(kr, pr, weight, omega, np.zeros(3))
    zm, _ = evolve(hm, inject(fm))
    em = zm[-2:].conj().T@zm[-2:]
    central_gap = norm(e0-em)
    assert central_gap > .0017
    return {"unpolarized_scalar_effect_residual": scalar_residual,
            "unpolarized_antipodal_effect_difference": antipodal,
            "scalar_effect_nonidentity_conditional_channel_residual": channel_residual,
            "conditional_output_distance_from_input_operator_norm": norm(post-rho),
            "fixed_source_max_traceless_projection": max(fixed),
            "same_zero_endpoint_different_source_effect_gap": central_gap,
            "scope": "controls exclude automatic source-position identity, not all possible direction measurements"}


def run():
    return {
        "candidate": "aligned_shell_after1062",
        "date": "2026-10-09",
        "passed": True,
        "counting": {"formal_round_assigned": False, "new_cognitive_axioms": 0,
                     "global_scientific_count_unchanged": True},
        "scope": {"new_joint_alignment_contract_explicit": True,
                  "configuration_rotation_permission_adopted_not_derived": True,
                  "full_normal_preparation_domain_not_all_motion_states": True,
                  "adopted_three_dimensional_geometry": True,
                  "continuum_all_direction_proof_in_markdown": True,
                  "no_386_or_425_full_operation_contract_claim": True,
                  "no_full_QED_or_physical_controller_certificate": True},
        "historical_sha256": {name: sha(ROOT/name) for name in HISTORY},
        "rational": rational(),
        "finite_covariance": finite_covariance(),
        "countercontrols": countercontrols(),
    }


def compare(actual, saved, path="result"):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys(), path
        for key in actual:
            compare(actual[key], saved[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(actual) == len(saved), path
        for i, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isfinite(saved), path
        assert abs(actual-saved) < 3e-13*max(1, abs(actual), abs(saved)), path
    else:
        assert actual == saved, (path, actual, saved)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-exclusive", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.save_exclusive:
        with RESULT.open("x", encoding="utf-8") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"candidate": result["candidate"], "passed": True,
                      "mode": "exclusive_first_save" if args.save_exclusive else "read_only",
                      "history": len(HISTORY),
                      "gap_lower": result["rational"]["exact"]["aligned_antipodal_gap_lower"]}))


if __name__ == "__main__":
    main()

