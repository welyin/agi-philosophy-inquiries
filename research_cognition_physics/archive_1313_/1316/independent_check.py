"""Independent finite-matrix checks for round 1316; author code is not imported.

Default and --check are read-only. --write-new exclusively creates this reviewer's
result. All samples are diagnostics of the stated finite tree-level contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ADMISSION = HERE.parent / "_admission" / "round1316_admission.md"
RESULT = HERE / "independent_results.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def herm(a):
    return (a + a.conj().T) / 2


def norm(a):
    return float(np.linalg.norm(a, 2))


def relative(a, b):
    return norm(np.atleast_2d(a - b)) / max(1.0, norm(np.atleast_2d(b)))


def invsqrt(a):
    w, u = np.linalg.eigh(herm(a))
    if w[0] <= 0:
        raise ValueError("Whitening requires positive definite input")
    return (u * w**-0.5) @ u.conj().T


def generalized(a, b):
    s = invsqrt(b)
    return np.linalg.eigvalsh(herm(s @ a @ s))


def unitary(rng, n):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return np.linalg.qr(z)[0]


def portal(rng, n):
    return unitary(rng, n) @ np.diag(np.linspace(0.7, 1.4, n)) @ unitary(rng, n)


def chain_matrix(length, q):
    a = np.zeros((length, length + 1))
    for j in range(length):
        a[j, j] = 1
        a[j, j + 1] = -q
    return a


def assembly(lengths, q, mass):
    n, h = len(lengths), sum(lengths)
    t = np.zeros((h, h + n))
    p = np.zeros((n, h + n))
    z = np.zeros((h + n, n))
    fs, ss, row, col = [], [], 0, 0
    for a, length in enumerate(lengths):
        qq = chain_matrix(length, q)
        t[row:row + length, col:col + length + 1] = mass * qq
        p[a, col + length] = 1
        vec = q**(-np.arange(length + 1, dtype=float))
        vec /= np.linalg.norm(vec)
        z[col:col + length + 1, a] = vec
        fs.append(float(vec[-1]))
        ki = np.linalg.inv(qq @ qq.T)
        ss.append(float(q*q / (mass*mass) * (ki @ ki)[-1, -1]))
        row += length
        col += length + 1
    # Derive heavy subspace from the whole right mass matrix, not closed modes.
    vals, vecs = np.linalg.eigh(t.T @ t)
    heavy_vals, heavy_vecs = vals[n:], vecs[:, n:]
    return t, p, z, np.array(fs), np.array(ss), heavy_vals, heavy_vecs


class Checks:
    def __init__(self):
        self.values = {}
        self.count = 0

    def close(self, name, a, b, tol=2e-9):
        error = relative(np.asarray(a), np.asarray(b))
        self.values[name] = max(self.values.get(name, 0.0), error)
        self.count += 1
        if error > tol:
            raise AssertionError((name, error, tol))

    def positive(self, name, a, tol=2e-9):
        violation = max(0.0, -float(np.linalg.eigvalsh(herm(a))[0])) / max(1., norm(a))
        self.values[name] = max(self.values.get(name, 0.0), violation)
        self.count += 1
        if violation > tol:
            raise AssertionError((name, violation, tol))

    def truth(self, name, value):
        self.count += 1
        if not bool(value):
            raise AssertionError(name)
        self.values.setdefault(name, 0.0)


def exercise(rng, lengths, q, mass, target, checks):
    n = len(lengths)
    t, p, z0, fs, ss, hv, vv = assembly(lengths, q, mass)
    y0 = portal(rng, n)
    y = y0 @ np.diag(fs)
    d = herm(y @ y.conj().T)
    gap = float(hv[0])
    rho = float(np.sqrt(target * gap / norm(d)))
    checks.close("zero_mode_kernel", t @ z0, np.zeros_like(t @ z0))
    checks.close("zero_mode_orthonormal", z0.T @ z0, np.eye(n))
    checks.close("endpoint_overlap", p @ z0, np.diag(fs))
    checks.close("heavy_complement", vv @ vv.T + z0 @ z0.T, np.eye(t.shape[1]))
    analytic_s = []
    for length in lengths:
        if q == 1:
            analytic_s.append(length*(2*length+1)/(6*mass*mass*(length+1)))
        else:
            powers = np.arange(1, length+1, dtype=float)
            analytic_s.append(np.sum((q**powers-q**(-powers))**2) /
                              (mass*mass*(q**(length+1)-q**(-(length+1)))**2))
    checks.close("closed_form_heavy_coefficient_including_q1", ss, np.array(analytic_s))

    # Independent whole-matrix Moore-Penrose inverse response.
    right_pinv = np.linalg.pinv(t.T @ t, rcond=1e-11, hermitian=True)
    c0 = herm(y0 @ p @ right_pinv @ p.T @ y0.conj().T)
    checks.close("matched_heavy_response", c0, y0 @ np.diag(ss) @ y0.conj().T)
    expected = np.sort(ss / fs**2)
    checks.close("portal_independent_generalized_spectrum", generalized(c0, d), expected)
    vleft = unitary(rng, n)
    checks.close("left_chart_covariance", generalized(vleft @ c0 @ vleft.conj().T,
                                                       vleft @ d @ vleft.conj().T), expected)
    another = portal(rng, n)
    dy = another @ np.diag(fs**2) @ another.conj().T
    cy = another @ np.diag(ss) @ another.conj().T
    checks.close("different_portal_same_spectrum", generalized(cy, dy), expected)

    def c_at(x):
        return herm(y0 @ p @ ((vv / (hv - x)) @ vv.T) @ p.T @ y0.conj().T)

    full = np.vstack((rho * y0 @ p, t))
    eigs = np.linalg.eigvalsh(herm(full @ full.conj().T))
    low = eigs[:n]
    base = t.T @ t
    update = full.conj().T @ full - base
    checks.close("same_parent_positive_update", update, rho*rho*p.T @ y0.conj().T @ y0 @ p)
    checks.positive("positive_update", update)
    checks.truth("exactly_n_positive_low_poles", low[0] > 0 and low[-1] < gap)
    checks.truth("minmax_low_bound", low[-1] <= target*gap + 2e-9 * max(1., gap))
    checks.truth("minmax_heavy_bound", eigs[n] >= gap - 2e-9 * max(1., gap))
    checks.close("left_right_mass_spectra", eigs, np.linalg.eigvalsh(herm(full.conj().T @ full)))
    lower = generalized(rho*rho*d, np.eye(n) + rho*rho*c0/(1-target))
    upper = generalized(rho*rho*d, np.eye(n) + rho*rho*c0)
    checks.truth("ordered_pole_lower_bound", np.all(low >= lower - 2e-9 * max(1., gap)))
    checks.truth("ordered_pole_upper_bound", np.all(low <= upper + 2e-9 * max(1., gap)))
    checks.truth("bare_Y_is_only_upper_bound", np.all(low <= np.linalg.eigvalsh(rho*rho*d) + 2e-9*max(1.,gap)))

    for x in np.linspace(0, target*gap, 6):
        cx = c_at(x)
        checks.positive("resolvent_lower_Loewner", cx-c0)
        checks.positive("resolvent_upper_Loewner", c0/(1-target)-cx)
        # Direct elimination from the full left block, a route independent of the
        # author's claimed resolvent expression.
        top, off, bottom = full[:n] @ full[:n].conj().T, full[:n] @ full[n:].conj().T, t @ t.T
        schur = top-x*np.eye(n)-off @ np.linalg.solve(bottom-x*np.eye(len(bottom)), off.conj().T)
        claimed = rho*rho*d-x*(np.eye(n)+rho*rho*cx)
        checks.close("exact_Schur_equation", schur, claimed, 1e-8)
    for x in low:
        k = rho*rho*d-x*(np.eye(n)+rho*rho*c_at(x))
        residual = np.linalg.svd(k, compute_uv=False)[-1]/max(1., norm(rho*rho*d))
        checks.close("all_exact_low_poles_solve_equation", np.array([residual]), np.zeros(1), 2e-8)
    for x in [-0.8*gap, -0.2*gap, 0.2*gap, 0.8*gap]:
        cwhiten = invsqrt(c0)
        remainder = norm(cwhiten @ (c_at(x)-c0) @ cwhiten)
        checks.truth("two_sided_normalized_resolvent_bound", remainder <= abs(x)/(gap-abs(x))+2e-9)

    desired = portal(rng, n)
    reconstructed_portal = desired @ np.diag(1/fs)
    checks.close("arbitrary_invertible_light_Y_still_free", reconstructed_portal @ np.diag(fs), desired)
    t2, p2, z2, f2, s2, _, _ = assembly(lengths, q, mass*2.3)
    checks.close("mass_scale_leaves_light_Y", f2, fs)
    checks.close("mass_scale_changes_heavy_response", s2, ss/(2.3**2))

    # A physical assumption removal: move first portal to the opposite endpoint.
    pp = p.copy()
    pp[0] = 0
    pp[0, 0] = 1
    new_d = herm(y0 @ pp @ z0 @ z0.T @ pp.T @ y0.conj().T)
    new_c = herm(y0 @ pp @ right_pinv @ pp.T @ y0.conj().T)
    endpoint_violation = relative(generalized(new_c, new_d), expected)
    if q == 1:
        checks.close("q1_opposite_endpoint_reflection_symmetry", generalized(new_c,new_d), expected)
    else:
        checks.truth("changed_endpoint_not_old_relation", endpoint_violation > 1e-5)
    raw = rng.normal(size=(n, t.shape[1])) + 1j*rng.normal(size=(n, t.shape[1]))
    extra_portal = raw @ (np.eye(t.shape[1])-z0@z0.T)
    full_portal = y0@p + extra_portal
    checks.close("multiport_preserves_same_light_Y", full_portal@z0, y)
    multi_c = herm(full_portal @ right_pinv @ full_portal.conj().T)
    checks.truth("multiport_changes_response_at_same_light_Y", relative(generalized(multi_c,d),expected) > 1e-5)
    extra = c0 + 0.7*norm(c0)*np.eye(n)
    checks.truth("extra_same_order_response_not_old_relation", relative(generalized(extra,d),expected) > 1e-5)
    singular_portal = y0.copy()
    singular_portal[:, -1] = 0
    checks.truth("singular_portal_excluded", np.linalg.matrix_rank(singular_portal @ np.diag(fs)) < n)

    # Parent-source derivatives: differentiate the whole mass-squared matrix.
    delta = 1e-5
    d_full_dm = np.vstack((np.zeros_like(full[:n]), t/mass))
    def source_mat(mm):
        tt = np.vstack((rho*y0@p, (mm/mass)*t))
        return tt.conj().T @ tt
    numerical = (source_mat(mass+delta)-source_mat(mass-delta))/(2*delta)
    exact = d_full_dm.conj().T @ full + full.conj().T @ d_full_dm
    checks.close("complete_mass_source_derivative", numerical, exact, 2e-8)
    return {
        "n": n, "lengths": lengths, "q": q, "m": mass,
        "rho": rho, "zeta": target, "unmixed_heavy_gap_squared": gap,
        "rho_squared_C0_norm": rho*rho*norm(c0),
        "light_poles_squared": low.tolist(),
        "lower_bounds": lower.tolist(), "upper_bounds": upper.tolist(),
        "generalized_spectrum": expected.tolist(),
        "changed_endpoint_relative_violation": endpoint_violation,
    }


def run():
    rng = np.random.default_rng(13160891)
    checks = Checks()
    samples = []
    for lengths, q, mass in [([1, 2, 4], 1.4, 1.7),
                              ([1, 2, 3, 5], 1.7, 2.1),
                              ([1, 2, 3, 4, 6], 2.0, 2.3),
                              ([1, 3, 6], 1.0, 1.3),
                              ([4, 5, 6], 2.0, 1.8)]:
        for target in [0.01, 0.1, 0.45, 0.9]:
            samples.append(exercise(rng, lengths, q, mass, target, checks))
    # Four-component convention: phases of barred fields are conjugated.
    l, higgs, nr, sl = 1j, -1, -1j, -1j
    checks.close("portal_r_protection", np.array([l.conjugate()*higgs*nr]), np.ones(1))
    checks.close("chain_r_protection", np.array([sl.conjugate()*nr]), np.ones(1))
    checks.truth("Majorana_r_forbidden", nr**2 == -1 and sl**2 == -1)
    return {
        "status": "PASS_WITH_SCOPE", "round": 1316,
        "scope": "Finite canonical Dirac chains; exact quadratic spectrum and tree heavy response; no instrument, loop, generation-selection or cognitive-necessity claim.",
        "independent_of_author_code": True,
        "bindings": {"independent_check.py": sha(__file__), "admission": sha(ADMISSION)},
        "sample_count": len(samples), "assertion_count": checks.count,
        "check_categories": len(checks.values),
        "max_diagnostic_error": max(checks.values.values()),
        "diagnostic_errors": checks.values, "samples": samples,
        "strong_mixing_case_max_rho_squared_C0": max(s["rho_squared_C0_norm"] for s in samples),
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write_new:
        with RESULT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
    elif RESULT.exists():
        old = json.loads(RESULT.read_text(encoding="utf-8"))
        if old != result:
            raise AssertionError("Saved independent result differs from recomputation")
    print(json.dumps({k: result[k] for k in ["status", "sample_count", "assertion_count", "check_categories", "max_diagnostic_error", "strong_mixing_case_max_rho_squared_C0"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
