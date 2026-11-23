"""Independent finite diagnostics for round 1333. No author-code imports.

All-n statements are proved in the review/proof, not inferred from sampling.
--write-new refuses an existing result; --check is read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
TOL = 3e-10


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shift(n):
    ans = np.zeros((n, n), complex)
    for k in range(n - 1):
        ans[k + 1, k] = 1
    return ans


def spectra(y):
    # Independent Gram-eigenvalue route; only the two largest are used for bounds.
    return np.sqrt(np.maximum(np.linalg.eigvalsh(y.conj().T @ y), 0))[::-1]


def op(y):
    return float(spectra(y)[0])


def unitary(rng, n):
    q, r = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return q @ np.diag(np.diag(r) / abs(np.diag(r)))


def invariants(yu, yd):
    hu, hd = yu @ yu.conj().T, yd @ yd.conj().T
    n = len(yu)
    aa = float(np.linalg.det(hu).real ** (1 / n))
    cc = float(np.linalg.det(hd).real ** (1 / n))
    bb = float((np.trace(hu).real - n * aa) / (n - 1))
    dd = float((np.trace(hd).real - n * cc) / (n - 1))
    cross = np.trace(hu @ hd).real
    return np.array([aa, bb, cc, dd, cross])


def run():
    rng = np.random.default_rng(13330029)
    count = 0
    largest_residual = 0.0

    def close(a, b, label):
        nonlocal count, largest_residual
        aa, bb = np.asarray(a), np.asarray(b)
        err = float(np.max(abs(aa - bb))) / (1 + float(np.max(abs(bb))))
        count += 1
        largest_residual = max(largest_residual, err)
        if err > TOL:
            raise AssertionError((label, err))

    def positive(x, label):
        nonlocal count
        count += 1
        if x < -TOL:
            raise AssertionError((label, x))

    spectrum_cases = []
    coefficients = [(0j, 1j), (1j, 0j), (1 + .7j, .03 - .04j),
                    (.01j, 2 - .3j), (1 + 0j, 1j)]
    for _ in range(12):
        coefficients.append(tuple(rng.normal(size=2) + 1j * rng.normal(size=2)))
    for n in (3, 4, 6, 9):
        s = shift(n)
        for a, b in coefficients:
            for reverse in (False, True):
                y = a * np.eye(n) + b * (s.conj().T if reverse else s)
                ss = spectra(y)
                aa, bb = abs(a), abs(b)
                # Forward first two input coordinates; reverse last two.
                e = np.eye(n, dtype=complex)[:, -2:] if reverse else np.eye(n, dtype=complex)[:, :2]
                compressed = e.conj().T @ y.conj().T @ y @ e
                close(np.linalg.eigvalsh(compressed),
                      [aa * aa + bb * bb - aa * bb, aa * aa + bb * bb + aa * bb], "compression")
                lower = np.sqrt(aa * aa + bb * bb - aa * bb)
                positive((ss[1] - lower) / (1 + ss[0]), "minmax")
                positive((aa + bb - ss[0]) / (1 + ss[0]), "norm upper")
                positive((lower - (aa + bb) / 2) / (1 + ss[0]), "elementary lower")
                positive(ss[1] / ss[0] - .5, "ratio")
                close(np.linalg.det(y), a ** n, "triangular determinant")
                spectrum_cases.append({"n": n, "reverse": reverse, "ratio": float(ss[1] / ss[0])})

    # Explicit exceptions, including ratio undefined for the zero matrix.
    n2 = spectra(.01 * np.eye(2) + shift(2))
    positive(.001 - n2[1] / n2[0], "n2 disproves all-n extension")
    close(spectra(np.zeros((3, 3))), np.zeros(3), "zero matrix")
    close(spectra(shift(3)), [1, 1, 0], "zero diagonal")
    close(spectra(2 * np.eye(3)), [2, 2, 2], "zero step")

    perturbation_cases = []
    for n in (3, 4, 7):
        y = (.7 + .2j) * np.eye(n) + (1.1 - .3j) * shift(n)
        old = spectra(y)
        for eps in (0, .01, .1, .3, .49):
            for _ in range(3):
                raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
                e = raw * (eps * old[0] / op(raw))
                now = spectra(y + e)
                close(op(e), eps * old[0], "perturbation norm")
                positive(now[1] - old[1] + op(e), "singular lower Weyl")
                positive(old[0] + op(e) - now[0], "singular upper Weyl")
                bound = (.5 - eps) / (1 + eps)
                positive(now[1] / now[0] - bound, "perturbed ratio")
                perturbation_cases.append({"n": n, "epsilon": eps, "ratio": float(now[1] / now[0]), "bound": bound})

    phase_cases = []
    for n in (3, 4, 6):
        s = shift(n)
        for _ in range(6):
            mags = rng.uniform(.6, 1.6, size=4)
            phases = rng.uniform(-2, 2, size=4)
            a, b, c, d = mags * np.exp(1j * phases)
            yu, yd = a * np.eye(n) + b * s, c * np.eye(n) + d * s.conj().T
            aa, bb, cc, dd = mags ** 2
            z = a.conjugate() * b * c.conjugate() * d
            hu, hd = yu @ yu.conj().T, yd @ yd.conj().T
            comm = hu @ hd - hd @ hu
            expected = 6j * z.imag * abs(b * d) ** 2 * (abs(b * c) ** 2 - abs(a * d) ** 2)
            close(np.trace(comm @ comm @ comm), expected, "CP formula")
            cross = n * aa * cc + (n - 1) * (aa * dd + bb * cc) + (n - 2) * bb * dd + 2 * (n - 1) * z.real
            close(invariants(yu, yd), [aa, bb, cc, dd, cross], "five invariants")
            # Explicit common-left rephasing and independent scalar right rephasings.
            theta = phases[0] - phases[1]
            q = np.diag(np.exp(1j * theta * np.arange(n)))
            newu = np.exp(-1j * phases[0]) * q @ yu @ q.conj().T
            newd = np.exp(-1j * phases[2]) * q @ yd @ q.conj().T
            phi = phases[1] + phases[3] - phases[0] - phases[2]
            close(newu, mags[0] * np.eye(n) + mags[1] * s, "phase reduction u")
            close(newd, mags[2] * np.eye(n) + mags[3] * np.exp(1j * phi) * s.conj().T, "phase reduction d")
            ul, ur, dr = [unitary(rng, n) for __ in range(3)]
            uy, dy = ul @ yu @ ur.conj().T, ul @ yd @ dr.conj().T
            close(invariants(uy, dy), invariants(yu, yd), "full Dirac basis invariance")
            h1, h2 = uy @ uy.conj().T, dy @ dy.conj().T
            c2 = h1 @ h2 - h2 @ h1
            close(np.trace(c2 @ c2 @ c2), expected, "basis invariant CP")
            phase_cases.append({"n": n, "CP_imaginary": float(expected.imag)})

    # A concrete exact CP-odd witness and its n=2 exception.
    cp_examples = []
    for n in (2, 3, 4, 8):
        s = shift(n)
        yu, yd = np.eye(n) + s, 2 * np.eye(n) + 1j * s.conj().T
        hu, hd = yu @ yu.conj().T, yd @ yd.conj().T
        c = hu @ hd - hd @ hu
        val = np.trace(c @ c @ c)
        close(val, 0j if n == 2 else 36j, "exact CP witness")
        cp_examples.append({"n": n, "trace_commutator_cubed_imaginary": float(val.imag)})
    # Generic five-dimensional invariant image: analytic Jacobian determinant.
    n = 3
    aa, bb, cc, dd, phi = 1.2, .7, 1.6, .9, .8
    jac = np.eye(5)
    jac[4, :4] = 0  # Other row entries do not alter this triangular determinant.
    jac[4, 4] = -2 * (n - 1) * np.sqrt(aa * bb * cc * dd) * np.sin(phi)
    positive(abs(np.linalg.det(jac)) - .1, "generic fifth invariant varies")

    # With noncanonical kinetic metric this is a different model, not a counterexample.
    masses = np.array([.001, .02, 1.0])
    kinetic_left = np.diag(1 / masses ** 2)
    canonical_y = np.diag(1 / np.sqrt(np.diag(kinetic_left))) @ np.eye(3)
    close(spectra(canonical_y), masses[::-1], "scope: noncanonical metric")

    historical = {}
    for label, central, upper in (
        ("charm_top", .619 / 171.7, (.619 + .084) / (171.7 - 3)),
        ("strange_bottom", .055 / 2.89, (.055 + .016) / (2.89 - .09)),
    ):
        required = (.5 - upper) / (1 + upper)
        positive(.5 - upper, "historical same-scale separation")
        close((.5 - required) / (1 + required), upper, "required perturbation inversion")
        historical[label] = {"central_ratio": central, "rectangle_upper_ratio": upper,
                             "necessary_relative_change_lower": required}

    return {
        "status": "PASS", "checks": count, "max_normalized_equality_residual": largest_residual,
        "input_hashes": {"independent_check.py": sha(Path(__file__)),
                         "round1333_admission.md": sha(HERE.parent / "_admission/round1333_admission.md")},
        "scope": "Necessary-condition rejection only; no quantum apparatus or common formation certified.",
        "spectrum_cases": spectrum_cases, "n2_ratio": float(n2[1] / n2[0]),
        "zero_matrix_ratio": None, "perturbation_cases": perturbation_cases,
        "phase_cases": phase_cases, "CP_examples": cp_examples,
        "generic_invariant_jacobian_det": float(np.linalg.det(jac)),
        "historical_2008_MSbar_MZ": {"source": "https://arxiv.org/html/0712.1419", "table": "II, MZ row",
              "warning": "Historical input illustration, not a new/current fit or statistical confidence bound.", **historical},
    }


def compare(saved, fresh, route="root"):
    if isinstance(fresh, dict):
        if saved.keys() != fresh.keys():
            raise AssertionError((route, "keys"))
        for k in fresh:
            compare(saved[k], fresh[k], route + "." + k)
    elif isinstance(fresh, list):
        if len(saved) != len(fresh):
            raise AssertionError((route, "length"))
        for i, (a, b) in enumerate(zip(saved, fresh)):
            compare(a, b, route + f"[{i}]")
    elif isinstance(fresh, float):
        if not np.isclose(saved, fresh, rtol=3e-9, atol=3e-12):
            raise AssertionError((route, saved, fresh))
    elif saved != fresh:
        raise AssertionError((route, saved, fresh))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write_new:
        with OUT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        compare(json.loads(OUT.read_text(encoding="utf-8")), result)
    print(json.dumps({k: result[k] for k in ("status", "checks", "max_normalized_equality_residual")}, ensure_ascii=False))
