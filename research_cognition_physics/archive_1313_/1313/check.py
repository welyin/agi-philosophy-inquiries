"""1313: full finite invariant sectors, actual vertices and all outcome branches.

No historical simulation is imported. --write-new creates evidence exclusively;
the default and --check recompute and verify the saved evidence without writes.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
ATOL = 7e-11


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def eps3(a, b, c):
    if len({a, b, c}) < 3:
        return 0
    return (-1) ** sum(x > y for x, y in itertools.combinations((a, b, c), 2))


def su(d, rng):
    q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    q[:, -1] /= np.linalg.det(q)
    return q


def make(role, d, m, energies=None):
    """Create the FULL particle-number sector from microscopic vertex indices."""
    n = len(m)
    if role == "C":
        assert d == 3
        left = list(itertools.product(range(d), range(n), range(d), range(d)))
    else:
        left = list(itertools.product(range(d), range(n), range(d)))
    right = left if role == "B" else [(i,) for i in range(n)]
    li, ri = {x: i for i, x in enumerate(left)}, {x: i for i, x in enumerate(right)}
    nl, nr = len(left), len(right)
    t = np.zeros((nr, nl), complex)
    w0, w1 = np.zeros((nl + nr, n), complex), np.zeros((nl + nr, n), complex)
    for idx, key in enumerate(left):
        a, j, b = key[:3]
        if role in ("A", "B"):
            if a == b:
                w0[idx, j] = 1 / np.sqrt(d)
        else:
            w0[idx, j] = eps3(a, b, key[3]) / np.sqrt(6)
        for i in range(n):
            if role == "A":
                t[ri[(i,)], idx] = m[i, j] * (a == b) / np.sqrt(d)
            elif role == "B":
                t[ri[(a, i, b)], idx] = m[i, j]
            else:
                t[ri[(i,)], idx] = m[i, j] * eps3(a, b, key[3]) / np.sqrt(6)
    for idx, key in enumerate(right):
        if role == "B":
            a, i, b = key
            if a == b:
                w1[nl + idx, i] = 1 / np.sqrt(d)
        else:
            w1[nl + idx, key[0]] = 1
    v = np.zeros((nl + nr, nl + nr), complex)
    v[nl:, :nl], v[:nl, nl:] = t, t.conj().T
    ql = np.diag(np.r_[np.ones(nl), np.zeros(nr)])
    qr = np.eye(nl + nr) - ql
    qb = np.eye(nl + nr) if role == "B" else ql.copy()
    q2 = ql.copy() if role == "C" else np.zeros_like(ql)
    e = energies or {"L": 1.1, "b": 0.7, "b2": 0.4,
                     "R": {"A": 1.8, "B": 1.1, "C": 2.2}[role]}
    q = {"L": ql, "R": qr, "b": qb, "b2": q2}
    h0 = sum(e[k] * q[k] for k in q)
    return {"v": v, "h0": h0, "w0": w0, "w1": w1, "q": q,
            "nl": nl, "nr": nr, "t": t, "energy": e}


def representation(role, g, n):
    gl = np.kron(np.kron(g, np.eye(n)), g.conj() if role != "C" else g)
    if role == "C":
        gl = np.kron(gl, g)
    gr = gl if role == "B" else np.eye(n)
    z = np.zeros((len(gl) + len(gr), len(gl) + len(gr)), complex)
    z[:len(gl), :len(gl)], z[len(gl):, len(gl):] = gl, gr
    return z


def evolution(v, t):
    eig, basis = np.linalg.eigh(v)
    return (basis * np.exp(-1j * t * eig)) @ basis.conj().T


def logical(m, t):
    u, s, vh = np.linalg.svd(m)
    a = (vh.conj().T * np.cos(t * s)) @ vh
    b = -1j * (u * np.sin(t * s)) @ vh
    return a, b, s


def run():
    rng = np.random.default_rng(1313)
    maxima, cases, source_cases = {}, [], []

    def eq(name, a, b, tol=ATOL):
        err = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
        maxima[name] = max(maxima.get(name, 0.0), err)
        if err > tol:
            raise AssertionError((name, err, tol))

    configurations = [(r, d, n) for r in ("A", "B")
                      for d, n in ((2, 3), (3, 4), (4, 5))]
    configurations += [("C", 3, n) for n in (2, 3, 4, 5)]
    for role, d, n in configurations:
        m = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2 * n)
        x = make(role, d, m)
        v, h0, w0, w1 = (x[k] for k in ("v", "h0", "w0", "w1"))
        w = np.hstack([w0, w1])
        h = np.block([[np.zeros((n, n)), m.conj().T], [m, np.zeros((n, n))]])
        eq("isometries", w.conj().T @ w, np.eye(2 * n))
        eq("vertex_intertwining", v @ w, w @ h)
        eq("hermitian", v, v.conj().T)
        eq("free_energy_conservation", h0 @ v, v @ h0)
        eq("fermion_number", x["q"]["L"] + x["q"]["R"], np.eye(len(v)))
        if role != "B":
            eq("reference_number_constraint", x["q"]["b"] + x["q"]["R"], np.eye(len(v)))
        if role == "C":
            eq("second_reference_constraint", x["q"]["b2"] + x["q"]["R"], np.eye(len(v)))
        gr = representation(role, su(d, rng), n)
        eq("group_invariance", gr @ v, v @ gr)
        eq("invariant_encoding", gr @ w, w)
        rank = int(np.linalg.matrix_rank(v, tol=1e-10))
        expected_rank = 2 * n * (d * d if role == "B" else 1)
        assert rank == expected_rank
        psi = rng.normal(size=(n, n + 1)) + 1j * rng.normal(size=(n, n + 1))
        psi /= np.linalg.norm(psi)
        for time in (0.0, 0.31, 0.87, 1.57):
            u = evolution(v, time)
            a, b, singular = logical(m, time)
            eq("full_unitarity", u.conj().T @ u, np.eye(len(v)))
            eq("full_evolution", u @ w0, w0 @ a + w1 @ b)
            eq("instrument_complete", a.conj().T @ a + b.conj().T @ b, np.eye(n))
            full = u @ w0 @ psi
            fail, succ = x["q"]["L"] @ full, x["q"]["R"] @ full
            eq("partner_failure", w0.conj().T @ fail, a @ psi)
            eq("partner_success", w1.conj().T @ succ, b @ psi)
            eq("branches_sum", np.vdot(fail, fail) + np.vdot(succ, succ), 1)
            eq("failure_partner", fail.T @ fail.conj(), (a @ psi).T @ (a @ psi).conj())
            p = float(np.vdot(succ, succ).real)
            for key, expect in (("L", 1 - p), ("R", p),
                                ("b", 1 if role == "B" else 1 - p),
                                ("b2", 1 - p if role == "C" else 0)):
                eq("resource_" + key, np.vdot(full, x["q"][key] @ full), expect)
            energy = h0 + v
            initial = w0 @ psi
            eq("full_energy", np.vdot(full, energy @ full), np.vdot(initial, energy @ initial))
            eq("interaction_energy", np.vdot(full, v @ full), 0)
            eq("readout_mean_energy", np.vdot(fail, energy @ fail) + np.vdot(succ, energy @ succ),
               np.vdot(full, energy @ full))
            if role == "B":
                right = succ[x["nl"]:].reshape(d, n, d, n + 1)
                # Trace the ACTUAL right activity and reference indices.
                reduced = np.einsum("aibz,akbw->izkw", right, right.conj())
                target = np.einsum("iz,kw->izkw", b @ psi, (b @ psi).conj())
                eq("receiver_partial_trace", reduced, target)
            else:
                eq("receiver_single_particle", succ[x["nl"]:], b @ psi)
            eq("worst_success", np.linalg.eigvalsh(b.conj().T @ b).min(),
               np.sin(time * singular).min() ** 2 if np.all(time * singular <= np.pi / 2)
               else np.min(np.sin(time * singular) ** 2))
        cases.append({"role": role, "d": d, "n": n, "full_dimension": len(v),
                      "vertex_rank": rank, "times": 4})

    for role in ("A", "B", "C"):
        d, n = 3, 2
        m = np.array([[1.0 + .2j, .3j], [.4, .8 - .1j]])
        x = make(role, d, m)
        # Independent source coordinates BEFORE imposing the resonance equation.
        delta = 1e-5
        for key in ("L", "R", "b") + (("b2",) if role == "C" else ()):
            ep, em = dict(x["energy"]), dict(x["energy"])
            ep[key] += delta
            em[key] -= delta
            derivative = (make(role, d, m, ep)["h0"] - make(role, d, m, em)["h0"]) / (2 * delta)
            eq("independent_energy_sources", derivative, x["q"][key])
        basis = np.zeros_like(m)
        basis[0, 1] = 1
        real_source = make(role, d, basis)["v"]
        imag_source = make(role, d, 1j * basis)["v"]
        eq("real_coupling_source", (make(role, d, m + delta * basis)["v"] -
           make(role, d, m - delta * basis)["v"]) / (2 * delta), real_source)
        eq("imag_coupling_source", (make(role, d, m + 1j * delta * basis)["v"] -
           make(role, d, m - 1j * delta * basis)["v"]) / (2 * delta), imag_source)
        source_cases.append(role)

    perfect = []
    for n in (3, 4, 5):
        u0 = su(n, rng)
        m = .8 * u0
        for role in ("A", "B", "C"):
            x = make(role, 3, m)
            evolution_map = evolution(x["v"], np.pi / 1.6) @ x["w0"]
            eq("perfect_all_input", evolution_map, -1j * x["w1"] @ u0)
            choi = evolution_map / np.sqrt(n)
            eq("perfect_partner", x["w1"].conj().T @ choi, -1j * u0 / np.sqrt(n))
            perfect.append({"role": role, "n": n, "p_min": 1.0})

    spectra = []
    for values in ([1., 2.], [1., 3.], [0., 1.], [1., 1.]):
        m, t = np.diag(values).astype(complex), np.pi / 2
        a, b, _ = logical(m, t)
        for role in ("A", "B", "C"):
            x = make(role, 3, m)
            eq("spectrum_instrument", evolution(x["v"], t) @ x["w0"],
               x["w0"] @ a + x["w1"] @ b)
        spectra.append({"singular_values": values, "rank": int(np.linalg.matrix_rank(m)),
                        "p_min": float(np.min(np.linalg.eigvalsh(b.conj().T @ b)))})

    timing = []
    for error in (0.0, .01, .1, .4):
        n = 3
        a, b, _ = logical(np.eye(n), np.pi / 2 + error)
        # Choi input lower bound: orthogonal erasure block saturates analytic upper bound.
        psi = np.eye(n).reshape(-1) / np.sqrt(n)
        ideal = np.outer(psi, psi.conj())
        output = np.sin(np.pi / 2 + error) ** 2 * ideal
        vacuum = (a / np.sqrt(n)).T @ (a / np.sqrt(n)).conj()
        distance = (np.abs(np.linalg.eigvalsh(output - ideal)).sum() +
                    np.abs(np.linalg.eigvalsh(vacuum)).sum()) / 2
        eq("timing_half_diamond_certificate", distance, np.sin(error) ** 2)
        timing.append({"angle_error": error, "half_diamond": float(distance)})

    return {"schema": "round1313_full_sector_v1", "status": "PASS", "seed": 1313,
            "tolerance": ATOL, "full_sector_cases": cases, "source_cases": source_cases,
            "perfect_handoffs": perfect, "spectral_contrasts": spectra, "timing": timing,
            "max_errors": maxima, "new_axioms": 0,
            "scope": "finite logical K handoff, not arbitrary bare activity nor complete SM",
            "bindings": {"check.py": sha(HERE / "check.py"), "proof.md": sha(HERE / "proof.md"),
                         "../_admission/round1313_admission.md": sha(HERE.parent / "_admission/round1313_admission.md")}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    current = run()
    dest = HERE / "results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(current, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        saved = json.loads(dest.read_text(encoding="utf-8"))
        assert saved["status"] == current["status"] == "PASS"
        assert saved["bindings"] == current["bindings"]
        for k in ("full_sector_cases", "source_cases", "perfect_handoffs", "scope", "new_axioms"):
            assert saved[k] == current[k], k
        for item in saved["max_errors"].values():
            assert item <= ATOL
        for a, b in zip(saved["spectral_contrasts"], current["spectral_contrasts"]):
            assert a["singular_values"] == b["singular_values"] and a["rank"] == b["rank"]
            assert abs(a["p_min"] - b["p_min"]) < ATOL
    print(json.dumps({"status": "PASS", "full_sector_cases": len(current["full_sector_cases"]),
                      "perfect_handoffs": len(current["perfect_handoffs"]),
                      "max_error": max(current["max_errors"].values()), "saved": args.write_new}))


if __name__ == "__main__":
    main()
