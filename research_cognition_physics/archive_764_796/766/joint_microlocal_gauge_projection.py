"""766: finite calibrations of the exact, frequency-compatible gauge repair.

Pseudodifferential order and wavefront properties are analytic in note766.
No matrix sample proves a continuum Hadamard state or physical positivity.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_covariant_gauge_complex as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/"joint_microlocal_gauge_projection_results.json"


def maximum(x):
    return float(np.max(np.abs(x)))


def block(a, b):
    return np.block([[a, np.zeros((a.shape[0], b.shape[1]))],
                     [np.zeros((b.shape[0], a.shape[1])), b]])


def adj(a, q_out, q_in):
    return np.linalg.solve(q_in, a.conj().T@q_out)


def auxiliary(c):
    minus = np.eye(c.shape[0])-c
    return c.conj().T@c+minus.conj().T@minus


def make_projection(K, c1, c0):
    h1, h0 = auxiliary(c1), auxiliary(c0)
    sharp = np.linalg.solve(h0, K.conj().T@h1)
    gram = sharp@K
    B = np.linalg.solve(gram, sharp)
    Pi = np.eye(K.shape[0])-K@B
    return B, Pi, h1, h0


def repair(c1, c0, K, B, Pi, q1, q0):
    kd = adj(K, q1, q0); bd = adj(B, q0, q1); pd = adj(Pi, q1, q1)
    R = c1@K-K@c0
    ct = pd@c1@Pi+bd@c0@kd+K@c0@B
    correction = bd@adj(R, q1, q0)+pd@R@B
    return ct, correction


def original_principal_setup():
    phi = np.array([.13, .64, -.11, .08, .37])
    weights = np.repeat([.73, 1.17, .41], [8, 3, 1])
    direction = np.array([1., 2., -3.]); direction /= np.linalg.norm(direction)
    H, qg, kp, *_ = old.principal(np.r_[1., direction], phi, weights)
    _, _, km, *_ = old.principal(np.r_[-1., direction], phi, weights)
    K = block(kp, km).astype(complex)
    q1, q0 = block(H, -H).astype(complex), block(qg, -qg).astype(complex)
    c1 = np.diag(np.r_[np.ones(63), np.zeros(63)]).astype(complex)
    c0 = np.diag(np.r_[np.ones(16), np.zeros(16)]).astype(complex)
    # A nonorthogonal change of Cauchy frame ensures that adjoints matter.
    rng = np.random.default_rng(766)
    m1 = np.eye(126, dtype=complex)
    m1[:63, 63:] = .025*rng.normal(size=(63, 63))
    m0 = np.eye(32, dtype=complex)
    m0[:16, 16:] = .06*rng.normal(size=(16, 16))
    i1, i0 = np.linalg.inv(m1), np.linalg.inv(m0)
    return (m1@K@i0, m1@c1@i1, m0@c0@i0,
            i1.conj().T@q1@i1, i0.conj().T@q0@i0)


def projection_check():
    K, c1, c0, q1, q0 = original_principal_setup()
    B, Pi, h1, h0 = make_projection(K, c1, c0)
    kd = adj(K, q1, q0)
    _, _, vh = np.linalg.svd(kd)
    N = vh.conj().T[:, K.shape[1]:]
    errors = dict(left_inverse=maximum(B@K-np.eye(32)),
                  projection=maximum(Pi@Pi-Pi), gauge_removed=maximum(Pi@K),
                  frequency_commutator=maximum(Pi@c1-c1@Pi),
                  auxiliary_self_adjoint=maximum(Pi.conj().T@h1-h1@Pi),
                  constraint_preserved=maximum(kd@Pi@N),
                  physical_charge_preserved=maximum(N.conj().T@(Pi.conj().T@q1@Pi-q1)@N),
                  isotropic_gauge=maximum(K.conj().T@q1@K))
    assert max(errors.values()) < 2e-14
    assert np.linalg.eigvalsh(h1)[0] > .5 and np.linalg.eigvalsh(h0)[0] > .5
    result = dict(errors=errors, auxiliary_min_eigenvalues=[
        float(np.linalg.eigvalsh(h1)[0]), float(np.linalg.eigvalsh(h0)[0])],
        nonorthogonal_frequency_projector_defect=maximum(c1-c1.conj().T),
        scope="One full 126-by-32 Cauchy principal-symbol calibration, not a PDE inverse or global zero-mode scan.")
    return result, (K, c1, c0, q1, q0, B, Pi)


def exact_repair_check(data):
    K, c1, c0, q1, q0, B, Pi = data
    rng = np.random.default_rng(1766)
    a = rng.normal(size=(126, 126))*.002
    skew = a-a.T
    generator = np.linalg.solve(q1, skew)
    eye = np.eye(126)
    U = np.linalg.solve(eye-generator/2, eye+generator/2)
    cd = U@c1@np.linalg.inv(U)
    ct, corr = repair(cd, c0, K, B, Pi, q1, q0)
    cm, _ = repair(eye-cd, np.eye(32)-c0, K, B, Pi, q1, q0)
    kd = adj(K, q1, q0)
    errors = dict(cayley_pseudo_unitary=maximum(U.conj().T@q1@U-q1),
                  original_charge_self_adjoint=maximum(q1@cd-cd.conj().T@q1),
                  repair_identity=maximum(ct+corr-cd),
                  repaired_charge_self_adjoint=maximum(q1@ct-ct.conj().T@q1),
                  sum_identity=maximum(ct+cm-eye),
                  exact_intertwining=maximum(ct@K-K@c0),
                  constraint_intertwining=maximum(kd@ct-c0@kd))
    initial_defect=maximum(cd@K-K@c0)
    assert initial_defect > .01 and max(errors.values()) < 2e-13
    return dict(initial_gauge_defect=initial_defect, correction_norm=maximum(corr),
                errors=errors, repaired_projector_defect=maximum(ct@ct-ct),
                scope="Finite smooth-error analogue. Repaired c need not be an idempotent; no positivity follows.")


def positivity_not_automatic():
    # Four-dimensional algebraic block: gauge null pair and physical +/- pair.
    q1 = np.array([[0, 1, 0, 0], [1, 0, 0, 0],
                   [0, 0, 1, 0], [0, 0, 0, -1.]], complex)
    q0 = np.ones((1, 1), complex)
    K = np.eye(4, dtype=complex)[:, :1]
    c1 = np.diag([1., 1., 1., 0.]).astype(complex)
    c0 = np.ones((1, 1), complex)
    B, Pi, _, _ = make_projection(K, c1, c0)
    v = np.array([1., -1., 0., 0.])/np.sqrt(2)
    w = np.array([0., 0., 0., 1.])
    theta = .37
    U = (np.eye(4)+(np.cos(theta)-1)*(np.outer(v, v)+np.outer(w, w))
         +np.sin(theta)*(np.outer(w, v)-np.outer(v, w)))
    assert maximum(U.conj().T@q1@U-q1)<1e-14
    cd = U@c1@U.T
    ct, correction = repair(cd, c0, K, B, Pi, q1, q0)
    value = float((w@q1@ct@w).real)
    expected = -float(np.sin(theta)**2)
    assert abs(value-expected)<1e-14 and value < -.1
    assert maximum(ct@K-K@c0)<1e-14
    assert maximum(K.conj().T@q1@w)==0
    return dict(physical_negative_variance=value, analytic_value=expected,
                exact_gauge_error=maximum(ct@K-K@c0),
                physical_test_is_constraint_compatible=True,
                scope="A finite algebraic failure witness for positivity inferred from gauge repair alone; not a counterexample to physical Hadamard existence.")


def run():
    p, data = projection_check()
    result = dict(round=766, tests_run=3, failures=0, errors=0,
                  frequency_adapted_projection=p,
                  exact_gauge_and_CCR_repair=exact_repair_check(data),
                  positivity_scope_witness=positivity_not_automatic())
    deps = ("research_note_753.md", "research_note_754.md", "research_note_764.md",
            "research_note_765.md", "joint_covariant_gauge_complex.py",
            "round766_drafts/physical_hadamard_entry.md")
    result["dependency_hashes"] = {n: hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps}
    result["scope"] = ("On the inherited irreducible on-shell background: elliptic Cauchy gauge "
                       "inverse, exact frequency-compatible projection, and a smooth repair of "
                       "Hadamard two-point candidates preserving exact gauge and CCR. "
                       "Physical positivity and all-sector renormalized sources remain unproved. "
                       "Numerics calibrate matrix identities only, not the full pseudodifferential construction.")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--write-results", action="store_true")
    args = p.parse_args(); result = run()
    if args.write_results:
        with TARGET.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text("utf8"))
    print(json.dumps({k: result[k] for k in ("round", "tests_run", "failures", "errors")}))
