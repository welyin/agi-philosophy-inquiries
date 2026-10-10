"""Independent 1060 review: Pati--Salam matrices, no author-code imports.

Finite checks diagnose the representation implementation.  The all-background
claim is the analytic group-structure pullback, not a numerical bordism result.
Default execution is read-only; --record creates the review result once.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import itertools
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
AUTHOR = [HERE.parent / "research_note_1060.md"] + [HERE / name for name in
    ("proof.md", "selection.md", "sources.md", "dependency_update.md",
     "global_anomaly_pullback.py", "results.json")]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block(*parts):
    out = np.zeros((sum(len(x) for x in parts),) * 2, complex)
    start = 0
    for part in parts:
        end = start + len(part)
        out[start:end, start:end] = part
        start = end
    return out


def wedge2(matrix):
    pairs = list(itertools.combinations(range(len(matrix)), 2))
    return np.array([[matrix[i, k] * matrix[j, l] - matrix[i, l] * matrix[j, k]
                      for k, l in pairs] for i, j in pairs])


def ps16(A, L, R):
    return block(np.kron(A, L), np.kron(A.conj(), R))


def ps10(A, L, R):
    return block(wedge2(A), np.kron(L, R))


def ps_map(k, C, W, theta):
    a = 1 + 6 * k
    A = block(np.exp(1j * a * theta) * C,
              np.array([[np.exp(-3j * a * theta)]]))
    R = np.diag(np.exp(1j * np.array([-3., 3.]) * theta))
    return A, W, R


def special_unitary(rng, n):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    q[:, 0] /= np.linalg.det(q)
    return q


def run():
    # Independent affine integer weights from (4,2,1)+(bar4,1,2).
    b = [1, 1, 1, -3]
    first = [(x, 6 * x) for x in b for _ in range(2)]
    second = [(-x + t, -6 * x) for x in b for t in (-3, 3)]
    expected = {
        "Q": [[1, 6]] * 6, "L": [[-3, -18]] * 2,
        "uc": [[-4, -6]] * 3, "dc": [[2, -6]] * 3,
        "N": [[0, 18]], "ec": [[6, 18]],
    }
    actual = {
        "Q": [list(x) for x in first[:6]], "L": [list(x) for x in first[6:]],
        "uc": [list(second[i]) for i in (0, 2, 4)],
        "dc": [list(second[i]) for i in (1, 3, 5)],
        "N": [list(second[6])], "ec": [list(second[7])],
    }
    assert actual == expected
    # Kernel descent in rational phases: all k-dependent terms are integers.
    # At theta=j*pi/3 the three PS factors equal the same (-1)^j.
    rational_cases = 0
    for j in range(6):
        for k in (-7, 0, 5):
            a = 1 + 6 * k
            target = Fraction(j, 2)
            phases = [Fraction((a + 2) * j, 6), Fraction(-3 * a * j, 6),
                      Fraction(j, 2), Fraction(-3 * j, 6), Fraction(3 * j, 6)]
            assert all((phase - target).denominator == 1 for phase in phases)
            rational_cases += 1
    # The all-k part is the two affine identities below, not the sample list.
    assert (1 + 2, 6) == (3, 6)
    assert (-3 * 1 - 3, -3 * 6) == (-6, -18)

    I4, I2 = np.eye(4), np.eye(2)
    c16 = ps16(1j * I4, I2, -I2)
    c10 = ps10(1j * I4, I2, -I2)
    assert np.array_equal(c16, 1j * np.eye(16))
    assert np.array_equal(c10, -np.eye(10))
    assert np.array_equal(c16 @ c16, -np.eye(16))
    assert np.array_equal(-(c16 @ c16), np.eye(16))
    assert np.array_equal(ps16(-I4, -I2, -I2), np.eye(16))
    assert np.array_equal(ps10(-I4, -I2, -I2), np.eye(10))

    rng = np.random.default_rng(1060007)
    errors = {key: 0.0 for key in ("group_law_16", "group_law_10", "periodicity",
                                 "center_descent", "higgs_restriction", "unitarity")}
    count = 0
    for k in (-2, 0, 3):
        for _ in range(4):
            C, D = special_unitary(rng, 3), special_unitary(rng, 3)
            W, V = special_unitary(rng, 2), special_unitary(rng, 2)
            theta, phi = rng.uniform(-1, 1, 2)
            p = ps_map(k, C, W, theta)
            q = ps_map(k, D, V, phi)
            pq = ps_map(k, C @ D, W @ V, theta + phi)
            for rep, label in ((ps16, "group_law_16"), (ps10, "group_law_10")):
                err = np.linalg.norm(rep(*p) @ rep(*q) - rep(*pq), 2)
                errors[label] = max(errors[label], float(err))
            R = ps16(*p)
            errors["unitarity"] = max(errors["unitarity"], float(np.linalg.norm(R.conj().T @ R - np.eye(16), 2)))
            errors["periodicity"] = max(errors["periodicity"], float(np.linalg.norm(R - ps16(*ps_map(k, C, W, theta + 2 * np.pi)), 2)))
            H = np.kron(W, p[2])
            # In L tensor R ordering, R weight +3 uses column indices 1,3.
            errors["higgs_restriction"] = max(errors["higgs_restriction"], float(np.linalg.norm(H[np.ix_([1, 3], [1, 3])] - np.exp(3j * theta) * W, 2)))
            count += 1
        for j in range(6):
            theta = j * np.pi / 3
            R = ps16(*ps_map(k, np.exp(2j * theta) * np.eye(3), np.exp(-3j * theta) * I2, theta))
            errors["center_descent"] = max(errors["center_descent"], float(np.linalg.norm(R - np.eye(16), 2)))
    assert max(errors.values()) < 1e-11, errors
    return {
        "review_of_round": 1060, "passed": True,
        "method": "independent Pati-Salam direct matrices and rational center phases",
        "author_code_imported": False,
        "author_sha256": {str(p.relative_to(HERE.parent)).replace('\\', '/'): sha(p) for p in AUTHOR},
        "charge_affine_coefficients": actual,
        "rational_center_diagnostics": rational_cases,
        "random_group_diagnostics": count,
        "operator_error_maxima": errors,
        "central_action_on_16": "i I16", "central_action_on_10": "-I10",
        "analytic_kernel": "exactly Z6; reverse inclusion proved in independent_review.md",
        "all_background_claim": "analytic natural pullback of adopted WWW section 5.1.3 theorem",
        "numerical_bordism_classification": False,
        "new_scientific_calibration_groups": 0,
        "whole_cognition_nonuniqueness_claimed": False,
    }


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=2e-12)
    return a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE / "independent_checks.json"
    if args.record:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert same(result, json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps(result, ensure_ascii=False))
