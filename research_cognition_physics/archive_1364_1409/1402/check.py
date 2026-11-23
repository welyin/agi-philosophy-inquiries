"""Round 1402: shared two-sided filters, full matrices and singular-value data."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def abs_det_ratio(v):
    return math.prod(v[:3]) / math.prod(v[3:])


def volume_interval(pairs, width):
    v, dx = np.asarray(pairs).T
    lo, hi = v - width * dx, v + width * dx
    assert np.all(lo > 0)
    lower = math.prod(lo[:3]) / math.prod(hi[3:])
    upper = math.prod(hi[:3]) / math.prod(lo[3:])
    return {"R": [lower, upper], "gain_ratio": [lower**(1/3), upper**(1/3)],
            "tau_if_equal_gains": [math.prod(lo[:3]) / math.prod(hi[3:5]),
                                   math.prod(hi[:3]) / math.prod(lo[3:5])],
            "matched_tau": [lo[5], hi[5]]}


def pos_sqrt(h):
    w, v = np.linalg.eigh(h)
    assert np.min(w) > 0
    return (v * np.sqrt(w)) @ v.conj().T


def recover(d, e, g):
    ident = np.eye(len(g))
    q = np.linalg.solve(e.T, d.T).T
    assert np.allclose(q.conj().T @ g @ q, g, rtol=1e-9, atol=1e-9)
    f = np.linalg.inv(pos_sqrt(g))
    b = pos_sqrt(e.conj().T @ g @ e)
    for item in (d, e):
        core = np.linalg.solve(f, item) @ np.linalg.inv(b)
        assert np.allclose(core.conj().T @ core, ident, atol=1e-9)
        assert np.allclose(f @ core @ b, item, atol=1e-11)


def spectral_witness(down, charged):
    kd = (math.prod(down) / math.prod(charged))**(1/3)
    d = np.asarray(down) / kd
    e = np.asarray(charged)
    # Choose b[0] arbitrarily; the final cyclic equation closes by equal volumes.
    b = np.ones(len(d))
    for i in range(len(d) - 1):
        b[i+1] = b[i] * e[i] / d[i]
    f = d / b
    assert math.isclose(f[-1] * b[0], e[-1], rel_tol=5e-14)
    balance = math.sqrt(max(f) / max(b))
    b = balance * b
    f = f / balance
    shift = np.roll(np.eye(len(d)), 1, axis=1)
    down_map = np.diag(f) @ np.diag(b)
    lepton_map = np.diag(f) @ shift @ np.diag(b)
    assert np.allclose(np.linalg.svd(down_map * kd, compute_uv=False), np.sort(down)[::-1], atol=1e-15)
    assert np.allclose(np.linalg.svd(lepton_map, compute_uv=False), sorted(charged, reverse=True), atol=1e-15)
    assert max(f) <= 1 and max(b) <= 1
    recover(down_map, lepton_map, np.diag(1 / f**2))
    return {"k_down_at_k_lepton_one": kd, "F_diag": f.tolist(), "B_diag": b.tolist(),
            "larger_filter_norm": max(max(f), max(b))}


def certificates():
    config = json.loads((HERE/"inputs.json").read_bytes())
    old = json.loads((HERE.parent/"1386"/"inputs.json").read_bytes())["yukawa"]
    high_down = json.loads((HERE.parent/"1397"/"inputs.json").read_bytes())["high_scale"]["down"]
    assert high_down == config["scales"]["1e16"][:3]
    scales = {"MZ": old["down"] + old["charged_leptons"], **config["scales"]}
    records = {}
    for label, pairs in scales.items():
        values = np.array(pairs)[:, 0]
        records[label] = {
            "volume_ratio": abs_det_ratio(values),
            "gain_ratio": abs_det_ratio(values)**(1/3),
            "one_width": volume_interval(pairs, 1), "two_width": volume_interval(pairs, 2),
            "center_spectral_construction": spectral_witness(values[:3], values[3:])
        }
    assert records["MZ"]["two_width"]["R"][0] > 1 > records["1e16"]["two_width"]["R"][1]
    assert records["1e5"]["two_width"]["R"][0] > 1 > records["1e9"]["two_width"]["R"][1]
    assert records["1e7"]["one_width"]["R"][0] < 1 < records["1e7"]["one_width"]["R"][1]
    # A nonnormal but unitarizable example, and a Jordan matrix with the same
    # unit-circle spectrum that is not diagonalizable and cannot preserve a metric.
    sim = np.array([[1, 1j], [.3, 2]], dtype=complex)
    diag = np.diag([1j, -1])
    q = sim @ diag @ np.linalg.inv(sim)
    inv = np.linalg.inv(sim)
    metric = inv.conj().T @ inv
    e = np.array([[.5, .2j], [.1+.2j, .7]], dtype=complex)
    recover(q @ e, e, metric)
    jordan = np.array([[1., 1.], [0., 1.]])
    assert np.linalg.matrix_rank(jordan - np.eye(2)) == 1
    same_det = np.allclose(np.linalg.det(.5*np.eye(3)), np.linalg.det(np.diag([1.,1.,.125])))
    assert same_det
    return {"scales": records, "nonunitarizable_jordan_norm_at_power_32":
            np.linalg.norm(np.linalg.matrix_power(jordan, 32), 2),
            "incompatible_contraction_example_equal_determinant": .125}


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for k in actual:
            compare(actual[k], saved[k])
    else:
        assert np.allclose(actual, saved, atol=1e-14, rtol=5e-10)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = certificates()
    out = HERE/"results.json"
    if args.write:
        with out.open("x", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        compare(result, json.loads(out.read_bytes()))
    print(json.dumps({k: [x["volume_ratio"], x["gain_ratio"], x["two_width"]["R"]]
                      for k, x in result["scales"].items()}))
