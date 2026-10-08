"""Independent review of old 929/956 formulae; no new scientific group.

No author modules are imported. Default replay is read-only. --write creates
the review result once. Universal statements are proved in the review text.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULT = HERE / "review_science_results.json"


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def scaled(a, denominator):
    return [[F(x, denominator) for x in row] for row in a]


def trace(a):
    return sum(a[j][j] for j in range(len(a)))


def review_929():
    # M_x=D_x/sqrt(5); U_pi/2=R/sqrt(2). Both two-step
    # products have rational entries, so this computation has no rounding.
    ds = [[[2, 0], [0, 1]], [[1, 0], [0, 2]]]
    r = [[1, -1], [1, 1]]
    out = {}
    traces = []
    distributions = []
    for name, rotated in [("theta_0", False), ("theta_pi_over_2", True)]:
        effects = []
        for x in range(2):
            for y in range(2):
                a = (scaled(mul(mul(mul(ds[y], r), ds[x]), r), 10)
                     if rotated else scaled(mul(ds[y], ds[x]), 5))
                effects.append(mul(transpose(a), a))
        total = [[sum(e[i][j] for e in effects) for j in range(2)]
                 for i in range(2)]
        assert total == [[1, 0], [0, 1]]
        ts = sorted(trace(e) for e in effects)
        ps = [e[0][0] for e in effects]
        traces.append(ts)
        distributions.append(ps)
        out[name] = {
            "effect_order": ["00", "01", "10", "11"],
            "effects_exact": [[[str(x) for x in row] for row in e]
                              for e in effects],
            "trace_multiset_exact": list(map(str, ts)),
            "input_zero_probabilities_exact": list(map(str, ps)),
            "effect_sum_is_identity": True,
        }
    assert traces == [[F(8, 25)]*2 + [F(17, 25)]*2, [F(1, 2)]*4]
    assert distributions == [[F(64,100),F(16,100),F(16,100),F(4,100)],
                             [F(13,100),F(37,100),F(13,100),F(37,100)]]
    gap = distributions[0][0] - distributions[1][0]
    margin = gap - 2*F(3, 1000)
    assert gap == F(51, 100) and margin == F(63, 125)
    out.update(event_00_gap_exact=str(gap),
               two_window_error_margin_exact=str(margin),
               inequivalent_under_common_input_unitary_and_event_permutation=True)
    return out


def basis_and_commutant(mats):
    n = mats[0].shape[0]
    array = np.array([m.reshape(-1) for m in mats])
    _, s, vh = np.linalg.svd(array, full_matrices=False)
    basis = vh[s > 1e-10].reshape(-1, n, n)
    constraints = np.vstack([np.kron(g.T, np.eye(n)) -
                              np.kron(np.eye(n), g) for g in basis])
    _, cs, cvh = np.linalg.svd(constraints, full_matrices=True)
    rank = int(np.sum(cs > 1e-10))
    commutant = cvh[rank:].reshape(-1, n, n)
    return len(basis), commutant


def channel_dimensions(kraus):
    products = [a.conj().T @ b for a in kraus for b in kraus]
    product_rank, algebra = basis_and_commutant(products)
    _, algebra_commutant = basis_and_commutant(list(algebra))
    return {"kraus_product_span_dimension": product_rank,
            "maximal_correctable_algebra_dimension": len(algebra),
            "its_commutant_dimension": len(algebra_commutant)}


def entropy(rho):
    eig = np.linalg.eigvalsh(rho)
    eig = eig[eig > 1e-14]
    return float(-np.sum(eig*np.log(eig)))


def binary_entropy(p):
    return float(-p*np.log(p) - (1-p)*np.log(1-p))


def review_956():
    u, v = 1., .1
    omega = np.sqrt(u*u + 16*v*v)
    j, d = (omega-u)/2, (1-u/omega)/2
    h = np.diag([0., 0., 0., 0., u, u])
    h[0,4] = h[4,0] = -2*v
    charge = np.diag([0., 0., 0., 0., 1., 1.])
    p = np.diag([1., 0., 0., 0.])
    w = np.zeros((6, 4))
    w[0,0], w[4,0] = np.sqrt(1-d), np.sqrt(d)
    w[1,1] = w[2,2] = w[3,3] = 1.
    ks = [(np.eye(6)-charge)@w, charge@w]
    encoding = np.stack(ks, axis=1).reshape(12,4)
    pointer_ks = list(encoding.reshape(6,2,4))
    identities = {
        "isometry": float(np.linalg.norm(encoding.T@encoding-np.eye(4))),
        "energy_intertwining": float(np.linalg.norm(h@w-w@(-j*p))),
        "charge_compression": float(np.linalg.norm(w.T@charge@w-d*p)),
        "material_product_00": float(np.linalg.norm(ks[0].T@ks[0]-(np.eye(4)-d*p))),
        "material_product_11": float(np.linalg.norm(ks[1].T@ks[1]-d*p)),
        "material_cross_product": float(np.linalg.norm(ks[0].T@ks[1])),
    }
    assert max(identities.values()) < 1e-13
    md, pd = channel_dimensions(ks), channel_dimensions(pointer_ks)
    assert list(md.values()) == [2,10,2]
    assert list(pd.values()) == [16,1,16]
    actual = sum(k@p@k.T for k in ks)
    compressed = [np.diag([np.sqrt(1-d),1.,1.,1.]), np.sqrt(d)*p]
    wrong = w@sum(k@p@k.T for k in compressed)@w.T
    before = w@p@w.T
    leakage = float(np.trace((np.eye(6)-w@w.T)@actual))
    actual_energy = float(np.trace(actual@h))
    wrong_energy = float(np.trace(wrong@h))
    energy_delta = float(np.trace((actual-wrong)@h))
    assert abs(leakage-2*d*(1-d)) < 1e-13
    assert abs(actual_energy-u*d) < 1e-13
    assert abs(wrong_energy+j) < 1e-13
    assert abs(energy_delta-(u*d+j)) < 1e-13
    pointer_s = sum(k@p@k.T for k in pointer_ks)
    triplet = np.diag([0.,1.,0.,0.])
    pointer_t = sum(k@triplet@k.T for k in pointer_ks)
    pointer_distance = float(np.sum(np.linalg.svd(pointer_s-pointer_t,
                                                compute_uv=False))/2)
    assert abs(pointer_distance-d) < 1e-13
    psi = np.array([1.,1.,0.,0.])/np.sqrt(2)
    rho = np.outer(psi,psi)
    coherent_material = sum(k@rho@k.T for k in ks)
    coherent_pointer = sum(k@rho@k.T for k in pointer_ks)
    ent = entropy(coherent_material)
    algebra_plus_center = float(np.log(2) + binary_entropy(d)/2)
    assert abs(ent-binary_entropy(d/2)) < 1e-13
    assert abs(ent-entropy(coherent_pointer)) < 1e-13
    assert algebra_plus_center > np.log(2) > ent
    return {
        "physical_basis": ["s","t1","t2","t3","d_plus","d_minus"],
        "logical_basis": ["s","t1","t2","t3"],
        "U": u, "v": v, "d": float(d), "J": float(j),
        "identity_residuals": identities,
        "material": md, "pointer": pd,
        "pointer_state_s": pointer_s.tolist(),
        "pointer_state_triplet": pointer_t.tolist(),
        "pointer_half_trace_distance": pointer_distance,
        "equal_prior_minimum_average_error": (1-pointer_distance)/2,
        "physical_poststate_leakage": leakage,
        "material_energy_before": float(np.trace(before@h)),
        "material_energy_actual_after": actual_energy,
        "material_energy_compressed_substitute_after": wrong_energy,
        "material_internal_energy_difference": energy_delta,
        "coherent_input_material_entropy": ent,
        "coherent_input_pointer_entropy": entropy(coherent_pointer),
        "inapplicable_algebra_entropy_plus_center": algebra_plus_center,
        "complementary_recovery_for_this_dilation": False,
        "instrument_manufacture_or_unread_pointer_permission_certified": False,
    }


def build():
    sources = ["archive_923_934/research_note_929.md",
               "archive_956_989/research_note_956.md",
               "archive_1009_/1042/common_parent_coverage_audit.md"]
    return {
        "status": "pass",
        "purpose": "Independent old-formula review for 1043; zero new science groups",
        "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in sources},
        "review_929_exact": review_929(),
        "review_956_chosen_physical_dilation": review_956(),
        "limitations": [
            "929 quotient is common input unitary and reversible event relabeling only",
            "956 numerical checks are diagnostic at U=1,v=.1; analytic proofs give scope",
            "956 preserves physical six-dimensional poststates",
            "956 minimal pure pointer is a selected extension, not supplied access",
            "No full E_nat or complete A1_D-A7 model is asserted",
        ],
    }


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for key in a: compare(a[key], b[key], path+"."+key)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x,y) in enumerate(zip(a,b)): compare(x,y,path+f"[{i}]")
    elif isinstance(a, float):
        assert np.isclose(a,b,rtol=1e-11,atol=1e-13), (path,a,b)
    else:
        assert a == b, (path,a,b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    computed = build()
    if args.write:
        with RESULT.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(computed,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
        print("PASS: independently checked; review result created")
    else:
        compare(computed,json.loads(RESULT.read_text(encoding="utf-8")))
        print("PASS: independent 929 exact / 956 physical-poststate read-only replay")
