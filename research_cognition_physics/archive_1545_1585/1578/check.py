"""Direct Weyl rank at the common vacuum and its equivariant block maps."""
import argparse
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def data(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def parent_blocks(cfg):
    m = cfg["block_dimension"]
    a = np.repeat(cfg["adjoint_eigenvalues"], m)
    s = np.zeros((len(a), len(a)))
    for left, right in cfg["symmetric_edges"]:
        ix = np.arange(m)
        s[left*m+ix, right*m+ix] = 1
        s[right*m+ix, left*m+ix] = 1
    return a, s


def direct_tensor_rank(a, y):
    antisym = list(combinations(range(len(a)), 2))
    sym = list(combinations_with_replacement(range(len(a)), 2))
    col = {pair: j for j, pair in enumerate(sym)}
    matrix = np.zeros((2*len(antisym), len(sym)), dtype=complex)
    for field, strength in enumerate(y):
        for i, pair in enumerate(antisym):
            matrix[field*len(antisym)+i, col[pair]] = strength*(a[pair[0]]-a[pair[1]])
    return 2*int(np.linalg.matrix_rank(matrix))


def mapped_blocks(a, s):
    ident = np.eye(len(a))
    center = ident-a@a
    outer = [(a@a+sign*a)/2 for sign in (-1,1)]
    return center, [p@s@center.T for p in outer]


def covariance(a, s):
    # Use a reproducible full parent basis change that mixes the three subspaces.
    rng = np.random.default_rng(1578)
    m = len(a)
    q, _ = np.linalg.qr(rng.normal(size=(m,m)) + 1j*rng.normal(size=(m,m)))
    q[:,0] /= np.linalg.det(q)
    A = np.diag(a)
    p, edges = mapped_blocks(A, s)
    P, moved = mapped_blocks(q@A@q.conj().T, q@s@q.T)
    assert np.allclose(P, q@p@q.conj().T)
    for edge, transformed in zip(edges, moved):
        assert np.allclose(transformed, q@edge@q.T)
        assert np.allclose(edge.conj().T@edge, p.T)
    assert all(np.linalg.matrix_rank(e) == 5 for e in edges)
    return len(moved)


def bilinear_channels():
    charges = {"X":2, "Y":-2, "L":-1}
    insertions = {"bare":0, "A":0, "S":2, "S*":-2}
    return [left+right+":"+key for (left,right),(key,s) in product(
        combinations_with_replacement(charges,2), insertions.items())
        if (charges[left]+charges[right]+s) % 15 == 0]


def run():
    cfg = data("inputs.json")
    a, s = parent_blocks(cfg)
    channels = bilinear_channels()
    assert set(channels) == {"XY:bare","XY:A","LL:S"}
    tensor_rank = direct_tensor_rank(a, [1+2j, 3-1j])
    assert tensor_rank == 150 and direct_tensor_rank(a, [0,0]) == 0
    allowed = []
    for k in range(4):
        h = np.diag([1j]*(k)+[0]*(3-k))
        r = int(np.linalg.matrix_rank(np.kron(h,s)))
        assert r == 10*k
        allowed.append(tensor_rank+r)
    # Independently count the maximal rank in the restricted representation.
    # Multiplicities of positive, dual irreps and dimension.
    complex_pairs = [(6,3,10),(3,3,15),(3,6,5)]
    real = [(6,24),(6,1)]
    full = sum((p+r)*m for p,r,m in complex_pairs) + sum(k*m for k,m in real)
    eft_rank = sum(2*min(p,r)*m for p,r,m in complex_pairs) + sum(k*m for k,m in real)
    assert (full, full-max(allowed), full-eft_rank) == (375,195,45)
    idx = 2*3 + 2*7 + 2*5
    delta_b = 2*idx/3
    base = data(cfg["comparison"])["SM_b"]
    return {"center_balanced_channels_before_index_symmetry":channels,
            "rank_tensor_offdiagonal":tensor_rank,"total_rank_at_nonzero_y_by_rank_h":allowed,
            "Weyl_count":full,"renormalizable_minimum_zero_modes":full-max(allowed),
            "EFT_maximum_mass_rank":eft_rank,"EFT_minimum_zero_modes":full-eft_rank,
            "covariant_full_rank_bridges":covariance(a,s),
            "unlifted_additional_Dynkin_index":idx,"all_open_delta_b":delta_b,
            "b_GUT_Y_weak_color_if_all_extra_zero_modes_active":[x+delta_b for x in base]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check", action="store_true")
    flags.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE / "results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result == json.loads(target.read_text())
    print(json.dumps(result))
