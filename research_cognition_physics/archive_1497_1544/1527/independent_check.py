"""Exact Weierstrass product ranks and a two-source bound using a nuclear weak current."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def product(left, right):
    """Monomials x^i y^j on y^2=x^3-x, allowing y/x in a translated bundle."""
    i, j = (a + b for a, b in zip(left, right))
    if j == 2:
        return {(i + 3, 0): 1, (i + 1, 0): -1}
    return {(i, j): 1}


def exact_rank(rows):
    mat = [[Fraction(v) for v in row] for row in rows]
    rank = 0
    for col in range(len(mat[0])):
        pivot = next((k for k in range(rank, len(mat)) if mat[k][col]), None)
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        lead = mat[rank][col]
        mat[rank] = [v / lead for v in mat[rank]]
        for row in mat[rank + 1:]:
            multiple = row[col]
            for k, v in enumerate(mat[rank]):
                row[k] -= multiple * v
        rank += 1
    return rank


def product_ranks():
    one, x, y, w = (0, 0), (1, 0), (0, 1), (-1, 1)
    degree_two = [one, x]
    degree_three = [one, x, y]
    cases = {
        "symmetric_3": (itertools.combinations_with_replacement(degree_three, 2), 6),
        "ordered_2_same": (itertools.product(degree_two, repeat=2), 4),
        "ordered_2_distinct": (itertools.product(degree_two, [one, w]), 4),
        "ordered_3_distinct": (itertools.product(degree_three, [one, x, w]), 6),
    }
    results = {}
    for name, (pairs, output_dim) in cases.items():
        polys = [product(a, b) for a, b in pairs]
        monomials = sorted(set().union(*polys))
        rank = exact_rank([[p.get(m, 0) for m in monomials] for p in polys])
        results[name] = dict(rank=rank, kernel=len(polys)-rank, cokernel=output_dim-rank)
    assert [item["rank"] for item in results.values()] == [6, 3, 4, 6]
    return results


def allowed_angle(cost, mass_ratio):
    return math.asin(math.sqrt(max(0.0, (cost**2 - mass_ratio**2) / (1 - mass_ratio**2))))


def joint_floor(overlap, ru, rd):
    angle = math.acos(overlap)
    lower, upper = max(ru, rd), 1.0
    for _ in range(80):
        middle = (lower + upper) / 2
        if allowed_angle(middle, ru) + allowed_angle(middle, rd) < angle:
            lower = middle
        else:
            upper = middle
    return upper


def physical_overlaps(t, grid):
    """Integrate the three Landau-gauge profiles on a torus of physical area one."""
    x = np.arange(grid)[None, :] / grid
    y = np.arange(grid)[:, None] / grid

    def modes(n):
        states = []
        for j in range(n):
            state = np.zeros((grid, grid), dtype=complex)
            for ell in range(-6, 7):
                k = j + n * ell
                state += np.exp(-math.pi * n * t * (y + k/n)**2 + 2j * math.pi * k * x)
            states.append((2*n*t)**.25 * state.ravel())
        states = np.array(states).T
        np.testing.assert_allclose(states.conj().T @ states / grid**2, np.eye(n), atol=2e-13)
        return states

    external, emitted = modes(3), modes(6)
    pairs = itertools.combinations_with_replacement(range(3), 2)
    products = np.column_stack([external[:,i]*external[:,j]*(1 if i == j else math.sqrt(2))
                                for i,j in pairs])
    matrix = emitted.conj().T @ products / grid**2
    np.testing.assert_allclose(matrix.imag, 0, atol=2e-13)
    np.testing.assert_allclose(emitted @ matrix, products, atol=2e-12)
    return matrix.real


def square_vertex():
    matrix = physical_overlaps(1., 32)
    np.testing.assert_allclose(matrix, physical_overlaps(1., 48), atol=1e-13)
    inv_gram = np.linalg.inv(matrix.T @ matrix)
    a, b, c = inv_gram[0,0], inv_gram[4,4], inv_gram[0,4]
    assert c < 0
    template = np.zeros_like(inv_gram)
    for i,j in [(0,4),(3,2),(5,1)]:
        template[i,i], template[j,j], template[i,j], template[j,i] = a,b,c,c
    np.testing.assert_allclose(inv_gram, template, atol=2e-13)
    equal_weight = (a + 2*b - 2*math.sqrt(2)*abs(c))/3
    pure_floor = min(a, equal_weight)
    for w, target in [(np.array([1.,0.,0.]), a), (np.ones(3)/math.sqrt(3), equal_weight)]:
        v = [w[i]*w[j]*(1 if i == j else math.sqrt(2))
             for i,j in itertools.combinations_with_replacement(range(3), 2)]
        np.testing.assert_allclose(np.linalg.norm(np.linalg.solve(matrix.T, v))**2, target)
    up = np.asarray(json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))["one_TeV_DR"]["up"])
    g, dg = json.loads((HERE.parent/"1523/inputs.json").read_text(encoding="utf-8"))["sm_dr_at_one_TeV"]["g2"]
    lo, hi = up[:,0]-2*up[:,1], up[:,0]+2*up[:,1]
    low, high = np.linalg.eigvalsh(inv_gram)[[0,-1]]
    light_high = np.linalg.norm(hi[:2])
    triangle = lo[2]*math.sqrt(pure_floor) - math.sqrt(high)*light_high
    quadratic = (lo[2]**2*pure_floor - (high-low)*hi[2]*light_high + low*(lo[:2] @ lo[:2]))
    floor = math.sqrt(quadratic)
    assert floor > triangle > g + 2*dg
    return dict(physical_gram_inverse_block=[a,b,c], pure_floor=pure_floor,
                rank1_values=[a,equal_weight], old_triangle_crosscheck=triangle,
                full_spectrum_floor=floor, weak_ceiling=g+2*dg, ratio=floor/(g+2*dg))


def rectangle_witness():
    # Evaluate the author's rounded fitted shape once through the physical integral.
    shape = 1.10951463
    matrix = physical_overlaps(shape, 40)
    raw = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))["one_TeV_DR"]
    down = json.loads((HERE/raw["down_from"]).read_text(encoding="utf-8"))["down_yukawa_dr_1TeV"]
    up = np.array(raw["up"])[:,0]
    down = np.asarray(down)[:,0]
    a,b,c,delta = raw["angles_order_12_13_23_delta"]
    def rotation(i,j,theta,phase):
        r = np.eye(3, dtype=complex)
        r[i,i] = r[j,j] = math.cos(theta)
        r[i,j] = math.sin(theta)*np.exp(-1j*phase)
        r[j,i] = -r[i,j].conjugate()
        return r
    v = rotation(1,2,c,0) @ rotation(0,2,b,delta) @ rotation(0,1,a,0)
    def norm_of_profile(y):
        vec = [y[i,j]*(1 if i == j else math.sqrt(2)) for i,j in itertools.combinations_with_replacement(range(3),2)]
        return np.linalg.norm(np.linalg.solve(matrix.T, vec))
    ru, rd = [norm_of_profile(y) for y in (np.diag(up), v@np.diag(down)@v.T)]
    target = json.loads((HERE.parent/"1523/inputs.json").read_text(encoding="utf-8"))["sm_dr_at_one_TeV"]["g2"][0]
    gain = math.hypot(ru,rd)
    assert abs(gain-target) < 1e-8
    return dict(shape=shape, gain=gain, tan_beta=ru/rd)


def data_readout():
    def archived(number):
        return json.loads((HERE.parent / number / "inputs.json").read_text(encoding="utf-8-sig"))

    overlap, uncertainty = archived("1381")["ud_beta"]
    mass = archived("1386")["yukawa"]
    results = {}
    for widths, transfer in [(0, 0), (2, 0), (2, .005)]:
        ceiling = min(1., overlap + widths * uncertainty + transfer)
        ratios = [(mass[sector][0][0] - widths * mass[sector][0][1]) /
                  (mass[sector][1][0] + widths * mass[sector][1][1])
                  for sector in ("up", "down")]
        floor = joint_floor(ceiling, *ratios)
        results[f"widths_{widths}_transfer_{transfer}"] = dict(
            overlap_upper=ceiling, ratios=ratios, angle_only=math.sqrt((1.-ceiling)/2),
            mass_and_angle=floor)

    # An independent 2-plane construction attains the bound with these overlap
    # and mass-ratio inputs. The third singular values are arbitrary larger ones.
    ru, rd = results["widths_0_transfer_0"]["ratios"]
    floor = joint_floor(overlap, ru, rd)
    a = allowed_angle(floor, ru)
    theta = math.acos(overlap)
    n = np.array([math.cos(a), math.sin(a), 0])
    plane = np.eye(3) - np.outer(n, n)
    c, s = math.cos(theta), math.sin(theta)
    rotation = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    for matrix in [np.diag([ru, 1, 9]), rotation @ np.diag([rd, 1, 7])]:
        assert np.isclose(np.linalg.norm(matrix-plane@matrix, 2), floor, atol=1e-12)
    return results


def certificate():
    points = [[1., -1., 0.], [1., 0., 0.], [1., 1., 0.], [1., 2., math.sqrt(6.)]]
    assert all(abs(y*y-(x**3-x)) < 1e-12 for _, x, y in points)
    ranks = [int(np.linalg.matrix_rank([points[i] for i in ids])) for ids in [(0, 1, 2),(0, 1, 3)]]
    assert ranks == [2, 3]
    return {"weierstrass": product_ranks(), "three_evaluations_collinear_and_general_rank": ranks,
            "square_torus_unit_vertex": square_vertex(),
            "rectangular_point_direct_integral": rectangle_witness(),
            "two_common_source_cost_over_second_yukawa": data_readout()}


def compare(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(compare(v, b[k]) for k, v in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(compare(x, y) for x, y in zip(a, b))
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-14)


if __name__ == "__main__":
    flags = argparse.ArgumentParser()
    group = flags.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = flags.parse_args()
    result = certificate()
    out = HERE / "independent_results.json"
    if args.write:
        with out.open("x", encoding="utf-8") as f:
            f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert compare(result, json.loads(out.read_text(encoding="utf-8")))
    print("1527 independent calculations passed")
