"""Finite odd scalar-mixing response of the fitted lower-row Nelson--Barr example."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
TREE = runpy.run_path(str(ARCHIVE/"1810/check.py"))
KERNEL = runpy.run_path(str(ARCHIVE/"1797/check.py"))


def read(path):
    return json.loads((ARCHIVE/path).read_bytes())


def matter(heavy, ordinary_vertex, bare_vertex):
    data = read("1386/inputs.json")
    weak = TREE["empirical_matrix"](data)
    masses = np.array(data["yukawa"]["down"])[:, 0]*read("1810/inputs.json")["vev_GeV"]/math.sqrt(2)
    gram = (weak*masses**2)@weak.conj().T
    val, vec = np.linalg.eigh(1j*gram.imag)
    z = vec[:, 0]*math.sqrt(-2*val[0])
    m = np.linalg.cholesky((gram+np.outer(z, z.conj())).real)
    xi = np.linalg.solve(m, z)
    row = heavy*xi.conj()
    axis = int(np.argmax(abs(row)))
    # Tree matter is invariant under a common vectorlike rephasing of T.
    # Fix the physical orientation relative to the added *real* couplings once.
    row *= np.exp(1j*(math.pi/4-np.angle(row[axis])))
    bare = heavy*math.sqrt(1-float(np.vdot(xi, xi).real))
    matrix = np.zeros((4, 4), complex)
    matrix[:3, :3] = m
    matrix[3] = np.r_[row, bare]
    left, values, right_h = np.linalg.svd(matrix)
    order = np.argsort(values)
    values, left, right = values[order], left[:, order], right_h.conj().T[:, order]
    np.testing.assert_allclose(left.conj().T@matrix@right, np.diag(values), atol=2e-11)
    assert abs(np.angle(np.linalg.det(matrix))) < 2e-12
    v = left[:3, :3]
    j = float((v[0, 1]*v[1, 2]*v[0, 2].conj()*v[1, 1].conj()).imag)
    assert abs(j-read("1810/results.json")["empirical_lower_block_competitor"]["J_at_L_3000_GeV"]) < 2e-13
    y, k = np.zeros((4, 4)), np.zeros((4, 4))
    y[3, axis], k[3, 3] = ordinary_vertex, bare_vertex
    vertices = (left.conj().T@y@right, left.conj().T@k@right)
    return values, vertices, row[axis], bare, j, axis


def full_kernel(masses, light, heavy, order, extent):
    # Three overlapping calls cover the six off-diagonal and four diagonal entries.
    result = np.full((4, 4), np.nan)
    splitting = heavy**2-light**2
    for indices in ([0, 1, 2], [0, 1, 3], [2, 3]):
        block = KERNEL["full_kernel"](masses[indices], light, splitting, order, extent)
        result[np.ix_(indices, indices)] = block[:len(indices), :len(indices)]
    assert np.isfinite(result).all()
    return result


def partial_widths(matrix, masses, scalar_mass):
    x = (masses/scalar_mass)**2
    widths = np.empty_like(matrix.real)
    for i, xi in enumerate(x):
        for j, xj in enumerate(x):
            if math.sqrt(xi)+math.sqrt(xj) >= 1:
                widths[i, j] = 0
                continue
            phase = math.sqrt((1-xi-xj)**2-4*xi*xj)
            mix = (abs(matrix[i, j])**2+abs(matrix[j, i])**2)*(1-xi-xj)
            mix -= 4*math.sqrt(xi*xj)*(matrix[i, j]*matrix[j, i]).real
            widths[i, j] = 3/(16*math.pi)*phase*mix
    return widths


def response():
    pole, y, k = 3000., .02, .01
    phi = math.pi/4
    s1, s2 = 3*pole, 5*pole
    m, (gy, gk), selected, bare, j, axis = matter(pole, y, k)
    a = full_kernel(m, s1, s2, 14, 96)
    b = full_kernel(m, s1, s2, 10, 84)
    cross = gy*gk.T+gk*gy.T
    coefficient = float(np.sum(a*cross).imag)
    difference = -math.sin(2*phi)*coefficient
    numeric_shift = abs(float(np.sum((a-b)*cross).imag))
    tail = KERNEL["kernel_tail"](m, s1, s2**2-s1**2, 96)
    tail_upper = float(np.sum(tail*abs(cross.imag)))
    assert numeric_shift < 1e-9*abs(coefficient)
    assert tail_upper < 1e-15*abs(coefficient)
    leading = math.log((s2/s1)**2)/(32*math.pi**2)*(m[:, None]/m[None, :]+m[None, :]/m[:, None])
    direct = float(np.sum(leading*cross).imag)
    weak_basis = k*y*selected.conjugate().imag*math.log((s2/s1)**2)/(16*math.pi**2*bare)
    assert math.isclose(direct, weak_basis, rel_tol=2e-10)
    distant = full_kernel(m, 100*s1, 100*s2, 12, 104)
    distant_difference = -math.sin(2*phi)*float(np.sum(distant*cross).imag)
    assert abs(distant_difference/(-weak_basis)-1) < .001
    positive = {}
    for tag, angle in [("plus", phi), ("minus", -phi)]:
        weights = []
        for mass, gg in [(s1, math.cos(angle)*gy+math.sin(angle)*gk),
                         (s2, -math.sin(angle)*gy+math.cos(angle)*gk)]:
            fractions = partial_widths(gg, m, mass)
            weights.append({"scalar_GeV": mass, "all_quark_width_over_scalar_mass": float(fractions.sum()),
                            "heavy_pair_width_over_scalar_mass": float(fractions[3, 3]),
                            "tagged_pair_fraction_in_all_quark_decays": float(fractions[3, 3]/fractions.sum())})
        positive[tag] = weights
    budget = read("1797/inputs.json")["theta_budget"]
    pair_plus = positive["plus"][0]["heavy_pair_width_over_scalar_mass"]
    pair_minus = positive["minus"][0]["heavy_pair_width_over_scalar_mass"]
    return {"scalar_poles_GeV": [s1, s2], "heavy_row_scale_GeV": pole,
            "common_quark_tree_poles_GeV": m.tolist(), "J_at_finite_tree_point": j,
            "ordinary_right_axis_zero_based_in_the_real_Cholesky_basis": axis,
            "ordinary_portal": y, "bare_heavy_portal": k, "bare_heavy_mass_GeV": bare,
            "selected_soft_component_GeV": [float(selected.real), float(selected.imag)],
            "strong_angle_plus_minus_difference": difference,
            "hierarchical_approximation_to_same_difference": -weak_basis*math.sin(2*phi),
            "exact_difference_when_both_scalar_masses_multiplied_by_100": distant_difference,
            "threshold_quadrature_change_upper_rounded_up": 1e-12 if numeric_shift < 1e-12 else numeric_shift,
            "absolute_integral_tail_upper_rounded_up": 1e-18 if tail_upper < 1e-18 else tail_upper,
            "at_least_one_absolute_theta_for_any_common_addition": abs(difference)/2,
            "theta_only_EDM_budget": budget,
            "largest_abs_sin_two_phi_for_both_signs_with_one_common_theta": 2*budget/abs(coefficient),
            "both_portals_scale_to_raise_halfdifference_to_budget": math.sqrt(2*budget/abs(difference)),
            "physical_two_body_reads": positive,
            "light_scalar_tagged_pair_width_plus_over_minus": pair_plus/pair_minus}


def round_data(obj):
    if isinstance(obj, dict):
        return {key: round_data(v) for key, v in obj.items()}
    if isinstance(obj, list):
        return [round_data(v) for v in obj]
    return float(f"{obj:.8g}") if isinstance(obj, (float, np.floating)) else obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    filename = HERE/"independent_results.json"
    value = round_data(response())
    serialized = (json.dumps(value, indent=2)+"\n").encode()
    if parser.parse_args().check:
        assert json.loads(filename.read_bytes()) == value
    else:
        filename.write_bytes(serialized)
    print(hashlib.sha256(serialized).hexdigest()[:12])
