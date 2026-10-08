"""1027: exact leading stress-source selection for the stated SM matter class.

Only Python/NumPy already present in the workspace are used.  Gaussian rational
arithmetic builds the actual representation equations; Fraction elimination gives
their full kernel.  This is a finite calibration of the proof, not a simulation
or construction of the interacting quantum field theory.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "charged_source_selection_results.json"


@dataclass(frozen=True)
class C:
    r: F = F(0)
    i: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "r", F(self.r))
        object.__setattr__(self, "i", F(self.i))

    def __add__(self, other):
        other = as_c(other)
        return C(self.r + other.r, self.i + other.i)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.r, -self.i)

    def __sub__(self, other):
        return self + -as_c(other)

    def __mul__(self, other):
        other = as_c(other)
        return C(self.r * other.r - self.i * other.i,
                 self.r * other.i + self.i * other.r)

    __rmul__ = __mul__

    def conj(self):
        return C(self.r, -self.i)

    def __bool__(self):
        return bool(self.r or self.i)


def as_c(value):
    return value if isinstance(value, C) else C(value)


def zeros(n, m=None):
    return [[C() for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    z = zeros(n)
    for i in range(n):
        z[i][i] = C(1)
    return z


def scale(a, x):
    return [[as_c(x) * v for v in row] for row in a]


def add(a, b):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def mm(a, b):
    return [[sum((x * y for x, y in zip(row, col)), C())
             for col in zip(*b)] for row in a]


def comm(a, b):
    return add(mm(a, b), scale(mm(b, a), -1))


def tr(a):
    return sum((a[i][i] for i in range(len(a))), C())


def dagger(a):
    return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]


def kron(a, b):
    return [[a[i][j] * b[k][l]
             for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def is_zero(a):
    return not any(v for row in a for v in row)


def hermitian_basis(n):
    out = []
    for i in range(n):
        a = zeros(n)
        a[i][i] = C(1)
        out.append((f"d{i}", a))
    for i in range(n):
        for j in range(i + 1, n):
            a = zeros(n)
            a[i][j] = a[j][i] = C(1)
            out.append((f"re{i}{j}", a))
            b = zeros(n)
            b[i][j], b[j][i] = C(0, 1), C(0, -1)
            out.append((f"im{i}{j}", b))
    return out


def symmetric_basis(n):
    out = []
    for i in range(n):
        for j in range(i, n):
            a = zeros(n)
            a[i][j] = a[j][i] = C(1)
            out.append((f"s{i}{j}", a))
    return out


def su_generators(n):
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = zeros(n), zeros(n)
            a[i][j] = a[j][i] = C(F(1, 2))
            b[i][j], b[j][i] = C(0, F(-1, 2)), C(0, F(1, 2))
            out.extend([a, b])
    for k in range(1, n):
        a = zeros(n)
        for j in range(k):
            a[j][j] = C(F(1, 2))
        a[k][k] = C(F(-k, 2))
        out.append(a)
    return out


def real_generator(t, include_sigma):
    """Realification of delta z = i T z, ordered (Re z, Im z)."""
    n = len(t)
    out = zeros(2 * n + int(include_sigma))
    for i in range(n):
        for j in range(n):
            out[i][j] = C(-t[i][j].i)
            out[i][j+n] = C(-t[i][j].r)
            out[i+n][j] = C(t[i][j].r)
            out[i+n][j+n] = C(-t[i][j].i)
    return out


SU2, SU3 = su_generators(2), su_generators(3)
# Physical right-handed fields are conjugated to a uniform left-Weyl convention.
MODULES = [
    ("q", 3, 2, F(1, 6), False),
    ("uc", 3, 1, F(-2, 3), True),
    ("dc", 3, 1, F(1, 3), True),
    ("l", 1, 2, F(-1, 2), False),
    ("ec", 1, 1, F(1), False),
]


def module_generators(spec):
    _, nc, nw, hypercharge, conjugate_color = spec
    out = [("Y", scale(eye(nc * nw), hypercharge))]
    if nc == 3:
        for t in SU3:
            color = [[-v.conj() for v in row] for row in t] if conjugate_color else t
            out.append(("3", kron(color, eye(nw))))
    if nw == 2:
        out.extend(("2", kron(eye(nc), t)) for t in SU2)
    return out


def push_equation(rows, equation):
    for part in ("r", "i"):
        row = {k: getattr(v, part) for k, v in equation.items()
               if getattr(v, part)}
        if row:
            pivot = row[min(row)]
            key = tuple((k, v / pivot) for k, v in sorted(row.items()))
            rows.add(key)


def rref(rows, n):
    mat = [[F(dict(row).get(i, 0)) for i in range(n)] for row in sorted(rows)]
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(mat)) if mat[i][c]), None)
        if p is None:
            continue
        mat[r], mat[p] = mat[p], mat[r]
        v = mat[r][c]
        mat[r] = [x / v for x in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][c]:
                v = mat[i][c]
                mat[i] = [x - v * y for x, y in zip(mat[i], mat[r])]
        pivots.append(c)
        r += 1
        if r == len(mat):
            break
    free = [c for c in range(n) if c not in pivots]
    kernel = []
    for c in free:
        v = [F(0)] * n
        v[c] = F(1)
        for i, p in enumerate(pivots):
            v[p] = -mat[i][c]
        kernel.append(v)
    return pivots, kernel


def solve(generations=3, include_sigma=True, portal=False, groups=("Y", "2", "3"),
          sterile_generations=0):
    labels = ["k_Y", "k_2", "k_3"]
    fl_basis = hermitian_basis(generations)
    indices = {}
    for spec in MODULES:
        name = spec[0]
        indices[name] = []
        for suffix, _ in fl_basis:
            indices[name].append(len(labels))
            labels.append(f"R_{name}_{suffix}")
    neutral_basis = hermitian_basis(sterile_generations)
    for suffix, _ in neutral_basis:
        labels.append(f"R_sterile_{suffix}")
    scalar_basis = symmetric_basis(4 + int(include_sigma))
    scalar_indices = []
    for suffix, _ in scalar_basis:
        scalar_indices.append(len(labels))
        labels.append(f"Q_{suffix}")
    index = {label: i for i, label in enumerate(labels)}
    rows = set()
    for spec in MODULES:
        name = spec[0]
        for factor, t in module_generators(spec):
            if factor not in groups:
                continue
            dim = len(t)
            for f in range(generations):
                for g in range(generations):
                    for i in range(dim):
                        for j in range(dim):
                            if not t[i][j]:
                                continue
                            eq = {idx: b[f][g] * t[i][j]
                                  for idx, (_, b) in zip(indices[name], fl_basis)
                                  if b[f][g]}
                            if f == g:
                                eq[index[f"k_{factor}"]] = -t[i][j]
                            push_equation(rows, eq)
    scalar_gens = [("Y", scale(eye(2), F(1, 2)))]
    scalar_gens += [("2", t) for t in SU2]
    for factor, complex_t in scalar_gens:
        if factor not in groups:
            continue
        t = real_generator(complex_t, include_sigma)
        products = [mm(b, t) for _, b in scalar_basis]
        commutators = [comm(b, t) for _, b in scalar_basis]
        for i in range(len(t)):
            for j in range(len(t)):
                # Gauge compatibility is explicit, though redundant after force locking.
                push_equation(rows, {idx: c[i][j] for idx, c in
                                     zip(scalar_indices, commutators) if c[i][j]})
                eq = {idx: p[i][j] for idx, p in zip(scalar_indices, products) if p[i][j]}
                if t[i][j]:
                    eq[index[f"k_{factor}"]] = -t[i][j]
                push_equation(rows, eq)
    if portal:
        assert include_sigma
        push_equation(rows, {index["Q_s44"]: C(1), index["k_Y"]: C(-1)})
    pivots, kernel = rref(rows, len(labels))
    common = [F(0)] * len(labels)
    for i, label in enumerate(labels):
        if label.startswith("k_") or (label.startswith("R_") and label.rsplit("_", 1)[1].startswith("d")):
            common[i] = F(1)
        if label.startswith("Q_s") and label[-1] == label[-2]:
            common[i] = F(1)
    assert all(sum((v * common[k] for k, v in row), F(0)) == 0 for row in rows)
    assert all(sum((v * vec[k] for k, v in row), F(0)) == 0
               for row in rows for vec in kernel)
    return {
        "unknown_count": len(labels), "independent_equation_rank": len(pivots),
        "deduplicated_equation_count": len(rows), "kernel_dimension": len(kernel),
        "labels": labels,
        "exact_constraint_rows": [{labels[k]: str(v) for k, v in row}
                                  for row in sorted(rows)],
        "kernel_basis": [{labels[i]: str(v) for i, v in enumerate(vec) if v}
                         for vec in kernel],
        "common_source_exactly_in_kernel": True,
        "groups_checked": list(groups), "portal_nonzero": portal,
        "sterile_free_flavors": sterile_generations,
    }


def representation_checks():
    lie_residuals = 0
    for generators in (SU2, SU3):
        for t in generators:
            assert dagger(t) == t and not tr(t)
        for a in generators:
            for b in generators:
                x = scale(comm(a, b), C(0, -1))
                reconstructed = zeros(len(a))
                for t in generators:
                    numerator, denominator = tr(mm(t, x)), tr(mm(t, t))
                    assert numerator.i == denominator.i == 0
                    reconstructed = add(reconstructed, scale(t, numerator.r / denominator.r))
                assert reconstructed == x
                lie_residuals += 1
    scalar_y = real_generator(scale(eye(2), F(1, 2)), True)
    scalar_square = mm(scalar_y, scalar_y)
    assert all(scalar_square[i][i] == C(F(-1, 4)) for i in range(4))
    assert scalar_square[4][4] == C()
    for spec in MODULES:
        for _, t in module_generators(spec):
            assert dagger(t) == t
        gens = module_generators(spec)
        for a, x in gens:
            for b, y in gens:
                if a != b:
                    assert is_zero(comm(x, y))
    return {"exact_lie_closure_pairs": lie_residuals,
            "three_generation_left_weyl_components": 45,
            "hypercharges": {s[0]: str(s[3]) for s in MODULES},
            "higgs_real_hypercharge_square": "-I_4/4",
            "all_matter_hypercharges_nonzero": True,
            "color_cartan_note": "Last su(3) Cartan generator is rescaled by sqrt(3); span and force constraints are unchanged."}


def ward_and_gauss_jet():
    eta = [F(-1), F(1), F(1), F(1)]
    ff = [[F(0)] * 4 for _ in range(4)]
    ff[1][0], ff[0][1] = F(2), F(-2)
    # ec spin polarization w=(0,1), bar-sigma=(I,-sigma), y=1: j=(1,0,0,1).
    current = [F(1), F(0), F(0), F(1)]
    dff = [[[F(0)] * 4 for _ in range(4)] for _ in range(4)]
    # j is delta S_m / delta A; the Maxwell equation is partial F + j = 0.
    dff[1][1][0], dff[1][0][1] = current[0], -current[0]
    for i in range(1, 4):
        dff[0][0][i], dff[0][i][0] = current[i], -current[i]
    div_f = [sum((eta[mu] * eta[nu] * dff[mu][mu][nu]
                  for mu in range(4)), F(0)) for nu in range(4)]
    bianchi = [dff[a][b][c] + dff[b][c][a] + dff[c][a][b]
               for a in range(4) for b in range(4) for c in range(4)]
    assert div_f == [-x for x in current] and all(x == 0 for x in bianchi)
    force = [sum((ff[nu][mu] * current[mu] for mu in range(4)), F(0))
             for nu in range(4)]
    delta = F(1, 2)
    weighted = [delta * x for x in force]
    assert weighted == [F(0), F(1), F(0), F(0)]
    # Independent vertex variation of L_Aj=+A_mu j^mu, modulo a total derivative:
    # (k-r) (partial_v A_mu) j^mu.  Select partial_1 A_0=2, partial_0 A_1=0.
    vertex_residual = -delta * ff[1][0] * current[0]
    assert vertex_residual == -weighted[1] == -1
    return {"F_lower": [[str(v) for v in r] for r in ff],
            "current_upper": list(map(str, current)),
            "partial_mu_F_upper_mu_nu": list(map(str, div_f)),
            "maxwell_equation_convention": "partial_mu F^{mu nu} + j^nu = 0, j=delta S_m/delta A",
            "bianchi_all_zero": True,
            "bad_source_weight_difference": str(delta),
            "bad_source_divergence_covector": list(map(str, weighted)),
            "integrated_vertex_variation_coefficient": str(vertex_residual),
            "common_even_monomial": "j and partial F entries are coefficients of one independent Grassmann-even bilinear; F_10=2 is a commuting background value.",
            "meaning": "Exact local Ward coefficient and compatible Maxwell first-jet coefficient; +/-1 is not a quantum-state expectation or a global coupled solution."}


def basis_transport():
    # Deterministic changes of flavor and Higgs coordinates; not active permissions.
    rng = np.random.default_rng(1027)
    worst = 0.0
    cases = 0
    for spec in MODULES:
        ninternal = spec[1] * spec[2]
        x = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        u, _ = np.linalg.qr(x)
        r = np.diag([1.0, 1.4, -0.3])
        full_u = np.kron(u, np.eye(ninternal))
        full_r = np.kron(r, np.eye(ninternal))
        for _, t in module_generators(spec):
            tt = np.array([[complex(float(v.r), float(v.i)) for v in row] for row in t])
            full_t = np.kron(np.eye(3), tt)
            original = (full_r - 0.7 * np.eye(3*ninternal)) @ full_t
            transported = ((full_u @ full_r @ full_u.conj().T -
                            0.7 * np.eye(3*ninternal)) @ full_t)
            worst = max(worst, float(np.max(np.abs(
                transported - full_u @ original @ full_u.conj().T))))
            cases += 1
    x = rng.normal(size=(5, 5))
    o, _ = np.linalg.qr(x)
    q = np.diag([1., 1., 1., 1., 2.])
    ty = real_generator(scale(eye(2), F(1, 2)), True)
    t = np.array([[float(v.r) for v in row] for row in ty])
    worst = max(worst, float(np.max(np.abs(
        (o @ q @ o.T) @ (o @ t @ o.T) - o @ (q @ t) @ o.T))))
    assert worst < 1e-12
    return {"fermion_generator_cases": cases, "max_transport_residual": worst,
            "scope": "Simultaneously transports representations and source matrices; not a graph connectivity argument."}


def portal_certificate():
    # Phi=(h1+i h3,h2+i h4)/sqrt(2): U_portal=kappa*sigma^2*|h|^2/4.
    # Mixed Hessian H_(h_i,sigma)=kappa*h_i*sigma.  Hence [Q,H] coefficient
    # is (q_H-q_sigma)*kappa, independently of all nonmixed potential terms.
    x, sigma, kappa, qh, qs = F(2), F(3), F(2, 5), F(1), F(7, 4)
    mixed = kappa * x * sigma
    residual = (qh - qs) * mixed
    assert residual == F(-9, 5)
    return {"normalization": "U_portal=kappa*sigma^2*sum(h_i^2)/4",
            "mixed_hessian": str(mixed), "bad_commutator_entry": str(residual),
            "vanishes_for_equal_weights_or_zero_portal": True,
            "theorem_reused": "1017; this is the neutral extension, not a new scalar classification."}


def run():
    cases = {
        "three_generation_charged_only": solve(include_sigma=False),
        "hypercharge_force_only_diagnostic": solve(groups=("Y",)),
        "B993_zero_portal": solve(),
        "B993_nonzero_portal": solve(portal=True),
        "one_generation_nonzero_portal": solve(generations=1, portal=True),
        "three_free_neutral_fermions_extra": solve(portal=True, sterile_generations=3),
    }
    expected = {"three_generation_charged_only": 1,
                "hypercharge_force_only_diagnostic": 4,
                "B993_zero_portal": 2, "B993_nonzero_portal": 1,
                "one_generation_nonzero_portal": 1,
                "three_free_neutral_fermions_extra": 10}
    for label, wanted in expected.items():
        assert cases[label]["kernel_dimension"] == wanted, (label, cases[label])
    sources = [
        BASE / "archive_1009_/research_note_1026.md",
        BASE / "archive_1009_/1026/input_dependency_update_v0_15.md",
        BASE / "archive_1009_/1026/NEXT.md",
        BASE / "archive_1009_/research_note_1017.md",
        BASE / "archive_1009_/research_note_1018.md",
        BASE / "archive_1009_/research_note_1024.md",
        BASE / "archive_301_341/research_note_326.md",
        BASE / "archive_702_741/research_note_732.md",
        BASE / "archive_702_741/research_note_735.md",
        BASE / "archive_956_989/981/drafts/common_parent_contract_v1.md",
        BASE / "archive_990_1008/993/common_candidate_v1.md",
    ]
    return {
        "round": 1027, "status": "scientific_calibration_verified",
        "new_calibration_groups": 1, "cumulative_test_groups": 3804,
        "new_cognitive_axioms": 0, "goal_complete": False,
        "all_scientific_calibrations_passed": True,
        "scope": "Constant gauge-compatible leading kinetic stress weights; no-derivative interaction completion; classical necessary Ward conditions.",
        "inputs_retained": ["common flat Lorentz kinetic principal part", "given SM gauge representations",
                            "nonzero three gauge couplings", "specified leading source class",
                            "arbitrary local on-shell jets", "nonzero portal only in its labeled case"],
        "not_proved": ["cognitive derivation", "selection of background metric or species",
                       "all higher-derivative/nonminimal currents", "full interacting quantum theory",
                       "nonlinear multi-metric completion or absolute gravitational coupling"],
        "representation_checks": representation_checks(),
        "exact_source_kernels": cases,
        "ward_and_gauss_jet": ward_and_gauss_jet(),
        "basis_transport": basis_transport(),
        "portal_certificate": portal_certificate(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "historical_source_sha256": {
            p.relative_to(BASE).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sources},
    }


def compare(actual, expected, path="root"):
    """Strict structure/exact arithmetic, with roundoff tolerance for NumPy checks."""
    if isinstance(actual, dict):
        assert isinstance(expected, dict) and actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(actual, list):
        assert isinstance(expected, list) and len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, float):
        assert isinstance(expected, (int, float)) and np.isfinite(actual) and np.isfinite(expected), path
        assert np.isclose(actual, expected, atol=1e-12, rtol=1e-10), (path, actual, expected)
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(output, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(output, json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps({"status": "passed", "cases": {
        k: {s: v[s] for s in ("unknown_count", "independent_equation_rank", "kernel_dimension")}
        for k, v in output["exact_source_kernels"].items()}, "transport": output["basis_transport"],
        "mode": "exclusive_first_write" if args.write else "read_only_recompute_compare",
        "output": str(OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
