"""1010: conditional hypercharge selection by anomaly and mass compatibility.

All charges are left-handed Weyl charges; scalars do not enter anomaly sums.
The note supplies analytic completeness proofs. These checks are finite
calibrations and do not derive field content, a mass mechanism, or a QFT.
Run --write once to create an exclusive result; default verifies that result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "hypercharge_neutrino_selection_results.json"
H = F(1, 2)
POSITIONS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
T_VALUES = (F(-2), F(-1, 3), F(1, 7), F(1), F(3))
SOURCES = (
    "archive_531_553/research_note_531.md",
    "archive_531_553/531/unified_physics_condition_ledger.md",
    "archive_531_553/research_note_532.md",
    "archive_531_553/research_note_539.md",
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
)


def species(x, zeta=None):
    """Build actual Weyl species, with dimensions and 2T(fund)=1 indices."""
    rows = []
    for family in range(3):
        if zeta is None:
            q, u, d = H / 3, -4 * H / 3, 2 * H / 3
            ell, e = x[family] - H, 2 * H - x[family]
        else:
            q, u, d = H / 3 + zeta / 3, -4 * H / 3 - zeta / 3, 2 * H / 3 - zeta / 3
            ell, e = -H - zeta, 2 * H + zeta
        for name, nc, nw, charge in (
            ("Q", 3, 2, q), ("u^c", 3, 1, u), ("d^c", 3, 1, d),
            ("L", 1, 2, ell), ("e^c", 1, 1, e),
        ):
            rows.append({"family": family, "species": name, "nc": nc,
                         "nw": nw, "y": charge})
        if zeta is not None:
            rows.append({"family": family, "species": "nu^c", "nc": 1,
                         "nw": 1, "y": zeta})
    return rows


def traces(rows):
    """Independent representation traces, without using reduced cubic form."""
    out = {key: F(0) for key in ("SU3_squared_U1", "SU2_squared_U1", "gravity_U1", "U1_cubed")}
    for row in rows:
        nc, nw, y = row["nc"], row["nw"], row["y"]
        if nc == 3:
            out["SU3_squared_U1"] += nw * y
        if nw == 2:
            out["SU2_squared_U1"] += nc * y
        out["gravity_U1"] += nc * nw * y
        out["U1_cubed"] += nc * nw * y ** 3
    return out


def rational_dict(values):
    return {key: str(value) for key, value in values.items()}


def branch(zero, t=F(1)):
    x = [F(0)] * 3
    nonzero = [i for i in range(3) if i != zero]
    x[nonzero[0]], x[nonzero[1]] = t, -t
    return x


def mass_matrix(zero, a, b):
    out = np.zeros((3, 3), dtype=complex)
    i, j = [k for k in range(3) if k != zero]
    out[i, j] = out[j, i] = a
    out[zero, zero] = b
    return out


def frobenius(x):
    return float(np.linalg.norm(x, "fro"))


def spectrum(c):
    return np.linalg.eigvalsh(c.conj().T @ c)


def constraint_matrix(c):
    """Real matrix for diagonal real X: XC+CX=0, together with Tr X=0."""
    rows = [np.ones(3)]
    for i, j in POSITIONS:
        coeff = np.zeros(3, dtype=complex)
        coeff[i] += c[i, j]
        coeff[j] += c[i, j]
        rows.extend((coeff.real, coeff.imag))
    return np.asarray(rows)


def hermitian_constraint_singular_values(c):
    """All nine real Hermitian directions, not only a chosen diagonal basis."""
    basis = []
    for i in range(3):
        e = np.zeros((3, 3), dtype=complex)
        e[i, i] = 1
        basis.append(e)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        e = np.zeros((3, 3), dtype=complex)
        e[i, j] = e[j, i] = 1
        basis.append(e)
        a = np.zeros((3, 3), dtype=complex)
        a[i, j], a[j, i] = 1j, -1j
        basis.append(a)
    columns = []
    for x in basis:
        value = x.T @ c + c @ x
        columns.append(np.concatenate((value.real.ravel(), value.imag.ravel(), [np.trace(x).real])))
    return np.linalg.svd(np.asarray(columns).T, compute_uv=False)


def main_result():
    exact = []
    for zero in range(3):
        for t in T_VALUES:
            x = branch(zero, t)
            rows = species(x)
            anomalies = traces(rows)
            assert all(value == 0 for value in anomalies.values())
            assert sum(x) == 0 and sum(a ** 3 for a in x) == 0
            # Same H and actual nonzero charged Yukawa selection rules.
            for family in range(3):
                y = {r["species"]: r["y"] for r in rows if r["family"] == family}
                assert y["Q"] + H + y["u^c"] == 0
                assert y["Q"] - H + y["d^c"] == 0
                assert y["L"] - H + y["e^c"] == 0
            exact.append({"zero_family": zero, "t": str(t), "x": list(map(str, x)),
                          "anomaly_traces": rational_dict(anomalies),
                          "weyl_components": sum(r["nc"] * r["nw"] for r in rows)})
    # Nonzero cubic must be detected; sum x=0 alone is insufficient.
    wrong = traces(species([F(1), F(1), F(-2)]))
    assert wrong["gravity_U1"] == 0 and wrong["U1_cubed"] == -6
    assert all(v == 0 for v in traces(species([F(0)] * 3)).values())

    supports = []
    for bits in itertools.product((0, 1), repeat=6):
        allowed = []
        for zero in range(3):
            x = branch(zero)
            if all(not used or x[i] + x[j] == 0 for used, (i, j) in zip(bits, POSITIONS)):
                allowed.append(zero)
        supports.append({"bits": "".join(map(str, bits)), "nonzero_branches": allowed})
    survivors = [s for s in supports if s["nonzero_branches"]]
    assert len(supports) == 64 and len(survivors) == 10
    assert sum(not s["nonzero_branches"] for s in supports) == 54
    assert [sum(z in s["nonzero_branches"] for s in supports) for z in range(3)] == [4, 4, 4]

    # Complex phases, paired mass degeneracy, and both zero boundaries.
    params = ((0j, 0j), (0j, 2 + 3j), (2 - 1j, 0j), (2 + 3j, 1 - 4j),
              (1j, 1 + 0j), (0.25 - 0.75j, -2 + 0.3j))
    max_rule_error = max_spectrum_error = max_gram_error = 0.0
    spectral_checks = []
    for zero in range(3):
        for a, b in params:
            c = mass_matrix(zero, a, b)
            x = np.diag([float(v) for v in branch(zero)])
            expected = np.array([abs(a) ** 2] * 3)
            expected[zero] = abs(b) ** 2
            max_rule_error = max(max_rule_error, frobenius(x @ c + c @ x))
            max_gram_error = max(max_gram_error, frobenius(c.conj().T @ c - np.diag(expected)))
            residual = float(np.max(np.abs(spectrum(c) - np.sort(expected))))
            max_spectrum_error = max(max_spectrum_error, residual)
            spectral_checks.append({"zero_family": zero, "a": [a.real, a.imag],
                                    "b": [b.real, b.imag], "mass_squared": spectrum(c).tolist()})
    assert max(max_rule_error, max_spectrum_error, max_gram_error) < 1e-12

    # Every surviving support can be populated with nonzero complex entries.
    max_support_gap = 0.0
    for support in survivors:
        for zero in support["nonzero_branches"]:
            c = np.zeros((3, 3), dtype=complex)
            for index, (bit, (i, j)) in enumerate(zip(support["bits"], POSITIONS)):
                if bit == "1":
                    c[i, j] = c[j, i] = (index + 1) + 0.7j
            eigs = spectrum(c)
            gap = float(np.min(np.diff(eigs)))
            max_support_gap = max(max_support_gap, abs(gap))
    assert max_support_gap < 1e-12

    diag = np.diag([0.0, 2.0, 5.0]).astype(complex)
    f = np.exp(2j * np.pi * np.outer(np.arange(3), np.arange(3)) / 3) / np.sqrt(3)
    dense = f @ diag @ f.T
    assert frobenius(f.conj().T @ f - np.eye(3)) < 1e-12
    positive = []
    for name, c in (("one_massless", diag), ("dense_unitary_congruence_one_massless", dense)):
        assert frobenius(c - c.T) < 1e-12
        sv = np.linalg.svd(constraint_matrix(c), compute_uv=False)
        assert sv[-1] > 0.1  # No diagonal real charge, including the massless direction.
        pair_residuals = [frobenius(np.diag(np.asarray(branch(z, 1), dtype=float)) @ c
                                   + c @ np.diag(np.asarray(branch(z, 1), dtype=float))) for z in range(3)]
        assert min(pair_residuals) > 1
        eigs = spectrum(c)
        assert np.max(np.abs(eigs - [0, 4, 25])) < 1e-12
        positive.append({"case": name, "mass_squared": eigs.tolist(),
                         "trace_plus_selection_smallest_singular_value": float(sv[-1]),
                         "nonzero_branch_rule_residuals_at_t_1": pair_residuals})

    hermitian_checks = []
    for name, c, expected_nullity in (
        ("one_massless", diag, 0),
        ("dense_unitary_congruence_one_massless", dense, 0),
        ("positive_degenerate_pair", np.diag([2., 2., 5.]).astype(complex), 1),
        ("two_zero_masses", np.diag([0., 0., 5.]).astype(complex), 3),
    ):
        sv = hermitian_constraint_singular_values(c)
        nullity = int(np.sum(sv < 1e-10))
        assert nullity == expected_nullity
        hermitian_checks.append({"case": name, "real_parameter_count": 9,
                                 "constraint_singular_values": sv.tolist(),
                                 "kernel_dimension": nullity})
    pair_generator = np.array([[0, 1j, 0], [-1j, 0, 0], [0, 0, 0]], dtype=complex)
    degenerate = np.diag([2., 2., 5.])
    assert frobenius(pair_generator.T @ degenerate + degenerate @ pair_generator) == 0
    assert np.trace(pair_generator) == 0

    dirac = []
    for zeta in (F(-2), F(-1, 3), F(0), F(1, 7), F(1), F(3)):
        rows = species(None, zeta)
        anomalies = traces(rows)
        assert all(value == 0 for value in anomalies.values())
        y = {r["species"]: r["y"] for r in rows if r["family"] == 0}
        assert y["L"] + H + y["nu^c"] == 0
        assert 2 * y["nu^c"] == 2 * zeta
        # All three generations carry the same charge: every Dirac Y entry is allowed.
        d = f @ np.diag([1.0, 2.0, 5.0]) @ f.conj().T
        assert np.max(np.abs(spectrum(d) - [1, 4, 25])) < 1e-12
        dirac.append({"zeta": str(zeta), "anomaly_traces": rational_dict(anomalies),
                      "neutral_yukawa_total_charge": str(y["L"] + H + y["nu^c"]),
                      "bare_majorana_total_charge": str(2 * zeta),
                      "bare_majorana_allowed": zeta == 0,
                      "dirac_mass_squared": [1, 4, 25]})

    source_hashes = {rel: hashlib.sha256((BASE / rel).read_bytes()).hexdigest() for rel in SOURCES}
    return {
        "round": 1010, "cumulative_research_groups": 3788, "new_calibration_groups": 1,
        "scope": "Given three-generation SM non-Abelian representations, one nonzero-charge Higgs, connected quark Yukawas, full-rank charged-lepton Yukawa, local anomalies, and a dimension-five neutrino mass operator; no right-handed neutrinos in the main branch.",
        "charge_convention": "All fermions left-handed; H=1/2; SU(N)^2 U(1) traces use 2T(fund)=1.",
        "global_scope": "Classification is at Lie-algebra charge level. Rational sample charges can be integrated after choosing a common charge lattice; they are not asserted to share a previously fixed U(1) period, quotient, or all global anomaly conditions.",
        "analytic_obligations": "The note proves branch completeness, exact singular-value pairing and the implication from pairwise distinct mass squares. Finite checks below are not substitutes for those proofs.",
        "nonzero_branches": exact,
        "negative_cubic_test": rational_dict(wrong),
        "support_order": [list(p) for p in POSITIONS],
        "support_classification": {"total": 64, "admitting_nonzero_charge_branch": 10,
                                   "excluding_all_nonzero_branches": 54, "rows": supports},
        "mass_pairing": {"samples": spectral_checks, "max_gauge_rule_error": max_rule_error,
                         "max_gram_diagonal_error": max_gram_error,
                         "max_spectrum_error": max_spectrum_error,
                         "max_nearest_mass_squared_gap_for_surviving_supports": max_support_gap},
        "nondegenerate_positive_examples": positive,
        "full_hermitian_generator_checks": hermitian_checks,
        "right_handed_dirac_sensitivity": dirac,
        "interpretation": "Without the neutrino spectral requirement, anomaly-free family-charge branches survive. Pairwise distinct mass squares, including one zero mass, exclude them under the stated main-branch inputs. Adding nu_R with purely Dirac masses changes those inputs and preserves a B-L deformation.",
        "claim_boundaries": {"derived_majorana": False, "SM_group_generated": False,
                             "generation_count_generated": False, "gauge_coupling_prediction": False,
                             "full_quantum_realization": False, "cognitive_principles_alone_select_hypercharge": False,
                             "arbitrary_global_group_certified": False, "neutrino_mass_values_predicted": False},
        "source_sha256": source_hashes,
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "all_checks_passed": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = main_result()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
        action = "written_exclusive"
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert saved == result, "Saved result differs from fresh calibration or source hashes."
        action = "verified_saved_result"
    print(json.dumps({"round": 1010, "status": action, "all_checks_passed": True,
                      "anomaly_branches": 15, "supports": 64, "surviving_supports": 10,
                      "pairing_samples": 18, "dirac_sensitivity_samples": 6}, ensure_ascii=False))


if __name__ == "__main__":
    main()
