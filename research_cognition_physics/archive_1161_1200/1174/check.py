"""1174 finite CAR-sector process. Default: recompute and compare, never write.

--write creates results.json once. Matrices are exact invariant N_S=N_P=1
blocks, not a Dirac-sea truncation. Hardware-error bounds are explicit inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
SOURCES = [
    "archive_653_701/research_note_662.md",
    "archive_342_369/research_note_349.md",
    "archive_342_369/research_note_350.md",
    "archive_1120_1158/1129/proof.md",
    "archive_1120_1158/1146/proof.md",
    "archive_1120_1158/1155/proof.md",
    "archive_1161_/1173/proof.md",
    "archive_1161_/_shared/goal.md",
]
EPS, OMEGA, KAPPA, T, BASE = .3, .6, .2, 2.5, 53.
MASSES = (2., 2.4)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
I6, I8, I36 = (np.eye(n, dtype=complex) for n in (6, 8, 36))


def kron(*xs):
    out = np.array([[1.]], complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def norm(x):
    return float(np.linalg.norm(x, 2))


def cjson(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def expval(v, a):
    return float(np.vdot(v, a @ v).real)


GX = [kron(*[X if i == j else I2 for j in range(3)]) for i in range(3)]
GY = [kron(*[Y if i == j else I2 for j in range(3)]) for i in range(3)]
GZ = [kron(*[Z if i == j else I2 for j in range(3)]) for i in range(3)]
PP = [kron(np.diag([float(i == j) for j in range(3)]), I2) for i in range(3)]
SOURCE_CHIR = np.array([1., 1.]) / np.sqrt(2)
PROBE = np.kron(np.array([1., 1., 2.]) * np.exp(2j * np.pi * np.arange(3) / 3) / np.sqrt(6),
                [np.cos(np.pi / 8), np.sin(np.pi / 8)])
G0 = np.ones(8) / np.sqrt(8)
E_CLOCK = (I2 - Y) / 2


def local_energy(m):
    hs = [m * kron(np.diag([float(i == j) for j in range(3)]), X) for i in range(3)]
    hops = []
    for i in range(3):
        j = (i + 1) % 3
        edge = np.zeros((3, 3), complex)
        edge[i, j], edge[j, i] = -.5j, .5j
        hop = kron(edge, Z)
        hs[i] += hop / 2
        hs[j] += hop / 2
        hops.append(hop)
    return hs, hops


def model(ms, mp, omega):
    hs, _ = local_energy(ms)
    hp, _ = local_energy(mp)
    local = [kron(hs[i], I6) + kron(I6, hp[i]) for i in range(3)]
    terms = {
        "source": kron(sum(hs), I6, I8),
        "probe": kron(I6, sum(hp), I8),
        "geometry": kron(I36, BASE * I8 - omega * sum(GX)
                         - KAPPA * sum(GZ[i] @ GZ[(i + 1) % 3] for i in range(3))),
        "interaction": EPS * sum(kron(local[i], GZ[i]) for i in range(3)),
    }
    return sum(terms.values()), terms, local


def prepared(r):
    src = np.kron(np.eye(3)[r], SOURCE_CHIR)
    return np.kron(np.kron(src, PROBE), G0)


def spectral_run(h):
    w, v = np.linalg.eigh(h)
    residual = norm(h @ v - v * w)
    orth = norm(v.conj().T @ v - np.eye(len(w)))
    polar = orth / (1 + math.sqrt(1 - orth))
    # Use the analytic full-parent bound, not an uncorrected approximate eigenvalue norm.
    hnorm = BASE + 3 * (8.8 + 4 + EPS * 8.8 + EPS * 4 + OMEGA + KAPPA)
    # A posteriori polar/Duhamel bound, plus declared floating arithmetic allowance.
    uniform = 2 * polar + polar * polar + T * (residual + (hnorm + max(abs(w))) * polar) + 3e-12
    out = {r: v @ (np.exp(-1j * T * w) * (v.conj().T @ prepared(r))) for r in (0, 1)}
    return out, {"eigen_residual": residual, "orthogonality_defect": orth,
                 "uniform_propagator_bound": float(uniform), "min_eigenvalue": float(w[0])}


def full_effect(outputs, op, internal_effect=None):
    """Compressed 4x4 effect on unknown D_S x D_P; preserves cross-mass phases."""
    if internal_effect is None:
        internal_effect = I2
    f = np.zeros((4, 4), complex)
    for ds in range(2):
        for a in range(2):
            for b in range(2):
                f[2 * ds + a, 2 * ds + b] = internal_effect[a, b] * np.vdot(
                    outputs[ds, a], op @ outputs[ds, b])
    return f


def witness_probability(f):
    w = np.array([1., 1., 0., 0.]) / np.sqrt(2)
    return expval(w, f)


def compute():
    source_hashes = {p: hashlib.sha256((RESEARCH / p).read_bytes()).hexdigest() for p in SOURCES}
    states, hams, terms, specs = {}, {}, {}, {}
    for label, omega, masses in (("active", OMEGA, MASSES), ("frozen", 0., MASSES),
                                  ("massless", OMEGA, (0.,)), ("massless_frozen", 0., (0.,))):
        for ds, ms in enumerate(MASSES):
            for dp, mp in enumerate(masses):
                key = (label, ds, dp)
                hams[key], terms[key], _ = model(ms, mp, omega)
                states[key], specs[str(key)] = spectral_run(hams[key])
    max_uniform = max(v["uniform_propagator_bound"] for v in specs.values())
    numerical_probability_bound = 2 * max_uniform + max_uniform * max_uniform
    assert numerical_probability_bound < 1e-9

    # Exact source equation, checked independently as double commutator.
    h = hams["active", 0, 0]
    _, _, loc = model(2., 2., OMEGA)
    eq_defects = []
    for i in range(3):
        z = kron(I36, GZ[i])
        dd = -(h @ (h @ z - z @ h) - (h @ z - z @ h) @ h)
        rhs = -4 * OMEGA**2 * z - 4 * OMEGA * EPS * kron(loc[i], GX[i])
        rhs += 4 * OMEGA * KAPPA * kron(I36, GX[i] @ (GZ[(i - 1) % 3] + GZ[(i + 1) % 3]))
        eq_defects.append(norm(dd - rhs))
    assert max(eq_defects) < 1e-10

    free_spectrum = {}
    for ds, mass in enumerate(MASSES):
        hs, _ = local_energy(mass)
        wf, vf = np.linalg.eigh(sum(hs))
        probs = []
        # Individual degenerate eigenvectors need not respect translations; sum by eigenspace.
        groups = []
        for index, value in enumerate(wf):
            if not groups or abs(value - wf[groups[-1][0]]) > 1e-10:
                groups.append([index])
            else:
                groups[-1].append(index)
        for r in (0, 1):
            source = np.kron(np.eye(3)[r], SOURCE_CHIR)
            q = abs(vf.conj().T @ source)**2
            probs.append([float(q[g].sum()) for g in groups])
        assert np.max(abs(np.array(probs[0]) - probs[1])) < 1e-12
        free_spectrum[str(ds)] = {"energies": [float(wf[g[0]]) for g in groups], "probabilities": probs}

    energy_records, energy_defects, total_variance = {}, [], []
    for r in (0, 1):
        v0 = prepared(r)
        record = {}
        for name in ("source", "probe", "geometry", "interaction"):
            before, after = 0., 0.
            for dp in range(2):
                key = ("active", 0, dp)
                before += expval(v0, terms[key][name]) / 2
                after += expval(states[key][r], terms[key][name]) / 2
            record[name] = {"before": before, "after": after, "change": after - before}
        defect = abs(sum(v["change"] for v in record.values()))
        energy_defects.append(defect)
        energy_records[str(r)] = record
        mean = sum(v["before"] for v in record.values())
        second = sum(float(np.vdot(hams["active", 0, dp] @ v0,
                                  hams["active", 0, dp] @ v0).real) / 2 for dp in range(2))
        total_variance.append(second - mean * mean)
    assert max(energy_defects) < 1e-10
    compressed_mean_diffs = []
    for ds in range(2):
        for dp in range(2):
            hh = hams["active", ds, dp]
            compressed_mean_diffs.append(abs(expval(prepared(0), hh) - expval(prepared(1), hh)))
    assert max(compressed_mean_diffs) < 1e-12

    clocks, g_pops, ports, effect_diffs = {}, {}, {}, {}
    for label in ("active", "frozen"):
        clocks[label], ports[label], g_pops[label] = {}, {}, {}
        full_clock = {}
        for r in (0, 1):
            oo = {(ds, dp): states[label, ds, dp][r] for ds in range(2) for dp in range(2)}
            clocks[label][str(r)], ports[label][str(r)], g_pops[label][str(r)] = [], [], []
            for j in range(3):
                op = kron(I6, PP[j], I8)
                fc = full_effect(oo, op, E_CLOCK)
                clocks[label][str(r)].append(witness_probability(fc))
                ports[label][str(r)].append(witness_probability(full_effect(oo, op)))
                if j == 0:
                    full_clock[r] = fc
                gp = full_effect(oo, kron(I36, (I8 + GZ[j]) / 2))
                g_pops[label][str(r)].append(witness_probability(gp))
        effect_diffs[label] = {"source0": cjson(full_clock[0]), "source1": cjson(full_clock[1]),
                               "operator_difference": norm(full_clock[0] - full_clock[1])}
    frozen_gap = max(abs(np.array(clocks["frozen"]["0"]) - clocks["frozen"]["1"]))
    assert frozen_gap < 1e-12 and effect_diffs["frozen"]["operator_difference"] < 1e-12

    signal = {}
    for label in ("massless", "massless_frozen"):
        signal[label] = {str(r): [expval(states[label, 0, 0][r], kron(I6, PP[j], I8))
                                   for j in range(3)] for r in (0, 1)}
    assert max(abs(np.array(signal["massless_frozen"]["0"]) - signal["massless_frozen"]["1"])) < 1e-12
    clock_gap = abs(clocks["active"]["0"][0] - clocks["active"]["1"][0])
    population_gap = abs(g_pops["active"]["0"][0] - g_pops["active"]["1"][0])
    signal_gap = abs(signal["massless"]["0"][1] - signal["massless"]["1"][1])
    assert clock_gap > .007 and population_gap > .007 and signal_gap > .028

    bound_sum = 3 * (8.8 + 4 + EPS * 8.8 + EPS * 4 + OMEGA + KAPPA)
    r_slices = math.ceil(T*T*bound_sum*bound_sum/(2e-4))
    trotter = T*T*bound_sum*bound_sum/(2*r_slices)
    eta = 4e-4 + 1e-9
    lowers = {"clock_site0": clock_gap - 2 * eta, "G_population_site0": population_gap - 2 * eta,
              "massless_site1": signal_gap - 2 * eta}
    assert min(lowers.values()) > .0065
    result = {
        "round": 1174, "date": "2026-10-13", "status": "PASS",
        "parameters": {"epsilon": EPS, "omega": OMEGA, "kappa": KAPPA, "T": T,
                       "masses": list(MASSES), "constant_energy": BASE, "sites": 3},
        "sources_sha256": source_hashes,
        "groups": {
            "complete_source_and_energy": {"status": "PASS", "parent_dimension": 2**27,
                "invariant_dimension": 1152, "mass_block_dimension": 288,
                "source_equation_defects": eq_defects, "free_source_energy_spectral_probabilities": free_spectrum,
                "full_H_initial_mean_differences_all_four_internal_basis": compressed_mean_diffs,
                "full_H_initial_variances_witness_source0_source1": total_variance,
                "full_H_spectral_distributions_claimed_equal": False,
                "energy_exchange": energy_records, "energy_total_defects": energy_defects},
            "same_clock_position_and_active_register": {"status": "PASS", "clock_minusY_site_probabilities": clocks,
                "ordinary_site_probabilities": ports, "G_plusZ_population_probabilities": g_pops,
                "compressed_unknown_4D_clock_site0_effects": effect_diffs,
                "clock_site0_gap": clock_gap, "G_population_site0_gap": population_gap,
                "frozen_nonselective_clock_max_gap": float(frozen_gap)},
            "massless_same_family": {"status": "PASS", "site_probabilities": signal,
                "site1_gap": signal_gap, "separate_experiment": True},
            "full_input_finite_effective_contract": {"status": "PASS", "local_term_count": 18,
                "full_CAR_term_norm_sum": bound_sum, "full_parent_positive_lower_bound": BASE-bound_sum,
                "trotter_slices": r_slices, "local_gate_calls": 18*r_slices,
                "trotter_uniform_error": trotter, "prep_isometry_error_input": 1e-4,
                "sum_actual_gate_error_input": 1e-4, "full_record_instrument_halfdiamond_input": 1e-4,
                "numerical_allowance": 1e-9, "max_computed_spectral_propagator_bound": max_uniform,
                "numerical_probability_bound": numerical_probability_bound,
                "spectral_checks": specs, "single_task_probability_error_budget": eta,
                "two_preparation_gap_lower_bounds": lowers,
                "actual_hardware_error_certified": False, "GR_or_classical_geometry_certified": False,
                "prepared_unknown_reference_preserved": True},
        },
    }
    return result


def compare(new, old, path="root"):
    if isinstance(new, dict):
        assert new.keys() == old.keys(), path
        for k in new:
            compare(new[k], old[k], path + "." + str(k))
    elif isinstance(new, list):
        assert len(new) == len(old), path
        for i, (a, b) in enumerate(zip(new, old)):
            compare(a, b, path + "." + str(i))
    elif isinstance(new, float):
        assert math.isclose(new, old, rel_tol=2e-8, abs_tol=2e-10), (path, new, old)
    else:
        assert new == old, (path, new, old)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="create results once; never overwrite")
    parser.add_argument("--preview", action="store_true", help="compute and print small summary without writing")
    args = parser.parse_args()
    result = compute()
    target = HERE / "results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
        status = "WRITE_ONCE"
    elif args.preview:
        status = "UNSAVED_PREVIEW"
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        status = "READ_ONLY_MATCH"
    g = result["groups"]
    print(json.dumps({"status": status, "groups": len(g), "sources": len(SOURCES),
                      "clock_gap": g["same_clock_position_and_active_register"]["clock_site0_gap"],
                      "massless_gap": g["massless_same_family"]["site1_gap"],
                      "lower_bounds": g["full_input_finite_effective_contract"]["two_preparation_gap_lower_bounds"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
