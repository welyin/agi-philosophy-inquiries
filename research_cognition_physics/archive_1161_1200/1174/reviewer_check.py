"""Independent 1174 audit. Default is read-only; --write creates a new report.

This implementation was written from proof.md/admission.md before author check.py
was available. It constructs the full 1152-dimensional N_S=N_P=1 sector, uses
full-matrix eigendecomposition, and checks state evolution with scaled Taylor
polynomials. Standard floating point diagnostics are not interval arithmetic.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
I12 = np.eye(12, dtype=complex)
I8 = np.eye(8, dtype=complex)
T, EPS, OMEGA, KAPPA, BASE = 2.5, .3, .6, .2, 53.
ETA = .000400001
CORE = ("admission.md", "proof.md", "check.py", "results.json", "NEXT.md", "../research_note_1174.md")


def kron(*args):
    out = np.array([[1.]], dtype=complex)
    for a in args:
        out = np.kron(out, a)
    return out


def gop(op, i):
    return kron(*[op if j == i else I2 for j in range(3)])


GX = [gop(X, i) for i in range(3)]
GY = [gop(Y, i) for i in range(3)]
GZ = [gop(Z, i) for i in range(3)]


def material(masses):
    onsite, edges = [], []
    for i in range(3):
        position = np.zeros((3, 3), complex)
        position[i, i] = 1
        onsite.append(kron(position, X, np.diag(masses)))
        edge = np.zeros((3, 3), complex)
        edge[i, (i + 1) % 3] = -.5j
        edge[(i + 1) % 3, i] = .5j
        edges.append(kron(edge, Z, I2))
    local = [onsite[i] + .5 * (edges[i - 1] + edges[i]) for i in range(3)]
    return sum(local), local, onsite, edges


def construct(probe_masses=(2., 2.4), omega=OMEGA):
    hs, ks, os_, es = material((2., 2.4))
    hp, kp, op_, ep = material(probe_masses)
    hg = -omega * sum(GX) - KAPPA * sum(GZ[i] @ GZ[(i + 1) % 3] for i in range(3))
    parts = {
        "source": kron(hs, I12, I8),
        "probe": kron(I12, hp, I8),
        "register_without_common_baseline": kron(I12, I12, hg),
        "interaction": EPS * sum(kron(ks[i], I12, GZ[i]) + kron(I12, kp[i], GZ[i]) for i in range(3)),
    }
    h = sum(parts.values())
    local_norms = [float(np.max(np.abs(np.linalg.eigvalsh(a)))) for a in ks + kp]
    bound = float(max(abs(np.linalg.eigvalsh(hs))) + max(abs(np.linalg.eigvalsh(hp)))
                  + 3 * abs(omega) + 3 * KAPPA + EPS * sum(local_norms))
    return h, parts, (hs, hp, ks, kp, os_, es), bound


def preparations():
    chi_s = np.array([1., 1.]) / np.sqrt(2)
    chi_p = np.array([np.cos(np.pi / 8), np.sin(np.pi / 8)])
    p = np.array([1., 1., 2.]) * np.exp(2j * np.pi * np.arange(3) / 3) / np.sqrt(6)
    g = np.ones(8) / np.sqrt(8)
    js = []
    for r in (0, 1):
        columns = []
        for ds in range(2):
            for dp in range(2):
                sv = kron(np.eye(3)[r].reshape(-1, 1), chi_s.reshape(-1, 1), np.eye(2)[ds].reshape(-1, 1)).ravel()
                pv = kron(p.reshape(-1, 1), chi_p.reshape(-1, 1), np.eye(2)[dp].reshape(-1, 1)).ravel()
                columns.append(np.kron(np.kron(sv, pv), g))
        js.append(np.column_stack(columns))
    witness = np.array([1., 1., 0., 0.]) / np.sqrt(2)
    return js, witness, g


def scaled_taylor(h, vectors, bound, degree=18):
    n = int(math.ceil(T * bound / .4))
    dt = T / n
    state = vectors.copy()
    for _ in range(n):
        term = state.copy()
        out = state.copy()
        for j in range(1, degree + 1):
            term = (-1j * dt / j) * (h @ term)
            out += term
        state = out
    x = dt * bound
    rem = math.exp(x) * x ** (degree + 1) / math.factorial(degree + 1)
    universal = n * rem * math.exp(n * rem)
    return state, {"steps": n, "degree": degree, "sector_norm_upper_bound": bound,
                   "exact_arithmetic_all_input_operator_error_upper_bound": universal}


def spectral(h, vectors):
    w, q = np.linalg.eigh(h)
    out = q @ (np.exp(-1j * T * w)[:, None] * (q.conj().T @ vectors))
    residual = np.linalg.norm(h @ q - q * w[None, :], "fro")
    orth = np.linalg.norm(q.conj().T @ q - np.eye(len(w)), "fro")
    # Frobenius norms upper-bound operator norms. Polar correction of q yields
    # a unitary diagonalization with a transparent conservative Duhamel bound.
    d = orth / (1 + math.sqrt(1 - orth))
    delta_h = residual * math.sqrt(1 + orth) + np.linalg.norm(h, "fro") * orth + max(abs(w)) * d * (2 + d)
    bound = T * delta_h + d * (2 + d)
    return out, (w, q), {"eigen_residual_frobenius": float(residual),
                         "orthogonality_residual_frobenius": float(orth),
                         "conservative_all_input_propagation_diagnostic": float(bound),
                         "conservative_all_input_probability_diagnostic": float(2 * bound + bound * bound),
                         "spectrum_without_common_baseline": [float(w[0]), float(w[-1])],
                         "residuals_are_float_diagnostics_not_interval_certificates": True}


def expval(state, op):
    return float(np.vdot(state, op @ state).real)


def rho_probe(state):
    a = state.reshape(12, 12, 8)
    return np.einsum("spg,sqg->pq", a, a.conj())


def observe(state):
    a = state.reshape(12, 3, 2, 2, 8)
    rp = rho_probe(state)
    eclock = (I2 - Y) / 2
    ey = kron(np.diag([1., 0., 0.]), I2, eclock)
    site = [float(np.trace(rp @ kron(np.diag(np.eye(3)[i]), I2, I2)).real) for i in range(3)]
    rz = np.einsum("spg,sph->gh", state.reshape(12, 12, 8), state.reshape(12, 12, 8).conj())
    gz = [float(np.trace(rz @ z).real) for z in GZ]
    # Complete terminal projective record: S computational basis, P site/chirality,
    # P internal Y basis, and all three G Z bits. No failure outcome is discarded.
    _, ybasis = np.linalg.eigh(Y)
    amplitudes = np.einsum("da,sjcag->sjcdg", ybasis.conj().T, a)
    probabilities = np.abs(amplitudes) ** 2
    joint = probabilities.sum(axis=(0, 2))  # P site, internal Y outcome, full G bits.
    return {"clock_minus_y_at_site0": float(np.trace(rp @ ey).real),
            "clock_minus_y_probability": float(np.trace(rp @ kron(np.eye(3), I2, eclock)).real),
            "clock_minus_y_site_probabilities": joint[:, 0].sum(axis=1).tolist(),
            "probe_site_probabilities": site, "register_z_expectations": gz,
            "register_plus_z_probabilities": [(1 + v) / 2 for v in gz],
            "complete_record_outcomes": int(probabilities.size),
            "complete_record_sum": float(probabilities.sum()),
            "complete_record_minimum": float(probabilities.min()),
            "clock_site0_record_marginal": float(joint[0, 0].sum()),
            "site_clock_register_record": joint.tolist()}


def energy(state, h, parts):
    vals = {(name if name != "register_without_common_baseline" else "register_including_common_baseline"): expval(state, value) + (BASE if name == "register_without_common_baseline" else 0.)
            for name, value in parts.items()}
    hv = h @ state
    mean = float(np.vdot(state, hv).real)
    vals["total_including_common_baseline"] = mean + BASE
    vals["total_variance"] = float(np.vdot(hv, hv).real - mean * mean)
    return vals


def trace_distance(a, b):
    return float(.5 * np.abs(np.linalg.eigvalsh((a - b + (a - b).conj().T) / 2)).sum())


def frozen_entangled_reference_check(eigensystem, g):
    rng = np.random.default_rng(1174001)
    psi = rng.normal(size=(12, 12, 3)) + 1j * rng.normal(size=(12, 12, 3))
    psi /= np.linalg.norm(psi)
    u, _ = np.linalg.qr(rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12)))
    alternate = np.einsum("st,tpr->spr", u, psi)
    w, q = eigensystem
    out = []
    for v in (psi, alternate):
        x = np.einsum("spr,g->spgr", v, g).reshape(1152, 3)
        final = (q @ (np.exp(-1j * T * w)[:, None] * (q.conj().T @ x))).reshape(12, 12, 8, 3)
        out.append(np.einsum("spgr,sqgt->prqt", final, final.conj()).reshape(36, 36))
    return {"reference_dimension": 3, "source_probe_reference_initially_entangled": True,
            "nonselective_probe_reference_trace_distance_after_local_source_unitary": trace_distance(*out)}


def run():
    js, witness, g = preparations()
    vectors = np.column_stack(js)
    report = {"reviewer_implementation": "full 1152 matrix, full spectral evolution, independent scaled Taylor action",
              "author_script_not_used_or_imported": True,
              "parameters": {"T": T, "epsilon": EPS, "omega": OMEGA, "kappa": KAPPA, "baseline": BASE},
              "preparation_isometry_residuals": [float(np.linalg.norm(j.conj().T @ j - np.eye(4), "fro")) for j in js],
              "models": {}}
    for name, masses, omega in (("active_massive", (2., 2.4), OMEGA),
                                ("frozen_massive", (2., 2.4), 0.),
                                ("active_massless_probe", (0., 0.), OMEGA),
                                ("frozen_massless_probe", (0., 0.), 0.)):
        h, parts, mats, hbound = construct(masses, omega)
        evolved, eigensystem, diagnostics = spectral(h, vectors)
        taylor, tb = scaled_taylor(h, vectors, hbound)
        # Restore the same common baseline phase to every flavor and source branch.
        common_phase = np.exp(-1j * BASE * T)
        evolved *= common_phase
        taylor *= common_phase
        observations, energies = [], []
        for r in (0, 1):
            initial = js[r] @ witness
            final = evolved[:, r * 4:(r + 1) * 4] @ witness
            observations.append(observe(final))
            energies.append({"initial": energy(initial, h, parts), "final": energy(final, h, parts)})
        differences = {
            "clock_minus_y_at_site0_r1_minus_r0": observations[1]["clock_minus_y_at_site0"] - observations[0]["clock_minus_y_at_site0"],
            "site1_r1_minus_r0": observations[1]["probe_site_probabilities"][1] - observations[0]["probe_site_probabilities"][1],
            "register_plus_z_r1_minus_r0": [b - a for a, b in zip(observations[0]["register_plus_z_probabilities"], observations[1]["register_plus_z_probabilities"])],
        }
        compressed_mean_difference = js[1].conj().T @ h @ js[1] - js[0].conj().T @ h @ js[0]
        clock_effect = kron(I12, np.diag([1., 0., 0.]), I2, (I2 - Y) / 2, I8)
        compressed_effects = []
        for r in (0, 1):
            jfinal = evolved[:, r * 4:(r + 1) * 4]
            compressed_effects.append(jfinal.conj().T @ clock_effect @ jfinal)
        item = {"hermiticity_frobenius": float(np.linalg.norm(h - h.conj().T, "fro")),
                "spectral": diagnostics, "taylor": tb,
                "spectral_vs_taylor_preparation_matrix_frobenius": float(np.linalg.norm(evolved - taylor, "fro")),
                "evolved_preparation_isometry_residuals": [float(np.linalg.norm(evolved[:, r * 4:(r + 1) * 4].conj().T @ evolved[:, r * 4:(r + 1) * 4] - np.eye(4), "fro")) for r in (0, 1)],
                "compressed_total_mean_difference_frobenius": float(np.linalg.norm(compressed_mean_difference, "fro")),
                "compressed_unknown_4d_clock_site0_effects": [{"real": f.real.tolist(), "imag": f.imag.tolist()} for f in compressed_effects],
                "compressed_clock_effect_difference_operator_norm": float(max(abs(np.linalg.eigvalsh(compressed_effects[1] - compressed_effects[0])))),
                "observations": observations, "energy_accounts": energies, "differences": differences}
        if omega == 0:
            item["entangled_reference_no_signal"] = frozen_entangled_reference_check(eigensystem, g)
        if name == "active_massive":
            hs, hp, ks, kp, onsite, edges = mats
            shift = np.roll(np.eye(3), 1, axis=0)
            shift12 = kron(shift, I2, I2)
            report["free_source_translation_commutator_frobenius"] = float(np.linalg.norm(shift12 @ hs - hs @ shift12, "fro"))
            report["car_one_species_onsite_positive_spectral_sum"] = float(sum(x for x in np.linalg.eigvalsh(onsite[0]) if x > 0))
            report["car_one_species_edge_positive_spectral_sum"] = float(sum(x for x in np.linalg.eigvalsh(edges[0]) if x > 0))
            reflection = [1, 0, 2]
            perm = []
            for s in range(12):
                for p in range(12):
                    for gi in range(8):
                        sr = 4 * reflection[s // 4] + s % 4
                        pr = 4 * reflection[p // 4] + p % 4
                        gr = ((gi & 4) >> 1) | ((gi & 2) << 1) | (gi & 1)
                        perm.append((sr * 12 + pr) * 8 + gr)
            perm = np.array(perm)
            anti_phase = np.exp(-2j * np.pi / 3)
            report["special_preparation_antiunitary_symmetry"] = {
                "definition": "Simultaneously reflect sites 0 and 1 for S, P, and G, then complex conjugate.",
                "hamiltonian_invariance_frobenius": float(np.linalg.norm(h[np.ix_(perm, perm)].conj() - h, "fro")),
                "prepared_isometry_reflection_conjugation_residual": float(np.linalg.norm(js[0][perm].conj() - anti_phase * js[1], "fro")),
                "consequence": "The fixed r0/r1 witnesses have equal full total-energy spectral probabilities; antiunitarity reverses time, so their forward-time outputs need not agree.",
                "not_inferred_from_source_translation_alone": True,
            }
            report["initial_register_z_second_derivatives"] = []
            for r in (0, 1):
                v = js[r] @ witness
                hv, hhv = h @ v, h @ (h @ v)
                direct, formula = [], []
                for i in range(3):
                    zi = kron(I12, I12, GZ[i])
                    direct.append(float((-2 * np.vdot(v, zi @ hhv).real + 2 * np.vdot(hv, zi @ hv).real)))
                    hi = kron(ks[i], I12, I8) + kron(I12, kp[i], I8)
                    formula.append(-4 * OMEGA * EPS * expval(v, hi))
                report["initial_register_z_second_derivatives"].append({"double_commutator": direct, "local_energy_formula": formula,
                    "maximum_difference": float(max(abs(a - b) for a, b in zip(direct, formula)))})
        report["models"][name] = item
        print(name, "clock/site1 delta:", differences["clock_minus_y_at_site0_r1_minus_r0"], differences["site1_r1_minus_r0"], flush=True)
    b = 3 * (8.8 + 4 + .3 * 8.8 + .3 * 4 + .6 + .2)
    slices = math.ceil(T * T * b * b / (2e-4))
    report["full_parent_contract"] = {"fock_and_register_dimension": 2 ** 27, "invariant_sector_dimension": 1152,
        "eigh_inputs_cover_entire_sector": True, "local_term_count": 18, "norm_sum_bound": b,
        "hamiltonian_positive_lower_bound": BASE - b,
        "product_formula_slices": slices, "local_gate_count": 18 * slices,
        "product_formula_operator_error_upper_bound": T * T * b * b / (2 * slices),
        "task_probability_error_budget": ETA, "two_run_probability_difference_penalty": 2 * ETA,
        "max_total_hamiltonian_norm_bound": BASE + b}
    act = report["models"]["active_massive"]["differences"]
    massless = report["models"]["active_massless_probe"]["differences"]
    report["robust_lower_bounds"] = {
        "original_moving_clock_site0_difference": abs(act["clock_minus_y_at_site0_r1_minus_r0"]) - 2 * ETA,
        "massless_probe_site1_difference": abs(massless["site1_r1_minus_r0"]) - 2 * ETA,
        "register_population_differences": [abs(x) - 2 * ETA for x in act["register_plus_z_r1_minus_r0"]]}
    for name, item in report["models"].items():
        assert item["spectral"]["conservative_all_input_propagation_diagnostic"] < 1e-9
        assert item["spectral"]["conservative_all_input_probability_diagnostic"] < 1e-9
        assert item["spectral_vs_taylor_preparation_matrix_frobenius"] < 1e-10
        assert item["compressed_total_mean_difference_frobenius"] < 1e-12
        for obs, en in zip(item["observations"], item["energy_accounts"]):
            assert abs(obs["complete_record_sum"] - 1) < 1e-12
            assert abs(obs["clock_minus_y_at_site0"] - obs["clock_site0_record_marginal"]) < 1e-12
            assert abs(en["initial"]["total_including_common_baseline"] - en["final"]["total_including_common_baseline"]) < 1e-11
        if name.startswith("frozen"):
            assert item["entangled_reference_no_signal"]["nonselective_probe_reference_trace_distance_after_local_source_unitary"] < 1e-11
            assert abs(item["differences"]["clock_minus_y_at_site0_r1_minus_r0"]) < 1e-11
            assert abs(item["differences"]["site1_r1_minus_r0"]) < 1e-11
            assert item["compressed_clock_effect_difference_operator_norm"] < 1e-11
    assert all(x["maximum_difference"] < 1e-11 for x in report["initial_register_z_second_derivatives"])
    assert report["special_preparation_antiunitary_symmetry"]["hamiltonian_invariance_frobenius"] < 1e-12
    assert report["special_preparation_antiunitary_symmetry"]["prepared_isometry_reflection_conjugation_residual"] < 1e-12
    assert report["robust_lower_bounds"]["original_moving_clock_site0_difference"] > 0
    assert report["robust_lower_bounds"]["massless_probe_site1_difference"] > 0
    report["checks_pass"] = True
    return report


def compare_author(report):
    """Read saved author results only after completing the independent calculation."""
    path = HERE / "results.json"
    if not path.exists():
        return {"status": "author results not present"}
    saved = json.loads(path.read_text(encoding="utf-8"))["groups"]
    all_deltas = []
    g = saved["same_clock_position_and_active_register"]
    for label, model in (("active", "active_massive"), ("frozen", "frozen_massive")):
        own = report["models"][model]
        for r in (0, 1):
            obs = own["observations"][r]
            for key, author_key in (("clock_minus_y_site_probabilities", "clock_minusY_site_probabilities"),
                                    ("probe_site_probabilities", "ordinary_site_probabilities"),
                                    ("register_plus_z_probabilities", "G_plusZ_population_probabilities")):
                all_deltas.extend(abs(np.asarray(obs[key]) - g[author_key][label][str(r)]).ravel().tolist())
            author_effect = g["compressed_unknown_4D_clock_site0_effects"][label]["source" + str(r)]
            own_effect = own["compressed_unknown_4d_clock_site0_effects"][r]
            for component in ("real", "imag"):
                all_deltas.extend(abs(np.asarray(own_effect[component]) - author_effect[component]).ravel().tolist())
    for label, model in (("massless", "active_massless_probe"), ("massless_frozen", "frozen_massless_probe")):
        for r in (0, 1):
            ours = report["models"][model]["observations"][r]["probe_site_probabilities"]
            theirs = saved["massless_same_family"]["site_probabilities"][label][str(r)]
            all_deltas.extend(abs(np.asarray(ours) - theirs).tolist())
    maximum = max(all_deltas)
    assert maximum < 1e-11, maximum
    return {"status": "PASS", "real_components_compared": len(all_deltas),
            "maximum_absolute_difference": maximum,
            "author_results_used_only_after_independent_calculation": True}


def compare_frozen(new, old, path="review"):
    if isinstance(new, dict):
        assert new.keys() == old.keys(), path
        for key in new:
            compare_frozen(new[key], old[key], path + "." + key)
    elif isinstance(new, list):
        assert len(new) == len(old), path
        for index, (a, b) in enumerate(zip(new, old)):
            compare_frozen(a, b, path + "." + str(index))
    elif isinstance(new, float):
        assert math.isclose(new, old, rel_tol=2e-8, abs_tol=2e-12), (path, new, old)
    else:
        assert new == old, (path, new, old)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Create review_checks.json after all six author cores exist; never overwrite.")
    args = parser.parse_args()
    if args.write:
        for rel in CORE:
            if not (HERE / rel).is_file():
                raise SystemExit("Author cores not complete; no report written: " + rel)
        if (HERE / "review_checks.json").exists():
            raise SystemExit("Refusing to overwrite frozen review_checks.json")
    hashes_before = {rel: hashlib.sha256((HERE / rel).read_bytes()).hexdigest() for rel in CORE if (HERE / rel).exists()}
    report = run()
    report["comparison_with_author_saved_results"] = compare_author(report)
    hashes_after = {rel: hashlib.sha256((HERE / rel).read_bytes()).hexdigest() for rel in CORE if (HERE / rel).exists()}
    report["author_six_core_sha256"] = hashes_before
    report["author_cores_unchanged_during_review"] = hashes_before == hashes_after
    report["reviewer_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.write:
        if len(hashes_before) != 6 or hashes_before != hashes_after:
            raise SystemExit("Author files changed during review; no report written")
        with (HERE / "review_checks.json").open("x", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
            f.write("\n")
        print("Created review_checks.json")
    else:
        saved_path = HERE / "review_checks.json"
        if saved_path.exists():
            compare_frozen(report, json.loads(saved_path.read_text(encoding="utf-8")))
            status = "READ_ONLY_MATCH"
        else:
            status = "UNSAVED_REVIEW"
        print(json.dumps({"status": status, "checks_pass": report["checks_pass"],
                          "author_match": report["comparison_with_author_saved_results"],
                          "robust_lower_bounds": report["robust_lower_bounds"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
