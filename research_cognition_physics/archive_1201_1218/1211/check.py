"""1211: finite native cooling resources; default replay is read-only.

Exact fractions/permutations certify native sector capacities and routing.
NumPy diagnoses general spectral capacities, the complete unknown-reference
SWAP, and standard Nernst examples. No old code is imported or executed.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fstr(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def entropy(a):
    ev = np.linalg.eigvalsh((a + a.conj().T) / 2)
    assert ev.min() > -1e-11
    return float(-sum(x * math.log(x) for x in ev if x > 1e-15))


def words(n):
    return list(product(range(4), repeat=n))


def counts(s):
    return tuple(s.count(i) for i in range(4))


def native_capacities(distributions):
    """Exact diagonal product populations; all joint energy-block unitaries."""
    n = len(distributions)
    ws = words(n)
    sectors = defaultdict(list)
    pops = []
    no_ground = F(0)
    for s in ws:
        p = math.prod(distributions[j][s[j]] for j in range(n))
        pops.append(p)
        sectors[counts(s)].append((s, p))
        if 0 not in s:
            no_ground += p
    selected = F(0)
    for block in sectors.values():
        cap = sum(s[0] == 0 for s, _ in block)
        selected += sum(sorted((p for _, p in block), reverse=True)[:cap], F(0))
    unrestricted = sum(sorted(pops, reverse=True)[:4 ** (n - 1)], F(0))
    routed = 1 - no_ground
    assert F(0) <= selected <= routed <= 1
    assert selected <= unrestricted <= 1
    assert no_ground == math.prod(1 - d[0] for d in distributions)
    # A blank n+1-state zero-H record is a new purity resource.
    unconstrained_record_cap = min(len(pops), (n + 1) * 4 ** (n - 1))
    all_with_record = sum(sorted(pops, reverse=True)[:unconstrained_record_cap], F(0))
    assert routed <= all_with_record
    return {
        "n_material": n,
        "energy_sectors": len(sectors),
        "ground_capacity_without_record_exact": fstr(selected),
        "unrestricted_unitary_without_record_exact": fstr(unrestricted),
        "no_ground_weight_exact": fstr(no_ground),
        "ground_capacity_with_pure_record_exact": fstr(routed),
        "unrestricted_with_record_exact": fstr(all_with_record),
        "record_dimension": n + 1,
        "values_for_reading": {
            "without_record": float(selected),
            "unrestricted_without_record": float(unrestricted),
            "with_record": float(routed),
            "unavoidable_failure": float(no_ground),
        },
    }


def routing_permutation(n):
    """Blank label 0; label j+1 records first ground site; 0 means no ground.

    Pair each initial (word,blank) with (swapped_word,first+1), leaving all
    other vectors unchanged. The destination never has blank label.
    """
    ws = words(n)
    index = {s: i for i, s in enumerate(ws)}
    dk = n + 1
    perm = list(range(len(ws) * dk))
    used = set()
    pairs = 0
    for i, s in enumerate(ws):
        if 0 not in s:
            continue
        j = s.index(0)
        t = list(s)
        t[0], t[j] = t[j], t[0]
        a, b = i * dk, index[tuple(t)] * dk + j + 1
        assert a != b and a not in used and b not in used
        used.update((a, b))
        perm[a], perm[b] = b, a
        pairs += 1
    assert sorted(perm) == list(range(len(perm)))
    assert all(perm[perm[i]] == i for i in range(len(perm)))
    assert all(counts(ws[i // dk]) == counts(ws[perm[i] // dk]) for i in range(len(perm)))
    # These two full input-basis identities prove the effect for every Omega/R.
    for i, s in enumerate(ws):
        t = ws[perm[i * dk] // dk]
        label = perm[i * dk] % dk
        assert (t[0] == 0) == (0 in s)
        assert label == (s.index(0) + 1 if 0 in s else 0)
        if 0 not in s:
            assert t == s
    return ws, dk, perm, pairs


def spectral_capacity_diagnostic():
    # Includes coherence between different total-energy sectors.
    rng = np.random.default_rng(1211)
    z = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    rho = z @ z.conj().T
    rho /= np.trace(rho)
    energies = np.array([0, 1, 1, 2, 2, 3])
    target = np.diag([1, 1, 0, 0, 0, 0])
    u = np.zeros((6, 6), complex)
    bound = 0.0
    for en in range(4):
        ix = np.flatnonzero(energies == en)
        block = rho[np.ix_(ix, ix)]
        ev, vec = np.linalg.eigh(block)
        order = np.argsort(ev)[::-1]
        g = [j for j in range(len(ix)) if target[ix[j], ix[j]] == 1]
        rest = [j for j in range(len(ix)) if j not in g]
        out_basis = np.eye(len(ix))[:, g + rest]
        u[np.ix_(ix, ix)] = out_basis @ vec[:, order].conj().T
        bound += float(sum(ev[order[:len(g)]]))
    achieved = float(np.trace(target @ u @ rho @ u.conj().T).real)
    err = float(np.max(np.abs(u.conj().T @ u - np.eye(6))))
    comm = float(np.max(np.abs(u @ np.diag(energies) - np.diag(energies) @ u)))
    all_bound = float(np.linalg.eigvalsh(rho)[-2:].sum())
    assert abs(bound - achieved) < 1e-12 and err < 1e-12 and comm < 1e-12
    assert bound <= all_bound + 1e-12 < 1
    return {"block_optimum": bound, "achieved": achieved,
            "unrestricted_optimum": all_bound, "unitarity_error": err,
            "energy_commutator_error": comm,
            "scope": "NumPy diagnostic; proof gives all finite states and sectors"}


def routing_diagnostics(native_e):
    rows = []
    for n in (1, 2, 3):
        ws, dk, perm, pairs = routing_permutation(n)
        d = len(perm)
        w = np.eye(d, dtype=int)[np.array(perm), :]
        # Construction is an involution; integer checks do not depend on eigenvalues.
        assert np.array_equal(w @ w, np.eye(d, dtype=int))
        v = w[:, np.arange(len(ws)) * dk]
        pg = np.array([int(s[0] == 0) for s in ws for _ in range(dk)])
        effect = v.T @ (pg[:, None] * v)
        desired = np.diag([int(0 in s) for s in ws])
        assert np.array_equal(effect, desired)
        if n == 2:
            # Fixed bounded positive parent, not the old P^2 clock.
            qdiag = np.repeat([sum(native_e[j] - native_e[0] for j in s) for s in ws], dk)
            g = math.pi / 2 * (np.eye(d) - w)
            h = np.diag(qdiag) + g
            actual = np.exp(-1j * qdiag)[:, None] * w
            assert np.max(np.abs(h @ actual - actual @ h)) < 1e-11
            assert np.linalg.eigvalsh(h).min() > -1e-11
            av = actual[:, np.arange(len(ws)) * dk]
            assert np.max(np.abs(av.conj().T @ (pg[:, None] * av) - desired)) < 1e-12
        rows.append({"n": n, "record_dimension": dk, "full_dimension": d,
                     "disjoint_transpositions": pairs,
                     "count_conservation_exact": True,
                     "complete_ground_effect_identity_exact": True,
                     "unknown_joint_state_and_passive_R_supported": True})
    return rows


def swap_unknown_reference(native_e):
    rng = np.random.default_rng(4211211)
    psi = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
    psi /= np.linalg.norm(psi)
    initial = np.zeros((4, 4, 3), complex)
    initial[:, 0, :] = psi
    after = initial.transpose(1, 0, 2).copy()
    phases = np.exp(-1j * (native_e - native_e[0]))
    after *= phases[:, None, None] * phases[None, :, None]
    expected = np.zeros_like(after)
    expected[0, :, :] = phases[:, None] * psi
    diff = float(np.linalg.norm(after - expected))
    assert diff < 1e-12
    rho_s = psi @ psi.conj().T
    ea_in = float(np.trace(np.diag(native_e) @ rho_s).real)
    stock_gain = ea_in - native_e[0]
    assert stock_gain >= -1e-12
    return {"complete_SR_isometry_error": diff,
            "target_ground_probability": float(np.linalg.norm(after[0]) ** 2),
            "stock_energy_gain": stock_gain,
            "target_energy_loss": -stock_gain,
            "stock_entropy_out_nats": entropy(rho_s),
            "source_information_transferred_not_erased": True,
            "stock_returned": False,
            "old_autonomous_clock_exactness_claimed": False}


def fixed_native_and_deletions(cold, hot):
    rows = {}
    for n in (1, 2, 3, 4):
        rows[f"independent_cold_{n}"] = native_capacities([cold] * n)
        assert F(rows[f"independent_cold_{n}"]["ground_capacity_without_record_exact"]) == cold[0]
    for m in (1, 2):
        rows[f"hot_plus_{m}_cold"] = native_capacities([hot] + [cold] * m)
    ground = [F(1), F(0), F(0), F(0)]
    excited = [F(0), F(1), F(0), F(0)]
    rows["hot_plus_pure_ground"] = native_capacities([hot, ground])
    rows["hot_plus_pure_excited"] = native_capacities([hot, excited])
    assert F(rows["hot_plus_pure_ground"]["ground_capacity_without_record_exact"]) == 1
    assert F(rows["hot_plus_pure_excited"]["ground_capacity_without_record_exact"]) == hot[0]
    # Six one-ground words: both marginals full rank, joint not full rank.
    corr = {s: F(1, 6) for j in (1, 2, 3) for s in ((0, j), (j, 0))}
    marginal = [sum((p for s, p in corr.items() if s[0] == i), F(0)) for i in range(4)]
    assert marginal == [F(1, 2), F(1, 6), F(1, 6), F(1, 6)]
    assert all(p > 0 for p in marginal)
    assert sum((p for s, p in corr.items() if 0 not in s), F(0)) == 0
    # One state from each two-dimensional sector fits without a record.
    corr_without_record = 3 * F(1, 6)
    corr_with_record = F(1)
    correlated_case = {
        "diagonal_words": {"".join(map(str, s)): fstr(p) for s, p in corr.items()},
        "both_marginals_exact": [fstr(p) for p in marginal],
        "without_record_ground_capacity_exact": fstr(corr_without_record),
        "with_record_ground_capacity_exact": fstr(corr_with_record),
        "joint_full_rank": False,
        "product_of_same_marginals_routed_probability_exact": fstr(1 - (1 - marginal[0]) ** 2),
    }
    # Pure zero-energy blank alone has enough spectrum capacity, but no energy capacity.
    blank_only = {"target_dimension": 4, "pure_record_dimension": 4,
                  "record_hamiltonian": "0",
                  "unrestricted_unitary_ground_capacity_exact": "1",
                  "exact_energy_preserving_ground_capacity_exact": fstr(hot[0]),
                  "nondegenerate_native_h_used": True}
    thresholds = []
    for exponent in (2, 6, 12):
        tol = F(1, 10 ** exponent)
        n = 1
        while (1 - cold[0]) ** n > tol:
            n += 1
        assert (1 - cold[0]) ** (n - 1) > tol
        thresholds.append({"failure_target": fstr(tol), "minimal_n_in_this_routing_family": n,
                           "record_dimension": n + 1,
                           "failure_exact": fstr((1 - cold[0]) ** n),
                           "initial_independent_pools": n,
                           "time_law_inferred": False})
    return rows, correlated_case, blank_only, thresholds


def gibbs_entropy(energies, beta):
    e = np.array(energies, float)
    e -= e.min()
    raw = np.exp(-beta * e)
    ground_rank = int(np.count_nonzero(e == 0))
    excitation_weight = float(raw[e > 0].sum())
    z = ground_rank + excitation_weight
    # Avoid subtracting nearly equal numbers or losing a tiny ground entropy.
    logz = math.log(ground_rank) + math.log1p(excitation_weight / ground_rank)
    s = logz + beta * float(np.dot(e, raw)) / z
    return s, excitation_weight / z


def nernst_diagnostics(native_e, gap_lower):
    rows = []
    for beta in (100, 300, 1000):
        s, q = gibbs_entropy(native_e, beta)
        r = 3 * math.exp(-beta * float(gap_lower))
        assert r < .5
        # h_2(r) <= r(1-log r); d-g=3, g=1.
        bound = r * (1 - math.log(r) + math.log(3))
        assert s <= bound + 1e-13 and q <= r + 1e-13
        rows.append({"beta": beta, "entropy": s, "excited_probability": q,
                     "uniform_entropy_bound_formula_value": bound})
    closing = []
    expected, _ = gibbs_entropy([0, 1], 1)
    for beta in (1, 10, 100, 1000):
        s, _ = gibbs_entropy([0, 1 / beta], beta)
        assert abs(s - expected) < 1e-13
        closing.append({"beta": beta, "gap": 1 / beta, "entropy": s})
    degenerate, _ = gibbs_entropy([0, 0, 1], 100)
    assert abs(degenerate - math.log(2)) < 1e-13
    return {"native_ground_rank": 1, "native_gap_strict_lower_exact": fstr(gap_lower),
            "fixed_native_spectrum": rows,
            "fixed_ground_rank_without_uniform_gap": closing,
            "degenerate_ground_entropy_limit": "log(2)",
            "degenerate_beta100_entropy": degenerate,
            "numerical_values_are_diagnostics": True,
            "no_dynamic_unattainability_inferred_from_nernst": True}


def actual_window_pool_sector_counts(p2, p4, cold, hot):
    """Only count N0=0 in the two old windows; do not redo full enumeration."""
    den = p2["preparation"]["root_denominator"]
    box = p2["preparation"]["root_boxes"]
    factorials = [math.factorial(i) for i in range(101)]
    rows = []
    cases = [("cold", F(39, 2), F(41, 2), int(p2["preparation"]["dimension"]), cold),
             ("hot", F(79, 2), F(81, 2), int(p4["new_reference"]["dimension"]), hot)]
    for name, low, high, dimension, marginal in cases:
        total, selected, examined, ambiguous = 0, 0, 0, 0
        for n1 in range(101):
            for n2 in range(101 - n1):
                ns = (0, n1, n2, 100 - n1 - n2)
                examined += 1
                lo = sum(n * interval[0] for n, interval in zip(ns, box))
                hi = sum(n * interval[1] for n, interval in zip(ns, box))
                if lo >= low * den and hi <= high * den:
                    selected += 1
                    total += factorials[100] // math.prod(factorials[n] for n in ns)
                elif not (hi < low * den or lo > high * den):
                    ambiguous += 1
        assert examined == 5151 and ambiguous == 0 and 0 < total < dimension
        failure = F(total, dimension)
        false_iid = (1 - marginal[0]) ** 100
        assert failure != false_iid
        assert failure < F(4, 10 ** 15)
        rows.append({"reference": name, "N": 100, "energy_window": [fstr(low), fstr(high)],
                     "only_N0_zero_types_examined": examined,
                     "selected_types": selected, "ambiguous_types": ambiguous,
                     "no_ground_word_count": str(total),
                     "full_window_dimension_reused": str(dimension),
                     "ideal_correlated_pool_failure_exact": fstr(failure),
                     "ideal_correlated_pool_failure": float(failure),
                     "incorrect_iid_substitution_failure": float(false_iid),
                     "pure_record_dimension_for_full_route": 101,
                     "below_old_1205_4e_minus15_preparation_transport_budget": True,
                     "positive_failure_for_actual_clock_prepared_state_claimed": False,
                     "scope": "ideal selected window pool; source labels/passive R retained"})
    return rows


def run():
    p2 = json.loads((STAGE / "1202/results.json").read_text(encoding="utf-8"))
    p4 = json.loads((STAGE / "1204/results.json").read_text(encoding="utf-8"))
    p6 = json.loads((STAGE / "1206/results.json").read_text(encoding="utf-8"))
    cold = list(map(F, p2["one_marginal_rational"]))
    hot = list(map(F, p4["new_reference"]["one_marginal_rational"]))
    assert sum(cold) == sum(hot) == 1 and min(cold + hot) > 0
    cert = p6["algebraic_certificate"]
    assert cert["galois_group"] == "S3" and cert["number_energy_2x2_determinant"] == "81/40"
    den = p2["preparation"]["root_denominator"]
    boxes = [[F(v, den) for v in pair] for pair in p2["preparation"]["root_boxes"]]
    gap_lower = boxes[1][0] - boxes[0][1]
    assert gap_lower > F(44, 1000)
    native_e = np.array([float(sum(box) / 2) for box in boxes])
    capacities, corr, blank, thresholds = fixed_native_and_deletions(cold, hot)
    deps = [
        STAGE / "1202/results.json", STAGE / "1204/results.json",
        STAGE / "1205/proof.md", STAGE / "1205/results.json",
        STAGE / "1206/proof.md", STAGE / "1206/results.json",
        STAGE / "1209/proof.md", STAGE / "1209/results.json",
        ROOT / "archive_370_428/research_note_421.md",
        ROOT / "archive_429_466/research_note_457.md",
    ]
    return {
        "schema": "native_cooling_resources_1211_v1",
        "status": "AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE",
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in deps},
        "scope": {
            "same_native_h_finitely_many_materials": True,
            "additional_record_hamiltonian_zero": True,
            "blank_record_initially_pure_independent": True,
            "exact_conservation_is_bare_sum_h": True,
            "arbitrary_joint_material_state_and_passive_R_in_theorem": True,
            "general_energy_block_class_is_mathematical_upper_menu": True,
            "explicit_route_not_claimed_as_old_QQ_or_1205_exact_clock": True,
            "no_general_third_law_refutation": True,
        },
        "spectral_block_reuse_diagnostic": spectral_capacity_diagnostic(),
        "native_capacities_exact": capacities,
        "complete_routing_checks": routing_diagnostics(native_e),
        "pure_ground_swap_unknown_SR": swap_unknown_reference(native_e),
        "correlations_delete_marginal_full_rank_argument": corr,
        "pure_zero_energy_blank_alone": blank,
        "finite_resource_growth_in_one_explicit_family": thresholds,
        "old_correlated_window_pools": actual_window_pool_sector_counts(p2, p4, cold, hot),
        "nernst_scope_diagnostics": nernst_diagnostics(native_e, gap_lower),
        "check_groups": 8,
    }


def compare(a, b, path="root"):
    if isinstance(a, float) and isinstance(b, (float, int)):
        assert math.isclose(a, b, rel_tol=1e-10, abs_tol=2e-12), (path, a, b)
    elif isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            compare(a[k], b[k], path + "." + k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path + f"[{i}]")
    else:
        assert a == b, (path, a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true", help="read-only default")
    args = parser.parse_args()
    computed = run()
    target = HERE / "results.json"
    if args.write:
        target.write_text(json.dumps(computed, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": "PASS", "all_checks_passed": True,
                          "results_match": True, "check_groups": 8,
                          "mode": "author_write", "independent_acceptance": "pending"}))
    else:
        compare(computed, json.loads(target.read_text(encoding="utf-8")))
        print(json.dumps({"status": "PASS", "all_checks_passed": True,
                          "results_match": True, "check_groups": 8,
                          "mode": "read_only_replay"}))
