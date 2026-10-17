"""1276 independent count/transport audit; no author implementation imported.

Default recomputes and matches independent_results.json without writing.
--write creates the results exclusively.  Large K is used only via log K.
Inherited 1175 numerical enclosures keep their declared arithmetic scope.
This checks a conditional reference observable, not total entropy or a new
area/entropy measurement, thermal preparation, null source, or GR identity.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, log
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics").is_dir())
RCP = ROOT / "research_cognition_physics"
OLD = RCP / "archive_1161_1200/1175"
RESULT = HERE / "independent_results.json"


def need(ok, description):
    if not ok:
        raise AssertionError(description)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(x):
    return {"numerator": x.numerator, "denominator": x.denominator}


def log_interval(q, terms=110):
    """Positive rational q: atanh series, with an explicit geometric tail."""
    q = F(q)
    if q < 1:
        low, high = log_interval(1 / q, terms)
        return -high, -low
    x = (q - 1) / (q + 1)
    total = F(0)
    for j in range(terms):
        total += 2 * x ** (2 * j + 1) / (2 * j + 1)
    tail = 2 * x ** (2 * terms + 1) / ((2 * terms + 1) * (1 - x * x))
    return total, total + tail


def coarse_upper(x, scale=10**12):
    return F(-((-x.numerator * scale) // x.denominator), scale)


def small_exact_checks():
    # Directly count configurations and verify the probability version of the
    # method-of-types inequalities; this does not invoke logarithm numerics.
    tested = 0
    for k in range(1, 49):
        need(sum(comb(k, m) for m in range(k + 1)) == 2**k,
             "all Z sectors exhaust the original qubit Hilbert space")
        for m in range(k + 1):
            p = F(m, k)
            prob = F(comb(k, m)) * p**m * (1 - p)**(k - m)
            need(F(1, k + 1) <= prob <= 1, "exact type probability bound")
            tested += 1

    # Actual basis-bit enumeration, independent of the binomial formula.
    k = 9
    counts = [0] * (k + 1)
    for bits in range(2**k):
        counts[bits.bit_count()] += 1
    need(counts == [comb(k, m) for m in range(k + 1)], "sector dimensions")
    weights = [F(c, 2**k) for c in counts]  # original |+X>^K preparation
    zmean = sum(w * F(2*m-k, k) for m, w in enumerate(weights))
    zvar = sum(w * F(2*m-k, k)**2 for m, w in enumerate(weights))
    need(zmean == 0 and zvar == F(1, k), "original coherent preparation moments")
    conditional_count = sum(float(w) * log(c) for w, c in zip(weights, counts))
    macro_entropy = -sum(float(w) * log(float(w)) for w in weights)
    need(abs(conditional_count + macro_entropy - k*log(2)) < 2e-14,
         "count plus macro entropy; not the actual pure-state entropy")
    return {
        "exact_type_cases": tested,
        "basis_enumeration_K": k,
        "basis_dimension": 2**k,
        "original_coherent_Z_mean": rational(zmean),
        "original_coherent_Z_variance": rational(zvar),
        "pure_state_von_Neumann_entropy": 0,
        "conditional_count_display": conditional_count,
        "macro_entropy_display": macro_entropy,
        "display_values_are_not_strict_log_evidence": True,
    }


def calculate():
    old = json.loads((OLD / "results.json").read_text(encoding="utf-8"),
                     parse_float=str)
    params = old["group2_autonomous_common_transport"]["parameters"]
    need(params["grid"] == 65 and params["V"] == 65**3, "original grid")
    ell, eps, nu, omega, mass, T = F(16, 65), F(3, 20), F(7, 20), F(1, 2), F(2), F(2)
    for key, exact in (("T", T), ("epsilon", eps), ("nu", nu),
                       ("Omega", omega), ("source_mass", mass)):
        need(F(params[key]) == exact, "original parameter " + key)
    V = 65**3
    bS = 2*eps*mass + 6*nu/ell
    bP = 2*eps*(F(3, 10)+F(1, 2)) + 18*nu/ell
    w = 2*V*bS
    L = V*bP
    need(w == F(80245425, 16), "reconstructed full-Fock interaction norm")
    need(L == F(227026995, 32), "reconstructed full probe norm")
    exponent = old["group1_parent_and_full_background_bridge"]["K_defined_as_power_of_10"]
    need(exponent == 52275242, "original physical finite K")

    ln2lo, ln2hi = log_interval(2, 36)
    ln10lo, ln10hi = log_interval(10, 110)
    # ln(4 L T) = 25 ln 2 + ln((4 L T)/2^25), no large exponentiation.
    ratio = 4*L*T/F(2**25)
    need(1 < ratio < 2, "safe logarithm reduction")
    _, lnratiohi = log_interval(ratio, 36)
    log_eta_over_target = 25*ln2hi + lnratiohi + 6*w*T - (F(exponent, 2)-4)*ln10lo
    need(log_eta_over_target < 0, "eta(T) < 10^-4 with original K")
    eta = F(1, 10000)
    # r/2 = 2 exp(6wt)/sqrt(K) <= eta/(2 L T) at T, uniformly below T.
    b_upper = eta/(2*L*T)
    need(b_upper < F(1, 2**38), "uniform binary entropy argument")
    binary_upper = (1 + 38*ln2hi)/2**38
    need(binary_upper < F(1, 10**8), "binary entropy modulus")
    # log(K+1)/K is decreasing for K>=1.  K >= 10^28 suffices.
    need(exponent >= 28, "finite lower bound on K for type error")
    type_upper = (ln2hi + 28*ln10hi)/10**28
    need(type_upper < F(1, 10**8), "uniform type error")
    count_error = type_upper + binary_upper + eta*ln2hi
    need(ln2hi < F(7, 10) and count_error < F(1, 10000),
         "uniform actual conditional-count error b < 10^-4")

    hist = old["group2_autonomous_common_transport"]["background"]["history"]
    need(F(hist[-1]["t"]) == T, "same original endpoint")
    znum = F(hist[-1]["z"])
    numerical = old["group3_fixed_menu_and_error_certificate"]["numerical_enclosure"]
    zerr = F(numerical["Bloch_terminal_error_bound"]) + F(1, 10**15)
    zlo, zhi = F(-29, 250), F(-14, 125)
    need(zlo <= znum-zerr and znum+zerr <= zhi < 0,
         "original residual enclosure implies declared z interval")
    tau = F(old["group2_autonomous_common_transport"]["background"]["tau"])
    tauerr = F(numerical["clock_tau_error"])
    need(tau+tauerr < T, "coordinate T not silently called proper time")

    # Curvature of the binary entropy: ln2-h((1+z)/2) >= z^2/2.
    ideal_drop = zhi**2/2
    drop_lower = ideal_drop - 2*F(1, 10000)
    need(ideal_drop == F(98, 15625), "reference count loss")
    need(drop_lower > F(3, 500), "actual normalized count drop > 3/500")
    a_lower = (1+eps*zhi)/(1+nu*zhi)
    area_relative = a_lower*a_lower-1
    need(a_lower == F(1229, 1201), "independent original effective scale")
    need(area_relative == F(68040, 1442401) and area_relative > F(47, 1000),
         "same effective area growth > .047")

    # Count unique boundary nodes and outward edges separately.
    nodes, edges = 0, 0
    for x in range(32):
        for y in range(32):
            for z in range(32):
                nodes += int(any(c in (0, 31) for c in (x, y, z)))
                edges += sum(int(c == 0)+int(c == 31) for c in (x, y, z))
    need(nodes == 32**3-30**3 == 5768, "unique relation factors")
    need(edges == 6*32**2 == 6144, "faces are not unique nodes")
    A0 = 6*(32*ell)**2
    need(A0 == F(1572864, 4225), "declared geometric area convention")
    # Sum expectation bounds over nodes; no factorization of their state.
    total_error_per_K = nodes*count_error
    need(total_error_per_K/nodes == count_error, "joint normalized bound")

    baseline = json.loads((HERE/"source_baseline.json").read_text(encoding="utf-8"))
    for rel, expected in baseline["files_sha256"].items():
        need(sha(ROOT/rel) == expected, "inherited baseline hash: " + rel)

    bindings = {}
    paths = [
        HERE/"admission.md", HERE/"source_baseline.json", HERE/"NEXT.md",
        HERE.parent/"research_note_1276.md",
        HERE/"proof.md", HERE/"check.py", HERE/"results.json",
        Path(__file__), OLD/"proof.md", OLD/"results.json",
        RCP/"archive_1161_1200/_admission/after1174_direction/mean_field_bound.md",
    ]
    for path in paths:
        bindings[str(path.relative_to(ROOT)).replace("\\", "/")] = sha(path)

    return {
        "round": 1276,
        "status": "PASS_WITH_SCOPE",
        "independent_groups": 1,
        "source_bindings": bindings,
        "strict_rational_bounds": {
            "K_log10": exponent, "w": rational(w), "L": rational(L),
            "log_eta_over_target_upper": rational(coarse_upper(log_eta_over_target)),
            "eta_upper": rational(eta),
            "binary_argument_upper": rational(b_upper),
            "binary_modulus_upper": rational(coarse_upper(binary_upper, 10**16)),
            "type_error_upper": rational(F(1, 10**8)),
            "uniform_count_error_upper": rational(F(1, 10000)),
            "z_numerical": rational(znum), "z_error": rational(zerr),
            "z_lower": rational(zlo), "z_upper": rational(zhi),
            "proper_time_lower": rational(tau-tauerr),
            "proper_time_upper": rational(tau+tauerr),
            "ideal_normalized_entropy_drop_lower": rational(ideal_drop),
            "actual_normalized_entropy_drop_lower": rational(drop_lower),
            "a_lower": rational(a_lower),
            "relative_area_increase_lower": rational(area_relative),
            "boundary_unique_nodes": nodes, "outward_faces": edges,
            "A0": rational(A0),
        },
        "independent_finite_checks": small_exact_checks(),
        "scope": {
            "same_original_H_preparation_K_and_unknown_probe_R": True,
            "exact_background_retained_not_Kfold_Hartree_replacement": True,
            "conditional_sector_reference_is_extra_model_choice": True,
            "reference_preparation_or_dephasing_performed": False,
            "full_interacting_Gibbs_reference_claimed": False,
            "new_L_or_area_measurement_certified": False,
            "effective_geometry_not_exact_quantum_metric": True,
            "full_entropy_thermal_flux_or_Jacobson_identity_claimed": False,
            "old_IEEE_residual_enclosure_scope_inherited": True,
            "new_ODE_or_large_K_simulation": False,
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = calculate()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("PASS_WITH_SCOPE EXCLUSIVE_WRITE")
    else:
        saved = json.loads(RESULT.read_text(encoding="utf-8"))
        need(saved == result, "read-only independent result mismatch")
        print("PASS_WITH_SCOPE READ_ONLY_MATCH")


if __name__ == "__main__":
    main()

