"""Non-numbered qualification diagnostics; no continuous-field simulation."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import math
import numpy as np

BASE = Path(__file__).resolve().parent
OUT = BASE / "after1147_spectral_boost_results.json"
SOURCE = BASE.parent / "1147" / "results.json"

def encoded(q):
    return {"fraction": str(q), "float": float(q)}

def norm(x):
    return float(np.linalg.norm(x, 2))

def car_group():
    dim = 16
    def ann(k):
        a = np.zeros((dim, dim), complex)
        for n in range(dim):
            if (n >> k) & 1:
                a[n ^ (1 << k), n] = (-1) ** ((n & ((1 << k)-1)).bit_count())
        return a
    b = [ann(k) for k in range(4)]
    ident = np.eye(dim)
    W = ident[:, [1, 2, 4, 8]]
    sx = sum((b[k].conj().T @ b[k+1] + b[k+1].conj().T @ b[k])/2 for k in (0, 2))
    sy = sum((-1j*b[k].conj().T @ b[k+1] + 1j*b[k+1].conj().T @ b[k])/2 for k in (0, 2))
    sz = sum((b[k].conj().T @ b[k] - b[k+1].conj().T @ b[k+1])/2 for k in (0, 2))
    rows = []
    for phi in [0.0, .713, math.pi/2]:
        av = [(b[s] + np.exp(1j*phi)*b[s+2])/math.sqrt(2) for s in range(2)]
        n = [x.conj().T @ x for x in av]
        E = ident - (ident-n[0]) @ (ident-n[1])
        v = np.array([1, np.exp(-1j*phi)]) / math.sqrt(2)
        e = np.kron(np.outer(v, v.conj()), np.eye(2))
        # Include an entangled, passive three-dimensional reference.
        q = np.arange(1, 13, dtype=float) + 1j*np.arange(12, 0, -1, dtype=float)
        q = q / np.linalg.norm(q)
        wr = np.kron(W, np.eye(3))
        full = np.kron(E, np.eye(3))
        expected = wr @ np.kron(e, np.eye(3)) @ q
        vals = {
            "phi": phi,
            "idempotence": norm(E@E-E),
            "spin_commutators": [norm(E@s-s@E) for s in (sx, sy, sz)],
            "code_compression": norm(W.conj().T@E@W-e),
            "code_leakage": norm(E@W-W@e),
            "entangled_reference_branch_residual": float(np.linalg.norm(full@wr@q-expected)),
            "single_spin_read_Sx_commutator": norm(n[0]@sx-sx@n[0]),
            "naive_number_sum_norm": norm(n[0]+n[1])
        }
        assert vals["idempotence"] < 1e-14
        assert max(vals["spin_commutators"]) < 1e-14
        assert vals["code_compression"] < 1e-14
        assert vals["code_leakage"] < 1e-14
        assert vals["entangled_reference_branch_residual"] < 1e-14
        assert abs(vals["single_spin_read_Sx_commutator"]-.5) < 1e-14
        assert abs(vals["naive_number_sum_norm"]-2) < 1e-14
        rows.append(vals)
    return {"pass": True, "kind": "exact_positive_mode_CAR_algebra_diagnostic",
            "full_positive_subalgebra_dimension": 16, "passive_reference_dimension": 3,
            "rows": rows, "not_a_full_sea_or_continuous_field_simulation": True}

def budget_group():
    old = json.loads(SOURCE.read_text(encoding="utf-8"))
    old_bound = Q(old["groups"][0]["with_phase_and_instrument_upper"]["fraction"])
    aw = Q(41, 100_000_000)
    rho_coarse = Q(1, 1_000_000)
    tail = 2*rho_coarse
    spectral = tail+aw
    common_mass = aw
    kin = Q(1, 1_000_000)
    time = Q(1, 1_000_000)
    extra_vector = kin+time+spectral+common_mass
    position = Q(573, 100_000_000)
    ry = Q(405, 1_000_000_000)
    ell = Q(1, 10)
    sqrt2_upper = Q(99, 70)
    assert sqrt2_upper**2 > 2
    assert sqrt2_upper*ry/ell < position
    total = old_bound+extra_vector+position
    assert total < Q(252, 100_000)
    # sin(5*pi/8)=sqrt(2+sqrt(2))/2 < 0.924.
    root2_upper = Q(283, 200)
    outer_upper = Q(231, 125)
    assert root2_upper**2 > 2
    assert outer_upper**2 > 2+root2_upper
    gap_lower = (1-outer_upper/2)/2
    assert gap_lower == Q(19, 500)
    frequency_penalty = Q(3, 400_000)*Q(22, 7)
    two_sided_gap = gap_lower-frequency_penalty-2*total
    assert two_sided_gap > Q(32, 1000)
    return {"pass": True, "kind": "conditional_exact_rational_budget",
            "conditional_on_uniform_graph_residual_and_domains": True,
            "formal_round_increment": 0,
            "h": 1,
            "angular_and_position_weight": encoded(aw),
            "spectral_tail_coarse": encoded(tail),
            "spectral_boost_error": encoded(spectral),
            "common_mass_boost_error": encoded(common_mass),
            "mass_evolution_error": encoded(time),
            "external_dispersion_error": encoded(kin),
            "additional_full_vector_error": encoded(extra_vector),
            "additional_joint_position_instrument_error": encoded(position),
            "old_1147_with_phase_and_instrument_error": encoded(old_bound),
            "new_joint_complete_instrument_error": encoded(total),
            "wrong_lab_time_ideal_probability_gap_lower": encoded(gap_lower),
            "actual_clock_frequency_two_phase_penalty": encoded(frequency_penalty),
            "conservative_two_sided_gap_after_errors": encoded(two_sided_gap)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    data = {"qualification_only": True, "formal_round": None, "formal_check_increment": 0,
            "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "all_passed": True, "groups": [car_group(), budget_group()]}
    value = json.dumps(data, ensure_ascii=False, indent=2) + "\n"
    if args.write_results:
        if OUT.exists() and OUT.read_text(encoding="utf-8") != value:
            raise SystemExit("refuse overwrite different qualification results")
        OUT.write_text(value, encoding="utf-8")
        print("saved", str(OUT))
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert saved == data, "qualification replay differs"
        print("PASS: two non-formal qualification groups; formal increment 0")
if __name__ == "__main__":
    main()
