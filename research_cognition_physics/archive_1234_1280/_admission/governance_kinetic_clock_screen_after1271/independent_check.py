"""Independent exact finite-query checks. No author-code imports.
Default: read-only compare. --write: create results exclusively.
Finite algebra certificates do not certify physical ports or model errors.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
STATES = [(n, N) for n in range(2) for N in range(2)]

def poly_add(a, b):
    return tuple((a[i] if i < len(a) else F(0)) +
                 (b[i] if i < len(b) else F(0))
                 for i in range(max(len(a), len(b))))

def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)

def neg(a):
    return tuple(-x for x in a)

def exact_material(rd):
    # z=exp(J/kBT)>1 is kept symbolic. Divide all rates by a shared
    # positive thermal factor: high->low gives z/R, reverse gives 1/R.
    weights = [(F(1),) if n == N else (F(0), F(1))
               for n, N in STATES]
    residual = [(F(0),) for _ in STATES]
    count = 0
    rate_coeffs = []
    for i, (n, N) in enumerate(STATES):
        coeffs = []
        for label, resist in (("D", rd), ("L", F(580)), ("R", F(580))):
            dst = (n, 1-N) if label == "D" else (1-n, N)
            j = STATES.index(dst)
            coeff = 1/resist
            rate = (F(0), coeff) if n == N else (coeff,)
            flux = poly_mul(weights[i], rate)
            residual[i] = poly_add(residual[i], neg(flux))
            residual[j] = poly_add(residual[j], flux)
            coeffs.append(coeff)
            count += 1
        rate_coeffs.append(coeffs)
    assert count == 12
    assert all(all(x == 0 for x in p) for p in residual)
    d, l, r = rate_coeffs[0]
    q = d/(d+l+r)
    assert all(c[0]/sum(c, F(0)) == q for c in rate_coeffs)
    return q, sum(rate_coeffs[0], F(0)), residual

def evaluate():
    q0, c0, residual0 = exact_material(F(43))
    q1, c1, residual1 = exact_material(F(43, 2))
    assert q0 == F(290,333) and q1 == F(580,623)
    gap = q1-q0
    assert gap == F(12470,207459)
    rate_ratio = c1/c0
    assert rate_ratio == F(623,333) > 1
    # W on the full {D,S,no event} record: no conditioning/deletion.
    w = (1-q0, -q0, F(0))
    old_y_coefficient = w[0]*q0+w[1]*(1-q0)
    new_y_coefficient = w[0]*q1+w[1]*(1-q1)
    assert old_y_coefficient == 0
    assert new_y_coefficient == gap
    # These exact coefficient equalities certify every y in [0,1].
    # The following finite classifications are additional diagnostics.
    classified = 0
    for y in (F(0), F(1,100), F(1,3), F(1,2), F(999,1000), F(1)):
        for q in (q0,q1):
            p = (y*q, y*(1-q), 1-y)
            assert sum(p,F(0)) == 1 and min(p) >= 0
            val = sum((x*z for x,z in zip(w,p)),F(0))
            assert val == y*(q-q0)
            # Retain L and R labels in the underlying four-result record.
            p4 = (p[0],p[1]/2,p[1]/2,p[2])
            w4 = (w[0],w[1],w[1],w[2])
            assert sum((x*z for x,z in zip(w4,p4)),F(0)) == val
            classified += 1
    # A different positive Taylor truncation from the author program.
    lower_exp7 = F(0)
    term = F(1)
    for n in range(16):
        if n:
            term *= F(7,n)
        lower_exp7 += term
    assert lower_exp7 > 1000
    # Changed survival is exp[-7*c1/c0] < exp(-7) < 1/1000.
    event_lower = F(999,1000)
    lower_W = event_lower*gap
    assert lower_W == F(3741,62300) > F(3,50)
    # TV dual bound uses oscillation, not twice the sup norm.
    oscillation = max(w)-min(w)
    assert oscillation == 1
    e = F(1,1000)
    observed_lower = lower_W-2*e*oscillation
    assert observed_lower == F(9041,155750) > F(29,500)
    # The constant one is attained by moving e probability from S to D.
    before = (F(0),F(1),F(0))
    after = (e,1-e,F(0))
    tv = sum((abs(x-y) for x,y in zip(before,after)),F(0))/2
    effect_shift = sum((z*(y-x) for z,x,y in zip(w,before,after)),F(0))
    assert tv == e and effect_shift == e
    # Deleted diagnostics: all three channel rates scaled together leave q
    # unchanged; dropping D/S labels leaves only waiting-time information.
    ratios = []
    for speed in (F(1,11),F(2,3),F(7,4),F(19)):
        d,s = F(1,43)*speed,F(2,580)*speed
        assert d/(d+s) == q0
        ratios.append(str(d/(d+s)))
    return {
        "verdict": "PASS_CONDITIONAL_KINETIC_CLOCK_DIAGNOSTIC_ONLY",
        "independent_groups": 1,
        "numbered_rounds": 0,
        "exact": {
            "D_branch_old": str(q0), "D_branch_changed": str(q1),
            "branch_difference": str(gap),
            "total_rate_ratio_at_same_state": str(rate_ratio),
            "old_W_y_coefficient": str(old_y_coefficient),
            "changed_W_y_coefficient": str(new_y_coefficient),
            "exp7_positive_sum_through_degree_15": str(lower_exp7),
            "finite_W_strict_lower": str(lower_W),
            "record_score_oscillation": str(oscillation),
            "TV_per_history_contract": str(e),
            "observed_W_difference_strict_lower": str(observed_lower)
        },
        "algebra": {
            "directed_channels_per_model": 12,
            "equilibrium_symbol": "z=exp(J/kBT), weights=(1,z,z,1)",
            "generator_times_Gibbs_numerator_zero_coefficientwise": True,
            "finite_record_classifications": classified,
            "old_W_identity_for_all_event_weights": True,
            "L_R_labels_only_coarsened_for_W_not_for_heat": True,
            "positive_common_speed_examples": ratios
        },
        "scope": [
            "Same marked state, fixed labels, same homogeneous rate model.",
            "Equilibrium statement requires V=0 and common temperature.",
            "Time changes multiply all outgoing channels together.",
            "No-event is retained; no postselection of first-jump records.",
            "Actual readout precision and effective-model error are inputs.",
            "No whole-geometry, arbitrary quantum-input, or GR exclusion.",
            "External resistance comparison is not an autonomous resistance update."
        ],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=evaluate()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write("\n")
        mode="WRITE_ONCE_PASS"
    else:
        assert json.loads(target.read_text(encoding="utf-8")) == result
        mode="READ_ONLY_MATCH"
    print(json.dumps({"mode":mode,"verdict":result["verdict"],
                      "strict_lower":result["exact"]["observed_W_difference_strict_lower"]}))
