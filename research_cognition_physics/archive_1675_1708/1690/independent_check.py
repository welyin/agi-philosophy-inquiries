"""R1690. Check the genuine Weyl weights and the order of a quotient-sensitive test."""
import argparse
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import json
from math import comb, gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
N = 12


def load(relative):
    return json.loads((HERE/relative).read_text(encoding="utf-8"))


def weights(rep, vector):
    if rep == "10":
        return [a+b for a,b in combinations(vector, 2)]
    if rep == "5":
        return list(vector)
    if rep == "5bar":
        return [-v for v in vector]
    return [0]


def flavor_test(charges):
    a1 = sum(charges)
    a3 = sum(q**3 for q in charges)
    return (Fraction(182*a3, 72).denominator == 1 and
            Fraction(2*a1, 12).denominator == 1)


def competitor_coefficients(m, t):
    # Alphabet is (c1(det E_2), c1(continuous F)).
    pairs = [(m,t),(-m,t),(0,-t),(0,-t)]
    return [Fraction(sum(comb(3,k)*y**k*f**(3-k) for y,f in pairs), 6)
            for k in range(4)]


def calculate():
    src = load("../1688/inputs.json")
    old = load("../1685/inputs.json")["flavor_charge"]
    reps = [("10",q) for q in src["su5_generation_charges_10"]]
    reps += [("5bar",q) for q in src["su5_generation_charges_5bar"]]
    for name in ["P", "T"]:
        reps += [("5",old[name]),("5bar",old[name+"bar"])]
    reps += [("1",old[key]) for key in src["neutral_charged_under_flavor"]]
    dim = [0]*5
    charges = [q for rep,q in reps for _ in weights(rep,dim)]
    singles = [q for q in range(N) if flavor_test([*charges,q])]
    pairs = [qs for qs in combinations_with_replacement(range(N),2)
             if flavor_test([*charges,*qs])]
    # Sum chiral weights instead of taking tabulated cubic/quadratic indices.
    probes = [[*v,-sum(v)] for v in product(range(-1,2), repeat=4)]
    for v in probes:
        assert sum(w**3 for rep,q in reps for w in weights(rep,v)) == 0
    residues = set()
    for a,b in product(probes, repeat=2):
        index = sum(x*y for x,y in zip(a,b))
        anomaly = sum(q*sum(x*y for x,y in zip(weights(rep,a),weights(rep,b)))
                      for rep,q in reps)
        assert anomaly == N*index
        residues.add(anomaly % N)
    # Realize each allowed evaluation a z^2 on S^1 x S^2 x S^2 via z=u+k v.
    even_values = sorted({(2*k) % N for k in range(N)})
    by_m = {}
    for m in range(1, 13):
        solutions = []
        for t in range(N):
            assert competitor_coefficients(m,t) == [0,0,t*m**2,0]
            assert competitor_coefficients(6*m,t) == [0,0,36*t*m**2,0]
            order = N//gcd(N,2*t*m*m)
            assert order == 6//gcd(6,t*m*m)
            passes = all((i*t*m*m) % N == 0 for i in even_values)
            assert passes == (order == 1)
            if passes:
                solutions.append(t)
        by_m[m] = solutions
    assert singles == [2,8] and len(pairs) == 6
    return {
        "restricted_field_table": {
            "flavor_charged_Weyl_linear_cubic_base": [sum(charges),sum(q**3 for q in charges)],
            "single_charge_choices": singles, "two_charge_choices": pairs,
            "traceless_integer_weight_probes": len(probes),
            "S2xS2_bundles_probed": len(probes)**2,
            "mixed_integer_index_per_fundamental": N,
            "mixed_residues_observed": sorted(residues)},
        "extension_by_conjugate_hypercharge_pair_and_two_neutrals": {
            "reachable_mod12_a_times_det_c1_squared": even_values,
            "flavor_charges_allowed_by_hypercharge_m": by_m,
            "m_t_1_1_phase_in_turns_quotient": str(Fraction(2, 12)),
            "same_probe_lifted_integer_circle_phase_in_turns_mod_one": str(Fraction(72,12)%1),
            "m_t_1_1_anomaly_order": 6}
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    target = HERE/"independent_results.json"
    result = json.dumps(calculate(), indent=2, ensure_ascii=False)+"\n"
    if args.check:
        assert json.loads(result) == load(target)
        print("R1690 independent OK")
    else:
        target.write_text(result, encoding="utf-8")
