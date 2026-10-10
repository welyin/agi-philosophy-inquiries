"""Independent finite checks for 1065; standard library only.

Default mode recomputes without writing.  --write creates this script's result
file exclusively; it never imports the author implementation.  Finite cyclic
quotients are algebraic controls, not simulations/proofs of solenoid topology.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def encoded(value: F) -> str:
    return str(value)


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def valuation_three(value: int) -> int:
    require(value > 0, "valuation is only used on nonzero residues")
    exponent = 0
    while value % 3 == 0:
        value //= 3
        exponent += 1
    return exponent


def finite_odd_quotients() -> dict:
    rows = []
    total_nonzero = 0
    total_nonzero_cycles = 0
    for depth in range(1, 9):
        modulus = 3**depth
        inverse_two = pow(2, -1, modulus)
        halve = [(inverse_two * x) % modulus for x in range(modulus)]
        require(sorted(halve) == list(range(modulus)), "halving is a permutation")
        for x in range(modulus):
            require((2 * halve[x]) % modulus == x, "right inverse")
            require(halve[(2 * x) % modulus] == x, "left inverse")
        for x in range(1, modulus):
            require(valuation_three(halve[x]) == valuation_three(x), "p-adic size invariant")

        visited = set()
        cycles = []
        for start in range(modulus):
            if start in visited:
                continue
            cycle = []
            current = start
            while current not in visited:
                visited.add(current)
                cycle.append(current)
                current = halve[current]
            require(current == start, "permutation cycle closes at its start")
            cycles.append(cycle)

        nonzero_cycles = [cycle for cycle in cycles if cycle != [0]]
        expected_lengths = [2 * 3**j for j in range(depth)]
        require(sorted(map(len, nonzero_cycles)) == expected_lengths, "finite cycle spectrum")
        # This is a symbolic certificate for arbitrary real costs: in each
        # cycle the coefficients of sum_x [C(halve(x)) - C(x)] all cancel.
        witnesses = []
        for cycle in nonzero_cycles:
            coefficients = Counter()
            for x in cycle:
                coefficients[halve[x]] += 1
                coefficients[x] -= 1
            require(all(v == 0 for v in coefficients.values()), "cycle cost differences telescope")
            # An unrelated positive rational cost supplies a concrete maximum
            # witness; it is not assumed to descend from a continuous cost.
            costs = {x: F(1 + (x*x + 3*x + 5) % (modulus + 7), modulus + 7) for x in cycle}
            maximizer = max(cycle, key=costs.__getitem__)
            double = (2 * maximizer) % modulus
            require(double in costs and halve[double] == maximizer, "same cycle contains double")
            require(costs[halve[double]] >= costs[double], "strict decrease fails at double of maximizer")
            require(sum((costs[halve[x]] - costs[x] for x in cycle), F(0)) == 0, "numeric telescope")
            witnesses.append({"cycle_length": len(cycle), "maximizer": maximizer, "double": double})
        rows.append({
            "depth": depth,
            "modulus": modulus,
            "nonzero_cycles": len(nonzero_cycles),
            "nonzero_cycle_lengths": expected_lengths,
            "maximum_witnesses": witnesses,
        })
        total_nonzero += modulus - 1
        total_nonzero_cycles += len(nonzero_cycles)
    return {
        "passed": True,
        "nonzero_residues_checked": total_nonzero,
        "nonzero_cycles_checked": total_nonzero_cycles,
        "rows": rows,
        "scope": "Finite Z/(3^k) only; cycle balance excludes strict decrease for every nonzero point of a cycle. Connectedness, small subgroups and non-Lie status of the solenoid are analytic claims, not numerical outputs.",
    }


def log_two_enclosure(terms: int = 40) -> tuple[F, F]:
    # log 2 = 2 atanh(1/3).  Replace each remaining odd denominator
    # by its first value to bound the positive geometric tail.
    lower = 2 * sum((F(1, (2*j + 1) * 3**(2*j + 1)) for j in range(terms)), F(0))
    tail_upper = F(9, 4 * (2*terms + 1) * 3**(2*terms + 1))
    upper = lower + tail_upper
    require(F(69, 100) < lower < upper < F(7, 10), "rational log-two enclosure")
    require(upper - lower < F(1, 10**39), "enclosure sufficiently narrow")
    return lower, upper


def slow_log_cost() -> dict:
    ell_lo, ell_hi = log_two_enclosure()
    maximum = F(1, 3)
    rows = []
    for log_start in (4, 5, 8, 16, 64, 256, 1024):
        initial = F(1, log_start)
        require(F(log_start) >= 3 + ell_hi, "nonempty D_a certified")
        # Work in the inverse-cost coordinate rather than evaluating tiny x:
        # |x|=exp(1-1/c), so doubling subtracts ell from 1/c.
        margin_lo = initial**2 * ell_lo / (1 - initial * ell_lo)
        margin_hi = initial**2 * ell_hi / (1 - initial * ell_hi)
        denominator_floor = 1 - maximum * ell_hi
        require(denominator_floor > 0 and margin_lo > 0, "whole rational cost-domain certificate")
        # For every c in [initial, maximum] and ell in [ell_lo, ell_hi],
        # numerator c^2 ell >= initial^2 ell_lo and
        # 0 < 1-c ell <= 1-initial ell_lo.  Thus Delta(c)>=margin_lo.
        for c in (initial, (initial + maximum) / 2, maximum):
            for ell in (ell_lo, (ell_lo + ell_hi) / 2, ell_hi):
                next_cost = 1 / (1/c - ell)
                delta = next_cost - c
                require(next_cost == c / (1 - c*ell), "inverse-coordinate identity")
                require(delta == c*c*ell / (1 - c*ell), "increment identity")
                require(delta >= margin_lo, "finite rational cross-check of domain bound")

        lower_exit_ratio = (log_start - 3) / ell_hi
        upper_exit_ratio = (log_start - 3) / ell_lo
        require(floor_fraction(lower_exit_ratio) == floor_fraction(upper_exit_ratio), "exact exit integer isolated")
        exact_exit = floor_fraction(lower_exit_ratio) + 1
        require(log_start - (exact_exit - 1)*ell_hi >= 3, "previous endpoint still in K")
        require(log_start - exact_exit*ell_lo < 3, "certified first exit from K")

        conservative_exit = floor_fraction((maximum - initial) / margin_lo) + 1
        require(initial + conservative_exit * margin_lo > maximum, "finite escape by cost budget")
        require(exact_exit <= conservative_exit, "independent exact exit respects conservative bound")
        epsilon = margin_lo / 4
        robust_margin = margin_lo - epsilon
        robust_exit = floor_fraction((maximum - initial) / robust_margin) + 1
        require(initial + robust_exit * robust_margin > maximum, "explicit true-increment slack certificate")
        rows.append({
            "L0": log_start,
            "a": encoded(initial),
            "M": encoded(maximum),
            "delta_lower": encoded(margin_lo),
            "delta_upper": encoded(margin_hi),
            "denominator_lower_on_rational_cost_domain": encoded(denominator_floor),
            "exact_first_exit": exact_exit,
            "conservative_exit_using_delta_lower": conservative_exit,
            "epsilon_true_increment_slack": encoded(epsilon),
            "conservative_exit_with_slack": robust_exit,
        })

    empty_log_start = F(7, 2)
    require(empty_log_start > 3 and empty_log_start - 3 < ell_lo, "empty D_a boundary control")
    no_uniform_examples = []
    for ratio in (F(1, 2), F(9, 10), F(99, 100), F(999, 1000)):
        log_coordinate = max(4, floor_fraction(ratio * ell_hi / (1 - ratio)) + 1)
        certified_ratio_lower = F(log_coordinate) / (log_coordinate + ell_hi)
        require(certified_ratio_lower > ratio, "proposed fixed q is violated")
        no_uniform_examples.append({
            "proposed_q": encoded(ratio),
            "L": log_coordinate,
            "halving_cost_ratio_lower": encoded(certified_ratio_lower),
        })
    return {
        "passed": True,
        "log2_lower": encoded(ell_lo),
        "log2_upper": encoded(ell_hi),
        "log2_interval_width": encoded(ell_hi - ell_lo),
        "finite_windows": rows,
        "empty_window_control_L0": encoded(empty_log_start),
        "fixed_q_counterchecks": no_uniform_examples,
        "scope": "Exact rational certificates using inverse cost and a proved log2 enclosure. No floating-point evaluation of exp(1-L), no author-code import. Failure of every q<1 on every identity neighborhood is established analytically by L/(L+log2)->1; finitely many q checks do not prove that quantifier. Epsilon is slack in the true increment, not automatically a one-endpoint measurement error.",
    }


def run() -> dict:
    return {
        "round": 1065,
        "checker": "independent inverse-cost and finite-cycle verification",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest().upper(),
        "passed": True,
        "odd_quotients": finite_odd_quotients(),
        "slow_log": slow_log_cost(),
        "analytic_only": [
            "Strict doubling cost increase excludes a compact subgroup by its positive cost maximum.",
            "Local compactness gives a compact closure of a sufficiently small identity neighborhood.",
            "Connected compact odd-p solenoid has unique continuous local roots but arbitrarily small compact subgroups.",
            "NSS-to-Lie is a mature theorem; no dimension or complete cognitive necessity follows from this script.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusively create independent_results.json")
    args = parser.parse_args()
    payload = run()
    destination = Path(__file__).with_name("independent_results.json")
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        mode = "created"
    elif destination.exists():
        require(json.loads(destination.read_text(encoding="utf-8")) == payload, "saved independent result matches recomputation")
        mode = "read_only_recomputed_and_matched"
    else:
        mode = "read_only_recomputed_no_saved_result"
    print(json.dumps({
        "passed": True,
        "mode": mode,
        "nonzero_residues": payload["odd_quotients"]["nonzero_residues_checked"],
        "nonzero_cycles": payload["odd_quotients"]["nonzero_cycles_checked"],
        "finite_windows": len(payload["slow_log"]["finite_windows"]),
        "exact_exit_layers": [row["exact_first_exit"] for row in payload["slow_log"]["finite_windows"]],
        "conservative_exit_layers": [row["conservative_exit_using_delta_lower"] for row in payload["slow_log"]["finite_windows"]],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
