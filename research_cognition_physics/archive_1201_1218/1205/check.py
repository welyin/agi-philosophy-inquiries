"""1205 fixed preparation circuit and normal-clock resource certificate.

No writes unless --write; --write creates results.json only if absent.
The bit circuit tests are finite diagnostics; proof.md gives uniform gate bounds.
The fixed word program uses one backend for exact semantics and macro counting.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext, ROUND_CEILING, ROUND_FLOOR
from collections import Counter
import argparse
import hashlib
import json
import math
import random

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


class Bits:
    """Each Boolean output is a fresh zero wire; all old wires stay."""
    def __init__(self):
        self.values = [0, 1]
        self.initial = [0, 1]
        self.gates = []
        self.booleans = 0

    def input(self, x, w):
        out = []
        for k in range(w):
            out.append(len(self.values))
            self.values.append((x >> k) & 1)
            self.initial.append((x >> k) & 1)
        return out

    def const(self, x, w):
        return [(x >> k) & 1 for k in range(w)]

    def gate(self, controls, target):
        controls = tuple(dict.fromkeys(controls))
        assert target not in controls and len(controls) <= 2
        self.values[target] ^= int(all(self.values[c] for c in controls))
        self.gates.append((controls, target))

    def boolean(self, kind, a, b=None):
        t = len(self.values)
        self.values.append(0)
        self.initial.append(0)
        self.booleans += 1
        if kind == "and":
            self.gate((a, b), t)
        elif kind == "xor":
            self.gate((a,), t)
            self.gate((b,), t)
        elif kind == "not":
            self.gate((), t)
            self.gate((a,), t)
        else:
            raise ValueError(kind)
        return t

    def add(self, a, b, carry=0):
        out = []
        for ai, bi in zip(a, b):
            x = self.boolean("xor", ai, bi)
            out.append(self.boolean("xor", x, carry))
            p = self.boolean("and", ai, bi)
            q = self.boolean("and", x, carry)
            carry = self.boolean("xor", p, q)
        return out, carry

    def sub(self, a, b):
        return self.add(a, [self.boolean("not", x) for x in b], 1)

    def mux(self, s, a, b):
        out = []
        for ai, bi in zip(a, b):
            x = self.boolean("xor", ai, bi)
            y = self.boolean("and", s, x)
            out.append(self.boolean("xor", ai, y))
        return out

    def mul(self, a, b):
        w = len(a)
        out = [0]*w
        for j in range(w):
            row = [0]*j + [self.boolean("and", a[k], b[j])
                           for k in range(w-j)]
            out, _ = self.add(out, row)
        return out

    def div_const(self, a, m):
        assert 1 <= m < 1 << len(a)
        w = len(a)
        rem = [0]*(w+1)
        divisor = self.const(m, w+1)
        quotient = [0]*w
        for k in range(w-1, -1, -1):
            shifted = [a[k]] + rem[:w]
            reduced, ge = self.sub(shifted, divisor)
            rem = self.mux(ge, shifted, reduced)
            quotient[k] = ge
        return quotient, rem

    def value(self, wires):
        return sum(self.values[x] << i for i, x in enumerate(wires))

    def reverse_check(self):
        values = self.values.copy()
        for controls, target in reversed(self.gates):
            values[target] ^= int(all(values[c] for c in controls))
        assert values == self.initial


def check_primitives():
    tests = 0
    maximum_ratio = F(0)
    for w in range(1, 5):
        for a in range(1 << w):
            for b in range(1 << w):
                for kind in ("add", "sub", "mul", "mux"):
                    c = Bits()
                    aa, bb = c.input(a, w), c.input(b, w)
                    if kind == "add":
                        z, carry = c.add(aa, bb)
                        assert c.value(z) == (a+b) % (1 << w)
                        assert c.values[carry] == (a+b) >> w
                    elif kind == "sub":
                        z, ge = c.sub(aa, bb)
                        assert c.value(z) == (a-b) % (1 << w)
                        assert c.values[ge] == int(a >= b)
                    elif kind == "mul":
                        z = c.mul(aa, bb)
                        assert c.value(z) == a*b % (1 << w)
                    else:
                        s = (a+b) % 2
                        z = c.mux(s, aa, bb)
                        assert c.value(z) == (b if s else a)
                    assert c.booleans <= 256*w*w
                    assert len(c.gates) <= 3*c.booleans
                    c.reverse_check()
                    maximum_ratio = max(maximum_ratio, F(c.booleans, w*w))
                    tests += 1
            for m in range(1, 1 << w):
                c = Bits()
                q, r = c.div_const(c.input(a, w), m)
                assert c.value(q) == a//m and c.value(r) == a % m
                assert c.booleans <= 256*w*w
                c.reverse_check()
                maximum_ratio = max(maximum_ratio, F(c.booleans, w*w))
                tests += 1
    return {"exhaustive_w1_to_w4_cases": tests,
            "maximum_boolean_count_over_w_squared": str(maximum_ratio),
            "reverse_restores_all_inputs_and_zero_ancillas": True,
            "uniform_bound_is_analytic_not_inferred_from_tests": True}


class Word:
    """Fixed-width macro interpreter; no data-dependent loop or gate count."""
    def __init__(self, w=256):
        self.w, self.mask = w, (1 << w)-1
        self.count = Counter()
        self.maximum_unsigned_input = 0

    def op(self, kind, a, b=0, s=0):
        assert 0 <= a <= self.mask and 0 <= b <= self.mask
        self.maximum_unsigned_input = max(self.maximum_unsigned_input, a, b)
        self.count[kind] += 1
        if kind == "and": return a & b
        if kind == "or": return a | b
        if kind == "notbit":
            assert a in (0, 1)
            return a ^ 1
        if kind == "lt": return int(a < b)
        if kind == "eq": return int(a == b)
        if kind == "add": return (a+b) & self.mask
        if kind == "sub": return (a-b) & self.mask
        if kind == "mul": return (a*b) & self.mask
        if kind == "div":
            assert 1 <= b <= 100
            return a//b
        if kind == "mux":
            assert s in (0, 1)
            return b if s else a
        raise ValueError(kind)


def fixed_program(machine, tables, dimensions, tag, seeds):
    o = machine.op
    mask = o("mux", (1 << 173)-1, (1 << 195)-1, tag)
    dim = o("mux", dimensions[0], dimensions[1], tag)
    selected, seen = 0, 0
    for j in range(10):
        x = o("and", seeds[j], mask)
        valid = o("lt", x, dim)
        take = o("and", valid, o("notbit", seen))
        selected = o("mux", selected, x, take)
        seen = o("or", seen, valid)
    counts, ways, residual = [0]*4, 0, 0
    for table_tag in range(2):
        for n, weight, start in tables[table_tag]:
            end = start + weight  # readonly adjacent prefix, not runtime add
            match = o("and", o("eq", tag, table_tag),
                      o("and", o("notbit", o("lt", selected, start)),
                        o("lt", selected, end)))
            residual = o("mux", residual, o("sub", selected, start), match)
            ways = o("mux", ways, o("sub", end, start), match)
            counts = [o("mux", counts[i], n[i], match) for i in range(4)]
    word = []
    for remaining in range(100, 0, -1):
        old_k, old_ways = residual, ways
        block = []
        for i in range(4):
            # The next assertion is part of the valid-input certificate only.
            assert old_ways*counts[i] < (1 << 207)
            product = o("mul", old_ways, counts[i])
            assert product % remaining == 0
            block.append(o("div", product, remaining))
        prefix = [0]
        for i in range(4):
            prefix.append(o("add", prefix[-1], block[i]))
        flags = [o("and", o("notbit", o("lt", old_k, prefix[i])),
                   o("lt", old_k, prefix[i+1])) for i in range(4)]
        assert sum(flags) == 1
        residual, ways = 0, 0
        for i in range(4):
            residual = o("mux", residual, o("sub", old_k, prefix[i]), flags[i])
            ways = o("mux", ways, block[i], flags[i])
        counts = [o("sub", counts[i], flags[i]) for i in range(4)]
        lo = o("or", flags[1], flags[3])
        hi = o("or", flags[2], flags[3])
        word.append(lo + 2*hi)  # two output wires, not a runtime arithmetic gate
    assert counts == [0]*4 and residual == 0 and ways == 1
    routes = [[o("and", seen, o("eq", x, j)) for j in range(4)] for x in word]
    assert all(sum(r) == seen for r in routes)
    return {"success": seen, "selected_index": selected,
            "word": word, "routes": routes}


def tables_from_saved_boxes(old):
    boxes = old["preparation"]["root_boxes"]
    den = old["preparation"]["root_denominator"]
    fac = [math.factorial(i) for i in range(101)]
    tables, dimensions = [[], []], [0, 0]
    for a in range(101):
        for b in range(101-a):
            for c in range(101-a-b):
                n = (a, b, c, 100-a-b-c)
                low = sum(n[i]*boxes[i][0] for i in range(4))
                high = sum(n[i]*boxes[i][1] for i in range(4))
                for tag, ends in enumerate(((39, 41), (79, 81))):
                    inside = 2*low >= ends[0]*den and 2*high <= ends[1]*den
                    outside = 2*high < ends[0]*den or 2*low > ends[1]*den
                    assert inside or outside
                    if not inside:
                        continue
                    weight = fac[100]//math.prod(fac[k] for k in n)
                    tables[tag].append((n, weight, dimensions[tag]))
                    dimensions[tag] += weight
    return tables, dimensions


def rank_word(word, rows):
    counts = [word.count(i) for i in range(4)]
    row = next(row for row in rows if tuple(counts) == row[0])
    index = row[2]
    for letter in word:
        for earlier in range(letter):
            if counts[earlier]:
                counts[earlier] -= 1
                index += math.factorial(sum(counts))//math.prod(
                    math.factorial(v) for v in counts)
                counts[earlier] += 1
        counts[letter] -= 1
    return index


def upward_fraction(x, places=30):
    with localcontext() as c:
        c.prec = 90
        z = D(x.numerator)/D(x.denominator)
        return str((z+D("1e-80")).quantize(D(10)**(-places), ROUND_CEILING))


def downward_fraction(x, places=30):
    with localcontext() as c:
        c.prec = 90
        z = D(x.numerator)/D(x.denominator)
        return str((z-D("1e-80")).quantize(D(10)**(-places), ROUND_FLOOR))


def entropy_upper(q, d):
    """Analytic - (1-q)log(1-q)<=q, log(d)<d; no float certificate."""
    # For q<1e-12, log(1/q) bounded by a supplied power-of-two bracket.
    k = 0
    while q < F(1, 2**k):
        k += 1
    # log 2 < 1, hence h2(q) <= q*(1+k).
    return q*(1+k+d)  # log(dimension-1) <= d; d is a log-dimension bound


def run():
    primitive = check_primitives()
    old = json.loads((STAGE/"1202/results.json").read_text(encoding="utf-8"))
    flow = json.loads((STAGE/"1204/results.json").read_text(encoding="utf-8"))
    tables, dimensions = tables_from_saved_boxes(old)
    assert list(map(len, tables)) == [1728, 2391]
    assert dimensions[1] == int(flow["new_reference"]["dimension"])
    counts = None
    cases = []
    for tag in range(2):
        dim = dimensions[tag]
        examples = [[dim]*10, [0]+[dim]*9, [dim]*9+[dim-1],
                    [dim, dim//2]+[0]*8]
        for seeds in examples:
            # Cold unused upper bits remain in the actual input register.
            if tag == 0:
                seeds = [x + (((1 << 22)-1) << 173) for x in seeds]
            machine = Word()
            out = fixed_program(machine, tables, dimensions, tag, seeds)
            assert rank_word(out["word"], tables[tag]) == out["selected_index"]
            expected = next((x & ((1 << (173 if tag == 0 else 195))-1)
                             for x in seeds
                             if (x & ((1 << (173 if tag == 0 else 195))-1)) < dim), None)
            assert out["success"] == int(expected is not None)
            assert out["selected_index"] == (0 if expected is None else expected)
            if counts is None:
                counts = machine.count
            assert machine.count == counts
            cases.append({"window": tag, "success": out["success"],
                          "rank": str(out["selected_index"])})
    word_bound = 2*(32*4119+8*10+64*100+8*100+50)
    same_h_menus = []
    for tags in ((0, 1), (1, 0)):
        machine = Word()
        for tag in tags:
            dim = dimensions[tag]
            seeds = [dim]*9+[dim-1]
            out = fixed_program(machine, tables, dimensions, tag, seeds)
            assert out["success"] == 1
            assert rank_word(out["word"], tables[tag]) == dim-1
        assert machine.count == Counter({k: 2*v for k,v in counts.items()})
        same_h_menus.append({"window_tags": list(tags),
                             "word_count": sum(machine.count.values()),
                             "same_readonly_label_controlled_program": True})
    assert 2*sum(counts.values()) <= word_bound
    nb = word_bound*256*256**2
    ng = 3*nb+800
    L = 2*10**13
    assert nb < 10**13 and ng < L
    T = 3*L+606
    B = F(66, 7)*(L+1)+150
    sigma, M = F(1, 10**43), 10**112
    correction = F(T*T, 8)/(M*M*sigma**4)
    eps = B*T*sigma*(1+correction)
    eps_cap = F(4, 10**15)
    assert eps < eps_cap
    tail = B*(sigma*(1+correction)+B*T*T/F(2*M))
    assert tail < F(2, 10**28)
    failures = [(1-F(dimensions[tag], 2**(173 if tag == 0 else 195)))**10
                for tag in range(2)]
    assert max(failures) < F(1, 1000)
    pstar = (1-max(failures))**2
    assert pstar > F(998, 1000)
    eta = 2*eps_cap/(F(998, 1000)-eps_cap)
    assert eta < F(9, 10**15)
    qclock = eps_cap**2/(F(998, 1000)-eps_cap)
    # logrank <3e13 log2 <3e13; h2 bound uses log2<1.
    clock_entropy = entropy_upper(qclock, 3*10**13)
    assert clock_entropy < F(5, 10**16)
    # Single contact AK<=16 or Y=K<=4: 256 is a conservative common dimension.
    fannes = entropy_upper(F(9, 10**15), 8)  # log255<8
    thermotransport = 2*fannes+2*(F(6, 5)+F(1, 3))*F(1094,1000)*F(9,10**15)
    assert thermotransport < F(2,10**12)
    # Preserve the already-certified two-reference smooth entropy inequality,
    # comparing finite A,B,K at both contact sections (not resetting C).
    d0 = F(flow["multi_reference_thermodynamics"]
           ["cold_hot_initial_relative_entropy_defect_upper"])
    clausius_defect = d0+thermotransport
    transport = {}
    menus = flow["actual_flow_predictions"]
    for name in ("t600_lambda1", "t300_lambda1"):
        row = menus[name]
        v = 600//row["t"]
        shape = F(1, 10000*v)
        final_error = F(9,10**15)+shape
        # Difference of energies at initial and final contact sections.
        eradius = F(1094,1000)*(2*F(9,10**15)+shape)+F(1,10**10)
        q = F(row["Q_dir_center"])
        assert q-eradius > 0
        transport[name] = {
            "Q_dir_center": str(q),
            "new_actual_energy_radius_upper": upward_fraction(eradius),
            "Q_dir_lower": downward_fraction(q-eradius),
            "rectangular_shape_error": str(shape),
            "pointer_probability_error_upper": upward_fraction(final_error+F(1,10**10)),
            "finite_proxy_positive_margin": downward_fraction(q-eradius-F(".01094")),
            "proxy_margin_requires_extra_success_flag_readout_contract": True,
            "old_K_systematic_budget_does_not_bound_success_flag_instrument": True,
            "full_old_statistical_protocol_not_automatically_certified": True,
            "actual_B_energy_still_prediction_not_a_new_pointer": True}
    return {
        "schema": "fixed_preparation_autonomous_bridge_1205_v1",
        "status": "AUTHOR_CANDIDATE_NOT_ACCEPTANCE",
        "primitive_checks": primitive,
        "fixed_program": {
            "N": 100, "word_width": 256, "types": list(map(len,tables)),
            "dimensions": list(map(str, dimensions)),
            "word_counts_each_pool": dict(sorted(counts.items())),
            "actual_word_count_two_pools": 2*sum(counts.values()),
            "word_operation_bound": word_bound, "boolean_output_bound": nb,
            "reversible_gate_bound": ng, "padded_program_slots": L,
            "working_bit_bound": 10**13, "packed_table_bits": 880889,
            "unpacked_constant_words_bits_upper": 10**7,
            "uniform_seed_bits": 3900, "fixed_tags_extra_bits": 3,
            "material_four_level_units": 1000, "pointer_dimension": 4,
            "program_hardware_description_bits_upper": 256*L,
            "all_seed_garbage_failure_stock_and_reference_kept": True,
            "test_cases": cases,
            "forward_reverse_same_H_menus": same_h_menus,
            "worst_joint_success_probability_lower": ".998",
            "recorded_initial_resources_are_inputs": True},
        "normal_clock": {
            "M": str(M), "sigma": str(sigma), "common_T": T,
            "force_bound": str(B),
            "epsilon_clock_upper": upward_fraction(eps),
            "epsilon_clock_contract": str(eps_cap),
            "endpoint_total_interaction_upper": upward_fraction(tail,45),
            "conditional_joint_trace_error_upper": upward_fraction(eta),
            "conditional_clock_entropy_upper": upward_fraction(clock_entropy),
            "clock_mean_energy_upper": "2e112 + 1.25e-27",
            "clock_energy_variance_upper": "2e86",
            "success_conditioned_reference_clock_not_reset": True,
            "same_conserved_success_event_at_both_contact_sections": True,
            "conserved_seed_event_is_not_a_physical_oracle_gate": True,
            "computed_flag_used_for_gates_and_terminal_data_has_separate_error": True,
            "shape_error_not_used_in_clock_entropy": True},
        "thermodynamic_transport": {
            "finite_entropy_and_energy_extra_upper": upward_fraction(thermotransport),
            "two_reference_clausius_defect_upper": upward_fraction(clausius_defect),
            "exact_initial_independence_not_assumed_after_real_preparation": True,
            "whole_preparation_energy_not_mislabeled_heat": True},
        "native_flow_transport": transport,
        "source_sha256": {f"{n}/results.json": hashlib.sha256(
            (STAGE/str(n)/"results.json").read_bytes()).hexdigest() for n in (1202,1203,1204)},
        "not_claimed": ["native QQ generates Boolean/inventory gates",
                        "free pure seeds, blanks, stocks, program or clock",
                        "universal controller or minimum resource cost",
                        "physical realization or stage completion",
                        "3D positioning, area entropy, Unruh or GR"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    out = HERE/"results.json"
    if args.write:
        if out.exists():
            raise FileExistsError("Refusing to overwrite results.json")
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        assert result == json.loads(out.read_text(encoding="utf-8"))
    print(json.dumps({"status": "AUTHOR_CHECKS_PASS", "results_match": True,
                      "results_file": str(out)}, ensure_ascii=False))
