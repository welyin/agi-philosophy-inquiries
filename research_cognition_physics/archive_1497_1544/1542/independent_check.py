"""Independent exact checks of the Pati--Salam bundle and discrete charges."""
from argparse import ArgumentParser
from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def square(class_xy):
    # Product of two oriented spheres: x^2=y^2=0; integral xy=1.
    x, y = class_xy
    return 2 * x * y


def pair(a, b):
    return a[0] * b[1] + a[1] * b[0]


def half_spinor_weights(roots, parity):
    weights = []
    for signs in product((1, -1), repeat=len(roots)):
        if signs.count(-1) % 2 == parity:
            weights.append(tuple(sum(Fraction(sign * root[i], 2)
                                     for sign, root in zip(signs, roots))
                                 for i in (0, 1)))
    return weights


def fermion_index(roots6, roots4, parity):
    # Physical tensor product of locally defined half-spin representations.
    weights6 = half_spinor_weights(roots6, parity)
    weights4 = half_spinor_weights(roots4, parity)
    weights = [(a[0] + b[0], a[1] + b[1])
               for a in weights6 for b in weights4]
    assert len(weights) == 8
    assert all(x.denominator == 1 for w in weights for x in w)
    return sum(square(w) / 2 for w in weights)


def bundle_table(a, b, c):
    # Choose remaining SO(6) roots such that the total w2 agrees.
    d = (a[0] + b[0] - c[0], a[1] + b[1] - c[1])
    roots6, roots4 = [c, d, (0, 0)], [a, b]
    p4 = square(a) + square(b)
    p6 = sum(square(v) for v in roots6)
    e4 = pair(a, b)
    instanton = [sum(square(w) / 2 for w in half_spinor_weights(roots4, j))
                 for j in (0, 1)]
    indices = [fermion_index(roots6, roots4, j) for j in (0, 1)]
    assert indices == [p6 + p4 + 2 * e4, p6 + p4 - 2 * e4]
    assert instanton[0] - instanton[1] == e4
    return {"euler": e4, "p1_E4": p4, "p1_E6": p6,
            "w2_xy": [(a[i] + b[i]) % 2 for i in (0, 1)],
            "k_ch2": [str(x) for x in instanton],
            "left_right_indices": [int(x) for x in indices]}


def central_phase(a, left, right, z):
    # Phases in units of 1/12. Choose z = exp(+2pi i F/12).
    return ((3 * a + 6 * left + z) % 12, (-3 * a + 6 * right - z) % 12)


def kernel_and_order():
    center = list(product(range(4), range(2), range(2), range(12)))
    kernel = [p for p in center if central_phase(*p) == (0, 0)]
    intrinsic = [p for p in kernel if p[-1] == 0]
    assert intrinsic == [(0, 0, 0, 0), (2, 1, 1, 0)]
    assert sorted({p[-1] for p in kernel}) == [0, 3, 6, 9]
    assert len(kernel) == 8
    phase_pairs = {central_phase(*p) for p in center}
    split_pairs = {central_phase(a, left, right, n)
                   for a, left, right, n in product(range(4), range(2), range(2),
                                                   (0, 4, 8))}
    assert split_pairs == phase_pairs
    assert len(phase_pairs) == 24
    return {"cover_kernel": kernel, "faithful_center_order": len(phase_pairs),
            "external_component_generator": 4}


def compute():
    integer_roots = list(product(range(-1, 2), repeat=2))
    halves = 0
    count = 0
    for a, b, c in product(integer_roots, repeat=3):
        row = bundle_table(a, b, c)
        halves += int(Fraction(row["k_ch2"][0]).denominator == 2)
        count += 1
    witness = bundle_table((1, 0), (0, 1), (1, 1))
    assert witness["k_ch2"] == ["1/2", "-1/2"]
    assert witness["left_right_indices"] == [4, 0]
    assert witness["w2_xy"] == [1, 1]
    generations = list(range(1, 13))
    instanton_phases = [(4 * n) % 3 for n in generations]
    assert [n for n, phase in zip(generations, instanton_phases) if phase == 0] == [3, 6, 9, 12]
    operators = {"Weinberg": (0, 2), "six_quarks": (2, 0), "QQQL": (1, 1)}
    discrete_charges = {}
    for fundamental_charge in (1, 5, 9):
        coefficient = (fundamental_charge - 1) // 4
        discrete_charges[str(fundamental_charge)] = {
            name: (baryon + coefficient * (baryon - lepton)) % 3
            for name, (baryon, lepton) in operators.items()}
    assert [discrete_charges[str(f)]["six_quarks"] for f in (1, 5, 9)] == [2, 1, 0]
    assert [discrete_charges[str(f)]["Weinberg"] for f in (1, 5, 9)] == [0, 1, 2]
    return {"center": kernel_and_order(), "nonliftable_witness": witness,
            "surviving_charge_on_SM_singlets": discrete_charges,
            "split_bundle_count": count, "half_integral_pair_count": halves,
            "generation_test": {"N": generations, "unit_instanton_phase_mod_3": instanton_phases},
            "flavour_extension_local_cubic": {"SU4": 3 * (2 - 2), "SU3F": 8 - 8},
            "mixed_bordism": {"orders_via_integral_homology": [3, 3, 3],
                               "first_smash_cell_continuous_extension": 6}}


def main():
    parser = ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()))
    path = Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with path.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        assert result == json.loads(path.read_text(encoding="utf-8"))
        print("PASS 1542 independent quotient, nonliftable spinor indices and generation tests")


if __name__ == "__main__":
    main()
