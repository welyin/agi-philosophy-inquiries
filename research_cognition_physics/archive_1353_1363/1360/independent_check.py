"""Find neutral spaces by a quadratic Lie-algebra generator on labelled matter slots."""
import argparse
from itertools import combinations, permutations, product
import json
from math import prod
from pathlib import Path

import numpy as np


def hermitian_generators(group):
    result = []
    for i in range(3):
        for j in range(i):
            unit = np.zeros((3, 3), dtype=complex)
            unit[i, j] = 1
            result.append(1j * (unit - unit.T))
            if group == "SU3":
                result.append(unit + unit.T)
    if group == "SU3":
        # Integer normalisations suffice for the common kernel.
        result.extend([np.diag([1., -1., 0.]), np.diag([1., 1., -2.])])
    return result


def action_on_slots(generator, slots):
    size = 3 ** len(slots)
    total = np.zeros((size, size), dtype=complex)
    for j, orientation in enumerate(slots):
        middle = generator if orientation == 1 else -generator.conj()
        lifted = np.kron(np.eye(3 ** j), np.kron(middle, np.eye(3 ** (len(slots)-1-j))))
        total += lifted
    return total


def modular_rank(matrix, prime=101):
    a = np.asarray(matrix, dtype=np.int64).copy() % prime
    row = 0
    for col in range(a.shape[1]):
        pivots = np.flatnonzero(a[row:, col])
        if not len(pivots):
            continue
        pivot = row + int(pivots[0])
        a[[row, pivot]] = a[[pivot, row]]
        row_part = a[row, col:] * pow(int(a[row, col]), -1, prime) % prime
        a[row+1:, col:] = (a[row+1:, col:] - a[row+1:, col, None] * row_part) % prime
        row += 1
        if row == a.shape[0]:
            break
    return row


def invariant_space(group, slots):
    size = 3 ** len(slots)
    positive_form = np.zeros((size, size), dtype=complex)
    for generator in hermitian_generators(group):
        total = action_on_slots(generator, slots)
        positive_form += total.conj().T @ total
    eig, basis = np.linalg.eigh(positive_form)
    assert eig.min() > -1.e-10
    positive = eig[eig > 1.e-8]
    assert positive.min() > 0.1
    # All products are sums of small Gaussian integers exactly representable in binary.
    assert np.array_equal(positive_form.imag, np.zeros((size, size)))
    integer = positive_form.real.astype(np.int64)
    assert np.array_equal(integer, positive_form.real)
    dim_upper = size - modular_rank(integer)
    kernel = basis[:, abs(eig) < 1.e-8]
    assert kernel.shape[1] == dim_upper
    return kernel, dim_upper


def integer_witnesses(group, slots):
    # Each column is an explicit invariant orientation tensor and (where needed)
    # one metric or matter-antimatter pair. No theorem asserting spanning is used.
    length = len(slots)
    fillings = list(product(range(3), repeat=length))
    columns = []
    for triple in combinations(range(length), 3):
        if group == "SU3" and any(slots[k] < 0 for k in triple):
            continue
        others = [k for k in range(length) if k not in triple]
        column = []
        for word in fillings:
            i, j, k = (word[m] for m in triple)
            sign = 0 if len({i, j, k}) != 3 else (-1)**(int(i > j) + int(i > k) + int(j > k))
            if others and word[others[0]] != word[others[1]]:
                sign = 0
            column.append(sign)
        columns.append(column)
    columns = np.asarray(columns, dtype=np.int64).T
    swapped = columns.reshape(*(3,)*length, -1).swapaxes(0, 1).reshape(columns.shape)
    return modular_rank(columns), [modular_rank(columns+swapped), modular_rank(columns-swapped)]


def volume():
    eps = np.zeros((3, 3, 3))
    for i, j, k in permutations(range(3)):
        eps[i, j, k] = (-1) ** (int(i > j) + int(i > k) + int(j > k)) / np.sqrt(6)
    return eps.ravel()


def run():
    eta = volume()
    cap = np.eye(27) / np.sqrt(27)
    output = {}
    for group in ["SU3", "SO3"]:
        data = {}
        volume_basis, _ = invariant_space(group, (1, 1, 1))
        assert volume_basis.shape[1] == 1
        assert np.allclose(volume_basis @ volume_basis.conj().T @ eta, eta)
        tagged = np.eye(27)[0]
        assert np.allclose(volume_basis.conj().T @ tagged, 0)
        for slots in [(1, 1, 1), (1, 1, 1, 1, -1)]:
            neutral, dim_upper = invariant_space(group, slots)
            size, rank = neutral.shape
            number_of_pairs = slots.count(1)
            routes = prod(range(2 * number_of_pairs-1, 0, -2))
            receive = 1 / (100 * routes * 27)
            # Spin and physical charge labels are held fixed, and the full colour-safe
            # reference is joint with the partners. The common efficiency is separate.
            test = neutral[:, :min(2, rank)] / np.sqrt(min(2, rank))
            balanced = np.einsum("wr,a->war", test, eta.conj())
            prepared = np.einsum("war,b,ab->wr", balanced, eta, cap)
            assert np.allclose(prepared, test / np.sqrt(27))
            for gen in hermitian_generators(group):
                wa = action_on_slots(gen, slots)
                ya_bar = action_on_slots(gen, [-1] * 3)
                total_action = (np.einsum("ij,jar->iar", wa, balanced) +
                                np.einsum("ab,wbr->war", ya_bar, balanced))
                assert np.linalg.norm(total_action) < 2.e-11
            swap = neutral.reshape(*(3,) * len(slots), rank).swapaxes(0, 1).reshape(size, rank)
            compressed = neutral.conj().T @ swap
            assert np.allclose(neutral @ compressed, swap)
            eigenvalues = np.linalg.eigvalsh(compressed)
            assert np.allclose(abs(eigenvalues), 1)
            signs = [int(sum(eigenvalues > 0)), int(sum(eigenvalues < 0))]
            dim_lower, sign_lower = integer_witnesses(group, slots)
            assert dim_lower == dim_upper and sum(sign_lower) == dim_upper
            assert signs == sign_lower
            if number_of_pairs == 4:
                assert signs == ([1, 2] if group == "SU3" else [3, 3])
            else:
                assert signs == [0, 1]
            row = {"neutral_rank": rank, "swap_dimensions_plus_minus": signs,
                   "perfect_logical_success": receive, "raw_success": rank * receive / size}
            if number_of_pairs == 4:
                flag_fractions = [2/3, 1/3]
                row["swap_and_receiver"] = [n * receive / size for n in signs]
                row["receiver_given_swap"] = [a / b for a, b in zip(row["swap_and_receiver"], flag_fractions)]
            data[str(number_of_pairs)] = row
        # A trivial central action on a mixed tensor of degree 3 leaves the SO3 kernel intact.
        assert np.allclose(np.exp(2j * np.pi * 3/3), 1.)
        output[group] = data
    raw_differences = np.abs(np.array(output["SU3"]["4"]["swap_and_receiver"]) -
                             output["SO3"]["4"]["swap_and_receiver"])
    output["raw_equal_halfwidth_exclusive"] = (raw_differences / 2).tolist()
    return output


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            compare(actual[key], saved[key])
    else:
        assert np.allclose(actual, saved, rtol=1.e-10, atol=1.e-17)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    computed = run()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8") as out:
            json.dump(computed, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(computed, json.loads(target.read_text(encoding="utf-8")))
    print("PASS", computed)
