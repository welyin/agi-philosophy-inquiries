"""Independent 1064 checks; imports no author or earlier research code.

The default only prints results. --write exclusively creates this script's own
independent_results.json and refuses to replace an existing file. Integer
matrix-unit identities and modular rank certificates are combined with explicit
integer epsilon tensors. No external density matrix or matrix exponential is
formed.
"""

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


PRIMES = (1009, 1013)


def modular_pivots(matrix, prime):
    """Exact row elimination over F_p; return pivot columns."""
    a = np.asarray(matrix, dtype=np.int64).copy() % prime
    row = 0
    pivots = []
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = a[row] * pow(int(a[row, col]), -1, prime) % prime
        for other in range(row + 1, a.shape[0]):
            if a[other, col]:
                a[other] = (a[other] - a[other, col] * a[row]) % prime
        pivots.append(col)
        row += 1
        if row == a.shape[0]:
            break
    return pivots


def permutation_sign(word):
    return (-1) ** sum(word[i] > word[j]
                      for i in range(len(word))
                      for j in range(i + 1, len(word)))


def word_index(word, d):
    answer = 0
    for entry in word:
        answer = d * answer + entry
    return answer


def matrix_unit_separation(d, rest_dimension=2):
    """Test the full linear identity on every matrix unit, not random P."""
    r = rest_dimension
    n = d * r
    total = n * d
    indices = np.arange(total).reshape(d, r, d)
    swap = np.eye(total, dtype=np.int64)[
        np.swapaxes(indices, 0, 2).reshape(-1)
    ]
    columns = []
    checked_blocks = 0
    for u, v in itertools.product(range(n), repeat=2):
        x = np.zeros((n, n), dtype=np.int64)
        x[u, v] = 1
        lifted = np.kron(x, np.eye(d, dtype=np.int64))
        commutator = lifted @ swap - swap @ lifted
        blocks = commutator.reshape(n, d, n, d)
        for q, qp in itertools.product(range(d), repeat=2):
            # S_aQ = sum_ij E_ij(a) tensor E_ji(Q).
            e = np.zeros((d, d), dtype=np.int64)
            e[qp, q] = 1
            on_a = np.kron(e, np.eye(r, dtype=np.int64))
            assert np.array_equal(blocks[:, q, :, qp], x @ on_a - on_a @ x)
            checked_blocks += 1
        columns.append(commutator.reshape(-1))
    constraints = np.stack(columns, axis=1)
    ranks = {str(p): len(modular_pivots(constraints, p)) for p in PRIMES}
    assert all(rank == n * n - r * r for rank in ranks.values())
    # An explicit r^2-dimensional kernel gives the matching characteristic-zero
    # upper bound; the modular ranks give the lower bound.
    for u, v in itertools.product(range(r), repeat=2):
        e = np.zeros((r, r), dtype=np.int64)
        e[u, v] = 1
        x = np.kron(np.eye(d, dtype=np.int64), e)
        assert not np.any(constraints @ x.reshape(-1))
    return dict(d=d, rest_dimension=r, matrix_units=n * n,
                exact_external_blocks_checked=checked_blocks,
                modular_ranks=ranks, exact_kernel_dimension=r * r,
                expected_kernel="I_a tensor M_rest")


def epsilon_seeds(d, members):
    """Integer epsilon, or every unordered product of two epsilons."""
    if members == 0:
        return np.ones((1, 1), dtype=np.int64)
    permutations = [(p, permutation_sign(p))
                    for p in itertools.permutations(range(d))]
    if members == d:
        result = np.zeros((d ** d, 1), dtype=np.int64)
        for word, sign in permutations:
            result[word_index(word, d), 0] = sign
        return result
    assert members == 2 * d
    first_blocks = [subset for subset in itertools.combinations(range(members), d)
                    if 0 in subset]
    result = np.zeros((d ** members, len(first_blocks)), dtype=np.int64)
    for column, first in enumerate(first_blocks):
        second = tuple(i for i in range(members) if i not in first)
        for (p, ps), (q, qs) in itertools.product(permutations, repeat=2):
            word = [0] * members
            for position, value in zip(first, p):
                word[position] = value
            for position, value in zip(second, q):
                word[position] = value
            result[word_index(word, d), column] = ps * qs
    return result


def collective_unit(vectors, d, members, to_level, from_level):
    if members == 0:
        return np.zeros_like(vectors)
    tensor = vectors.reshape((d,) * members + (vectors.shape[1],))
    answer = np.zeros_like(tensor)
    for position in range(members):
        source = [slice(None)] * (members + 1)
        target = source.copy()
        source[position] = from_level
        target[position] = to_level
        answer[tuple(target)] += tensor[tuple(source)]
    return answer.reshape(vectors.shape)


def singlet_certificate(d, members):
    """A zero-weight upper bound and explicit invariant lower bound meet."""
    assert members % d == 0
    counts = members // d
    words = [word for word in itertools.product(range(d), repeat=members)
             if all(word.count(level) == counts for level in range(d))]
    row_keys = {}
    entries = []
    for column, word in enumerate(words):
        for root in range(d - 1):
            for position, level in enumerate(word):
                if level != root + 1:
                    continue
                moved = list(word)
                moved[position] = root
                key = (root, tuple(moved))
                row = row_keys.setdefault(key, len(row_keys))
                entries.append((row, column))
    raising = np.zeros((len(row_keys), len(words)), dtype=np.int64)
    for row, column in entries:
        raising[row, column] += 1
    seeds = epsilon_seeds(d, members)
    coordinates = [word_index(word, d) for word in words]
    restricted = seeds[coordinates]
    assert not np.any(raising @ restricted)
    # Check every off-diagonal collective matrix unit and every simple Cartan,
    # not merely the positive-root constraints used for the upper bound.
    generator_checks = 0
    for a, b in itertools.permutations(range(d), 2):
        assert not np.any(collective_unit(seeds, d, members, a, b))
        generator_checks += 1
    for a in range(d - 1):
        diagonal = (collective_unit(seeds, d, members, a, a)
                    - collective_unit(seeds, d, members, a + 1, a + 1))
        assert not np.any(diagonal)
        generator_checks += 1
    ranks = {}
    for p in PRIMES:
        root_rank = len(modular_pivots(raising, p))
        seed_rank = len(modular_pivots(restricted, p))
        assert len(words) - root_rank == seed_rank
        ranks[str(p)] = dict(raising_rank=root_rank, seed_rank=seed_rank)
    dimension = ranks[str(PRIMES[0])]["seed_rank"]
    expected = 1 if members < 2 * d else math.comb(2 * d, d) // (d + 1)
    assert dimension == expected
    pivots = modular_pivots(restricted, PRIMES[0])
    basis, _ = np.linalg.qr(seeds[:, pivots].astype(float), mode="reduced")
    assert float(np.max(np.abs(basis.T @ basis - np.eye(dimension)))) < 2e-12
    record = dict(d=d, rest_members=members, zero_weight_dimension=len(words),
                  raising_shape=list(raising.shape), integer_seed_columns=seeds.shape[1],
                  collective_generator_checks=generator_checks,
                  modular_rank_certificates=ranks, exact_singlet_dimension=dimension)
    return record, basis


def literal_port_intertwiner(d, singlet_basis):
    """Thin isometry check; the largest array is 6561 by 45 for d=3."""
    rest_members = 2 * d
    members = rest_members + 1
    memory = singlet_basis.shape[1]
    wq = np.kron(np.eye(d), np.kron(singlet_basis, np.eye(d)))
    swapped = np.swapaxes(wq.reshape((d,) * (members + 1)
                                    + (d * memory * d,)),
                          0, members).reshape(wq.shape)
    expected = np.swapaxes(wq.reshape(wq.shape[0], d, memory, d),
                           1, 3).reshape(wq.shape)
    residual = float(np.max(np.abs(swapped - expected)))
    isometry = float(np.max(np.abs(wq.T @ wq - np.eye(d * memory * d))))
    assert max(residual, isometry) < 2e-12
    return dict(d=d, total_organization_members=members,
                full_singlet_memory_dimension=memory, thin_isometry_shape=list(wq.shape),
                physical_swap_vs_port_swap_max_residual=residual,
                isometry_max_residual=isometry)


def delocalized_code_leakage(d):
    """Exact integer projection, with a normalized full-SWAP witness."""
    members = d + 1
    physical = d ** members
    permutations = [(p, permutation_sign(p))
                    for p in itertools.permutations(range(d))]
    embeddings = np.zeros((physical, members * d), dtype=np.int64)
    for port, g in itertools.product(range(members), range(d)):
        for word, sign in permutations:
            moved = list(word)
            moved.insert(port, g)
            embeddings[word_index(moved, d), port * d + g] = (-1) ** port * sign
    # The full fundamental isotypic projector is EE^T / denominator.
    denominator = (d + 1) * math.factorial(d - 1)
    gram = embeddings.T @ embeddings
    assert np.array_equal(gram @ gram, denominator * gram)
    assert all(len(modular_pivots(embeddings, p)) == d * d for p in PRIMES)
    projector_numerator = embeddings @ embeddings.T
    initial = np.kron(embeddings, np.eye(d, dtype=np.int64))
    swapped = np.swapaxes(initial.reshape((d,) * (members + 1)
                                         + (initial.shape[1],)),
                          0, members).reshape(initial.shape)
    projected_numerator = np.einsum(
        "ij,jqc->iqc", projector_numerator,
        swapped.reshape(physical, d, initial.shape[1]), optimize=True
    ).reshape(initial.shape)
    leakage_numerator = denominator * swapped - projected_numerator
    squared = np.sum(leakage_numerator * leakage_numerator, axis=0)
    best = int(np.argmax(squared))
    probability = Fraction(int(squared[best]),
                           denominator ** 2 * math.factorial(d))
    assert 0 < probability <= 1
    embedding_column, external_level = divmod(best, d)
    embedded_port, input_level = divmod(embedding_column, d)
    return dict(d=d, old_organization_members=members,
                full_fundamental_code_dimension=d * d,
                raw_contact_member=0, input_epsilon_embedding_port=embedded_port,
                input_fundamental_level=input_level, external_input_level=external_level,
                normalized_full_swap_leakage_probability=str(probability),
                leakage_probability_float=float(probability), arithmetic="exact integers and Fraction")


def run():
    for p in PRIMES:
        assert all(p % divisor for divisor in range(2, math.isqrt(p) + 1))
    commutants = [matrix_unit_separation(d) for d in (2, 3)]
    singlets = []
    counts = []
    literal = []
    for d in (2, 3):
        dimensions = []
        largest_basis = None
        for members in range(2 * d + 1):
            if members % d:
                dimensions.append(0)  # central element exp(2 pi i / d)
                continue
            certificate, basis = singlet_certificate(d, members)
            singlets.append(certificate)
            dimensions.append(certificate["exact_singlet_dimension"])
            largest_basis = basis
        assert all(x < 2 for x in dimensions[:-1]) and dimensions[-1] >= 2
        counts.append(dict(d=d, rest_member_range=[0, 2 * d],
                           exact_singlet_dimensions=dimensions,
                           first_nontrivial_memory_rest_members=2 * d,
                           minimum_literal_port_organization_members=2 * d + 1,
                           excluded_nonmultiple_reason="nontrivial SU(d) center character"))
        literal.append(literal_port_intertwiner(d, largest_basis))
    return dict(
        scope="Independent finite calibration of the fixed-code literal-member contract; no author imports",
        matrix_unit_commutants=commutants,
        singlet_certificates=singlets,
        minimal_member_counts=counts,
        literal_port_intertwiners=literal,
        old_delocalized_code_witnesses=[delocalized_code_leakage(d) for d in (2, 3)],
        passed=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    answer = run()
    if args.write:
        destination = Path(__file__).with_name("independent_results.json")
        with destination.open("x", encoding="utf-8") as handle:
            json.dump(answer, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    print(json.dumps(answer, ensure_ascii=False))
