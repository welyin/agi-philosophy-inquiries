"""712 entry: unique one-generation, undressed Q^3 L scalar tensor.

Exact integer characters give multiplicities; NumPy checks the original
710 CAR tensor, not a new physical dynamics or an instanton amplitude.
"""
import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from math import factorial
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_charge_changing_vertex as original

TARGET = HERE / 'minimal_left_vertex_entry_results.json'
OCC = [n for n in product(range(4), repeat=4) if sum(n) == 3]
INDEX = {n: i for i, n in enumerate(OCC)}
WEIGHTS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def character(j2, k2):
    return Counter({(j, k): 1 for j in range(-j2, j2+1, 2) for k in range(-k2, k2+1, 2)})


def sym_generator(a):
    result = np.zeros((20, 20), complex)
    for col, n in enumerate(OCC):
        for j in range(4):
            if not n[j]:
                continue
            result[col, col] += a[j, j]*n[j]
            for i in range(4):
                if i == j:
                    continue
                m = list(n)
                m[j] -= 1
                m[i] += 1
                result[INDEX[tuple(m)], col] += a[i, j]*np.sqrt(n[j]*(n[i]+1))
    return result


def run():
    sym_char = Counter(tuple(sum(n[i]*WEIGHTS[i][s] for i in range(4)) for s in range(2)) for n in OCC)
    assert sym_char == character(3, 3) + character(1, 1)
    full_char = Counter()
    for (a, b), multiplicity in sym_char.items():
        for x, y in WEIGHTS:
            full_char[a+x, b+y] += multiplicity
    decomposition = [(4, 4, 1), (4, 2, 1), (2, 4, 1), (2, 2, 2), (2, 0, 1), (0, 2, 1), (0, 0, 1)]
    expected_char = Counter()
    for j2, k2, multiplicity in decomposition:
        for weight in character(j2, k2):
            expected_char[weight] += multiplicity
    assert full_char == expected_char and sum(full_char.values()) == 80

    generators = []
    for weak in (True, False):
        for p in PAULI:
            a = np.kron(p/2, np.eye(2)) if weak else np.kron(np.eye(2), p/2)
            generators.append(np.kron(sym_generator(a), np.eye(4)) + np.kron(np.eye(20), a))
    casimir = sum(t@t for t in generators)
    eig, vectors = np.linalg.eigh(casimir)
    multiplicities = Counter(int(round(x)) for x in eig)
    expected = {0: 1, 2: 6, 4: 18, 8: 30, 12: 25}
    assert dict(multiplicities) == expected
    spectral_error = float(max(abs(eig-np.rint(eig))))
    assert spectral_error < 1e-12

    # One quark in each color: the unique color-singlet embedding is symmetric
    # in the weak x spin assignments because the CAR/color signs cancel.
    assignments = list(product(range(4), repeat=4))
    embedding = np.zeros((256, 80))
    original_vector = np.zeros(256)
    for row, (u0, u1, u2, ell) in enumerate(assignments):
        n = tuple((u0, u1, u2).count(i) for i in range(4))
        multiplicity = factorial(3)
        for count in n:
            multiplicity //= factorial(count)
        embedding[row, 4*INDEX[n]+ell] = 1/np.sqrt(multiplicity)
        ids = (u0, 4+u1, 8+u2, 24+ell)
        original_vector[row] = original.COEFF.get(ids, 0)/np.sqrt(original.NORM2)
    isometry_error = float(np.linalg.norm(embedding.T@embedding-np.eye(80)))
    projected = embedding.T@original_vector
    reconstruction_error = float(np.linalg.norm(embedding@projected-original_vector))
    generator_error = float(max(np.linalg.norm(t@projected) for t in generators))
    overlap = float(abs(np.vdot(vectors[:, 0], projected)))
    assert max(isometry_error, reconstruction_error, generator_error, abs(overlap-1)) < 1e-12
    assert abs(np.linalg.norm(projected)-1) < 1e-12 and np.count_nonzero(original_vector) == 36
    assert 3*1-3 == 0  # Original integer U(1) charges of Q and L.
    names = ('research_note_605.md', 'research_note_606.md', 'research_note_607.md',
             'research_note_608.md', 'research_note_609.md', 'research_note_629.md',
             'research_note_710.md', 'research_note_711.md', 'joint_charge_changing_vertex.py')
    return dict(entry_round=712, new_formal_round=False,
        exact_character_identity=True, symmetric_cube_dimension=20, color_singlet_sector_dimension=80,
        irreps_twice_spin=decomposition, invariant_multiplicity=1,
        casimir_multiplicities={str(k): v for k, v in sorted(multiplicities.items())},
        casimir_spectral_error=spectral_error, original_CAR_modes=32,
        original_tensor_monomials=int(np.count_nonzero(original_vector)), raw_norm_squared=original.NORM2,
        embedding_isometry_error=isometry_error, original_tensor_reconstruction_error=reconstruction_error,
        total_generator_error=generator_error, normalized_unique_singlet_overlap=overlap,
        scope='One generation, Q^3 L, no derivatives or Higgs insertions, holomorphic left-spin scalar. Coefficient, generations, dynamics and topological origin remain inputs/open.',
        dependencies={n: hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    # JSON canonicalization makes tuple-valued exact representation data reproducible.
    result = json.loads(json.dumps(result))
    if args.write_results:
        with TARGET.open('x', encoding='utf8') as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result, ensure_ascii=False, indent=2))
