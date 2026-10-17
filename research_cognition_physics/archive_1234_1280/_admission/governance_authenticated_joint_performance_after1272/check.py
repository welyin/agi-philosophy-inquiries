"""Exact finite diagnostics; default compares results without modifying files.

The coherent/adaptive security theorem is imported in proof.md, not established
by this Pauli enumeration. No physical energy or experimental error is inferred.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent

def fstr(x):
    return f'{x.numerator}/{x.denominator}'

def symplectic(a, b, q):
    mask = (1 << q) - 1
    ax, az = a & mask, a >> q
    bx, bz = b & mask, b >> q
    return ((ax & bz).bit_count() + (az & bx).bit_count()) % 2

def sp4_images():
    counts = [[0] * 16 for _ in range(16)]
    number = 0
    for x1 in range(1, 16):
        for z1 in range(1, 16):
            if symplectic(x1, z1, 2) != 1:
                continue
            perp = [v for v in range(1, 16)
                    if symplectic(x1, v, 2) == symplectic(z1, v, 2) == 0]
            for x2 in perp:
                for z2 in perp:
                    if symplectic(x2, z2, 2) != 1:
                        continue
                    number += 1
                    basis = [x1, x2, z1, z2]
                    images = []
                    for v in range(1, 16):
                        w = 0
                        for j, u in enumerate(basis):
                            if (v >> j) & 1:
                                w ^= u
                        images.append(w)
                        counts[v][w] += 1
                    assert len(set(images)) == 15 and 0 not in images
    assert number == 720
    assert all(counts[v][w] == 48 for v in range(1, 16) for w in range(1, 16))
    return {'symplectic_bases': number, 'nonzero_inputs': 15, 'each_image_count': 48,
            'clifford_classes_mod_phase': number * 16}

def classify(m, n):
    counts = {'reject': 0, 'benign_accept': 0, 'harmful_accept': 0}
    # I=0, X=1, Y=2, Z=3; all identities excluded.
    for word in product(range(4), repeat=m+n):
        if not any(word):
            continue
        if any(p in (1, 2) for p in word[m:]):
            counts['reject'] += 1
        elif any(word[:m]):
            counts['harmful_accept'] += 1
        else:
            counts['benign_accept'] += 1
    total = 4**(m+n) - 1
    assert sum(counts.values()) == total
    assert counts['benign_accept'] == 2**n - 1
    assert counts['harmful_accept'] == (4**m - 1) * 2**n
    accept = counts['benign_accept'] + counts['harmful_accept']
    assert Fraction(accept, total) <= Fraction(1, 2**n)
    return dict(m=m, n=n, total=total, **counts,
                wrong_accept_weight=fstr(Fraction(counts['harmful_accept'], total)),
                conditional_wrong_given_accept_for_bell_input=fstr(Fraction(counts['harmful_accept'], accept)))

def clifford_order(q):
    answer = 2**(q*q+2*q)
    for j in range(1, q+1):
        answer *= 4**j - 1
    return answer

def compute():
    m, n, r = 2, 12, 4
    q = m+n
    key_choices = clifford_order(q) * 4**q
    capacity = (key_choices - 1).bit_length()
    assert 2**(capacity-1) < key_choices <= 2**capacity
    pure_error = Fraction(r, 2**n)
    assert pure_error == Fraction(1, 1024)
    return {'status': 'PASS_EXACT_DIAGNOSTICS_AND_CONDITIONAL_BUDGET',
            'new_scientific_rounds': 0,
            'sp4_complete_enumeration': sp4_images(),
            'pauli_cases': [classify(1, 1), classify(1, 3), classify(2, 2)],
            'joint_contract_budget': {'message_qubits_per_slot': m, 'tag_qubits_per_slot': n,
                'slots': r, 'sent_code_qubits': r*q, 'tag_preparations': r*n,
                'key_choices_per_slot': str(key_choices), 'capacity_bits_per_key_register': capacity,
                'three_retained_key_register_capacity_bits': 3*r*capacity,
                'composable_error_ideal_model': fstr(pure_error),
                'error_if_total_preparation_and_implementation_distance_le_1_over_4096': fstr(pure_error+Fraction(1,4096))},
            'deleted_authentication': {'bell_payload_pauli_attack_wrong_accept': '1/1'},
            'scope': ['Pauli enumeration is not a coherent-attack security proof.',
                      'Full security uses Portmann Appendix D.1 in the declared attack/resource domain.',
                      'No key recycling, no empirical error, no physical energy ceiling.',
                      'Classical facts, source honesty, identity uniqueness and gravity are not proved.']}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    result = compute()
    destination = HERE / 'results.json'
    if args.record:
        with destination.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        assert json.loads(destination.read_text(encoding='utf-8')) == result
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()
