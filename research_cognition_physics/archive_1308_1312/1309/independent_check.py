"""Independent finite-window audit: divisor classification and extremal controls.

No imports from the author's checker. Exact arithmetic tests finite diagnostic
instances; the all-integer and measurable-control proofs are in the review.
"""

import argparse
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def beta(nc, ng):
    return F(8 * ng * ((nc + 1) // 4) - 43, 6)


def divisor_candidates(bound):
    """m=k*n; enumerate the allowed divisors of each product m."""
    top = (43 + 6 * bound) // 8
    answer = set()
    for m in range(3, top + 1):
        for ng in range(3, m + 1):
            if m % ng == 0:
                k = m // ng
                answer.add((4 * k - 1, ng))
    return sorted(answer)


def extrema(initial, b, controls, lengths):
    y = F(initial)
    nodes = [y]
    for u, dt in zip(controls, lengths):
        y += (u - b) * dt
        nodes.append(y)
    return min(nodes), nodes


def calculate():
    source = json.loads((HERE / 'results.json').read_text(encoding='utf-8'))
    classification = []
    for row in source['exact_budget_classification']:
        bound = F(row['B'])
        values = divisor_candidates(bound)
        assert [list(p) for p in values] == row['candidates']
        classification.append({'B': str(bound), 'divisor_pairs': [list(p) for p in values]})

    base = [(3, 3), (3, 4), (3, 5)]
    first = F(5, 6)
    assert divisor_candidates(first - F(1, 1000000)) == base
    assert divisor_candidates(first) == base + [(3, 6), (7, 3)]
    assert divisor_candidates(first + F(1, 1000000)) == base + [(3, 6), (7, 3)]
    assert beta(3, 6) == beta(7, 3) == first
    # The next factor-product threshold is m=7, hence b=13/6.
    assert divisor_candidates(F(13, 6)) == base + [(3, 6), (3, 7), (7, 3)]

    # Only corners are needed for each finite partition: every nodal y is
    # affine in all control values. The analytic envelope extends to all
    # bounded measurable errors and every time, not just these samples.
    partitions = [
        [F(1)],
        [F(1, 4), F(3, 4)],
        [F(1, 7), F(2, 7), F(4, 7)],
        [F(1, 10), F(1, 5), F(3, 10), F(2, 5)],
    ]
    cases = [
        (F(1, 3), F(1, 4)),
        (F(5, 6), F(1, 4)),
        (F(1), F(1, 12)),
        (F(1, 10), F(2)),
    ]
    controls_checked = 0
    records = []
    for budget, eta in cases:
        for nc, ng in base + [(3, 6), (7, 3), (11, 3)]:
            b = beta(nc, ng)
            for lengths in partitions:
                mins = []
                for signs in itertools.product((-1, 1), repeat=len(lengths)):
                    low, nodes = extrema(budget, b, [s * eta for s in signs], lengths)
                    mins.append(low)
                    controls_checked += 1
                worst_min = min(mins)
                best_min = max(mins)
                assert (worst_min >= 0) == (b <= budget - eta)
                assert (best_min >= 0) == (b <= budget + eta)
            records.append({'B': str(budget), 'eta': str(eta), 'Nc': nc, 'ng': ng,
                            'robust': b <= budget - eta,
                            'envelope_possible': b <= budget + eta})

    # Finite admissibility really is weaker than non-increasing g:
    # b=5/6>0 passes B=1 even with eta=1/12 against every allowed error.
    assert first > 0 and first <= F(1) - F(1, 12)

    # Endpoint admissibility alone is not sufficient for a particular error.
    interior_min, nodes = extrema(F(1, 10), first, [-F(2), F(2)],
                                  [F(1, 4), F(3, 4)])
    assert nodes[0] > 0 and nodes[-1] > 0 and interior_min < 0

    # The old three branches are all robust under exactly the advertised
    # two-sided margin; omitted boundary cases show why both sides matter.
    budget, eta = F(1, 3), F(1, 4)
    assert divisor_candidates(budget - eta) == base
    assert divisor_candidates(budget + eta) == base
    budget, eta = F(1, 10), F(3, 5)
    assert budget - eta == -F(1, 2)
    assert divisor_candidates(budget - eta) == base  # equality survives
    assert divisor_candidates(budget - eta - F(1, 1000)) == base[:2]
    # Equality at the upper exclusion margin already admits competitors.
    assert divisor_candidates(F(1, 3) + F(1, 2)) == base + [(3, 6), (7, 3)]

    no_floor = []
    for nc, ng in [(3, 3), (7, 3), (11, 17), (99, 1000)]:
        b, eta = beta(nc, ng), F(2, 7)
        y0 = 1 + max(b + eta, 0)
        lower = min(y0, y0 - b - eta)
        assert lower >= 1
        no_floor.append({'Nc': nc, 'ng': ng, 'y0': str(y0), 'lower': str(lower)})

    inherited = [
        HERE / 'proof.md', HERE / 'check.py', HERE / 'results.json',
        HERE.parent / '1308' / 'proof.md',
        ROOT / 'archive_531_553' / 'research_note_553.md',
        ROOT / 'archive_585_628' / 'research_note_604.md',
        ROOT / 'archive_1290_1307' / '1303' / 'proof.md',
    ]
    bindings = {str(p.relative_to(ROOT)).replace('\\', '/'):
                hashlib.sha256(p.read_bytes()).hexdigest() for p in inherited}
    return {
        'round': 1309, 'review_status': 'PASS_WITH_DECLARED_SCOPE',
        'method': 'divisor-product classification and exact control-polytope vertices',
        'classification': classification,
        'control_sequences_checked': controls_checked, 'quantifier_diagnostics': records,
        'positive_beta_robust_example': {'Nc': 7, 'ng': 3, 'B': '1', 'eta': '1/12'},
        'endpoint_only_counterexample': {'nodes': [str(z) for z in nodes],
                                         'interior_minimum': str(interior_min)},
        'no_uniform_initial_floor': no_floor,
        'source_sha256': bindings,
        'limits': ['analytic proof supplies general quantifiers',
                   'rho is an additional certified envelope, not a high-loop calculation',
                   'initial coupling floor is not a whole-window floor',
                   'no full interacting instrument or unique three generations'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE / 'independent_results.json'
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        assert json.loads(path.read_text(encoding='utf-8')) == result
    print(json.dumps({'status': result['review_status'],
                      'control_sequences': result['control_sequences_checked'],
                      'mode': 'first-write' if args.write else 'read-only'},
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
