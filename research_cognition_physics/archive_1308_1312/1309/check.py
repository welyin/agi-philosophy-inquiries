"""Round 1309: rational certificates; no history writes or numerical scans as proofs."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path


def beta(nc, ng):
    # C_A=2, each Weyl doublet has T=1/2; one complex Higgs doublet.
    return -F(11, 3) * 2 + F(2, 3) * ng * (nc + 1) * F(1, 2) + F(1, 3) * F(1, 2)


def candidates(budget):
    # Exact analytic upper bounds; not an arbitrary scan cutoff.
    product_cap = (43 + 6 * budget) / 2
    answer = []
    for k in range(1, int(product_cap // 12) + 1):
        for ng in range(3, int(product_cap // (4 * k)) + 1):
            answer.append([4 * k - 1, ng])
    return answer


def fr(x):
    return str(x)


def calculate():
    base = [[3, 3], [3, 4], [3, 5]]
    rows = []
    for b in [F(1, 100), F(1, 3), F(5, 6), F(7, 6), F(3), F(8)]:
        got = candidates(b)
        # Second route: each possible N and n constrained separately by n>=3,N>=3.
        cap = (43 + 6 * b) / 2
        independent = [[nc, ng] for nc in range(3, int(cap // 3))
                       for ng in range(3, int(cap // 4) + 1)
                       if nc % 4 == 3 and beta(nc, ng) <= b]
        assert got == independent
        for nc, ng in got:
            assert beta(nc, ng) == F(-43 + 2 * ng * (nc + 1), 6)
        rows.append({'B': fr(b), 'candidates': got})
    assert candidates(F(1, 3)) == base
    assert candidates(F(5, 6)) == base + [[3, 6], [7, 3]]
    assert beta(3, 6) == beta(7, 3) == F(5, 6)

    B, eta = F(1, 3), F(1, 4)
    robust, possible = candidates(B - eta), candidates(B + eta)
    assert robust == possible == base
    uncertainty = []
    for nc, ng in base + [[3, 6], [7, 3]]:
        b = beta(nc, ng)
        # s=t/L, y=8pi²(x-x_floor)/L; y(0)=B, y'= -b+u, |u|<=eta.
        low_end, high_end = B - b - eta, B - b + eta
        robust_pass = min(B, low_end) >= 0
        exists_pass = min(B, high_end) >= 0
        assert robust_pass == (b <= B - eta)
        assert exists_pass == (b <= B + eta)
        uncertainty.append({'Nc': nc, 'ng': ng, 'b': fr(b),
                            'worst_endpoint_y': fr(low_end),
                            'best_endpoint_y': fr(high_end),
                            'robust': robust_pass, 'possible': exists_pass})

    # Tight logical boundaries: endpoint zero is legal; the next rational perturbation is not.
    b = beta(7, 3)
    edge_B = b + eta
    assert min(edge_B, edge_B - b - eta) == 0
    assert min(edge_B - F(1, 1000), edge_B - F(1, 1000) - b - eta) < 0
    # A gap between necessary and robust sufficient is real.
    gap_B = b
    assert b <= gap_B + eta and not b <= gap_B - eta

    no_floor = []
    for nc, ng in [(3, 3), (7, 3), (11, 17), (99, 1000)]:
        b = beta(nc, ng)
        y0 = F(1) + max(b + eta, 0)
        worst_end = y0 - b - eta
        assert min(y0, worst_end) >= 1
        no_floor.append({'Nc': nc, 'ng': ng, 'b': fr(b), 'y0': fr(y0),
                         'uniform_lower_y': fr(min(y0, worst_end))})

    # Piecewise, nonconstant errors: exact affine segments attain minima at their endpoints.
    segments = [(F(1, 5), eta), (F(2, 5), -eta), (F(2, 5), eta / 2)]
    piecewise = []
    for nc, ng in base:
        y, minimum = B, B
        for length, u in segments:
            y += (-beta(nc, ng) + u) * length
            minimum = min(minimum, y)
        assert minimum >= 0
        piecewise.append({'Nc': nc, 'ng': ng, 'minimum_y': fr(minimum), 'end_y': fr(y)})
    return {
        'round': 1309, 'status': 'PASS', 'new_axioms': 0,
        'scope': 'conditional finite-window leading running with a separately supplied error envelope',
        'exact_budget_classification': rows,
        'error_budget': {'B': fr(B), 'eta': fr(eta), 'robust': robust,
                         'possible': possible, 'rows': uncertainty},
        'tight_boundary_and_quantifier_checks': True,
        'nonconstant_error_checks': piecewise,
        'without_uniform_initial_floor': no_floor,
        'general_claim_source': 'proof.md analytic inequalities, not finite samples',
        'not_certified': ['physical origin of coupling floor', 'actual higher-order rho',
                          'full interacting instrument', 'unique three generations', 'full SM generation'],
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = Path(__file__).with_name('results.json')
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        assert json.loads(path.read_text(encoding='utf-8')) == result, 'saved result mismatch'
    print(json.dumps({'round': 1309, 'status': 'PASS', 'cases': result['error_budget']['robust'],
                      'critical_B': '5/6', 'mode': 'first-write' if args.write else 'read-only'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
