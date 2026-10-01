"""587 entry: full kinetic/potential incidence audit, not a finished round.

The D_A W_Z values are generic input coefficients here, not simulated gauge
field values. This checks allocation and support, not physical continuum flux.
"""
from pathlib import Path
import argparse
import json
import numpy as np

TARGET = Path(__file__).with_name('local_current_entry_results.json')


def run():
    n = 4
    sites = list(np.ndindex((n, n, n)))
    index = {p: i for i, p in enumerate(sites)}
    def shifted(p, axis):
        q = list(p); q[axis] = (q[axis]+1) % n
        return tuple(q)
    edges = {}
    for p in sites:
        for axis in range(3):
            edges[p, axis] = (index[p], index[shifted(p, axis)])
    # Each entry contains allocation of a kinetic factor, a potential and D_A W_Z.
    terms = []
    rng = np.random.default_rng(587)
    for p in sites:
        i = index[p]
        terms.append(('onsite', {i: 1.}, {i: 1.}, float(rng.normal())))
        for axis in range(3):
            i, j = edges[p, axis]
            b = {i: .5, j: .5}
            terms.append(('node_edge', {i: 1.}, b, float(rng.normal())))
            terms.append(('node_edge', {j: 1.}, b, float(rng.normal())))
            terms.append(('electric_edge', b, b, float(rng.normal())))
        for mu in range(3):
            for nu in range(mu+1, 3):
                pm, pn = shifted(p, mu), shifted(p, nu)
                corners = (index[p], index[pm], index[pn], index[shifted(pm, nu)])
                b = {v: .25 for v in corners}
                for edge in ((p, mu), (pm, nu), (pn, mu), (p, nu)):
                    i, j = edges[edge]
                    terms.append(('electric_face', {i: .5, j: .5}, b, float(rng.normal())))
    currents = np.zeros((len(sites), len(sites)))
    direct_rate = np.zeros(len(sites))
    cancellation = dict(onsite=0., electric_edge=0.)
    max_distance = 0
    for kind, a, b, derivative in terms:
        support = set(a) | set(b)
        for i in support:
            direct_rate[i] += (b.get(i, 0)-a.get(i, 0))*derivative
            for j in support:
                term = (a.get(j, 0)*b.get(i, 0)-a.get(i, 0)*b.get(j, 0))*derivative
                currents[i, j] += term
                if kind in cancellation:
                    cancellation[kind] = max(cancellation[kind], abs(term))
                if abs(term) > 0:
                    displacement = abs(np.array(sites[i])-np.array(sites[j]))
                    max_distance = max(max_distance, int(np.minimum(displacement, n-displacement).sum()))
    antisymmetry = float(abs(currents+currents.T).max())
    balance = float(abs(currents.sum(axis=1)-direct_rate).max())
    assert antisymmetry < 1e-14 and balance < 1e-13
    assert max_distance == 2 and all(v == 0 for v in cancellation.values())
    # Conservation alone does not eliminate circulation around an elementary face.
    ring = [index[p] for p in ((0,0,0), (1,0,0), (1,1,0), (0,1,0))]
    curl = np.zeros_like(currents)
    for i, j in zip(ring, ring[1:]+ring[:1]):
        curl[j, i] += .17
        curl[i, j] -= .17
    assert np.max(abs(curl.sum(axis=1))) == 0
    return dict(status='587 entry only; not a completed scientific round',
                size=n, vertices=len(sites), derivative_incidence_terms=len(terms),
                antisymmetry_error=antisymmetry, full_allocation_balance_error=balance,
                max_nonzero_current_taxicab_distance=max_distance,
                exact_pair_current_cancellations=cancellation,
                nonzero_divergence_free_cycle_norm=float(np.linalg.norm(curl)),
                gauge_differential_and_continuum_flux_not_tested=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result, ensure_ascii=False))
