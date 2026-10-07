"""Unnumbered candidate after 500: a locally readable returned label.

This changes the data carrier and adds a local identity-writing interaction.
It does NOT implement original-qubit H or inherit 500's retry guarantee.
Only exact finite identities and a conservative analytic budget are checked.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import current_neighbor_detection as frozen500


def report():
    trees, h, _, flip = frozen500.six_model()
    n = len(trees)
    colors = 6  # query plus the five non-root reply labels in this sector
    hh = np.kron(h, np.eye(colors, dtype=np.int64))
    for vertex in range(1, n):
        for graph in range(n):
            query = (vertex*n+graph)*colors
            hh[query, query+vertex] += 1
            hh[query+vertex, query] += 1
    ids = np.arange(n)*colors
    p = np.eye(len(hh), dtype=np.int64)[:, ids]
    reply = np.zeros(len(hh), bool)
    bad = reply.copy()
    for graph, tree in enumerate(trees):
        for vertex in range(1, n):
            k = graph*colors+vertex
            reply[k] = True
            bad[k] = (0, vertex) not in tree
    power = p.copy()
    rows = []
    for order in range(5):
        rows.append(dict(order=order, return_max=int(np.max(np.abs(power[reply]))),
                         bad_max=int(np.max(np.abs(power[bad])))))
        if order < 3:
            assert not np.any(power[reply])
        if order == 3:
            expected = np.zeros_like(power)
            for graph, tree in enumerate(trees):
                for vertex in range(1, n):
                    expected[graph*colors+vertex, graph] = int((0, vertex) in tree)
            assert np.array_equal(power*reply[:, None], expected)
            assert not np.any(power[bad])
        power = hh@power
    graph_h = np.kron(np.kron(np.eye(n, dtype=np.int64), flip),
                      np.eye(colors, dtype=np.int64))
    comm = graph_h@hh-hh@graph_h
    # Proposed general Dyson budget, distinct from finite matrix identities.
    c, j, g, kappa = 3, 1, 1, 1
    f = 2*c*(c-1)**2
    v = 2*abs(j)*c+abs(g)
    lip = 2*f*abs(kappa*j)
    r = F(v*v*lip, 4)+F(v**4, 12)
    a0 = abs(j*j*g)
    t = F(1, 131072)
    assert v*t <= 1 and 12*r*t <= a0
    success_lower = a0*a0*t**6/144
    error_upper = (12*r/a0+2*abs(kappa)*f)**2*t*t
    assert error_upper < F(1, 100)
    return dict(kind='unnumbered_candidate_diagnostic', frozen_baseline=500,
        dependency_result_sha256=hashlib.sha256(
            (HERE.parent/'current_neighbor_detection_results.json').read_bytes()).hexdigest(),
        dimension=len(hh), exact_low_order_blocks=rows,
        exact_commutator_row_sum=int(max(np.sum(abs(comm), axis=1))),
        proposed_uniform_budget=dict(C=c, v=v, lipschitz=lip, remainder_coefficient=str(r),
            time=str(t), success_lower=str(success_lower), conditional_bad_upper=str(error_upper)),
        scope=dict(new_multilevel_data_carrier=True, new_local_identity_coupling=True,
            readout_only_at_root=True, complete_analytic_note_pending=True,
            same_graph_retry_proved=False, original_qubit_generator_implemented=False,
            autonomously_generated_identity_reference=False, physical_positions_derived=False,
            formal_round_created=False, full_GR_goal_completed=False))


if __name__ == '__main__':
    result = report()
    if '--write-results' in sys.argv:
        with (HERE/'local_return_probe_candidate_results.json').open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
