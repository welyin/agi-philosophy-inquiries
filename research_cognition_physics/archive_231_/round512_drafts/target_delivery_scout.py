"""Unfrozen numerical scout for a possible round 512; no theorem certificate.

Compare an already existing passive drift with the already existing round 503
identity instrument. Eigenspace tolerances are diagnostic, not exact evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(BASE))
import coherent_graph_mean_obstruction as passive
import uniform_identity_receipt as identity


def diagnose(h, n, graphs, colours, target):
    values, vectors = np.linalg.eigh(h)
    blocks = []
    for j, value in enumerate(values):
        if not blocks or abs(value-values[blocks[-1][0]]) > 1e-8:
            blocks.append([])
        blocks[-1].append(j)
    dark = np.eye(len(h))
    lo, hi = target*graphs*colours, (target+1)*graphs*colours
    for block in blocks:
        eig = vectors[:, block]
        _, singular, vh = np.linalg.svd(eig[lo:hi], full_matrices=False)
        bright = eig @ vh[singular > 1e-8].T
        dark -= bright @ bright.T
    rows = []
    for start in range(n):
        source = (start*graphs+np.arange(graphs))*colours
        overlap = np.linalg.eigvalsh(dark[np.ix_(source, source)])
        rows.append(dict(start=start, maximum_query_source_dark_weight=float(overlap[-1])))
    return dict(target=target, numerical_dark_rank=float(np.trace(dark)), sources=rows,
                idempotence_residual=float(np.linalg.norm(dark @ dark-dark)),
                invariance_residual=float(np.linalg.norm(h @ dark-dark @ h)))


def run():
    trees, f, _, _, h, _, _, _ = passive.model(2)
    tagged = identity.model(trees, f, 6)
    return dict(status='exploratory_unfrozen_not_a_completed_scientific_round',
        complete_baseline_round=511, completed_check_count_unchanged=2500,
        problem='first arrival at another fixed target, not monitored return to the source',
        diagnostic_eigenvalue_and_singular_tolerance=1e-8,
        common_vertex_roles=dict(internal=[0,1], leaves=[2,3,4,5]),
        both_models_share_exact_tree_order_and_flip=True,
        existing_passive_model=[diagnose(h,6,len(trees),1,a) for a in (0,2)],
        existing_identity_model=[diagnose(tagged,6,len(trees),7,a) for a in (0,2)],
        exact_membership_in_observable_span_proved=False,
        finite_delivery_budget_proved=False, dimension_or_GR_proved=False,
        source_sha256={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py','uniform_identity_receipt.py',
            'research_note_502.md','research_note_503.md')})


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with (HERE/'target_delivery_aligned_scout_results.json').open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
