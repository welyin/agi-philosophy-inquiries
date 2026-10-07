"""Unnumbered scope check: a returned packet is not a persistent remote record.

No new research-round credit. This checks a direct consequence of the frozen
503 one-packet carrier and does not claim that all remote information vanishes.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import uniform_identity_receipt as frozen503
import current_neighbor_detection as evolution

HERE=Path(__file__).resolve().parent
TARGET=HERE/'bilateral_receipt_scope_probe_results.json'


def run():
    checks=[]
    for n in [2,3,4]:
        d=n+2
        words=list(itertools.product(range(d),repeat=n))
        packet=np.array([sum(a!=0 for a in word)==1 for word in words])
        for root in range(n):
            for label in range(n):
                if label==root:continue
                # Local root projection intersected with the whole one-packet code.
                kept=packet & np.array([word[root]==label+2 for word in words])
                expected=[0]*n;expected[root]=label+2
                assert np.flatnonzero(kept).tolist()==[words.index(tuple(expected))]
        checks.append(dict(N=n,local_packet_dimension=d,
            tested_root_label_pairs=n*(n-1),rank_one_carrier_factorization_exact=True))

    trees=[frozenset({(0,1)})]
    h=frozen503.model(trees,np.zeros((1,1),dtype=np.int64),2)
    source,success,_,_=frozen503.root_interface(trees,2,0)
    u=evolution.evolution(h,Fraction(1,8))
    state=u[:,source[0]]
    p=float(np.linalg.norm(state[success])**2)
    assert p>0
    selected=np.zeros_like(state);selected[success]=state[success]/np.sqrt(p)
    # Trace every root colour: at port zero all remote packet coordinates are vacuum.
    remote=np.zeros((4,4),complex)
    for colour in range(3):remote[0,0]+=abs(selected[colour])**2
    for a in range(3):
        for b in range(3):remote[a+1,b+1]+=selected[3+a]*selected[3+b].conjugate()
    idle=np.diag([1,0,0,0]).astype(complex)
    distance=float(np.linalg.svd(remote-idle,compute_uv=False).sum()/2)
    assert distance<1e-14
    remote_inflight=float(np.linalg.norm(state[3:])**2)
    assert remote_inflight>0.001
    return dict(numbered_round=None,scientific_baseline_round=503,
        diagnostic_checks=2,numbered_scientific_test_count_increment=0,
        operator_checks=checks,
        static_two_node=dict(time='1/8',success_probability=p,
            selected_remote_packet_vs_idle_trace_distance=distance,
            unconditioned_remote_nonvacuum_probability=remote_inflight),
        scope=dict(remote_packet_and_fixed_identity_only=True,
            general_remote_graph_or_old_memory_claim=False,
            zero_communication_claim=False,persistent_peer_receipt_automatic=False,
            complete_cognitive_countermodel=False),
        frozen_dependency_sha256=hashlib.sha256((HERE/'uniform_identity_receipt.py').read_bytes()).hexdigest(),
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
