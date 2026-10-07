"""Working876: compression cross terms and a graph-norm bridge.
Reuses706's original finite diagnostic, not a full gauge-field computation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'local_graph_norm_probe_results.json'
def run():
    with ResearchRuntime(Layout()).installed():
        import joint_local_history_source_limit as old
        d=old.data()
    T,W=d['parts'][:2];S=T+W+np.eye(len(T))
    rows=[]
    for k in range(1,d['r']+1):
        rank=4*k
        ix=np.array([a*d['nloc']+b for a in range(rank) for b in range(rank)])
        q=np.setdiff1d(np.arange(len(T)),ix)
        tr=T[np.ix_(ix,ix)];wr=W[np.ix_(ix,ix)]
        tail=W[np.ix_(q,ix)]
        correction=2*tail.conj().T@tail
        defect=(tr@wr+wr@tr)-(T@W+W@T)[np.ix_(ix,ix)]
        identity_error=float(np.linalg.norm(defect-correction))
        reference_leak=float(np.linalg.norm(S[np.ix_(q,ix)]))
        assert identity_error<1e-11 and reference_leak<1e-11
        assert np.linalg.eigvalsh(correction).min()>-1e-12
        graph=[];source=[]
        for sigma in (-.08,0.,.08):
            H=sum(np.exp(power*sigma)*part for power,part in zip(d['powers'],d['parts']))
            G=sum(power*np.exp(power*sigma)*part for power,part in zip(d['powers'],d['parts']))
            shift=max(1.,1-float(np.linalg.eigvalsh(H).min()))
            A=H[np.ix_(ix,ix)]+shift*np.eye(len(ix))
            ai=np.linalg.inv(A)
            graph.append(float(np.linalg.norm(S[np.ix_(ix,ix)]@ai,2)))
            source.append(float(np.linalg.norm(G[np.ix_(ix,ix)]@ai,2)))
        rows.append(dict(bands=k,dimension=len(ix),
            positive_cross_term_identity_error=identity_error,
            comparison_operator_cutoff_leakage=reference_leak,
            separate_potential_cutoff_leakage=float(np.linalg.norm(tail)),
            correction_norm=float(np.linalg.norm(correction)),
            graph_inverse_norms_at_three_geometries=graph,
            source_inverse_norms_at_three_geometries=source))
    assert max(x['separate_potential_cutoff_leakage'] for x in rows)>1e-5
    return dict(working_round=876,formal_reports=875,cumulative_numbered_groups=3660,
        fresh_numbered_groups=0,all_probe_checks_passed=True,
        declared_numerical_scope='Inherited706 two-node radial/neutral-CAR finite diagnostic.',
        product_projection_commutes_with_sum_not_separate_terms=True,
        rows=rows,full_original_graph_norm_argument_in_working_derivation=True,
        complete_second_source_history_match_signed_off=False,
        spatial_refinement_or_Q_E_match_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
