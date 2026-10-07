"""832 working: coherent success branches and graded-tensor bookkeeping.

Finite CP instruments check algebra only. They are not constructed native
clocks or autonomous decoders in the original model.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'conditional_record_probe_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code


def choi(kraus):
    d=kraus[0].shape[1]
    vec=lambda a:a.reshape(-1,order='F')/np.sqrt(d)
    return sum(np.outer(vec(a),vec(a).conj()) for a in kraus)


def trace_distance(a,b):return float(np.sum(abs(np.linalg.eigvalsh(a-b)))/2)


def run():
    gamma,_,_,_,comp,_,_=code.code_data()
    parity=(0,(1<<10)-1,0)
    assert comp(parity)==('scalar',2)
    internal=[code.I]+[code.mul((0,0,1),code.mul(g,parity)) for g in gamma]
    internal += [code.mul((0,0,1),code.mul(gamma[i],gamma[j]))
                 for i in range(20) for j in range(i+1,20)]
    assert len(internal)==211 and all(code.adj(x)==x for x in internal)
    counts={'zero':0,'scalar':0,'logical':0}
    for x in internal:
        for y in internal:counts[comp(code.mul(x,y))[0]]+=1
    assert counts==dict(zero=44280,scalar=241,logical=0)
    v=np.array([[1.,0.],[0.,0.],[0.,0.],[0.,1.]])
    x=np.array([[0.,1.],[1.,0.]])
    z=np.diag([1.,-1.])
    ideal=choi([v]);rows=[]
    for p in (.2,.5,.8):
        for epsilon in (0.,.03):
            successful=[np.sqrt(p*(1-epsilon))*v,np.sqrt(p*epsilon)*v@z]
            failed=[np.sqrt(1-p)*v@x]
            kraus=successful+failed
            assert np.max(abs(sum(k.conj().T@k for k in kraus)-np.eye(2)))<1e-14
            success=choi(successful);unconditional=choi(kraus)
            assert abs(np.trace(success)-p)<1e-14
            err=trace_distance(unconditional,ideal)
            assert abs(err-(1-p+p*epsilon))<1e-14
            rows.append(dict(success_probability=p,conditional_error=epsilon,
                             unconditional_Choi_trace_error=err,
                             upper_bound_saturated=True))
    # Input-dependent heralding preserves the two basis vectors but not an
    # unknown coherent superposition. Normalizing each basis test hides it.
    filt=np.diag(np.sqrt([.2,.8]));plus=np.ones(2)/np.sqrt(2)
    output=filt@plus;prob=float(output@output);output/=np.sqrt(prob)
    fidelity=float(abs(plus@output)**2)
    assert abs(fidelity-.9)<1e-14
    return dict(round=832,status='working_not_formal',formal_test_groups_added=0,
                all_working_checks_passed=True,
                graded_internal_factors=211,ordered_products=counts,
                CP_success_failure_examples=rows,
                input_dependent_filter_basis_success_probabilities=[.2,.8],
                input_dependent_filter_plus_fidelity=fidelity,
                native_clock_instrument_or_success_probability_constructed=False,
                original_positive_gap_numerically_computed=False,
                scope='Exact coherent successful encoding has an input-independent success weight. The full unconditional channel retains failed branches. Finite CP examples and the explicit carrier-parity factors are checked; no original clock or decoder has been implemented.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
