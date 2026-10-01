"""659 entry: actual nonadjacent-time source in the original auxiliary weight.

Excludes only a nearest-time scalar product on raw E slices, not a quantum
process with retained fermions, additional state, or a controlled limit.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_spatial_auxiliary_geometry as old


def run():
    f=old.frame(1,3,4);w=old.weights(f)
    rp=np.sqrt(4-3/np.sqrt(2));rm=np.sqrt(4+3/np.sqrt(2));delta=(1/rp-1/rm)**2
    expected=delta*np.array([4,1,1])/9;actual=w[0,6:9]
    assert max(abs(expected-actual))<3e-14 and expected[0]>0
    rows=[]
    for step in (.004,.002):
        logs={}
        for s in (-1,1):
            for t in (-1,1):
                angle=np.zeros(12);angle[0]=s*step;angle[6]=t*step
                e,_=old.plane(f,angle);value=old.ratio(f,e)
                assert abs(value.imag)<1e-13 and value.real>0
                logs[(s,t)]=float(np.log(value.real))
        fd=(logs[(1,1)]-logs[(1,-1)]-logs[(-1,1)]+logs[(-1,-1)])/(4*step**2)
        assert abs(fd-expected[0])<2e-6
        rows.append(dict(step=step,actual_signed_Pfaffian_log_mixed_derivative=fd,
            projection_prediction=float(expected[0]),difference=float(abs(fd-expected[0]))))
    zero=float(old.weights(old.frame(0,3,4))[0,6]);assert zero<1e-25
    return dict(date='2026-10-02',status='659 entry; incomplete',nx=3,nt=4,
        nonadjacent_time_weights=actual.tolist(),analytic_radical_weights=expected.tolist(),
        finite_difference_rows=rows,old_spatially_decoupled_lag_two_weight=zero,
        dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
            ('joint_spatial_auxiliary_geometry.py','research_note_656.md','research_note_658.md')},
        scope='Original finite raw auxiliary weight has a nonzero mixed log source between time0 and time2. This excludes only a twice-differentiable scalar nearest-time factorization in a positive neighborhood. No all-process or all-route no-go.')


if __name__=='__main__':
    r=run()
    with (HERE/'temporal_memory_probe_results.json').open('x',encoding='utf8') as stream:
        json.dump(r,stream,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=r['status'],actual_nonadjacent_response=r['nonadjacent_time_weights'][0])))
