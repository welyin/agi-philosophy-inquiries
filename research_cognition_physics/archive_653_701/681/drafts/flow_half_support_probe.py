"""681 executed mapping diagnostic, not a new physics round.
Abelian transverse linearized gauge flow on a finite periodic time regulator.
Tests preservation of the physical half-algebra, NOT reflection positivity itself.
"""
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'flow_half_support_probe_results.json'

def run():
    rows=[]
    for n in (8,16):
        shift=np.roll(np.eye(n),1,axis=1)
        lap=2*np.eye(n)-shift-shift.T
        vals,vec=np.linalg.eigh(lap)
        vals=np.maximum(vals,0)
        vals[vals<1e-13]=0
        reflection=np.eye(n)[(1-np.arange(n))%n]
        positive=np.arange(1,n//2+1)
        negative=np.array([i for i in range(n) if i not in positive])
        mask=np.zeros(n);mask[positive]=1;pi=np.diag(mask)
        assert np.max(abs(reflection@pi@reflection-(np.eye(n)-pi)))<1e-14
        for s in (.2,1.,2.):
            for kind,frequencies in [('gradient',vals),('half_space_EOM',np.sqrt(vals))]:
                flow=(vec*np.exp(-s*frequencies))@vec.T
                comm=float(np.max(abs(flow@reflection-reflection@flow)))
                cross=flow[np.ix_(positive,negative)]
                support=float(np.linalg.norm(cross,2))
                # A boundary impulse on the negative half changes a positive-half flowed field.
                i,j=np.unravel_index(np.argmax(abs(cross)),cross.shape)
                response=float(cross[i,j])
                assert comm<3e-14 and support>1e-3 and response>1e-3
                rows.append(dict(time_points=n,flow_distance=s,kind=kind,
                    reflection_commutator_error=comm,
                    opposite_half_operator_norm=support,
                    impulse_source_time=int(negative[j]),
                    positive_read_time=int(positive[i]),response=response))
    return dict(entry_round=681,not_formal_round=True,
        regulator='declared periodic one-dimensional Euclidean time Laplacian',
        rows=rows,original_full_model_not_replaced=True,
        gradient_linear_Abelian_transverse_test=True,
        EOM_only_decaying_half_space_branch_not_entire_paper_slab=True,
        reflection_covariance_does_not_imply_half_support=True,
        no_RP_failure_or_success_conclusion=True,
        no_new_spatial_dimension_or_physical_time_claim=True)

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=681,rows=len(result['rows']),all_checks_passed=True)))
