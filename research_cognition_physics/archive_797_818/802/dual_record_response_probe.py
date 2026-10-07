"""Working 802: retain both outgoing observables in a mixed response.

Finite exact diagnostic, not a bosonic CCR or original continuum simulation.
"""
from pathlib import Path
from fractions import Fraction as Q
import sys,json
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'798'))
import internal_register_probe as a
def run():
    # Same numerical Yukawa components as the old finite calibrations;
    # the state and two-dimensional first factor are diagnostic choices only.
    record=a.scale(a.I+a.Z,Q(1,2))
    current=a.scale(a.X,Q(31,100))+a.scale(a.Y,Q(9,100))
    rho_b=a.proj(2,0)
    rho_f=a.scale(a.I+a.scale(a.X,Q(1,3))+a.scale(a.Y,Q(1,5))+a.scale(a.Z,Q(1,7)),Q(1,2))
    assert Q(1,9)+Q(1,25)+Q(1,49)<1
    rho=np.kron(rho_b,rho_f)
    field=np.kron(a.X,a.I);pointer=np.kron(a.I,record)
    vertex=np.kron(a.Y,current)
    comm=lambda x,y:x@y-y@x
    delta=lambda x:a.scale(comm(x,vertex),a.C(0,1))
    mean=lambda x:np.trace(rho@x)
    partial=mean(delta(field)@pointer)-mean(delta(field))*mean(pointer)
    pointer_term=mean(field@delta(pointer))-mean(field)*mean(delta(pointer))
    joint=mean(delta(field@pointer))-mean(delta(field))*mean(pointer)-mean(field)*mean(delta(pointer))
    wb=np.trace(rho_b@a.X@a.Y);wrev=np.trace(rho_b@a.Y@a.X)
    cov_pj=np.trace(rho_f@record@current)-np.trace(rho_f@record)*np.trace(rho_f@current)
    cov_jp=np.trace(rho_f@current@record)-np.trace(rho_f@current)*np.trace(rho_f@record)
    formula=a.C(0,1)*(wb*cov_pj-wrev*cov_jp)
    assert joint==partial+pointer_term==formula==a.C(Q(13,750))
    assert partial==a.C(Q(13,750),Q(4,125))
    assert pointer_term==a.C(0,Q(-4,125))
    assert a.same(comm(field,pointer),a.scale(field,0))
    # Common evolution preserves their commutator also at first order.
    assert not any((comm(delta(field),pointer)+comm(field,delta(pointer))).flat)
    return dict(working_round=802,all_checks_passed=True,
                field_only_connected_response=str(partial),
                required_pointer_response=str(pointer_term),
                joint_connected_response=str(joint),
                direct_commutator_and_two_covariance_formula_agree=True,
                common_evolution_preserves_mutual_commutation=True,
                frozen_pointer_false_imaginary_part=True,
                original_W_or_full_continuum_vertex_used=False,
                continuous_physical_nonzero_signal_proven=False,
                formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    result=run()
    (HERE/'dual_record_response_probe_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
