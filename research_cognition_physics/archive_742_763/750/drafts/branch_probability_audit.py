"""750: original reference histories cannot be assembled outcome-by-outcome for free.

This is a necessary probability-normalization audit. The outcome-dependent
geometry rules below are explicit diagnostic choices, NOT Einstein solutions.
"""
import json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_dynamic_continuum_reference as native
import marked_source_entry as marked
TARGET=HERE/'branch_probability_audit_results.json'
def run():
    P,dP,B,dB=native.flow(.08,96,derivative=True)
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    p,branches,dp,_=marked.conditional(P,e,dP)
    energies=[float(np.trace(B@(Pr-P)).real/2) for Pr in branches]
    rules=[('opposite_history_deformations',(-1.,1.)),
           ('record_energy_proportional_history_deformations',tuple(energies))]
    rows=[]
    for label,alpha in rules:
        predicted=dp*(alpha[1]-alpha[0]);assert abs(predicted)>1e-9
        entries=[]
        for lam in (.04,.02,.01):
            sums=[]
            for sign in (1,-1):
                pp=[]
                for a in alpha:
                    Pr,*_=native.flow(.08,96,gamma=sign*lam*a)
                    pp.append(float(np.vdot(e,Pr@e).real))
                sums.append((1-pp[0])+pp[1])
            centered=(sums[0]-sums[1])/(2*lam)
            entries.append(dict(lambda_parameter=lam,naive_probability_sum=sums[0],
                normalization_defect=sums[0]-1,centered_derivative=centered,
                predicted_derivative=predicted,centered_derivative_error=abs(centered-predicted)))
        assert abs(entries[-1]['normalization_defect'])>1e-10
        assert entries[-1]['centered_derivative_error']<entries[0]['centered_derivative_error']/5
        rows.append(dict(rule=label,alpha=list(alpha),geometric_feedback_not_solved=True,rows=entries))
    return dict(entry_round=750,latest_completed_round=749,formal_tests_unchanged=3436,
        actual_reference_probability=p,actual_probability_derivative=dp,
        candidate_rule='Use outcome-specific deformed PREPARATION histories, then take outcome0 probability from history0 and outcome1 from history1.',
        rows=rows,
        exact_necessary_condition='For backgrounds B+lambda b_r, sum_r p_r(B+lambda b_r)=1 requires sum_r Dp_r[b_r]=0 at first order; positivity of separate branches does not imply this.',
        safe_joint_instrument_identity='A fixed common preparation followed by Kraus L_r and conditional trace-preserving future maps retains the original p_r and sum p_r=1.',
        renormalization_is_not_an_unknown_input_CPTP_completion=True,
        no_retrocausal_interpretation=True,no_full_continuum_or_Einstein_simulation=True)
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
