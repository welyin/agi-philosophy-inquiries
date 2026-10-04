"""Candidate 558: fixed matter background only, not full slow dynamics."""
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'fast_link_matching_probe_results.json'


def solve(alpha,b,cutoff):
    n=np.arange(cutoff+1,dtype=float)
    casimir=n*(n+2)/4
    cosine=np.diag(np.ones(cutoff)/2,1)+np.diag(np.ones(cutoff)/2,-1)
    h=b*np.diag(casimir)-alpha*cosine
    eigenvalues,vectors=np.linalg.eigh(h)
    psi=vectors[:,0]
    lam=float(eigenvalues[0])
    dc=float(psi@cosine@psi)
    electric=float(b*np.sum(casimir*psi*psi))
    return dict(lambda_fast=lam,cosine_mean=dc,
                edge_shift=-alpha*dc,electric_energy=electric,
                ground_first_omitted_amplitude=float(abs(psi[-1])),
                gap_truncated=float(eigenvalues[1]-eigenvalues[0]))


def run():
    rows=[]
    b=1.
    for alpha in (.025,.05,.1,.2,.4,1.,4.,8.):
        row=solve(alpha,b,40)
        other=solve(alpha,b,64)
        assert abs(row['lambda_fast']-other['lambda_fast'])<1e-11
        assert abs(row['lambda_fast']-row['edge_shift']-row['electric_energy'])<1e-11
        # r1=r2 gives F0=alpha: positivity belongs to the full fast link.
        assert alpha+row['lambda_fast']>=0
        approximation=-alpha*alpha/(3*b)
        row.update(alpha=alpha,b=b,quadratic_lambda=approximation,
                   matched_edge_shift_at_second_order=2*approximation,
                   equal_radius_exact_total=alpha+row['lambda_fast'],
                   equal_radius_quadratic_total=alpha+approximation,
                   independent_cutoff_difference=abs(row['lambda_fast']-other['lambda_fast']))
        delta=3*b/4
        if alpha<delta:
            # Schur complement bound, to be independently reviewed.
            bound=alpha**4/(16*delta*(delta-alpha))*(1/(delta-alpha)+1/(2*b-alpha))
            error=abs(row['lambda_fast']-approximation)
            assert error<=bound
            row['candidate_fourth_order_error_bound']=bound
            row['actual_quadratic_error']=error
        rows.append(row)
    assert rows[-1]['equal_radius_quadratic_total']<0
    return dict(status='preliminary_558_not_a_completed_round',rows=rows,
                exact_compression_identity='P0 F^2 P0 - (P0 F P0)^2 = alpha^2/4',
                inherited_identity_source='research_note_363.md Eq.(5), not a new theorem',
                scope=dict(frozen_matter_configuration=True,
                    graph_group_and_link_inherited_inputs=True,
                    finite_cutoff_convergence_not_a_full_spectral_proof=True,
                    no_uniform_all_amplitude_perturbation_claim=True,
                    no_Born_Oppenheimer_dynamics_or_autonomous_pointer_proved=True,
                    candidate_error_bound_pending_independent_review=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():
        assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],cases=len(result['rows']),
         example=result['rows'][2],
         outside_window=result['rows'][-1]),ensure_ascii=False))
