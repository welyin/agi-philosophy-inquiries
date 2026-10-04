"""Preliminary 550 check: independent Einstein jets versus Jordan trace identity."""
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_reference_gravity_constraints as model


def run():
    pars=model.bare_parameters(); rows=[]
    for state in ([1.,np.sqrt(1/3)],[.7,.4],[1.2,-.2]):
        for eps in (.5,.125,.03125):
            data,grad,hess,_,_,_=model.origin_jets(state,eps,pars)
            f,df=data['F'],data['dF']
            contraction=grad@model.ETA@grad.T
            boxphi=np.einsum('mn,imn->i',model.ETA,hess)
            dlog=df/f
            ddlog=-np.eye(2)/(3*f)-np.outer(df,df)/f**2
            boxlog=np.einsum('ij,ij',ddlog,contraction)+dlog@boxphi
            dlogsq=dlog@contraction@dlog
            re=np.einsum('ij,ij',data['G'],contraction)+4*data['U']
            rj=f*(re+3*boxlog-1.5*dlogsq)
            r_trace=(4*data['V']-np.array(state)@data['dV'])/float(pars['M02'])
            mu2=float(pars['C0'])/4
            chi=2.  # Gaussian cutoff; not asserted for general effective matching
            ratio=24*((chi-1)+f/float(pars['M02']))
            assert abs(rj-r_trace)<1e-12
            assert abs(rj/mu2-ratio)<1e-11
            rows.append(dict(fields=list(state),epsilon=eps,R_J_from_Einstein=rj,
                             R_J_from_trace=r_trace,R_J_over_moment_scale_squared=rj/mu2,
                             R_J_over_cutoff_squared=rj/float(pars['cutoff_squared'])))
    return dict(status='preliminary_not_completed_round',samples=rows,
                determinant_and_constraint_existence_inherited_from_549=True,
                full_spectral_control_proved=False,completed_unification=False)


if __name__=='__main__':
    result=run()
    with (HERE/'trace_scale_probe_results.json').open('x',encoding='utf8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['status'],samples=len(result['samples']),
                         first=result['samples'][0]),ensure_ascii=False))
