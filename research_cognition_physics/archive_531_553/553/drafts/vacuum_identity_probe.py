"""Unreviewed 553 probe: running-tree vacuum and the same matter supertrace."""
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_singlet_common_mass_rg as rg


def run():
    saved=json.loads((HERE.parent/'joint_singlet_common_mass_rg_results.json').read_text('utf8'))
    rows=[]
    for index in (2,4,5):
        row=saved['examples'][index]
        s=row['state']
        state=np.array([s[k] for k in ('q','S','z','lambda_H','p','lambda_s','x','y')])
        q,S,z,lh,p,ls,x,y=state
        gy,gw,gc=rg.gauge_squared(-row['u'],rg.XSTAR)
        beta=np.array(rg.beta_numerator(state,(gy,gw,gc)))
        C=.25*np.array([x,y])
        L=np.array([[lh,p],[p,ls]])
        vev=np.linalg.solve(L,C)
        assert min(vev)>0 and np.linalg.eigvalsh(L)[0]>0
        bC=.25*beta[6:]
        bL=np.array([[beta[3],beta[4]],[beta[4],beta[5]]])
        bV0=2*C[0]**2+.5*C[1]**2
        bDelta=bV0-.5*bC@vev+.25*vev@bL@vev
        D=np.diag(np.sqrt(vev))
        # Matter effective potential on a fixed Jordan background: canonical matter masses.
        # NOT the coupled Einstein-frame scalar-gravity mass eigenvalues of 552.
        radial=np.linalg.eigvalsh(2*D@L@D)
        h2,s2=vev
        mw2=gw*h2/4
        mz2=(gy+gw)*h2/4
        mt2=q*h2/2
        mn2=z*s2+S*h2/2
        masses=np.array([*radial,mw2,mz2,mt2,mn2])
        degeneracies=np.array([1,1,6,3,-12,-4])
        supertrace=float(degeneracies@masses**2)
        residual=abs(bDelta-supertrace/2)
        assert residual<1e-12
        constants=np.array([1.5,1.5,5/6,5/6,1.5,1.5])
        def one_loop(log_mu):
            return float((degeneracies*masses**2)@(np.log(masses)-2*log_mu-constants)/(64*math.pi**2))
        step=1e-4
        explicit_derivative=(one_loop(step)-one_loop(-step))/(2*step)
        cancellation=abs(bDelta/rg.LOOP+explicit_derivative)
        assert cancellation<1e-11
        rows.append(dict(source_index=index,q0=row['q0'],beta_vacuum_numerator=float(bV0),
            beta_tree_minimum_numerator=float(bDelta),half_matter_supertrace=supertrace/2,
            supertrace_identity_residual=residual,explicit_scale_cancellation_residual=cancellation,
            matter_masses_squared=masses.tolist(),degrees_of_freedom=degeneracies.tolist()))
    return dict(status='preliminary_unreviewed_not_completed_round',examples=rows,
        MS_one_loop_fixed_background_matter_only=True,graviton_loops_not_included=True,
        curved_space_counterterms_and_nonminimal_beta_not_yet_derived=True,
        running_tree_vacuum_is_not_physical_scale_dependence=True)


if __name__=='__main__':
    result=run()
    with (HERE/'vacuum_identity_probe_results.json').open('x',encoding='utf8',newline='\n') as stream:
        stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
