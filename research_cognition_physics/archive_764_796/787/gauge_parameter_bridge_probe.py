"""Unpublished 787 probe: physical free kernels through a gauge change.

Finite matrices only; no claim about a continuum singular endpoint or quantum
renormalization. Does not increment the signed-off research test count.
"""
from pathlib import Path
import json
import numpy as np


def run():
    n=5
    eye=np.eye(n)
    d=np.diag(np.linspace(.7,1.3,n))@(np.roll(eye,1,axis=1)-eye)
    k=np.vstack((d,eye))
    f=np.hstack((eye,-d))
    l=np.hstack((np.zeros((n,n)),eye))+.11*f
    p=f.T@(eye+d.T@d)@f
    pi=np.eye(2*n)-k@l
    gp=pi@np.linalg.inv(p-k@k.T)@pi.T
    cases=[]
    for t in (0.,.25,.5,.75,1.):
        gf=(1-t)*k.T+t*l
        a=gf@k
        ia=np.linalg.inv(a)
        tangent=np.eye(2*n)-k@ia@gf
        xi=.4
        top=tangent@gp@tangent.T-xi*k@ia@ia.T@k.T
        inv=np.block([[top,k@ia],[ia.T@k.T,np.zeros((n,n))]])
        op=np.block([[p,gf.T],[gf,xi*eye]])
        residual=max(np.max(np.abs(op@inv-np.eye(3*n))),np.max(np.abs(inv@op-np.eye(3*n))))
        physical=np.max(np.abs(pi@top@pi.T-gp))
        assert residual < 3e-12 and physical < 3e-12
        cases.append(dict(t=t,inverse_residual=float(residual),physical_projection_residual=float(physical)))
    return dict(round=787,status='working_probe_not_signed_off',cases=cases,
                physical_kernel_algebraic_cancellation=True,
                continuum_endpoint_proven=False,quantum_source_matching_proven=False)


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('gauge_parameter_bridge_probe_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
