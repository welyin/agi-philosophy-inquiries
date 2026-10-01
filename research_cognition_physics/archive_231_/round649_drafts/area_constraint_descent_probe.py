"""Unfinished649 probe: original590 geometry principal symbol vs area matching.

This is not a gravity Hamiltonian. Test a specifically added corner-matching
constraint on two copies of the original internal geometry controller.
"""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_full_spatial_metric as original

Q0=np.array([[.1,.04,0.],[.04,-.12,.03],[0.,.03,.02]])
RADIUS=.2
BASIS=[]
for i in range(3):
    a=np.zeros((3,3));a[i,i]=1;BASIS.append(a)
for i,j in ((0,1),(0,2),(1,2)):
    a=np.zeros((3,3));a[i,j]=a[j,i]=1/np.sqrt(2);BASIS.append(a)
BASIS=np.array(BASIS)
TANGENT=np.array([0,2])


def area(theta):
    g=original.shape_exp(Q0+RADIUS*np.einsum('a,aij->ij',np.sin(theta),BASIS))
    return float(np.sqrt(np.linalg.det(g[np.ix_(TANGENT,TANGENT)])))


def analytic_gradient():
    values,V=np.linalg.eigh(Q0)
    dd=(values[:,None]-values[None,:])/2
    ratio=np.ones_like(dd);np.divide(np.sinh(dd),dd,out=ratio,where=abs(dd)>1e-14)
    divided=np.exp((values[:,None]+values[None,:])/2)*ratio
    g=original.shape_exp(Q0);q=g[np.ix_(TANGENT,TANGENT)];c=np.sqrt(np.linalg.det(q))
    ans=[]
    for b in BASIS:
        dg=V@(divided*(V.T@(RADIUS*b)@V))@V.T
        ans.append(.5*c*np.trace(np.linalg.solve(q,dg[np.ix_(TANGENT,TANGENT)])))
    return np.array(ans)


def run():
    zero=np.zeros(6);c0=area(zero);grad=analytic_gradient()
    # Both replicas use the same declared positive inertia; hbar/inertia only
    # give units for this probe, not a new model parameter selection theorem.
    hbar=.7;inertia=1.7
    predicted=-hbar*hbar*2*np.sum(grad*grad)/inertia
    rows=[]
    for step in (.004,.002,.001):
        numerical_gradient=[];lap=0.
        for a in range(6):
            d=np.zeros(6);d[a]=step
            cp=area(d);cm=area(-d)
            numerical_gradient.append((cp-cm)/(2*step))
            # 12 directions: each replica contributes the same second difference
            # of (area_L-area_R)^2 at matching theta_L=theta_R=0.
            lap+=2*((cp-c0)**2+(cm-c0)**2)/step**2
        actual=-hbar*hbar*lap/(2*inertia)
        rows.append(dict(step=step,finite_difference_H_Q_squared=actual,
                         error=abs(actual-predicted),
                         gradient_error=float(np.max(abs(np.array(numerical_gradient)-grad)))))
    assert predicted<-.001 and rows[-1]['error']<rows[0]['error']/12
    return dict(candidate_round=649,complete_round=False,
        original590_metric_parameters_reused=True,area=c0,area_gradient=grad.tolist(),
        area_constraint_regular_at_test_point=bool(np.linalg.norm(grad)>0),
        zero_physical_trace_Q_squared=0.,nonzero_H_trace_analytic=float(predicted),
        hbar=hbar,replica_inertia=inertia,rows=rows,
        scope='Only the declared area-matching addition to two original590 controllers; no ADM constraint or full gravity state assumed.')


if __name__=='__main__':
    result=run()
    target=HERE/'area_constraint_descent_probe_results.json'
    with target.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
