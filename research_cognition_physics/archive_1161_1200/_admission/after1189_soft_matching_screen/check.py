"""Read-only diagnostic on original649 geometry; no new research round or GR model."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[2]
Q0=np.array([[.1,.04,0.],[.04,-.12,.03],[0.,.03,.02]])
BASIS=[]
for i in range(3):
    b=np.zeros((3,3));b[i,i]=1.;BASIS.append(b)
for i,j in ((0,1),(0,2),(1,2)):
    b=np.zeros((3,3));b[i,j]=b[j,i]=1/np.sqrt(2);BASIS.append(b)
BASIS=np.array(BASIS)

def metric_area_gradient(theta):
    q=Q0+.2*np.einsum('a,aij->ij',np.sin(theta),BASIS)
    w,v=np.linalg.eigh(q)
    metric=(v*np.exp(w))@v.T
    ix=np.ix_([0,2],[0,2]);cut=metric[ix]
    area=float(np.sqrt(np.linalg.det(cut)))
    d=(w[:,None]-w[None,:])/2
    ratio=np.ones_like(d);np.divide(np.sinh(d),d,out=ratio,where=abs(d)>1e-13)
    divided=np.exp((w[:,None]+w[None,:])/2)*ratio
    grad=[]
    for i,b in enumerate(BASIS):
        derivative=v@(divided*(v.T@(.2*np.cos(theta[i])*b)@v))@v.T
        grad.append(.5*area*np.trace(np.linalg.solve(cut,derivative[ix])))
    return area,np.array(grad)

def run():
    saved=json.loads((RESEARCH/'archive_629_652/649/drafts/area_constraint_descent_probe_results.json').read_text('utf8'))
    rows=[]
    hbar=.7;inertia=1.7;lam=3.
    for t in (-.3,0.,.3):
        theta=np.zeros(6);theta[0]=t
        c,g=metric_area_gradient(theta)
        n=np.r_[g,-g]
        penalty_hessian=2*lam*np.outer(n,n)
        common=np.vstack((np.eye(6),np.eye(6)))/np.sqrt(2)
        null_residual=float(np.linalg.norm(penalty_hessian@common,2))
        x=np.zeros(6);x[0]=g[2];x[2]=-g[0]
        df=np.zeros(6);df[0]=.2*np.cos(t);df[2]=-.2
        xf=float(x@df)
        a=float(n@n/inertia)
        fd=[];step=1e-5
        for j in range(6):
            dt=np.zeros(6);dt[j]=step
            fd.append((metric_area_gradient(theta+dt)[0]-metric_area_gradient(theta-dt)[0])/(2*step))
        error=float(np.max(abs(g-np.array(fd))))
        assert error<1e-8 and null_residual<1e-14
        assert abs(x@g)<1e-16 and abs(xf)>.03
        rows.append(dict(theta0=t,area=c,gradient=g.tolist(),penalty_hessian_rank=int(np.linalg.matrix_rank(penalty_hessian,tol=1e-12)),common_direction_residual=null_residual,XQ=float(x@g),Xf=xf,added_f_double_commutator_per_mu=2*hbar*hbar*xf*xf,normal_coefficient_a=a,local_frozen_normal_zero_point_per_sqrt_lambda=hbar*np.sqrt(a/2),finite_difference_gradient_residual=error))
    assert abs(rows[1]['area']-saved['area'])<1e-13
    assert np.max(abs(np.array(rows[1]['gradient'])-saved['area_gradient']))<1e-13
    assert all(r['penalty_hessian_rank']==1 for r in rows)
    return dict(status='PASS_SCOPE_ONLY',original649_parameters_reproduced=True,rows=rows,numbered_round=False,actual_finite_record_prediction=False,finite_time_equal_matching_claim=False,normal_zero_point_is_local_frozen_diagnostic_only=True,hard_limit_statement_is_analytic_not_from_truncation=True,scope='Finite penalties remain allowed; this screen does not derive or refute GR and does not claim a complete soft-constraint evolution simulation.')

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-9,(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        with (HERE/'results.json').open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:compare(result,json.loads((HERE/'results.json').read_text('utf8')))
    print(json.dumps({'status':result['status'],'diagnostic_points':3,'new_numbered_round':False}))
