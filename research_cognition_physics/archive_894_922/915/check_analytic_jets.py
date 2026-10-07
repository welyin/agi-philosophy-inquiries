"""915 independent checks of the new finite-interpolant coordinate jets."""
from pathlib import Path
import json, argparse, hashlib
import numpy as np
import analytic_material_jets as aj
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'analytic_jet_checks.json'

def maxabs(x):return float(np.max(np.abs(x)))
def implicit_solution(A,b):
    x=np.linalg.solve(A.v,b.v)
    d=np.stack([np.linalg.solve(A.v,b.d[i]-A.d[i]@x) for i in range(4)])
    dd=np.stack([np.stack([np.linalg.solve(A.v,b.dd[i,j]-A.dd[i,j]@x-A.d[i]@d[j]-A.d[j]@d[i]) for j in range(4)]) for i in range(4)])
    return x,d,dd

def run():
    rng=np.random.default_rng(2915);n=5
    a=rng.normal(size=(n,n));a=a@a.T+2*np.eye(n)
    da=rng.normal(size=(4,n,n));dda=rng.normal(size=(4,4,n,n));dda=(dda+dda.swapaxes(0,1))/2
    b=rng.normal(size=(n,2));db=rng.normal(size=(4,n,2));ddb=rng.normal(size=(4,4,n,2));ddb=(ddb+ddb.swapaxes(0,1))/2
    A=aj.Jet(a,da,dda);B=aj.Jet(b,db,ddb);x=np.linalg.solve(A,B);y=implicit_solution(A,B)
    matrix={k:maxabs(v-w) for k,v,w in zip(('value','first','second'),(x.v,x.d,x.dd),y)}
    assert max(matrix.values())<1e-12
    # Independently differentiated normal equations on the original reference.
    bg,field,_=aj.ex.previous.rebuilt.build_pair()
    src=aj.ex.loop.LoopSource(bg,2,8);points=src.points[::64]
    bc=aj.Cache(bg,points);vc=aj.Cache(field,points)
    bp=bc.reference_inputs();vp=vc.reference_inputs()
    J=np.stack([aj.ex.oldref.reference_delta(bp,bc.reference_inputs((i,))) for i in range(4)],axis=-1)
    dX=aj.ex.oldref.reference_delta(bp,vp)[...,None]
    xi=np.linalg.solve(J,dX);direct=implicit_solution(J,dX)
    actual={k:maxabs(v-w) for k,v,w in zip(('value','first','second'),(xi.v,xi.d,xi.dd),direct)}
    rel={k:actual[k]/max(1.,maxabs(w)) for k,w in zip(('value','first','second'),direct)}
    assert max(rel.values())<1e-10
    # Check scalar product/chain rules independently, including the mixed term.
    u=aj.Jet(2.,np.array([.2,-.3,.1,.4]),np.eye(4)*.03)
    v=aj.Jet(3.,np.array([-.1,.2,.5,.1]),np.eye(4)*.02)
    q=np.sqrt(u*v)
    val=6.;first=u.d*v.v+u.v*v.d
    second=u.dd*v.v+np.outer(u.d,v.d)+np.outer(v.d,u.d)+u.v*v.dd
    expected=.5/np.sqrt(val)*second-.25/val**1.5*np.outer(first,first)
    chain=maxabs(q.dd-expected);assert chain<1e-14
    return dict(round=915,all_checks_passed=True,manufactured_matrix_implicit_derivative_errors=matrix,
        original_reference_sample_count=len(points),original_reference_implicit_derivative_errors=actual,
        original_reference_relative_errors=rel,scalar_mixed_chain_error=chain,
        physical_A_error_certified=False,complex_step_roundoff_certified=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'analytic_material_jets.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    r=run()
    if a.write:TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
