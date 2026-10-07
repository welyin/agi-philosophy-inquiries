"""Independent internal covariance of916's Euler tuple on original source points."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import covariant_joint_residual as c
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;TARGET=HERE/'internal_covariance_results.json'
def run():
    bg,field,_=c.ex.previous.rebuilt.build_pair();src=c.ex.loop.LoopSource(bg,2,8);points=src.points[::128]
    p=c.aj.Cache(bg,points).jets();par=c.parameters();r=c.euler(p,par)
    alpha=np.random.default_rng(3916).normal(size=(len(points),12))*.03
    v={k:np.zeros_like(a) for k,a in p.items()}
    for k in ('phi','dphi','ddphi'):
        v[k]=np.einsum('ba,aij,b...j->b...i',alpha,c.REP,p[k])
    for k in ('A','dA','ddA'):
        aa=alpha.reshape((len(points),)+(1,)*(p[k].ndim-2)+(12,))
        v[k]=c.bracket(aa,p[k])
    h=1e-24;shifted=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    expected=dict(gravity=np.zeros_like(r['gravity']),scalar=np.einsum('ba,aij,bj->bi',alpha,c.REP,r['scalar']),YM=c.bracket(alpha[:,None,:],r['YM']))
    errors={k:c.maximum(shifted[k].imag/h-e) for k,e in expected.items()}
    assert max(errors.values())<1e-11,errors
    return dict(round=916,all_checks_passed=True,actual_source_points=len(points),constant_internal_gauge_Euler_covariance_errors=errors,
      full_goal_completed=False,source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),HERE/'covariant_joint_residual.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    r=run()
    if a.write:TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
