"""916 independent equations, original canonical stress, off-shell covariance."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
import covariant_joint_residual as c
aj=c.aj;ex=c.ex;HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'covariant_residual_validation_results.json'

def flatjets():
    return dict(g=np.diag([-1.,1,1,1])[None],dg=np.zeros((1,4,4,4)),ddg=np.zeros((1,4,4,4,4)),
        phi=np.zeros((1,6)),dphi=np.zeros((1,4,6)),ddphi=np.zeros((1,4,4,6)),
        A=np.zeros((1,4,12)),dA=np.zeros((1,4,4,12)),ddA=np.zeros((1,4,4,4,12)))

def manufactured(par):
    p=flatjets();H=.37;a=1.2
    p['g'][0,1:,1:]=a*a*np.eye(3);p['dg'][0,0,1:,1:]=2*H*a*a*np.eye(3);p['ddg'][0,0,0,1:,1:]=4*H*H*a*a*np.eye(3)
    z=c.geometry(p);expected=3*H*H*p['g']
    ric=c.maximum(z['Ric']-expected);scalar=c.maximum(z['R']-12*H*H)
    p=flatjets();value=.23;omega=.7;k=1.1
    p['phi'][0,5]=value;p['ddphi'][0,0,0,5]=-omega**2*value;p['ddphi'][0,1,1,5]=-k*k*value
    kg=c.maximum(c.euler(p,par)['scalar'][:,5]-(omega*omega-k*k-1)*value)
    p=flatjets();p['A'][0,2,11]=value;p['ddA'][0,0,0,2,11]=-omega**2*value;p['ddA'][0,1,1,2,11]=-k*k*value
    maxwell=c.maximum(c.euler(p,par)['YM'][:,2,11]-par.K[11]*(omega*omega-k*k)*value)
    result=dict(deSitter_Ricci_error=ric,deSitter_scalar_curvature_error=scalar,off_shell_Klein_Gordon_error=kg,off_shell_Maxwell_error=maxwell)
    assert max(result.values())<1e-12,result
    return result

def coord_derivatives(p,dp,par):
    r=c.euler(p,par);h=1e-24;d={k:[] for k in ('gravity','scalar','YM')}
    for mu in range(4):
        rr=c.euler({k:p[k].astype(complex)+1j*h*dp[mu][k] for k in p},par)
        for k in d:d[k].append(rr[k].imag/h)
    return r,{k:np.stack(v,axis=1) for k,v in d.items()}

def affine_field_lie(p,dp,xi,dx):
    out={};coordcount=dict(g=2,dg=3,ddg=4,phi=0,dphi=1,ddphi=2,A=1,dA=2,ddA=3)
    for key,a in p.items():
        v=np.einsum('br,br...->b...',xi,np.stack([q[key] for q in dp],axis=1))
        for axis in range(1,coordcount[key]+1):
            v+=np.moveaxis(np.einsum('bir,br...->bi...',dx,np.moveaxis(a,axis,1)),1,axis)
        out[key]=v
    return out

def euler_lie(r,d,xi,dx,alpha=None):
    out={k:np.einsum('br,br...->b...',xi,d[k]) for k in d}
    out['gravity']+=np.einsum('bmr,brn->bmn',dx,r['gravity'])+np.einsum('bnr,bmr->bmn',dx,r['gravity'])
    out['YM']-=np.einsum('brm,bra->bma',dx,r['YM'])
    if alpha is not None:
        out['scalar']+=np.einsum('ba,aAB,bB->bA',alpha,c.REP,r['scalar'])
        out['YM']+=c.bracket(alpha[:,None,:],r['YM'])
    return out

def canonical_stress(p,par,r):
    z=r['z'];G,_,_,_=c.target(p['phi'],par);D=r['D'];F=r['F']
    vel=(D[:,0]-np.einsum('bi,biA->bA',z['beta'],D[:,1:]))/z['alpha'][:,None]
    electric=(F[:,0,1:]-np.einsum('bj,bjia->bia',z['beta'],F[:,1:,1:]))/z['alpha'][:,None,None]
    pi=z['vol'][:,None]*np.einsum('bAB,bB->bA',G,vel)
    E=z['vol'][:,None,None]*par.K*np.einsum('bij,bja->bia',z['invgamma'],electric)
    y=dict(g=p['g'],phi=p['phi'],pi=pi,A=p['A'][:,1:],E=E)
    class LocalModel(c.ev.Model):
        def __init__(self):self.K=par.K;self.L=par.L;self.u=par.u
        def grad(self,a):
            if a is y['phi']:return p['dphi'][:,1:]
            if a is y['A']:return p['dA'][:,1:,1:]
            raise AssertionError('Unexpected canonical gradient')
    old=LocalModel().matter(y,z)
    return dict(full_canonical_stress_error=c.maximum(old['T']-r['stress']),canonical_scalar_velocity_error=c.maximum(old['v']-vel))

def run():
    par=c.parameters();known=manufactured(par)
    bg,field,_=ex.previous.rebuilt.build_pair();src=ex.loop.LoopSource(bg,2,8);x=src.points[::64]
    cache=aj.Cache(bg,x);p=cache.jets();dp=[cache.jets((mu,)) for mu in range(4)]
    r,dr=coord_derivatives(p,dp,par);canonical=canonical_stress(p,par,r)
    assert max(canonical.values())<1e-12,canonical
    rng=np.random.default_rng(916);B=len(x)
    xi=rng.normal(size=(B,4))*.03;dx=rng.normal(size=(B,4,4))*.02
    v=affine_field_lie(p,dp,xi,dx);h=1e-24
    shifted=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    expected=euler_lie(r,dr,xi,dx)
    err={k:c.maximum(shifted[k].imag/h-expected[k]) for k in expected}
    assert max(err.values())<1e-9,err
    # Original LA, not an invented gauge direction, applied to base residuals.
    weighted=aj.prior.weighted
    actual_xi,actual_dx,_,_,alpha,*_=weighted.extract_values(bg,field,x)
    correction=euler_lie(r,dr,actual_xi,actual_dx,alpha)
    # The residual derivative includes coordinate and internal covariance terms.
    # This is exact at the analytic level, but these arrays are diagnostic only.
    return dict(round=916,date='2026-10-06',samples=len(x),all_checks_passed=True,
        manufactured_equation_checks=known,canonical_dictionary_checks=canonical,
        independent_affine_diffeomorphism_covariance_errors=err,
        actual_extractor_background_residual_correction_max={k:c.maximum(v) for k,v in correction.items()},
        maximum_actual_xi=c.maximum(actual_xi),maximum_actual_dxi=c.maximum(actual_dx),maximum_actual_alpha=c.maximum(alpha),
        projected_response_on_approximate_background_requires_correction=True,
        physical_solution_error_certified=False,full_goal_completed=False,
        source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),HERE/'covariant_joint_residual.py',STAGE/'914/weighted_receiver_pairing.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
