"""1192: conditional finite-clock/weak-Ricci certificate. No apparatus certification.
Default recomputes and compares saved results; --write-results creates the first result.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def encoded(x):
    return {'exact': str(x), 'float': float(x)}

def inverse(a):
    n = len(a)
    m = [list(row) + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if m[i][j])
        m[j], m[p] = m[p], m[j]
        z = m[j][j]
        m[j] = [v/z for v in m[j]]
        for i in range(n):
            if i != j:
                z = m[i][j]
                m[i] = [v-z*w for v, w in zip(m[i], m[j])]
    return [r[n:] for r in m]

def design():
    u = [[F(1), F(0), F(0), F(0)]]
    for i in range(1, 4):
        for sign in [1, -1]:
            v = [F(5, 4), F(0), F(0), F(0)]
            v[i] = sign*F(3, 4)
            u.append(v)
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        v = [F(3, 2), F(0), F(0), F(0)]
        v[i], v[j] = F(1), F(1, 2)
        u.append(v)
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    assert all(v[0]**2-sum(z*z for z in v[1:]) == 1 for v in u)
    a = [[v[i]*v[j]*(1 if i == j else 2) for i, j in pairs] for v in u]
    ai = inverse(a)
    assert all(sum(a[i][k]*ai[k][j] for k in range(10)) == F(i == j)
               for i in range(10) for j in range(10))
    norm = max(sum(abs(x) for x in r) for r in ai)
    assert norm == F(24, 5)
    # A^-1 e_last is invisible to the other nine exact directions.
    missing = [r[-1] for r in ai]
    reads = [sum(x*y for x, y in zip(row, missing)) for row in a]
    assert reads == [F(0)]*9+[F(1)]
    return u, norm, {'eta_unit': True, 'inverse_row_norms': [str(sum(map(abs,r))) for r in ai],
                    'inverse_infinity_norm': str(norm), 'missing_last_direction_null': list(map(str, missing))}

def finite_budget(norm):
    B0=L=F(1, 100); V=F(2); tau=F(1, 1000)
    dx=du=F(1, 10**12); hb=F(1, 10000); h=F(1, 20)
    J=16*V**3*(L+32*B0**2)
    b=2*dx/tau**2+2*du/tau+32*B0*V*du+J*tau/3
    B=min(B0,norm*(b+16*V**2*L*hb)+L*h)
    Q=F(256,315); I0=Q**3; I1=8*Q**2
    curvature=lambda z: 192*z*I1+12288*z*z*I0
    source=F(75,16); es=eE=F(1,8); kappa=F(1)
    margin=kappa*(source-es)-curvature(B)-eE
    assert margin > 0 and kappa*(source-es)-curvature(B0)-eE < 0
    # Coupling remains an input: finite tolerances do not exclude every G>0.
    kcrit=(curvature(B)+eE)/(source-es)
    assert 0 < kcrit < 1
    assert F(1,10)*(source-es) < curvature(B)+eE
    return {k:encoded(v) for k,v in {
        'B0':B0,'L':L,'V':V,'tau':tau,'delta_x':dx,'delta_u':du,
        'h_beam':hb,'h_cover':h,'J':J,'b':b,'Gamma_updated':B,
        'I0':I0,'I1':I1,'weak_Ricci_prior_bound':curvature(B0),
        'weak_Ricci_updated_bound':curvature(B),'source_reference_lower':source,
        'source_matching_error':es,'Einstein_tolerance':eE,'kappa':kappa,
        'strict_margin':margin,'kappa_exclusion_threshold':kcrit}.items()}

# Exact even polynomial integrals on [-1,1]. Coefficients use ascending order.
def add(a,b):
    return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0))
            for i in range(max(len(a),len(b)))]
def scale(a,s): return [s*x for x in a]
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def diff(a): return [i*a[i] for i in range(1,len(a))]
def integ(a): return sum((2*x/F(i+1) for i,x in enumerate(a) if i%2==0),F(0))

def curved_diagnostic(u):
    # Nonflat prescribed metric ds²=a(x)²dt²-dx²-dy²-dz², a=1+eps*x².
    # Mathematical diagnostic only; its probe backreaction is not certified here.
    eps=F(1,1000); a=[F(1),F(0),eps]; ap=diff(a); app=diff(ap)
    q=[F(1),F(0),F(-4),F(0),F(6),F(0),F(-4),F(0),F(1)]
    Q=integ(q); assert Q==F(256,315)
    direct=F(1,2)*integ(mul(mul(mul(a,a),app),q))*Q**2
    ibp=add(scale(mul(mul(mul(a,a),ap),diff(q)),F(-1,2)),
            scale(mul(mul(mul(a,ap),ap),q),F(-1)))
    by_parts=integ(ibp)*Q**2
    assert direct == by_parts and direct > 0
    B=2*eps*(1+eps); L=2*eps+6*eps**2
    bound=192*B*(8*Q**2)+12288*B**2*Q**3
    assert direct < bound
    E=float(eps)
    def rhs(z):
        v=z[4:]; x=z[1]; av=1+E*x*x; dav=2*E*x
        dv=np.zeros(4)
        dv[0]=-2*dav/av*v[0]*v[1]; dv[1]=-av*dav*v[0]**2
        return np.r_[v,dv]
    def rk(z,T,n):
        z=z.copy(); dt=T/n
        for _ in range(n):
            k1=rhs(z); k2=rhs(z+dt*k1/2); k3=rhs(z+dt*k2/2); k4=rhs(z+dt*k3)
            z+=dt*(k1+2*k2+2*k3+k4)/6
        return z
    worst_ratio=0.0; max_discretization=0.0; max_unit_error=0.0
    V=2.; J=16*V**3*(float(L)+32*float(B)**2)
    for x0 in [-.25,0.,.25]:
        av=1+E*x0*x0
        for vec in u:
            v=np.array(vec,dtype=float); v[0]/=av
            z=np.r_[np.array([0.,x0,0.,0.]),v]
            accel=rhs(z)[4:]
            for T in [.05,.1,.2]:
                out=rk(z,T,128); finer=rk(z,T,256)
                rem=np.max(np.abs(out[:4]-z[:4]-T*v-T*T*accel/2))
                worst_ratio=max(worst_ratio,rem/(J*T**3/6))
                max_discretization=max(max_discretization,float(np.max(np.abs(out-finer))))
                unit=(1+E*out[1]**2)**2*out[4]**2-np.dot(out[5:],out[5:])
                max_unit_error=max(max_unit_error,abs(unit-1))
                assert np.max(np.abs(out[4:])) < V and np.max(np.abs(out[1:4])) < 1
    assert worst_ratio < 1 and max_discretization < 1e-10 and max_unit_error < 1e-10
    return {'metric':'diag((1+x^2/1000)^2,-1,-1,-1)',
            'exact_Ricci_pairing':encoded(direct),'exact_integrated_connection_pairing':encoded(by_parts),
            'analytic_weak_bound':encoded(bound),'taylor_cases':90,'max_remainder_bound_ratio':worst_ratio,
            'RK_refinement_difference':max_discretization,'unit_speed_error':max_unit_error,
            'interpretation':'prescribed-metric numerical diagnostic; no same-parent matter/geometry implementation'}

def regularity_deletion():
    # Compact bump in an unsampled interval: exact flat sample germs do not fix the gap.
    eps=F(1,1000); radius=F(1,10); y=F(1,2)
    q=(1-y*y)**4; qp=-8*y*(1-y*y)**3
    a=1+eps*q; ax=eps*qp/radius
    gamma=a*ax
    assert gamma != 0
    # Center at x=.5, radius=.1: samples x=0,1 have flat open neighborhoods.
    return {'samples':[0,1],'support_interval':['2/5','3/5'],
            'Gamma_at_samples':['0','0'],'Gamma_inside_gap':encoded(gamma),
            'scope':'invalidates extension from finite samples without regularity/coverage; not an Einstein countermodel'}

def baseline_check():
    baseline=json.loads((HERE/'baseline.json').read_text(encoding='utf-8-sig'))
    p=ROOT/baseline['prior_progress_path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==baseline['prior_progress_sha256']
    snap=json.loads(p.read_text(encoding='utf-8-sig'))
    for name,sha in snap['science_files_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    for n,sha in snap['acceptance_sha256'].items():
        p=HERE.parent/n/'acceptance.json'
        assert hashlib.sha256(p.read_bytes()).hexdigest()==sha,str(p)
    return {'old_scientific_files_unchanged':len(snap['science_files_sha256']),
            'old_acceptances_unchanged':len(snap['acceptance_sha256'])}

def run():
    u,norm,d=design()
    return {'status':'PASS_WITH_SCOPE','check_groups':3,'clock_design':d,
            'finite_certificate':finite_budget(norm),
            'fluid_tracefree_pairs':{name:str(F(3,4)*(F(1)+pressure)) for name,pressure in
                 [('dust',F(0)),('vacuum',F(-1)),('negative_pressure',F(-2,3))]},
            'curved_diagnostic':curved_diagnostic(u),
            'deletion_check':regularity_deletion(),'frozen_history':baseline_check(),
            'new_adopted_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables_completed':0}

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float): assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-13),(a,b)
    else: assert a==b,(a,b)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--write-results',action='store_true'); args=p.parse_args()
    result=run(); target=HERE/'results.json'
    if args.write_results:
        if target.exists(): raise SystemExit('results exists; refusing overwrite')
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: compare(result,json.loads(target.read_text(encoding='utf-8-sig')))
    print(json.dumps({'status':result['status'],'check_groups':3,
          'strict_margin':result['finite_certificate']['strict_margin'],
          'old_scientific_files_unchanged':result['frozen_history']['old_scientific_files_unchanged']},ensure_ascii=False))
