"""1193 finite stress-centroid bridge; conditional mathematics, not a GR model.
Default checks saved results, --write-results creates them once.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def enc(v): return {'exact':str(v),'float':float(v)}
def matvec(A,v): return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]

def positive_integral(A,g,T,N=22):
    # Sum n=0..N A^n g T^(n+1)/(n+1)!, with componentwise positive tail.
    term=[T*x for x in g]; partial=term[:]
    for n in range(1,N+1):
        term=[x*T/F(n+1) for x in matvec(A,term)]
        partial=[x+y for x,y in zip(partial,term)]
    mu=max(sum(row) for row in A); G=max(g)
    assert all(x>=0 for row in A for x in row) and T*mu < N+3
    tail=G*mu**(N+1)*T**(N+2)/math.factorial(N+2)/(1-mu*T/F(N+3))
    return partial,[x+tail for x in partial],tail

def budget():
    B=F(1,100); L=F(1,50); v=F(1,2); r=F(1,10000)
    rf=F(1,10**7); sigma=F(1,10000); eta=F(1,100000); T=F(1)
    S=1+3*v; Ax=(1+v)*L*S*S; Av=B*(S*S+6*(1+v)*S)
    D=B*r+rf; W=(1+v)*(B*sigma+L*r+eta)
    assert (Ax,Av,D,W)==(F(3,16),F(23,80),F(11,10**7),F(39,2*10**6))
    low,high,tail=positive_integral([[F(0),F(1)],[Ax,Av]],[D,W],T)
    assert high[0]<F(2,100000) and high[1]<F(3,100000)
    return {k:enc(z) for k,z in {'B':B,'L':L,'v':v,'r_tau':r,'r_f':rf,'sigma':sigma,
        'eta':eta,'T':T,'A_x':Ax,'A_v':Av,'D':D,'W':W,
        'position_lower':low[0],'position_upper':high[0],'velocity_lower':low[1],
        'velocity_upper':high[1],'Taylor_tail':tail}.items()}

def dust_diagnostic():
    # Prescribed genuine-curvature metric N=1+x²/100. Point dust is diagnostic
    # for moment identities only, not a C2 self-gravitating solution or a quantum model.
    ep=.01; B=.0125; L=.021; vbound=.4
    def lapse(x): return 1+ep*x*x
    def gam(x):
        N=lapse(x); Np=2*ep*x
        return Np/N,N*Np  # Gamma0_01 and Gamma1_00
    def rhs(y):
        x,z,E=y; a,b=gam(x)
        return np.array([z,-b+2*a*z*z,-2*a*z*E])
    def rk(y,dt):
        a=rhs(y); b=rhs(y+dt*a/2); c=rhs(y+dt*b/2); d=rhs(y+dt*c)
        return y+dt*(a+2*b+2*c+d)/6
    x=.2+np.linspace(-.005,.005,9)
    z=.2+np.linspace(-.002,.002,9)[::-1]
    masses=1+np.arange(9)/20
    E=masses/np.sqrt(lapse(x)**2-z*z)
    y=np.array([x,z,E]); worst_d=worst_e=identity_error=0.; sample_count=0
    E0=E.sum(); X0=np.dot(E,x)/E0; w0=np.dot(E,z)/E0
    ref=np.array([[X0],[w0],[1.]])
    max_position_difference=0.
    for step in range(201):
        x,z,E=y; total=E.sum(); X=np.dot(E,x)/total; w=np.dot(E,z)/total
        dy=rhs(y); dE=dy[2]; dP=np.dot(dE,z)+np.dot(E,dy[1])
        Xprime=(np.dot(dE,x)+np.dot(E,z)-X*dE.sum())/total
        wprime=(dP-w*dE.sum())/total
        r=np.dot(E,np.abs(x-X)*(1+np.abs(z))**2)/total
        sigma=np.dot(E,(z-w)**2)/total
        a,b=gam(x); ac,bc=gam(X)
        d= -np.dot((x-X)*E,2*a*z)/total
        # Integrated conservation P'= -int Gamma^1_ab tau^ab.
        identity_error=max(identity_error,abs(dP+np.dot(E,b)),abs(Xprime-w-d))
        e=wprime-(-bc+2*ac*w*w)
        Db=B*r; Wb=(1+vbound)*(B*sigma+L*r)
        worst_d=max(worst_d,abs(d)/Db); worst_e=max(worst_e,abs(e)/Wb)
        assert abs(d)<=Db+1e-13 and abs(e)<=Wb+1e-13
        assert max(np.max(np.abs(x)),abs(X),abs(ref[0,0]))<.6
        assert max(np.max(np.abs(z)),abs(w),abs(ref[1,0]))<vbound
        max_position_difference=max(max_position_difference,abs(X-ref[0,0]))
        sample_count+=1
        if step<200: y=rk(y,.005); ref=rk(ref,.005)
    assert identity_error<1e-13
    assert max_position_difference>1e-9
    return {'samples':sample_count,'moment_identity_error':identity_error,
        'largest_centroid_slip_ratio':worst_d,'largest_momentum_defect_ratio':worst_e,
        'centroid_geodesic_position_difference':max_position_difference,
        'scope':'prescribed-background dust-flow diagnostic; sampled ratios are not certified global moment budgets'}

def deletion():
    # Flat Minkowski in x-coordinate Y=x+alpha*x²/2, x=0 event.
    # Two streams have equal energy and coordinate velocities +/-v.
    alpha=F(1,5); v=F(1,2); E=[F(1),F(1)]; speed=[v,-v]
    w=sum(e*z for e,z in zip(E,speed))/sum(E)
    sigma=sum(e*z*z for e,z in zip(E,speed))/sum(E)-w*w
    wprime=-alpha*sum(e*z*z for e,z in zip(E,speed))/sum(E)
    assert w==0 and sigma==F(1,4) and wprime==F(-1,20)
    # A translated smooth flat packet accelerated by explicit exchange:
    # tau00=A(t)r(x-X), tau01=A v r, tau11=A v² r,
    # with constant A and v'=a: f0=0, f1=A a r, r_f=0.
    acceleration=F(1,100); eta=abs(acceleration)
    assert acceleration!=0 and eta>0
    return {'cold_closure':{'geometry':'flat, nonlinear spatial coordinate',
        'instant_mean_speed':str(w),'Sigma_11':str(sigma),'centroid_acceleration':str(wprime),
        'rest_geodesic_acceleration':'0','scope':'single-event limit; not a finite smooth self-gravitating body or an equivalence violation'},
        'exchange':{'flat_packet_acceleration':str(acceleration),'required_eta':str(eta),
                    'scope':'specified exchanged stress; the environment must carry the opposite exchange; no complete quantum model claimed'}}

def history():
    b=json.loads((HERE/'baseline.json').read_text(encoding='utf-8-sig'))
    p=ROOT/b['prior_progress_path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==b['prior_progress_sha256']
    d=json.loads(p.read_text(encoding='utf-8-sig'))
    for name,h in d['science_files_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    for n,h in d['acceptance_sha256'].items():
        assert hashlib.sha256((HERE.parent/n/'acceptance.json').read_bytes()).hexdigest()==h,n
    return {'science_files_unchanged':len(d['science_files_sha256']),'acceptances_unchanged':len(d['acceptance_sha256'])}

def run(): return {'status':'PASS_WITH_SCOPE','check_groups':3,'finite_budget':budget(),
    'dust_diagnostic':dust_diagnostic(),'deletion':deletion(),'frozen_history':history(),
    'new_adopted_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables_completed':0}

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-13),(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    d=run();target=HERE/'results.json'
    if args.write_results:
        with target.open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(d,json.loads(target.read_text(encoding='utf-8-sig')))
    print(json.dumps({'status':d['status'],'groups':3,'position_upper':d['finite_budget']['position_upper'],
                     'velocity_upper':d['finite_budget']['velocity_upper'],'history':d['frozen_history']},ensure_ascii=False))
