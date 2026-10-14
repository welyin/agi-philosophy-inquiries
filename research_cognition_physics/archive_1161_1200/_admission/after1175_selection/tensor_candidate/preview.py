"""Unnumbered tensor feasibility preview; no finite-K parent simulation.

Default read-only comparison; --write creates results_preview.json once.
RK4 convergence diagnostics are not certified error bounds.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
GRID,LENGTH,TIME,STEPS=65,16.0,2.0,512
EPS,NU,OMEGA,COUPLING,MS=.15,.35,.5,.2,2.0
ELL=LENGTH/GRID
VOL=GRID**3
MP=np.array([0.,.3,.5])
I2=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
BETA=np.kron(SZ,I2)
ALPHA=np.array([np.kron(SX,s) for s in (SX,SY,SZ)])
PAIRS=list(itertools.combinations(range(4),2))
MODES=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
PMOM=2*np.pi*MODES/LENGTH
LMOM=np.sin(ELL*PMOM)/ELL
CHI=np.array([1,0,0,1],complex)/np.sqrt(2)
CLOCK=np.zeros((3,3),complex); CLOCK[1:,1:]=(I2-SY)/2
MASSLESS=np.diag([1.,0.,0.]).astype(complex)
KABS=np.sin(2*np.pi*np.array(list(itertools.product(range(33),repeat=3)))/GRID)/ELL
MULT=np.prod(np.where(KABS==0,1,2),axis=1)
WEIGHTS=MULT.astype(float)/VOL
assert len(KABS)==33**3 and int(sum(MULT))==VOL

def source_observables(r):
    em=2*MS*np.dot(WEIGHTS,r[:,0])
    et=2*np.sum(WEIGHTS[:,None]*KABS*r[:,1:],axis=0)
    return em,et

def coeffs(g,coupling):
    n=1+EPS*g[0,2]
    v=np.ones(3)*(1+NU*g[0,2])
    v[:2]+=coupling*g[1,2]*np.array([1.,-1.])
    return n,v

def derivative(state,coupling):
    r,s,g,acc,u=state
    n,v=coeffs(g,coupling)
    d=np.empty_like(r); d[:,0]=n*MS; d[:,1:]=KABS*v
    dr=np.zeros_like(r); ds=np.empty_like(s)
    for j,(a,b) in enumerate(PAIRS):
        dr[:,a]-=2*s[:,j]*d[:,b]
        dr[:,b]+=2*s[:,j]*d[:,a]
        ds[:,j]=2*(d[:,b]*r[:,a]-d[:,a]*r[:,b])
    em,et=source_observables(r)
    fields=np.array([EPS*em+NU*sum(et),coupling*(et[0]-et[1])])
    dg=np.empty((2,3))
    dg[:,0]=-2*fields*g[:,1]
    dg[:,1]=2*(fields*g[:,0]+OMEGA*g[:,2])
    dg[:,2]=-2*OMEGA*g[:,1]
    # B is obtained by source rotation and G_t X-conjugation, not a new ODE fit.
    velocities=np.array([v,v[[1,0,2]]])
    kin=np.einsum('ab,pb,bij->apij',velocities,LMOM,ALPHA)
    h=kin[:,:,None,:,:]+n*MP[None,None,:,None,None]*BETA
    return dr,ds,dg,np.r_[n,v],-1j*(h@u)

def plus(x,y,a): return tuple(xx+a*yy for xx,yy in zip(x,y))

def evolve(coupling,steps):
    r=np.zeros((len(KABS),4)); r[:,1]=1
    zero=KABS[:,0]==0; r[zero,1]=0; r[zero,0]=1
    s=np.zeros((len(KABS),6))
    g=np.zeros((2,3)); g[:,0]=1
    u=np.broadcast_to(np.eye(4,dtype=complex),(2,4,3,4,4)).copy()
    state=(r,s,g,np.zeros(4),u)
    em0,et0=source_observables(r)
    energy0=em0+sum(et0)-2*OMEGA
    initial_fields=np.array([EPS*em0+NU*sum(et0),coupling*(et0[0]-et0[1])])
    zsecond=-4*OMEGA*initial_fields
    asecond=np.ones(3)*(EPS-NU)*zsecond[0]
    asecond[:2]-=coupling*zsecond[1]*np.array([1.,-1.])
    energy_error=0.; history=[]; dt=TIME/steps
    for k in range(steps+1):
        r,s,g,acc,u=state
        if k%(steps//16)==0 or k==steps:
            em,et=source_observables(r); n,v=coeffs(g,coupling)
            ene=n*em+np.dot(v,et)-OMEGA*sum(g[:,0])
            energy_error=max(energy_error,abs(ene-energy0))
            history.append({'t':k*dt,'zs':float(g[0,2]),'zt':float(g[1,2]),
                'N':float(n),'v':v.tolist(),'a':(n/v).tolist(),'Em':float(em),'Et':et.tolist()})
        if k==steps: break
        a=derivative(state,coupling)
        b=derivative(plus(state,a,dt/2),coupling)
        c=derivative(plus(state,b,dt/2),coupling)
        d=derivative(plus(state,c,dt),coupling)
        state=tuple(x+dt*(aa+2*bb+2*cc+dd)/6 for x,aa,bb,cc,dd in zip(state,a,b,c,d))
    return {'u':state[4],'integrals':state[3], 'history':history,
            'initial_source':[float(em0),*et0.tolist()],
            'initial_a_second':asecond.tolist(),
            'diagnostics':{'energy_per_site_error':float(energy_error),
                'source_quadratic_Casimir_error':float(np.max(np.abs(np.sum(state[0]**2,axis=1)+np.sum(state[1]**2,axis=1)-1))),
                'G_norm_error':float(np.max(abs(np.sum(state[2]**2,axis=1)-1))),
                'prepared_48code_gram_error':float(max(np.linalg.norm(x.conj().T@x-np.eye(4),2) for x in state[4].reshape(-1,4,4)))}}

def site_effect(axis):
    e=np.eye(4,dtype=complex)/2
    e[0,axis+1]=1j/4; e[axis+1,0]=-1j/4
    return e

def effect_matrix(u,axis,c):
    f=site_effect(axis); e=np.zeros((48,48),complex)
    for p,q,m,n in itertools.product(range(4),range(4),range(3),range(3)):
        e[12*p+4*m:12*p+4*m+4,12*q+4*n:12*q+4*n+4]=f[p,q]*c[m,n]*(u[p,m].conj().T@u[q,n])
    return (e+e.conj().T)/2

def witness(massless):
    c=np.array([1,0,0]) if massless else np.array([0,1,1])/np.sqrt(2)
    return np.tile(np.kron(c,CHI),4)/2

def readouts(u):
    out={}
    for b,name in enumerate('xyz'):
        for label,c,ml in [('clock_minusY_joint',CLOCK,False),('massive_position',np.diag([0,1,1]),False),('massless_position',MASSLESS,True)]:
            e0=effect_matrix(u[0],b,c); e1=effect_matrix(u[1],b,c)
            w=witness(ml)
            p0=float(np.vdot(w,e0@w).real); p1=float(np.vdot(w,e1@w).real)
            out[name+'_'+label]={'source_A':p0,'source_B':p1,'signed_gap':p0-p1,
                'complete_48code_effect_distance':float(np.linalg.norm(e0-e1,2))}
    return out

def calculate():
    active=evolve(COUPLING,STEPS)
    fine=evolve(COUPLING,2*STEPS)
    scalar=evolve(0.,STEPS)
    maxdu=max(np.linalg.norm(a-b,2) for a,b in zip(active['u'].reshape(-1,4,4),fine['u'].reshape(-1,4,4)))
    bs=2*EPS*MS+6*NU/ELL+4*COUPLING/ELL
    bp=2*EPS*sum(MP)+18*NU/ELL+12*COUPLING/ELL
    wbound=2*VOL*bs; lbound=VOL*bp; etak=1e-5
    logk=2*(math.log10(4*lbound*TIME/etak)+6*wbound*TIME/math.log(10))
    # For these four modes, each momentum lies on one coordinate axis.
    lattice_pref=ELL**2*(2*np.pi/LENGTH)**3/6
    lattice_each=lattice_pref*float(max(active['integrals'][1:]))
    summaries={}
    for name,run in [('active',active),('scalar_rule',scalar)]:
        summaries[name]={k:v for k,v in run.items() if k!='u'}
        summaries[name]['integrals']=summaries[name]['integrals'].tolist()
    return {'kind':'UNNUMBERED_DIAGNOSTIC_PREVIEW',
        'parameters':{'grid':GRID,'ell':ELL,'T':TIME,'steps':STEPS,'epsilon':EPS,'nu':NU,'g':COUPLING,'Omega':OMEGA,'M_source':MS,'M_probe':MP.tolist(),'source_classes':len(KABS)},
        'scope':{'finite_K_parent_simulated':False,'numerical_error_certified':False,'instruments_constructed':False,'Einstein_or_healthy_spin2_selected':False,'unknown_prepared_code_dimension':48,'all_frequency_branches_kept':True,'old_scalar_kernel_only_probe_R_menu':True},
        'background':summaries,'active_original_menu':readouts(active['u']),
        'scalar_original_menu':readouts(scalar['u']),
        'diagnostic_refinement':{'steps':2*STEPS,'whole_prepared_operator_difference':float(maxdu),
            'max_integral_difference':float(max(abs(active['integrals']-fine['integrals']))),
            'fine_original_menu':readouts(fine['u']),'fine_residuals':fine['diagnostics']},
        'future_budget_NOT_A_CERTIFICATE':{'new_bS':bs,'new_bP':bp,'w_bound':wbound,'L_bound':lbound,
            'etaK_target':etak,'sufficient_log10_K':logk,
            'lattice_pair_from_uncertified_integral':2*lattice_each,
            'generic_lattice_pair_bound':2*TIME*(1+NU+COUPLING)*lattice_pref,
            'needed':'Certified ODE residual; complete finite preparation/readout contract. Step doubling is diagnostic only.'}}

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{k}]')
    elif isinstance(a,float): assert math.isclose(a,b,rel_tol=3e-8,abs_tol=3e-10),(path,a,b)
    else: assert a==b,(path,a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate(); out=HERE/'results_preview.json'
    if args.write:
        with out.open('x',encoding='utf-8') as f: json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False)
        status='WRITE_ONCE'
    else:
        compare(result,json.loads(out.read_text(encoding='utf-8'))); status='READ_ONLY_MATCH'
    print(json.dumps({'status':status,'initial_source':result['background']['active']['initial_source'],
        'initial_a_second':result['background']['active']['initial_a_second'],
        'end':result['background']['active']['history'][-1],
        'menu':result['active_original_menu'],'refinement':{k:v for k,v in result['diagnostic_refinement'].items() if k!='fine_original_menu'},
        'budget':result['future_budget_NOT_A_CERTIFICATE']},indent=2,ensure_ascii=False))
