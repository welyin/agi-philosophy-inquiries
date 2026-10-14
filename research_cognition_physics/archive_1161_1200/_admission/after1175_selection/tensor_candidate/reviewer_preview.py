"""Independent power-series integration of the fixed tensor candidate.

Implemented from candidate.md and the checked Clifford equations. Does not import
author preview.py or numerical_bound.py. Default recomputes and compares read-only;
--write creates the reviewer JSON once. Numerical agreement is not a residual proof.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
GRID=65; LENGTH=16.; ELL=LENGTH/GRID; TIME=2.
EPS=.15; NU=.35; OMEGA=.5; MS=2.; GT=.2
MASSES=np.array([0.,.3,.5])
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1.,-1.])
BETA=np.kron(Z,I2); ALPHA=np.array([np.kron(X,a) for a in (X,Y,Z)])
GAMMA=np.array([BETA,*ALPHA])
PAIRS=[(a,b) for a in range(4) for b in range(a+1,4)]
INTS=np.indices((33,33,33)).reshape(3,-1).T
MULT=np.prod(1+(INTS!=0),axis=1).astype(np.int64)
WEIGHTS=MULT.astype(float)/GRID**3
MOM=np.sin(2*np.pi*INTS/GRID)/ELL
PROBE_K=math.sin(2*math.pi/GRID)/ELL
PMOM=np.vstack([np.zeros(3),np.eye(3)*PROBE_K])
P_MINUS=np.zeros((3,3),complex); P_MINUS[1:,1:]=(I2-Y)/2
P_MASSIVE=np.diag([0.,1.,1.]); P_ZERO=np.diag([1.,0.,0.])
EFFECTS={'clock_minusY_joint':P_MINUS,'massive_position':P_MASSIVE,
         'massless_position':P_ZERO}


def source_force(q):
    energy=2*np.sum(WEIGHTS[:,None]*q[:,:4]*np.column_stack([
        np.full(len(MOM),MS),MOM]),axis=0)
    return energy


def rotate_coefficients(q,d):
    """Linear A(d)q, without assembling a dense 10x10 matrix per mode."""
    out=np.zeros_like(q)
    for j,(a,b) in enumerate(PAIRS):
        out[:,a]-=2*d[:,b]*q[:,4+j]
        out[:,b]+=2*d[:,a]*q[:,4+j]
        out[:,4+j]=2*(d[:,b]*q[:,a]-d[:,a]*q[:,b])
    return out


def source_energy(q,g,coupling):
    e=source_force(q); n=1+EPS*g[0,2]
    v=1+NU*g[0,2]+coupling*g[1,2]*np.array([1.,-1.,0.])
    return float(n*e[0]+np.dot(v,e[1:])-OMEGA*np.sum(g[:,0]))


def integrate(steps,order,coupling,limit=None):
    dt=TIME/steps; count=steps if limit is None else limit
    d0=np.column_stack([np.full(len(MOM),MS),MOM])
    ds=d0*np.array([EPS,NU,NU,NU])
    dtensor=d0*np.array([0.,coupling,-coupling,0.])
    q=np.zeros((len(MOM),10)); q[:,1]=1
    zero=INTS[:,0]==0; q[zero,1]=0; q[zero,0]=1
    g=np.array([[1.,0.,0.],[1.,0.,0.]])
    u=np.broadcast_to(np.eye(4,dtype=complex),(2,4,3,4,4)).copy()
    integrals=np.zeros(4)
    kin=np.einsum('pb,bij->pij',PMOM,ALPHA)[:,None]
    mass=MASSES[None,:,None,None]*BETA
    hb=kin+mass; hs=NU*kin+EPS*mass
    hta=coupling*np.einsum('pb,b,bij->pij',PMOM,np.array([1.,-1.,0.]),ALPHA)
    ht=np.array([hta,-hta])[:,:,None]
    qp=np.zeros((order+1,len(MOM),10)); gp=np.zeros((order+1,2,3))
    up=np.zeros((order+1,2,4,3,4,4),complex)
    ip=np.zeros((order+1,4)); bp=np.zeros((order+1,2))
    initial_energy=source_energy(q,g,coupling)
    max_energy=0.; start=time.perf_counter()
    for step in range(count):
        qp[0]=q; gp[0]=g; up[0]=u; ip[0]=integrals
        for n in range(order):
            e=source_force(qp[n])
            bp[n]=[EPS*e[0]+NU*np.sum(e[1:]),coupling*(e[1]-e[2])]
            # Convolution in theta=(t-t0)/dt, then divide by n+1.
            mixtures=np.tensordot(gp[:n+1,:,2].T,qp[n::-1],axes=(1,0))
            dq=rotate_coefficients(qp[n],d0)
            dq+=rotate_coefficients(mixtures[0],ds)
            dq+=rotate_coefficients(mixtures[1],dtensor)
            qp[n+1]=dt/(n+1)*dq
            xy=np.sum(bp[:n+1,:,None]*gp[n::-1,:,:2],axis=0)
            dg=np.column_stack([-2*xy[:,1],2*xy[:,0]+2*OMEGA*gp[n,:,2],
                                -2*OMEGA*gp[n,:,1]])
            gp[n+1]=dt/(n+1)*dg
            us=np.tensordot(gp[:n+1,0,2],up[n::-1],axes=(0,0))
            ut=np.tensordot(gp[:n+1,1,2],up[n::-1],axes=(0,0))
            du=hb@up[n]+hs@us+ht@ut
            up[n+1]=-1j*dt/(n+1)*du
            base=1. if n==0 else 0.
            zs,zt=gp[n,:,2]
            ip[n+1]=dt/(n+1)*np.array([base+EPS*zs,base+NU*zs+coupling*zt,
                                     base+NU*zs-coupling*zt,base+NU*zs])
        q=np.sum(qp,axis=0);g=np.sum(gp,axis=0);u=np.sum(up,axis=0)
        integrals=np.sum(ip,axis=0)
        max_energy=max(max_energy,abs(source_energy(q,g,coupling)-initial_energy))
        if limit is None and (step+1)%32==0:
            print(json.dumps({'method':'independent_Taylor','step':step+1,'steps':steps,
                              'order':order,'g':coupling}),flush=True)
    residuals={'energy_max_drift':max_energy,
       'source_quadratic_invariant':float(np.max(abs(np.sum(q*q,axis=1)-1))),
       'G_norm':float(np.max(abs(np.sum(g*g,axis=1)-1))),
       'probe_unitarity':float(np.max(abs(u.conj().swapaxes(-1,-2)@u-np.eye(4))))}
    return {'u':u,'g':g,'q':q,'integrals':integrals,'residuals':residuals,
            'elapsed':time.perf_counter()-start}


def full_u(blocks):
    result=np.zeros((48,48),complex)
    for p in range(4):
        for m in range(3):
            start=12*p+4*m; result[start:start+4,start:start+4]=blocks[p,m]
    return result


def effects(axis,clock):
    f=np.eye(4,dtype=complex)/2
    f[0,axis+1]=1j/4; f[axis+1,0]=-1j/4
    return np.kron(f,np.kron(clock,np.eye(4)))


def witness(massless=False):
    flavor=np.array([1.,0.,0.]) if massless else np.array([0.,1.,1.])/math.sqrt(2)
    chir=np.array([1.,0.,0.,1.])/math.sqrt(2)
    return np.kron(np.ones(4)/2,np.kron(flavor,chir))


def menu(blocks):
    us=[full_u(x) for x in blocks]; result={}
    for axis,label in enumerate('xyz'):
        for name,c in EFFECTS.items():
            effect=effects(axis,c); es=[u.conj().T@effect@u for u in us]
            psi=witness(name=='massless_position')
            probs=[float(np.vdot(psi,e@psi).real) for e in es]
            result[label+'_'+name]={'source_A':probs[0],'source_B':probs[1],
               'signed_gap':probs[0]-probs[1],
               'complete_48code_effect_distance':float(np.linalg.norm(es[0]-es[1],2))}
    return result


def reflection_checks(blocks):
    ud=BETA@np.kron(I2,-1j*(X-Y)/math.sqrt(2))
    perm=np.eye(4)[[0,2,1,3]]
    reflection=np.kron(perm,np.kron(np.eye(3),ud))
    ua,ub=[full_u(x) for x in blocks]
    rng=np.random.default_rng(1176)
    psi=rng.normal(size=(48,3))+1j*rng.normal(size=(48,3));psi/=np.linalg.norm(psi)
    rotated=reflection@psi
    effect=effects(1,P_MINUS)
    ea=ua@psi;eb=ub@rotated
    pa=float(np.trace(ea.conj().T@effect@ea).real)
    pb=float(np.trace(eb.conj().T@(reflection@effect@reflection.conj().T)@eb).real)
    out=ua@psi
    rho0=psi.T@psi.conj()
    rhopre=out.T@out.conj()
    # Three mass outcomes: |0>, |-Y>, |+Y>. Keep site, Dirac and R amplitudes.
    columns=np.array([[1,0,0],[0,1/math.sqrt(2),1/math.sqrt(2)],
                      [0,-1j/math.sqrt(2),1j/math.sqrt(2)]],complex)
    rho_record=np.zeros((3,3),complex); norm=0.; pmin=1.
    wave=out.reshape(4,3,4,3)
    sums=np.zeros(3)
    for first in range(0,GRID**3,4096):
        flat=np.arange(first,min(first+4096,GRID**3))
        xyz=np.column_stack([flat//GRID**2,(flat//GRID)%GRID,flat%GRID])
        phase=np.column_stack([np.ones(len(flat)),np.exp(2j*np.pi*xyz/GRID)])
        position=(phase@wave.reshape(4,-1)).reshape(-1,3,4,3)/math.sqrt(GRID**3)
        amplitudes=np.einsum('mk,amdr->akdr',columns.conj(),position)
        probabilities=np.sum(abs(amplitudes)**2,axis=(2,3))
        norm+=float(np.sum(probabilities));pmin=min(pmin,float(probabilities.min()))
        sums+=np.sum(probabilities,axis=0)
        rho_record+=np.einsum('akdr,akds->rs',amplitudes,amplitudes.conj())
    density=out@out.conj().T
    mass0=np.arange(48).reshape(4,3,4)[:,0,:].ravel()
    mass1=np.arange(48).reshape(4,3,4)[:,1,:].ravel()
    return {'source_reflection_intertwiner':float(np.linalg.norm(ub-reflection@ua@reflection.conj().T,2)),
            'complete_passive_probability_difference':abs(pa-pb),
            'arbitrary_R_before_record_error':float(np.linalg.norm(rhopre-rho0,2)),
            'complete_site_clock_record_normalization_error':abs(norm-1),
            'record_R_error':float(np.linalg.norm(rho_record-rho0,2)),
            'record_min_probability':pmin,'fine_outcomes':3*GRID**3,
            'mass0_mass03_coherence_before_record':float(np.linalg.norm(density[np.ix_(mass0,mass1)],'fro')),
            'record_outcome_weights':sums.tolist(),
            'R_dimension_sample':3,
            'all_R_claim_from_operator_norm_not_sample':True}


def covariance_sample(q):
    inds=np.unique(np.r_[np.arange(0,len(q),701),np.flatnonzero(INTS[:,0]==0)[::57]])
    p=np.broadcast_to(np.eye(4),(len(inds),4,4)).astype(complex).copy()
    p+=np.einsum('na,aij->nij',q[inds,:4],GAMMA)
    bilinear=np.array([1j*GAMMA[a]@GAMMA[b] for a,b in PAIRS])
    p+=np.einsum('na,aij->nij',q[inds,4:],bilinear);p/=2
    return {'sampled_classes':len(inds),'projector_error':float(np.max(abs(p@p-p))),
            'trace2_error':float(np.max(abs(np.trace(p,axis1=-2,axis2=-1)-2)))}


def delta_menu(a,b):
    return max(abs(a[k][field]-b[k][field]) for k in a for field in a[k])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    coarse=integrate(64,16,GT)
    fine=integrate(128,12,GT)
    scalar=integrate(64,16,0.)
    active_menu=menu(fine['u']);scalar_menu=menu(scalar['u'])
    author=json.loads((HERE/'results_preview.json').read_text(encoding='utf-8'))
    comparison={'author_fine_menu_max_difference':delta_menu(active_menu,author['diagnostic_refinement']['fine_original_menu']),
       'author_scalar_menu_max_difference':delta_menu(scalar_menu,author['scalar_original_menu']),
       'author_active_integral_max_difference':float(np.max(abs(fine['integrals']-author['background']['active']['integrals'])))}
    convergence={'whole_prepared_operator_difference':float(np.max(np.linalg.norm(coarse['u']-fine['u'],ord=2,axis=(-2,-1)))),
        'integrals_max_difference':float(np.max(abs(coarse['integrals']-fine['integrals']))),
        'all_source_coefficients_max_difference':float(np.max(abs(coarse['q']-fine['q']))),
        'G_max_difference':float(np.max(abs(coarse['g']-fine['g'])))}
    checks=reflection_checks(fine['u'])
    result={'kind':'UNNUMBERED_INDEPENDENT_TAYLOR_REVIEW',
       'implementation_provenance':{'author_preview_imported':False,'author_source_implementation_read_before_writing':False,
          'candidate_equations_and_checked_10_component_bridge_read':True,'author_reported_results_seen_before_run':True,
          'external_blind_review':False,'interval_certificate_reimplemented':False},
       'method':{'full_source_classes':len(MOM),'active_runs':[[64,16],[128,12]],'scalar_run':[64,16],
          'fourier_coefficients':'analytic 1/2 and +/-i/4','finite_K_simulated':False},
       'weight_checks':{'integer_sum':int(MULT.sum()),'target':GRID**3,
          'zero_x_plane_integer_weight':int(MULT[INTS[:,0]==0].sum()),'zero_x_plane_target':GRID**2,
          'rank2_source_particles':2*GRID**3},
       'active_integrals':fine['integrals'].tolist(),'active_G':fine['g'].tolist(),
       'active_original_menu':active_menu,'scalar_original_menu':scalar_menu,
       'convergence':convergence,'author_comparison':comparison,
       'active_invariants':fine['residuals'],'coarse_invariants':coarse['residuals'],
       'scalar_invariants':scalar['residuals'],'covariance_projection':covariance_sample(fine['q']),
       'common_transport_and_reference':checks,
       'inputs_sha256':{name:sha(HERE/name) for name in ['candidate.md','results_preview.json']}}
    assert result['weight_checks']['integer_sum']==GRID**3
    assert result['weight_checks']['zero_x_plane_integer_weight']==GRID**2
    assert max(convergence.values())<2e-9
    assert max(comparison.values())<2e-9
    assert checks['source_reflection_intertwiner']<1e-11
    assert checks['record_R_error']<1e-11
    assert checks['complete_site_clock_record_normalization_error']<1e-11
    return result


def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,float): assert math.isclose(a,b,rel_tol=2e-10,abs_tol=5e-13),(path,a,b)
    else:assert a==b,(path,a,b)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true')
    ap.add_argument('--bench',action='store_true');args=ap.parse_args()
    if args.bench:
        r=integrate(128,12,GT,limit=4)
        print(json.dumps({'bench_seconds':r['elapsed'],'projected_128_steps':32*r['elapsed']}))
    else:
        r=run();dest=HERE/'reviewer_preview.json'
        if args.write:
            with dest.open('x',encoding='utf-8',newline='\n') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
            print('WRITE_ONCE_PASS')
        else:
            compare(r,json.loads(dest.read_text(encoding='utf-8')));print('READ_ONLY_MATCH')
