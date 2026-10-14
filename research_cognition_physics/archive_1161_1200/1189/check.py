"""1189: driven full-Fock comparison, exact bounds and finite diagnostics.

Bounds are rational certificates from proof.md. Finite occupation matrices only
check actual noncommuting Dirac drive and source identities, never certify the
full Fock channel. No historical results are rewritten.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, itertools, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OLD=BASE/'archive_935_955/940/native_exchange_bridge.py'
spec=importlib.util.spec_from_file_location('old940_read_only',OLD)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def up(q):
 q=F(q);x=float(q)
 while F.from_float(x)<q:x=math.nextafter(x,math.inf)
 return x

def cert(q):return {'exact':str(q),'upper_display':up(q)}

def bounds():
 eps=F(16,5)*F(3,5)**33+F(3,10**12)
 v=2*eps+eps*eps;d=5*v;r=6*v
 eta=d*F(347,150)
 z0=F(17,10);zplus=z0/(1-r/2)
 source=eta*(F(6,5)+4*(z0+zplus))
 assert eps<F(1529,10**10) and eta<F(4,10**6) and source<F(6,10**5)
 # One-particle stress bounds after hermitian cosine/sine lift, n=2 CAR.
 R,S,vx,vy,vz=F(8,5),F(3,2),F(13,10),F(4,5),F(21,20)
 V0=(R*R/2+R*S+2*vx*vx+2*vy*vy+F(3,2)*vz*vz)/4
 assert V0<3
 return {'epsilon':cert(eps),'frequency_squared_relative_error':cert(v),
  'd':cert(d),'r':cert(r),'full_Fock_prepared_distance':cert(eta),
  'same_kappa_source_insertion_error':cert(source),'constraint_V0_norm_bound':cert(V0),
  'source_norm_bound_per_mode':'1','Nmax':1,'unknown_input_code_dimension':30,
  'full_Fock_output_preserved':True,'arbitrary_reference_preserved':True,
  'T_max':'1','kappa':'1/10','modes':4,'m':'2','b_upper':'2/5',
  'u_upper':'3/5','w_upper':'4','proof_required_for_quantifiers':True}

def model():
 q,j,t,e0,e1=old.source(np.array([1.,0.,0.]),np.array([1.,0.,1.]),1.)
 masks,dg,normal=old.fermion_sector()
 hm=dg(np.diag([e0,e0,e1,e1]))
 def lift(a):
  b=np.zeros((4,4),complex);b[2:,0:2]=a
  return [(b+b.conj().T)/np.sqrt(2),(b-b.conj().T)/(1j*np.sqrt(2))]
 rho=lift(t[0,0]);sigma=lift(t[1,1]+t[2,2]+t[3,3])
 mom=[lift(t[0,i]) for i in (1,2,3)]
 cg=sum((.5*normal(rho[r],rho[r])+normal(rho[r],sigma[r])
  -2*normal(mom[0][r],mom[0][r])-2*normal(mom[1][r],mom[1][r])
  -1.5*normal(mom[2][r],mom[2][r])) for r in range(2))
 V0=-cg/8
 js=[dg(a) for pol in (lift((t[1,1]-t[2,2])/np.sqrt(2)),lift(np.sqrt(2)*t[1,2])) for a in pol]
 ward=max(float(np.max(abs(np.einsum('a,abij->bij',old.ETA@q,t)))),float(np.max(abs(np.einsum('a,aij->ij',old.ETA@q,j)))))
 comm=max(old.norm(hm@a-a@hm) for a in js)
 pairs=max(old.norm(a@b-b@a) for a in js for b in js)
 assert ward<2e-14 and comm>.01 and pairs>.01
 assert max(old.norm(a) for a in js)<=1 and old.norm(V0)<3
 return masks,hm,V0,js,{'Ward_diagnostic':ward,'Hm_source_commutator_norm':comm,'source_pair_commutator_norm':pairs,'source_norms':[old.norm(a) for a in js],'V0_norm':old.norm(V0),'energy_levels':[e0,e1]}

def diagnostic(hm,V0,js,masks):
 # Projection of full Q^2 is assembled directly, not square of projected Q.
 ns=list(itertools.product(range(5),repeat=4));ns=[n for n in ns if sum(n)<=4]
 index={n:i for i,n in enumerate(ns)};nb=len(ns);nf=6;I=np.eye(nb);If=np.eye(nf)
 Q=[];Q2=[]
 for mode in range(4):
  a=np.zeros((nb,nb));a2=np.zeros((nb,nb))
  for col,n in enumerate(ns):
   a2[col,col]=n[mode]+.5
   for shift in (-1,1):
    m=list(n);m[mode]+=shift;m=tuple(m)
    if min(m)>=0 and m in index:
     a[index[m],col]=math.sqrt(n[mode] if shift<0 else n[mode]+1)/math.sqrt(2)
   for shift in (-2,2):
    m=list(n);m[mode]+=shift;m=tuple(m)
    if min(m)>=0 and m in index:
     a2[index[m],col]=math.sqrt(n[mode]*(n[mode]-1) if shift<0 else (n[mode]+1)*(n[mode]+2))/2
  Q.append(a);Q2.append(a2)
 N=np.diag([sum(n) for n in ns]);hb=np.diag([sum(n)+2 for n in ns])
 kappa=.1;coupling=sum(np.kron(a,b) for a,b in zip(Q,js))
 delta_exact=F(2,5)*sum(F(math.comb(2*j,j),4**j)*F(3,5)**j for j in range(33))**2+F(1,10**12)-1
 delta=float(delta_exact)
 D=.5*delta*np.kron(sum(Q2),If)
 def H(k):return np.kron(I,hm+k*k*V0)+np.kron(hb,If)-k*coupling/2
 H0=H(kappa);Hp=H0+D;S=np.kron(I,2*kappa*V0)-coupling/2
 v=np.zeros(nb*nf,complex);v[index[(0,0,0,0)]*nf+masks.index(3)]=1/math.sqrt(2);v[index[(0,0,0,0)]*nf+masks.index(5)]=1/math.sqrt(2)
 def evolution(h):
  ev,U=np.linalg.eigh(h);return U@(np.exp(-1j*ev)*(U.conj().T@v))
 v0=evolution(H0);vp=evolution(Hp)
 ex=lambda a,x:float(np.vdot(x,a@x).real)
 outside=np.diag([float(sum(n)>1) for n in ns])
 initial_indices=[j*nf+i for j,n in enumerate(ns) if sum(n)<=1 for i in range(nf)]
 outside_indices=[j*nf+i for j,n in enumerate(ns) if sum(n)>1 for i in range(nf)]
 generator_leak=float(np.linalg.norm(H0[np.ix_(outside_indices,initial_indices)],2))
 assert generator_leak>1e-3
 source_error=abs(ex(S,vp)-ex(S,v0))
 step=1e-5;deriv=float(np.max(abs((H(kappa+step)-H(kappa-step))/(2*step)-S)))
 leak=ex(np.kron(outside,If),v0);occ=ex(np.kron(N,If),v0)
 assert occ>1e-5 and leak>0 and deriv<1e-9
 return {'finite_total_occupation_cutoff':4,'dimension':nb*nf,
  'initial_state':'boson vacuum times (two-fermion mask 3 + mask 5)/sqrt(2)',
  'actual_delta_frequency_squared':str(delta_exact),'driven_mean_number':occ,
  'population_outside_initial_code':leak,'initial_code_generator_leak_norm':generator_leak,'prepared_one_input_vector_difference':float(np.linalg.norm(vp-v0)),
  'same_source_expectation_difference':source_error,'full_source_derivative_residual':deriv,
  'energy_conservation_residual':abs(ex(H0,v0)-ex(H0,v)),
  'truncation_used_as_full_Fock_certificate':False,'physical_instrument_claimed':False}

def calculate():
 b=bounds();m,hm,V0,js,mat=model();num=diagnostic(hm,V0,js,m)
 return {'round':1189,'date':'2026-10-14','status':'PASS_CONDITIONAL_DRIVEN_FINITE_BAND',
  'groups':{'1_complete_noncommuting_Dirac_source':mat,'2_exact_full_Fock_and_source_bounds':b,'3_driven_finite_diagnostic_only':num},
  'old_source_sha256':hashlib.sha256(OLD.read_bytes()).hexdigest(),
  'new_axioms':0,'Einstein_generated':False,'actual_internal_clock_implemented':False,'full_goal_complete':False}

def same(a,b):
 if isinstance(a,dict):
  assert a.keys()==b.keys()
  for k in a:same(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):same(x,y)
 elif isinstance(a,float):assert abs(a-b)<=2e-11,(a,b)
 else:assert a==b,(a,b)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();result=calculate();target=HERE/'results.json'
 if args.write:
  with target.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
  label='WRITE_ONCE_PASS'
 else:same(result,json.loads(target.read_text(encoding='utf-8')));label='READ_ONLY_MATCH'
 q=result['groups']['2_exact_full_Fock_and_source_bounds']
 print(label,'1189: 3 groups; full_output_upper=',q['full_Fock_prepared_distance']['upper_display'],'source_upper=',q['same_kappa_source_insertion_error']['upper_display'])
