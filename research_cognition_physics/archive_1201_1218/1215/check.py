"""1215 author checks. No earlier program imports. --write creates results once.
Finite matrices diagnose the analytic finite-window theorem; not a local device.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
SOURCES=[
 'archive_370_428/research_note_421.md',
 'archive_429_466/research_note_457.md',
 'archive_1009_1043/research_note_1037.md',
 'archive_1201_/_shared/four_laws_dependency_after1211.md',
 'archive_1201_/1209/proof.md','archive_1201_/1209/results.json',
 'archive_1201_/1213/proof.md','archive_1201_/1214/proof.md',
 'archive_1201_/1214/results.json',
]
TOL=8e-11

def req(c,m):
 if not bool(c): raise AssertionError(m)

def near(a,b,m,tol=TOL):
 e=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
 req(e<tol,m+': '+str(e))
 return e

def num(x): return round(float(x),13)
def maxerr(*xs): return num(max(xs))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def trace_distance(a,b):
 z=(a-b+(a-b).conj().T)/2
 return float(np.abs(np.linalg.eigvalsh(z)).sum()/2)
def ent(x):
 p=np.linalg.eigvalsh((x+x.conj().T)/2)
 req(p.min()>-TOL,'entropy positive input')
 p=p[p>1e-14]
 return float(-np.dot(p,np.log(p)))
def pures(v): return np.outer(v,v.conj())
def marginal(x,ds,da,keep=0):
 t=x.reshape(ds,da,ds,da)
 return np.trace(t,axis1=1,axis2=3) if keep==0 else np.trace(t,axis1=0,axis2=2)
def channel_iso(v,z,ds,da): return marginal(v@z@v.conj().T,ds,da)
def prepared(w,ds,b): return w@np.kron(np.eye(ds),b.reshape(-1,1))
def kraus(v,ds,da):
 a=v.reshape(ds,da,ds)
 return [a[:,k,:] for k in range(da)]
def apply(ks,z): return sum((k@z@k.conj().T for k in ks),np.zeros_like(z,dtype=complex))
def tp_unital(v,ds,da):
 ks=kraus(v,ds,da)
 return maxerr(near(sum(k.conj().T@k for k in ks),np.eye(ds),'TP'),near(sum(k@k.conj().T for k in ks),np.eye(ds),'unital'))
def random_density(rng,d):
 z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
 z=z@z.conj().T
 return z/np.trace(z)
def f_value(b,balances,delta):
 return sum(b[a,c] for a,ea in enumerate(balances) for c,ec in enumerate(balances) if abs(ea-ec-delta)<1e-12)
def formula(u,energies,b,balances,rho):
 d=len(energies);out=np.zeros((d,d),complex)
 for i,l,j,k in itertools.product(range(d),repeat=4):
  out[i,l]+=u[i,j]*rho[j,k]*u[l,k].conjugate()*f_value(b,balances,energies[i]-energies[l]-energies[j]+energies[k])
 return out

def qblocks(n,rule):
 # C levels0..n, complete q=1..n, both extreme singleton sectors retained.
 dc=n+1;w=np.eye(2*dc,dtype=complex)
 for q in range(1,n+1):
  inds=[q,dc+q-1];w[np.ix_(inds,inds)]=rule(q)
 return w

def conservation(w,e,balances):
 q=np.array([x+y for x in e for y in balances],float)
 return maxerr(near(w.conj().T@w,np.eye(len(q)),'full unitarity'),near((q[:,None]-q[None,:])*w,0,'bare total energy'))

def copy_isometry(v,ds,dc,label,nk):
 # Flatten environment as C,K. Includes retained coherent record, not a measurement.
 out=np.zeros((ds*dc*nk,ds),complex)
 for s in range(ds):
  for c in range(dc):out[(s*dc+c)*nk+label[s],:]=v[s*dc+c,:]
 return out

def exp_h(h,t=1):
 e,u=np.linalg.eigh((h+h.conj().T)/2)
 return (u*np.exp(-1j*t*e))@u.conj().T

def group1(rng):
 dc=6;bal=list(range(dc));en=[0,1]
 had=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 b=np.zeros(dc,complex);b[1:5]=.5
 w=qblocks(5,lambda q:had);v=prepared(w,2,b)
 errors=[conservation(w,en,bal),tp_unital(v,2,dc)]
 # Test the matrix-unit identity and arbitrary correlated S,R input.
 for j,k in itertools.product(range(2),repeat=2):
  z=np.zeros((2,2),complex);z[j,k]=1
  errors.append(near(channel_iso(v,z,2,dc),formula(had,en,pures(b),bal,z),'all matrix units'))
 rho=random_density(rng,6)
 sr_actual=apply([np.kron(k,np.eye(3)) for k in kraus(v,2,dc)],rho)
 sr_formula=np.zeros((6,6),complex)
 rr=rho.reshape(2,3,2,3)
 for j,k in itertools.product(range(2),repeat=2):
  z=np.zeros((2,2),complex);z[j,k]=1
  sr_formula+=np.kron(formula(had,en,pures(b),bal,z),rr[j,:,k,:])
 errors.append(near(sr_actual,sr_formula,'full unknown SR formula'))
 # A genuinely mixed coherent carrier, not just a diagonal or pure input.
 b2=np.zeros(dc,complex);b2[2]=1/math.sqrt(2);b2[3]=1j/math.sqrt(2)
 bm=.4*pures(b)+.6*pures(b2)
 joint=w@np.kron(np.eye(2)/2,bm)@w.conj().T
 errors.append(near(marginal(joint,2,dc),np.eye(2)/2,'mixed carrier unital'))
 exchanges=[]
 for j in range(2):
  r=np.diag([1-j,j]);o=v@r@v.conj().T
  sx=marginal(o,2,dc);cx=marginal(o,2,dc,1)
  ds=float(sx[1,1].real-j)
  db=float(np.dot(bal,np.diag(cx).real)-np.dot(bal,abs(b)**2))
  errors.append(near(ds+db,0,'exchange ledger'))
  exchanges.append({'input_energy':j,'material_change':num(ds),'carrier_change':num(db)})
 return {'name':'complete finite amplitude calibration','max_error':maxerr(*errors),'material_dimension':2,'carrier_dimension':dc,'exchanges':exchanges,'full_SR_checked':True}

def group2(rng):
 dc=4;m=1;ident=np.eye(2);x=np.array([[0,1],[1,0]])
 w=qblocks(3,lambda q:x if q==m+1 else ident)
 b=np.eye(dc)[:,m];v=prepared(w,2,b)
 j=np.zeros((dc,2));j[m,0]=1;j[m+1,1]=1
 expected=np.zeros_like(v);expected[:dc,:]=j
 errors=[conservation(w,[0,1],range(dc)),near(v,expected,'exact reset isometry')]
 rho=random_density(rng,6);sr=apply([np.kron(k,np.eye(3)) for k in kraus(v,2,dc)],rho)
 rr=marginal(rho,2,3,1)
 errors.append(near(sr,np.kron(np.diag([1,0]),rr),'unknown R marginal reset'))
 r=marginal(rho,2,3);out=v@r@v.conj().T;cr=marginal(out,2,dc,1)
 gain=float(np.dot(range(dc),np.diag(cr).real)-m)
 errors.extend([near(gain,r[1,1].real,'address reset energy'),near(ent(cr),ent(r),'stored entropy')])
 defect=trace_distance(channel_iso(v,np.eye(2)/2,2,dc),np.eye(2)/2)
 errors.append(near(defect,.5,'reset defect'))
 return {'name':'balance address deletion','max_error':maxerr(*errors),'normalized_unital_defect':'1/2','unknown_information_map':'|i> -> |1+i>_C','full_SR_isometry_checked':True,'example_carrier_gain':num(gain)}

def exact_phase_gram(l):
 # Each prepared amplitude is sign/sqrt(2*l); mixed material contributes1/2.
 # Integer signs give exact output entries divided by4*l.
 dc=l+2;amps=np.zeros((2,2,dc),dtype=int)
 had=((1,1),(1,-1))
 for initial in range(2):
  for c in range(1,l+1):
   q=c+initial
   for out in range(2):
    phase=(-1)**q if initial==0 else 1
    amps[initial,out,c+initial-out]=had[out][initial]*phase
 return [[F(sum(int(amps[j,i,c])*int(amps[j,k,c]) for j in range(2) for c in range(dc)),4*l) for k in range(2)] for i in range(2)]

def group3():
 l=4;dc=l+2;had=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 w=qblocks(l+1,lambda q:had@np.diag([(-1)**q,1]))
 b=np.zeros(dc);b[1:l+1]=1/math.sqrt(l);v=prepared(w,2,b)
 exact=exact_phase_gram(l)
 req(exact==[[F(1,2),F(-3,8)],[F(-3,8),F(1,2)]],'exact signed Gram')
 output=channel_iso(v,np.eye(2)/2,2,dc)
 errors=[conservation(w,[0,1],range(dc)),near(output,np.array(exact,dtype=float),'phase exact output')]
 for q in range(1,dc):errors.append(near(abs(had@np.diag([(-1)**q,1]))**2,.5,'same probability table'))
 deph=np.diag(np.diag(output));errors.append(near(deph,np.eye(2)/2,'late energy copy erases this witness'))
 return {'name':'population calibration is insufficient','max_error':maxerr(*errors),'exact_output':[[str(x) for x in row] for row in exact],'output_eigenvalues':['1/8','7/8'],'normalized_unital_defect':'3/8','all_transition_probabilities':'1/2','late_energy_copy_defect':'0'}

def group4():
 dc=4;had=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 b=np.eye(dc)[:,1];work=qblocks(3,lambda q:had);v=prepared(work,2,b)
 copied=copy_isometry(v,2,dc,[0,1],2)
 errors=[tp_unital(copied,2,dc*2)]
 # S,K,C ordering, full copy followed by conditional energy-compensated flip.
 d=4*dc;copy=np.zeros((d,d),complex)
 ix=lambda s,k,c:(2*s+k)*dc+c
 for s,k,c in itertools.product(range(2),range(2),range(dc)):
  copy[ix(s,k^s,c),ix(s,k,c)]=1
 fb=np.eye(d,dtype=complex)
 for c in range(dc-1):
  a=ix(0,1,c+1);z=ix(1,1,c)
  fb[a,a]=fb[z,z]=0;fb[a,z]=fb[z,a]=1
 w=fb@copy
 errors.append(conservation(w,[0,0,1,1],range(dc)))
 vsk=prepared(w,4,b);errors.append(tp_unital(vsk,4,dc))
 # Prepare K=0 without silently changing which material was proved unital.
 select=np.zeros((4,2));select[0,0]=1;select[2,1]=1
 vs=vsk@select # S,K,C already correct S then environmentK,C.
 expected=np.zeros_like(vs);expected[ix(0,0,1),0]=1;expected[ix(0,1,2),1]=1
 errors.append(near(vs,expected,'QND feedback exact reset'))
 defect=trace_distance(channel_iso(vs,np.eye(2)/2,2,2*dc),np.eye(2)/2)
 errors.append(near(defect,.5,'feedback material defect'))
 pointer=marginal(copied@np.eye(2)/2@copied.conj().T,2,dc*2,1).reshape(dc,2,dc,2).trace(axis1=0,axis2=2)
 errors.append(near(pointer,np.eye(2)/2,'physical QND record'))
 return {'name':'terminal record and feedback distinction','max_error':maxerr(*errors),'terminal_QND_unital':True,'record_entropy_for_maximally_mixed_input':num(ent(pointer)),'record_blank_distance':'1/2','whole_SK_feedback_unital':True,'reduced_S_feedback_defect':'1/2'}

def first_hit(t):
 out=list(t)
 for i,a in enumerate(t):
  if a in (0,2):out[i]=2-a;return tuple(out)
 return tuple(out)
def first_x(t):
 for a in t:
  if a in (0,2):return int(a==0)
 return 0

def group5(rng):
 strings=list(itertools.product(range(3),repeat=2));ds=len(strings);lookup={s:i for i,s in enumerate(strings)}
 t=[lookup[first_hit(s)] for s in strings];en=np.array([sum(s) for s in strings]);l=5;dc=l+2;bal=2*np.arange(dc)
 w=np.eye(ds*dc,dtype=complex)
 for a in range(ds):
  z=t[a]
  if en[a]<en[z]:
   req(en[z]-en[a]==2,'first hit gap')
   for k in range(1,dc):
    i=a*dc+k;j=z*dc+k-1
    w[i,i]=w[j,j]=0;w[i,j]=w[j,i]=1
 b=np.zeros(dc);b[1:l+1]=1/math.sqrt(l)
 v=prepared(w,ds,b);u=np.eye(ds)[:,t] # t involution: columns map j to t[j]
 errors=[conservation(w,en,bal),tp_unital(v,ds,dc)]
 # All prepared columns, not only diagonal population tests.
 expected=np.zeros_like(v)
 for j in range(ds):
  for k in range(1,l+1):expected[t[j]*dc+k+int((en[j]-en[t[j]])/2),j]=1/math.sqrt(l)
 errors.append(near(v,expected,'1209 exact window amplitude membership'))
 cp=copy_isometry(v,ds,dc,[first_x(s) for s in strings],2)
 errors.append(tp_unital(cp,ds,2*dc))
 # Full S,R and Gram channel checked directly from complete output environment.
 a=[]
 for j in range(ds):a.append(cp.reshape(ds,dc*2,ds)[t[j],:,j])
 a=np.array(a);gram=a@a.conj().T # use explicit index definition below
 gram=np.array([[np.vdot(a[k],a[j]) for k in range(ds)] for j in range(ds)])
 rho=random_density(rng,2*ds)
 actual=apply([np.kron(k,np.eye(2)) for k in kraus(cp,ds,dc*2)],rho)
 rr=rho.reshape(ds,2,ds,2);expected_sr=np.zeros_like(rho)
 for j,k in itertools.product(range(ds),repeat=2):
  expected_sr[2*t[j]:2*t[j]+2,2*t[k]:2*t[k]+2]=gram[j,k]*rr[j,:,k,:]
 errors.append(near(actual,expected_sr,'first hit copied full SR Gram'))
 r=np.array([.5,1/3,1/6]);pn=np.kron(r,r);out=channel_iso(cp,np.diag(pn),ds,dc*2)
 errors.append(near(np.diag(out).real,pn[t],'exact population permutation'))
 return {'name':'1209 finite member with retained copy','max_error':maxerr(*errors),'material_dimension':ds,'carrier_dimension':dc,'pointer_dimension':2,'full_SR_Gram_checked':True,'population_drift':num(trace_distance(out,np.diag(pn))),'amplitude_rule_checked_for_all_safe_columns':True}

def group6(rng):
 dc=6;theta=.05;zeta=.01;had=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 w0=qblocks(5,lambda q:had)
 w=qblocks(5,lambda q:np.diag(np.exp(-1j*theta*np.array([1,-1])))@had if q==2 else had)
 b=np.zeros(dc,complex);b[1:5]=.5
 b_leak=math.sqrt(1-zeta)*b;b_leak[0]=math.sqrt(zeta)
 v0=prepared(w0,2,b);v=prepared(w,2,b);vl=prepared(w,2,b_leak)
 a=2*math.sin(theta/2)
 core_indices=[s*dc+c for s in range(2) for c in range(1,5)]
 norm=float(np.linalg.svd((w-w0)[:,core_indices],compute_uv=False)[0])
 req(norm<=a+TOL and a<theta,'global safe-column norm')
 rho=random_density(rng,6)
 o0=apply([np.kron(k,np.eye(3)) for k in kraus(v0,2,dc)],rho)
 ol=apply([np.kron(k,np.eye(3)) for k in kraus(vl,2,dc)],rho)
 drift=trace_distance(ol,o0)
 defect=trace_distance(channel_iso(vl,np.eye(2)/2,2,dc),np.eye(2)/2)
 bound=a+math.sqrt(zeta)
 req(max(drift,defect)<bound,'leak plus calibration bound')
 req(trace_distance(pures(b),pures(b_leak))<math.sqrt(zeta)+TOL,'projection distance')
 return {'name':'finite calibration and leakage transport','all_checks_passed':True,'angle':theta,'support_leakage':zeta,'safe_operator_error':num(norm),'analytic_half_diamond_upper':num(bound),'simple_upper':'3/20','tested_full_SR_distance':num(drift),'tested_normalized_unital_defect':num(defect),'diamond_bound_from_analytic_isometry_not_numerical_optimization':True}

def group7(rng):
 n=5;dc=n+1;had=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 w=qblocks(n,lambda q:had@np.diag([(-1)**q,1]))
 g=np.zeros_like(w)
 for q in range(1,dc):
  ids=[q,dc+q-1];u=w[np.ix_(ids,ids)]
  vals,vec=np.linalg.eig(u)
  near(vec.conj().T@vec,np.eye(2),'normal block eigenvectors')
  phases=np.mod(-np.angle(vals),2*math.pi)
  g[np.ix_(ids,ids)]=(vec*phases)@vec.conj().T
 g=(g+g.conj().T)/2
 qvals=np.array([s+c for s in range(2) for c in range(dc)],float);q=np.diag(qvals)
 h=q+g
 errors=[near(g@q,q@g,'resident interaction commutes bare energy'),near(exp_h(g),w,'positive logarithm endpoint'),near(exp_h(h),np.diag(np.exp(-1j*qvals))@w,'full free-phase endpoint')]
 req(np.linalg.eigvalsh(g).min()>-TOL and np.linalg.eigvalsh(g).max()<2*math.pi+TOL,'positive finite G')
 rho=random_density(rng,2*dc);out=exp_h(h)@rho@exp_h(h).conj().T
 errors.append(near(np.trace(g@out),np.trace(g@rho),'resident G expectation retained'))
 req(np.max(abs(w-w.conj().T))>.1,'noninvolution witness')
 return {'name':'positive autonomous endpoint with resident interaction','max_error':maxerr(*errors),'G_min':num(np.linalg.eigvalsh(g).min()),'G_max':num(np.linalg.eigvalsh(g).max()),'bare_energy_and_G_conserved':True,'free_phase_retained':True,'not_spatial_clock_or_QQ_implementation':True}

def compute():
 rng=np.random.default_rng(1215007)
 gs=[group1(rng),group2(rng),group3(),group4(),group5(rng),group6(rng),group7(rng)]
 return {'schema':'1215-finite-balance-window-author-v1','status':'PASS','all_checks_passed':True,'check_groups':len(gs),'groups':gs,'code_sha256':sha(Path(__file__)),'source_sha256':{s:sha(RESEARCH/s) for s in SOURCES},'scope':'Finite mathematical window contract and deletion checks. No adoption of new cognitive axioms; no native QQ, spatial apparatus, or four-law completion certificate.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
 dest=HERE/'results.json'
 if args.write and dest.exists():raise SystemExit('--write refused: results.json already exists')
 result=compute()
 if args.write:
  dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:
  req(dest.exists(),'results.json missing; initial --write required')
  old=json.loads(dest.read_text(encoding='utf-8-sig'))
  req(result==old,'saved results/hash mismatch')
 print(json.dumps({'status':'PASS','all_checks_passed':True,'results_match':True,'check_groups':result['check_groups']},ensure_ascii=False))

if __name__=='__main__':main()
