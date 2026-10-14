"""1218 author replay, finite complete matrices; no prior code imports.
--write creates results once. Default checks saved values and source hashes only.
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
TOL=9e-11
SOURCES=[
 'archive_370_428/research_note_421.md',
 'archive_429_466/research_note_457.md',
 'archive_1009_1043/research_note_1037.md',
 'archive_1201_/1209/proof.md',
 'archive_1201_/1214/proof.md',
 'archive_1201_/1215/proof.md','archive_1201_/1215/results.json',
 'archive_1201_/1216/proof.md','archive_1201_/1216/results.json',
]

def req(c,m):
 if not bool(c):raise AssertionError(m)
def near(a,b,m):
 z=float(np.max(abs(np.asarray(a)-np.asarray(b))))
 req(z<TOL,m+': '+str(z));return z
def number(x):return round(float(x),13)
def mx(*a):return number(max(a))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pure(v):return np.outer(v,v.conj())
def td(a,b):return float(np.abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)).sum()/2)
def entropy(r):
 v=np.linalg.eigvalsh((r+r.conj().T)/2);req(v.min()>-TOL,'positive density')
 v=v[v>1e-14];return float(-np.dot(v,np.log(v)))
def partial(r,ds,dc,side=0):
 t=r.reshape(ds,dc,ds,dc)
 return np.trace(t,axis1=1,axis2=3) if side==0 else np.trace(t,axis1=0,axis2=2)
def prep(w,ds,b):return w@np.kron(np.eye(ds),b.reshape(-1,1))
def kraus(v,ds,dc):
 t=v.reshape(ds,dc,ds);return [t[:,c,:] for c in range(dc)]
def apply(ks,r):return sum((k@r@k.conj().T for k in ks),np.zeros_like(r,dtype=complex))
def reduced(v,r,ds,dc):return partial(v@r@v.conj().T,ds,dc)
def unital(v,ds,dc):
 ks=kraus(v,ds,dc)
 return mx(near(sum(k.conj().T@k for k in ks),np.eye(ds),'TP'),near(sum(k@k.conj().T for k in ks),np.eye(ds),'unital'))
def density(rng,d):
 z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));r=z@z.conj().T
 return r/np.trace(r)
def safety_expected(u,energies,b):
 ds=len(energies);dc=len(b);v=np.zeros((ds*dc,ds),complex)
 for j,i,c in itertools.product(range(ds),range(ds),range(dc)):
  if abs(u[i,j]*b[c])>1e-15:
   k=c+int(energies[j]-energies[i]);req(0<=k<dc,'safe target exists')
   v[i*dc+k,j]+=u[i,j]*b[c]
 return v

def qubit_w(dc,rule):
 w=np.eye(2*dc,dtype=complex)
 for e in range(1,dc):
  ix=[e,dc+e-1];w[np.ix_(ix,ix)]=rule(e)
 return w

def contact_w(dc,target,rule):
 # Basis A,B,C. Both singleton energy boundaries remain part of W.
 w=np.eye(4*dc,dtype=complex)
 for other in range(2):
  for e in range(1,dc):
   if target==0:ix=[other*dc+e,(2+other)*dc+e-1]
   else:ix=[2*other*dc+e,(2*other+1)*dc+e-1]
   w[np.ix_(ix,ix)]=rule(e)
 return w

def energy_error(w,energies,dc):
 q=np.array([e+c for e in energies for c in range(dc)])
 return mx(near(w.conj().T@w,np.eye(len(q)),'complete unitary'),near((q[:,None]-q[None,:])*w,0,'complete energy'))

def main_setup():
 dc=7;en=[0,1,1,2]
 h=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 u=np.array([[1,1j],[1j,1]],complex)/math.sqrt(2)
 a=contact_w(dc,0,lambda e:h);b=contact_w(dc,1,lambda e:u)
 carrier=np.zeros(dc,complex);carrier[2]=1/math.sqrt(2);carrier[3]=1j/math.sqrt(2)
 return dc,en,h,u,a,b,carrier

def group1(rng):
 dc,en,ua,ub,wa,wb,b=main_setup();u=np.kron(ua,ub);w=wb@wa;v=prep(w,4,b)
 e=[energy_error(wa,en,dc),energy_error(wb,en,dc),unital(v,4,dc),near(v,safety_expected(u,en,b),'whole safe isometry')]
 rho=density(rng,12)
 actual=apply([np.kron(k,np.eye(3)) for k in kraus(v,4,dc)],rho)
 expected=apply([np.kron(k,np.eye(3)) for k in kraus(safety_expected(u,en,b),4,dc)],rho)
 e.append(near(actual,expected,'full unknown correlated ABR'))
 return {'name':'safe cross-material closure','max_error':mx(*e),'carrier_dimension':dc,'material_dimension':4,'unknown_R_dimension':3,'initial_materials_need_not_be_independent':True,'same_carrier_not_reset':True}

def group2(rng):
 dc,en,ua,ub,wa,wb,b=main_setup();u=np.kron(ua,ub);v=prep(wb@wa,4,b);n=19
 p=[];us=[]
 for k in range(n):
  theta=2*math.pi*k/n;f=np.exp(1j*np.arange(dc)*theta)
  pk=float(np.vdot(f,pure(b)@f).real/n);p.append(pk)
  d=np.diag(np.exp(1j*np.array(en)*theta));us.append(d@u@d.conj().T)
 e=[near(sum(p),1,'phase probability normalization')];req(min(p)>-TOL,'phase positivity from b')
 for j,k in itertools.product(range(4),repeat=2):
  r=np.zeros((4,4),complex);r[j,k]=1
  expected=sum((pk*uk@r@uk.conj().T for pk,uk in zip(p,us)),np.zeros_like(r))
  e.append(near(reduced(v,r,4,dc),expected,'phase all matrix units'))
 rho=density(rng,8)
 actual=apply([np.kron(k,np.eye(2)) for k in kraus(v,4,dc)],rho)
 expected=sum((pk*np.kron(uk,np.eye(2))@rho@np.kron(uk.conj().T,np.eye(2)) for pk,uk in zip(p,us)),np.zeros_like(rho))
 e.append(near(actual,expected,'phase full unknown R'))
 return {'name':'positive common-phase representation','max_error':mx(*e),'phase_terms':n,'minimum_weight':number(min(p)),'weight_sum':number(sum(p)),'same_phase_for_all_contacts':True,'not_a_new_physical_register':True}

def group3():
 dc=5;h=np.array([[1,1],[1,-1]],complex)/math.sqrt(2);w=qubit_w(dc,lambda e:h);b=np.eye(dc)[:,2]
 v1=prep(w,2,b);v2=prep(w@w,2,b);r=np.diag([1.,0.])
 first=reduced(v1,r,2,dc);actual=reduced(v2,r,2,dc);reset=apply(kraus(v1,2,dc),first)
 e=[near(w@w,np.eye(2*dc),'complete involution'),near(first,np.eye(2)/2,'first loses local record'),near(actual,r,'same carrier echo'),near(reset,np.eye(2)/2,'hypothetical reset channel'),unital(v2,2,dc)]
 return {'name':'shared carrier cannot be resampled','max_error':mx(*e),'same_carrier_final_zero_probability':'1','reset_comparison_final_zero_probability':'1/2','output_distance':'1/2','not_a_new_claim_about_1216_actual_A':True}

def copy_gate(ds,dc,labels,nk):
 u=np.zeros((ds*dc*nk,ds*dc*nk),complex)
 for s,c,k in itertools.product(range(ds),range(dc),range(nk)):
  out=(s*dc+c)*nk+(k^labels[s]);col=(s*dc+c)*nk+k;u[out,col]=1
 return u

def group4():
 dc=5;h=np.array([[1,1],[1,-1]],complex)/math.sqrt(2);w=qubit_w(dc,lambda e:h);copy=copy_gate(2,dc,[0,1],2)
 env=np.kron(np.eye(dc)[:,2],np.array([1.,0.]));lift=np.kron(w,np.eye(2))
 mid=lift@copy@lift;terminal=copy@lift@lift;vm=prep(mid,2,env);vt=prep(terminal,2,env)
 e=[unital(vm,2,dc*2),unital(vt,2,dc*2),near(mid.conj().T@mid,np.eye(4*dc),'record complete unitary')]
 q=np.array([s+c for s in range(2) for c in range(dc) for k in range(2)])
 e.append(near((q[:,None]-q[None,:])*mid,0,'record complete energy'))
 r=np.diag([1.,0.]);m=reduced(vm,r,2,dc*2);t=reduced(vt,r,2,dc*2)
 e.extend([near(m,np.eye(2)/2,'mid record destroys echo'),near(t,r,'late record retains zero event')])
 envout=partial(vm@r@vm.conj().T,2,dc*2,1)
 kp=envout.reshape(dc,2,dc,2).trace(axis1=0,axis2=2)
 e.append(near(kp,np.eye(2)/2,'pointer retained and used'))
 return {'name':'actual middle nonfeedback QND','max_error':mx(*e),'middle_record_zero_probability':'1/2','terminal_record_zero_probability':'1','both_material_channels_unital':True,'pointer_entropy':number(entropy(kp)),'pointer_blank_distance':'1/2','all_records_retained':True}

def exact_joint_output(w,dc,support):
 # All selected initial carrier amplitudes1/sqrt(len(support)); materialI/4.
 perm=np.argmax(w,axis=0);req(np.array_equal(w,np.eye(4*dc)[:,perm]),'integer complete permutation')
 raw=np.zeros((4*dc,4*dc),dtype=object);raw[:]=0
 for s in range(4):
  dest=[int(perm[s*dc+c]) for c in support]
  for a,b in itertools.product(dest,repeat=2):raw[a,b]+=1
 den=4*len(support)
 ab=[[F(sum(int(raw[s*dc+c,t*dc+c]) for c in range(dc)),den) for t in range(4)] for s in range(4)]
 car=[[F(sum(int(raw[s*dc+c,s*dc+d]) for s in range(4)),den) for d in range(dc)] for c in range(dc)]
 return ab,car,raw,den

def group5():
 dc=5;x=np.array([[0,1],[1,0]],complex);wa=contact_w(dc,0,lambda e:x);wb=contact_w(dc,1,lambda e:x);w=(wb@wa).real.astype(int)
 ab,car,raw,den=exact_joint_output(w,dc,[1,3]);expected=[F(1,8),F(3,8),F(3,8),F(1,8)]
 req(ab==[[expected[i] if i==j else F(0) for j in range(4)] for i in range(4)],'exact joint deletion diagonal')
 ce=[[F(0) for j in range(dc)] for i in range(dc)]
 ce[0][0]=ce[4][4]=F(1,8);ce[1][1]=ce[3][3]=F(3,8);ce[1][3]=ce[3][1]=F(1,4)
 req(car==ce,'exact carrier final state')
 req(sum(abs(expected[i]-F(1,4)) for i in range(4))/2==F(1,4),'exact unital defect')
 req(expected[0]+expected[1]==expected[2]+expected[3]==F(1,2),'A unchanged')
 req(expected[0]+expected[2]==expected[1]+expected[3]==F(1,2),'B unchanged')
 req(sum(F(c)*car[c][c] for c in range(dc))==2,'C mean energy unchanged')
 safe=np.eye(dc)[:,2];vs=prep(w,4,safe);e=[unital(vs,4,dc)]
 b=np.zeros(dc);b[1]=b[3]=1/math.sqrt(2)
 # Each contact is separately certified with the ORIGINAL preparation.
 e.extend([unital(prep(wa,4,b),4,dc),unital(prep(wb,4,b),4,dc),energy_error(w,[0,1,1,2],dc)])
 globalout=np.array(raw,dtype=float)/den;af=np.array(ab,dtype=float);cf=np.array(car,dtype=float)
 e.append(near(entropy(globalout),math.log(4),'complete entropy'))
 req(entropy(af)<math.log(4) and entropy(cf)>0,'purity and complete resource ledger')
 return {'name':'same finite carrier path deletion','max_error':mx(*e),'exact_joint_probabilities':[str(v) for v in expected],'normalized_unital_defect':'1/4','different_energy_labels_initial_probability':'1/2','different_energy_labels_final_probability':'3/4','A_B_means_each':'1/2','C_initial_and_final_mean':'2','carrier_exact_matrix':[[str(v) for v in row] for row in car],'carrier_exact_spectrum':['5/8','1/8','1/8','1/8','0'],'material_final_entropy':number(entropy(af)),'carrier_final_entropy':number(entropy(cf)),'total_final_entropy':number(entropy(globalout)),'AB_C_final_mutual_information':number(entropy(af)+entropy(cf)-entropy(globalout)),'fresh_preparation_each_stage_is_unital':True}

def group6(rng):
 dc=5;x=np.array([[0,1],[1,0]],complex);wa=contact_w(dc,0,lambda e:x);vb=contact_w(dc,1,lambda e:x)
 # P excludes precisely B,C singleton boundary columns, keeps A unrestricted.
 diagonal=[int(1<=b+c<=4) for a in range(2) for b in range(2) for c in range(dc)];p=np.diag(diagonal)
 safe=np.eye(dc)[:,2];edge=(np.eye(dc)[:,1]+np.eye(dc)[:,3])/math.sqrt(2)
 a_safe=prep(wa,4,safe);a_edge=prep(wa,4,edge)
 le_safe=np.linalg.svd((np.eye(4*dc)-p)@a_safe,compute_uv=False)[0]
 le_edge=np.linalg.svd((np.eye(4*dc)-p)@a_edge,compute_uv=False)[0]
 e=[near(le_safe,0,'safe whole reachable subspace'),near(le_edge**2,.5,'edge all-input leakage')]
 # Exact rational rotation, so analytic operator error^2=4/2501<1/625.
 rotation=np.array([[2499,-100],[100,2499]],float)/2501
 req(F(2499,2501)**2+F(100,2501)**2==1,'exact rotation unitary')
 wb=contact_w(dc,1,lambda q:rotation@x if q==2 else x)
 norm=np.linalg.svd((wb-vb)@p,compute_uv=False)[0]
 e.append(near(norm**2,4/2501,'calibration exact norm squared'))
 req(F(4,2501)<F(1,625),'strict one twenty-fifth bound')
 v=prep(wb@wa,4,safe);v0=prep(vb@wa,4,safe)
 err=np.linalg.svd(v-v0,compute_uv=False)[0]
 req(err<=norm+TOL,'telescoping full prepared isometry')
 rho=density(rng,12)
 actual=apply([np.kron(k,np.eye(3)) for k in kraus(v,4,dc)],rho)
 ideal=apply([np.kron(k,np.eye(3)) for k in kraus(v0,4,dc)],rho)
 drift=td(actual,ideal);defect=td(reduced(v,np.eye(4)/4,4,dc),np.eye(4)/4)
 req(max(drift,defect)<1/25,'full SR and normalized unital error')
 return {'name':'reachable operator calibration and path error','max_error':mx(*e),'safe_eta_squared':'0','deleted_path_eta_squared':'1/2','calibration_a_squared':'4/2501','strict_half_diamond_upper':'1/25','prepared_isometry_error':number(err),'tested_full_SR_distance':number(drift),'tested_normalized_unital_defect':number(defect),'uses_operator_supremum_not_average_path_probability':True}

def group7():
 # The theorem permits arbitrary initial correlations inside materials/R,
 # but not C with those unknown inputs. An explicit deletion of that distinction.
 dc=5;h=np.array([[1,1],[1,-1]],complex)/math.sqrt(2);w=qubit_w(dc,lambda e:h)
 psi=np.zeros(2*dc);psi[2]=psi[dc+1]=1/math.sqrt(2)
 initial=pure(psi);out=w@initial@w.conj().T
 si=partial(initial,2,dc);ci=partial(initial,2,dc,1);so=partial(out,2,dc)
 independent=w@np.kron(si,ci)@w.conj().T
 e=[near(si,np.eye(2)/2,'correlated input marginal'),near(so,np.diag([1,0]),'correlation-assisted reset'),near(partial(independent,2,dc),np.eye(2)/2,'same marginals independent')]
 req(abs(float(np.trace(initial@initial).real)-1)<TOL,'correlated complete pure input')
 return {'name':'initial carrier independence is not recreated midway','max_error':mx(*e),'initial_material_marginal':'I/2','correlated_final_material':'|0><0|','independent_same_marginals_final_material':'I/2','not_a_counterexample_to_safe_independent_protocol':True}

def group8(rng):
 dc=7;en=[0,1,1,2];h=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 last=np.array([[1,1j],[1j,1]],complex)/math.sqrt(2)
 ws=[contact_w(dc,0,lambda e:h),contact_w(dc,1,lambda e:h),contact_w(dc,0,lambda e:last)]
 w=ws[2]@ws[1]@ws[0];u=np.kron(last@h,h);b=np.eye(dc)[:,3];v=prep(w,4,b)
 e=[energy_error(stage,en,dc) for stage in ws]
 e.extend([near(v,safety_expected(u,en,b),'three segment telescoping'),unital(v,4,dc)])
 # An initially entangled AB state and an unknown 4-dimensional R are allowed.
 bell=np.zeros(16,complex)
 for j in range(4):bell[4*j+j]=.5
 rho=pure(bell);sr=apply([np.kron(k,np.eye(4)) for k in kraus(v,4,dc)],rho)
 e.append(near(partial(sr,4,4),np.eye(4)/4,'maximally mixed material image of entangled input'))
 for stage in ws:
  for a,bbit,c in itertools.product(range(2),range(2),range(dc)):
   req(abs(np.vdot(stage[:,(2*a+bbit)*dc+c],stage[:,(2*a+bbit)*dc+c])-1)<TOL,'all boundary columns retained')
 return {'name':'all boundaries and longer same-carrier protocol','max_error':mx(*e),'number_of_contacts':3,'all_material_dimension':4,'carrier_dimension':dc,'same_A_contacted_again':True,'all_complete_boundary_columns_checked':True,'arbitrary_material_reference_channel_checked':True}

def compute():
 rng=np.random.default_rng(1218008)
 groups=[group1(rng),group2(rng),group3(),group4(),group5(),group6(rng),group7(),group8(rng)]
 return {'schema':'1218-cross-segment-safe-window-v1','status':'PASS','all_checks_passed':True,'check_groups':len(groups),'groups':groups,'code_sha256':sha(Path(__file__)),'source_sha256':{s:sha(RESEARCH/s) for s in SOURCES},'scope':'Finite endpoint qualification and path deletion. Known lift homomorphism reused. No Markov, controller, four-law completion, or uniform time/resource assumption.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();dest=HERE/'results.json'
 if a.write and dest.exists():raise SystemExit('--write refused: results already exist')
 r=compute()
 if a.write:dest.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:
  req(dest.exists(),'no saved results; initial --write required')
  req(json.loads(dest.read_text(encoding='utf-8-sig'))==r,'saved result/source hash mismatch')
 print(json.dumps({'status':'PASS','all_checks_passed':True,'results_match':True,'check_groups':r['check_groups']}))
if __name__=='__main__':main()
