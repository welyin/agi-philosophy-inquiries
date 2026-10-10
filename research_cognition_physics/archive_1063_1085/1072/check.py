"""1072: same-source consensus revisions versus the fixed three-qubit relation code.
No writes on execution. Existing NumPy only. Integer commutators certify a finite window.
"""
import json,math
from fractions import Fraction as F
import numpy as np
I=np.eye(2,dtype=complex)
sig=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
k=np.kron
checks=0;residual=0.
def eq(a,b,tol=3e-12):
 global checks,residual
 e=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
 residual=max(residual,e);assert e<tol,e;checks+=1

def swap(a,b,count=3):
 m=np.zeros((2**count,2**count),dtype=np.int64)
 for i in range(2**count):
  j=i
  if ((i>>(count-1-a))&1)!=((i>>(count-1-b))&1):j=i^(1<<(count-1-a))^(1<<(count-1-b))
  m[j,i]=1
 return m
S=[swap(0,1),swap(0,2),swap(1,2)]
Ps=sum(S)/3
P=np.eye(8)-Ps
Sab=(np.eye(4)-swap(0,1,2))/2
Tab=np.eye(4)-Sab
Qs=k(Sab,I);Qt=P-Qs
Ds=[I,*[-1j*s for s in sig]]
Vs=[k(k(d,d),I) for d in Ds]
eq(Ps@Ps,Ps);eq(P@Qs,Qs);eq(Qt@Qt,Qt)
eq(sum(v.conj().T@Ps@v for v in Vs)/4,2*k(Tab,I)/3)
eq(sum(v.conj().T@Ps@v for v in Vs[1:]),8*k(Tab,I)/3-Ps)

def ptc(r):return np.trace(r.reshape(4,2,4,2),axis1=1,axis2=3)
def unitary(n,angle):return math.cos(angle/2)*I-1j*math.sin(angle/2)*np.einsum('i,ijk->jk',n,sig)
def leak(v,r):return float(np.trace(Ps@v@r@v.conj().T).real)
rows=[]
for w in [0.,.25,1/3,.7,1.]:
 ps=(1+3*w)/4;pt=3*(1-w)/4
 rho=ps*Qs/2+pt*Qt/2
 omega=ps*Sab+pt*Tab/3
 eq(np.trace(rho),1);eq(P@rho@P,rho);eq(ptc(rho),omega)
 assert np.linalg.eigvalsh(rho).min()>-1e-12;checks+=1
 values=[leak(v,rho) for v in Vs]
 eq(values[0],0);eq(values[1:],np.full(3,2*(1-w)/3))
 for n in [*np.eye(3),np.ones(3)/np.sqrt(3),np.array([1.,2.,-2.])/3]:
  for angle in [.17,.7,1.9,math.pi]:
   u=unitary(n,angle);v=k(k(u,u),I)
   eq(leak(v,rho),(1-w)*(1-math.cos(angle))/3)
   eq(P@v.conj().T@Ps@v@P,(8/9)*math.sin(angle/2)**2*Qt)
   eq(ptc(v@rho@v.conj().T),omega)
 rows.append({'w':w,'pauli_leakages':values,'sharp_worst_case':2*(1-w)/3})
# Original round433 model; no imported archived code and no new coupling settings.
def embed(q,j):
 return k(k(q if j==0 else I,q if j==1 else I),q if j==2 else I)
ns=[embed(np.diag([0,1]),j).real.astype(np.int64) for j in range(3)]
xs=[embed(sig[0],j).real.astype(np.int64) for j in range(3)]
H=sum(k(np.eye(8,dtype=np.int64),xx+nn)+k(ss,nn) for ss,xx,nn in zip(S,xs,ns))
N=sum(ns);O=k(np.eye(8,dtype=np.int64),N)
assert np.array_equal(H,H.T);checks+=1
# In int64, all intermediate products through order14 are bounded by 3*18**14<2**63.
assert 3*18**14<2**63
ad=O.copy(); coeff={};poly=F(0);t=F(1,8)
for order in range(15):
 a2=2*int(ad[0,0])
 b2=int(ad[16,16])+int(ad[32,32])-int(ad[16,32])-int(ad[32,16])
 if order%2:
  assert a2-b2==0;checks+=1
 else:
  c=F(((-1)**(order//2))*(a2-b2),6*math.factorial(order))
  coeff[order]=c;poly+=c*t**order
 if order<14:ad=H@ad-ad@H
assert coeff[0]==coeff[2]==0 and coeff[4]==F(-1,6);checks+=1
# Difference of two expectation remainders; exp(18t)<exp(3)<27.
tail=2*27*(18*t)**15/F(math.factorial(15))
lo,hi=poly-tail,poly+tail
assert hi<0;checks+=1
gap_lower=-hi/2
# Spectral calculation is a separate finite floating-point check of the exact certificate.
eig,U=np.linalg.eigh(H.astype(float));V=(U*np.exp(-1j*float(t)*eig))@U.T
E=(V.conj().T@(O/3)@V)[::8,::8]
a=float(np.trace(Ps@E).real/4);b=float(np.trace(P@E).real/4)
eq(E,a*Ps+b*P)
assert float(lo)<a-b<float(hi);checks+=1
w=.25;rho=(1+3*w)/4*Qs/2+3*(1-w)/4*Qt/2
p_before=float(np.trace(rho@E).real)
p_after=[]
for v in Vs[1:]:
 r=v@rho@v.conj().T
 pp=float(np.trace(r@E).real);p_after.append(pp)
 eq(pp,b+(a-b)/2)
 assert abs(pp-p_before)>float(gap_lower);checks+=1
# Weak initial-support version follows from the same exact identity (not a new optimisation).
arbitrary=k((1+3*w)/4*Sab+3*(1-w)/4*Tab/3,I/2)
l0=leak(Vs[0],arbitrary)
eq(sum(leak(v,arbitrary) for v in Vs[1:]),2*(1-w)-l0)
print(json.dumps({'round':1072,'checks':checks,'maximum_residual':residual,
 'operator_identity_all_inputs':True,'sharp_minimax_formula':'2*(1-w)/3','rows':rows,
 'finite_contact_readout':{'time':str(t),'coefficient_polynomial':{str(j):str(c) for j,c in coeff.items()},
 'polynomial_value':str(poly),'remainder_bound':str(tail),'contrast_interval':[str(lo),str(hi)],
 'contrast_interval_float':[float(lo),float(hi)],'certified_gap_for_w_1_4':str(gap_lower),
 'certified_gap_float':float(gap_lower),'numeric_a':a,'numeric_b':b,'numeric_contrast':a-b,
 'initial_probability':p_before,'post_pi_probabilities':p_after},
 'spatial_dimension_three_derived':False,'general_cognition_refuted':False},ensure_ascii=False,indent=2))
