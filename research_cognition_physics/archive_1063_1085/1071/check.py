"""Round1071 author finite certificates. Read-only, no geometry-generation claim."""
from fractions import Fraction as Q
import numpy as np
import json
checks=0;res=[]
def eq(a,b):
 global checks
 checks+=1
 assert a==b,(a,b)
def near(a,b):
 global checks
 checks+=1;v=float(np.max(np.abs(np.asarray(a)-np.asarray(b))));res.append(v)
 assert v<2e-12,v
# Exact radial iteration, no uniform q and no precise half.
radial=[]
for r in [Q(1,2),Q(1,7),Q(1,100),Q(1,10000)]:
 y=r
 for k in range(1,21):
  y=y/(1+y);eq(y,r/(1+k*r));eq(y<=Q(1,k+1),True)
 ratio=1/(1+r);eq(ratio<Q(1),True)
 radial.append({'initial':str(r),'first_ratio':str(ratio),'after20':str(y)})
eq(2*Q(1,2)/(1+Q(1,2))==Q(1,2),False)
def R(x):return x/(1+np.linalg.norm(x))
def B(x,s):return (1-s+s/(1+np.linalg.norm(x)))*x
def F(x,t):
 if t==1:return np.zeros_like(x)
 k=int(np.floor(-np.log2(1-t)))
 lo=1-2.**(-k);hi=1-2.**(-k-1)
 y=x/(1+k*np.linalg.norm(x))
 return B(y,(t-lo)/(hi-lo))
for dim in [1,2,3,5]:
 x=np.arange(1,dim+1,dtype=float);x=x/np.linalg.norm(x)*.7
 near(B(x,0),x);near(B(x,1),R(x))
 for k in [0,1,2,7,15]:
  y=x/(1+k*np.linalg.norm(x))
  near(B(y,1),x/(1+(k+1)*np.linalg.norm(x)))
  near(F(x,1-2.**(-k)),y)
  for s in [.25,.5,.75]:
   t=(1-2.**(-k))+s*2.**(-k-1)
   near(F(x,t),B(y,s))
 near(F(x,1),np.zeros(dim))
 # For a fixed prepared x,s, a translation implements B(x,s), but is not R globally.
 v=B(x,.4)-x; z=-x/2
 near(np.linalg.norm((x+v)-(z+v)),np.linalg.norm(x-z))
# Existing425 shift model: exact finite-support calculation of the infinite sum.
for ts in [[Q(1,3)],[Q(0),Q(1,2)],[Q(1,7),Q(2,5),Q(0),Q(1,11)]]:
 C=lambda vs:sum((Q(1,2**j)*min(v%1,1-v%1) for j,v in enumerate(vs,1)),Q(0))
 eq(C([Q(0)]+ts),C(ts)/2)
# A tail circle has degree1 before shifting and degree0 in its same projection after.
angles=np.linspace(0,2*np.pi,129)
z=np.exp(1j*angles)
degree=float(np.sum(np.angle(z[1:]/z[:-1]))/(2*np.pi));near(degree,1)
near(float(np.sum(np.angle(np.ones(128,dtype=complex)))/(2*np.pi)),0)
print(json.dumps({'round':1071,'passed':True,'checks':checks,'max_residual':max(res),
 'radial_cases':radial,'dimensions_checked':[1,2,3,5],'tail_projection_degrees':[1,0],
 'universal_topology_proved_by_sampling':False,'space_dimension_three_derived':False,
 'physical_finite_time_infinite_steps_claimed':False},ensure_ascii=False,indent=2,sort_keys=True))
