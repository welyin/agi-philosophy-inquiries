import sys, itertools, math
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, '.')
import rigid_leaf_reference_audit as base
P=list(itertools.permutations(range(6))); ix={p:i for i,p in enumerate(P)}; ident=ix[tuple(range(6))]
trees,_,fg=base.system()
edges=sorted(set().union(*trees))
masks=[np.array([int(e in g) for g in trees],dtype=np.int64) for e in edges]
left=[];right=[]
for a,b in edges:
 s=list(range(6));s[a],s[b]=s[b],s[a]
 left.append(np.array([ix[tuple(s[p[k]] for k in range(6))] for p in P]))
 right.append(np.array([ix[tuple(p[s[k]] for k in range(6))] for p in P]))
def ad(z):
 out=fg@z-z@fg
 for l,r,n in zip(left,right,masks):
  out[l]+=n[None,:,None]*z
  out[r]-=z*n[None,None,:]
 return out
I=np.eye(2,dtype=complex)
pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1])]
def etr(p,axis):
 a=[I,I if axis<0 else I+pauli[axis],I,I+pauli[0],I+pauli[1],I+pauli[2]]
 seen=set();total=1+0j
 for i in range(6):
  if i in seen:continue
  cyc=[];j=i
  while j not in seen:
   cyc.append(j);seen.add(j);j=p[j]
  mat=I.copy()
  for k in reversed(cyc):mat=mat@a[k]
  total*=np.trace(mat)
 return complex(round(total.real),round(total.imag))
tr=[[etr(p,a) for p in P] for a in [-1,0,1,2]]
o=3*base.distance(0,3)-8
z=np.zeros((720,6,6),dtype=np.int64);z[ident]=np.diag(o)
for k in range(11):
 sums=z.sum(axis=(1,2)).astype(object)
 vals=[]
 for a,b in [(-1,-1),(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]:
  num=sum(int(v)*tr[a+1][j]*tr[b+1][j] for j,v in enumerate(sums))*(1j**k)
  assert abs(num.imag)<1e-5
  val=F(int(num.real),3*4096*6*math.factorial(k))
  vals.append(str(val))
 print(k,vals,flush=True)
 z=ad(z)
