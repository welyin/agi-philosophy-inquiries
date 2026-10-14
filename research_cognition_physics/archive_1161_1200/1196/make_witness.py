from pathlib import Path
import numpy as np,json,argparse
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent
z=F(989835554078098,10**15);mu=2*F(16,65)*(1+F(3,20)*z)/(1+F(7,20)*z);M=32
B=np.zeros((2*M,2*M))
for j in range(M):
 B[2*j,2*j]=float(mu)*(j+.5);B[2*j+1,2*j+1]=-float(mu)*(j+.5)
 if j<M-1:
  for s in (0,1):B[2*j+s,2*(j+1)+1-s]=B[2*(j+1)+1-s,2*j+s]=(j+1)/2
e,u=np.linalg.eigh(B);j=np.flatnonzero(e>0)[0];w=u[:,j]
if w[np.argmax(abs(w))]<0:w=-w
out={'sites':M,'mu':str(mu),'eigenvalue':str(F(float(e[j]))),'vector':[str(F(float(t)))for t in w],
 'purpose':'Frozen rational test mode. Approximate eigenvalue/eigenvector are not trusted without certificate.py residual bounds.'}
p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
if a.write:
 with (HERE/'mode_witness.json').open('x',encoding='utf-8')as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps({'sites':M,'eigenvalue':float(e[j]),'mu':float(mu),'wrote':a.write}))
