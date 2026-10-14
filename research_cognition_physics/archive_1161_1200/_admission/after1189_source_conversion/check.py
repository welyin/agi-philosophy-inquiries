"""Unnumbered source-conversion admission; exact low-block identity, not a new GR result."""
from pathlib import Path
import argparse,importlib.util,json,math,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'1189'/'check.py'
spec=importlib.util.spec_from_file_location('accepted1189',OLD)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def calculate():
 masks,hm,v,js,_=old.model();k=.1;nf=len(hm);M=len(js)
 states=[tuple(n) for n in __import__('itertools').product(range(3),repeat=M) if sum(n)<=2]
 index={n:i for i,n in enumerate(states)};nb=len(states);qs=[]
 for a in range(M):
  q=np.zeros((nb,nb),complex)
  for col,n in enumerate(states):
   for step in [-1,1]:
    t=list(n);t[a]+=step;t=tuple(t)
    if min(t)>=0 and t in index:q[index[t],col]=math.sqrt(n[a]+(step==1))/math.sqrt(2)
  qs.append(q)
 hb=np.diag([sum(n)+M/2 for n in states]);B=sum(np.kron(q,j)/2 for q,j in zip(qs,js))
 A=np.kron(np.eye(nb),hm)+np.kron(hb,np.eye(nf));V=np.kron(np.eye(nb),v)
 H=A+k*k*V-k*B;S=2*k*V-B;C=H@S-S@H
 z=index[(0,)*M];cols=list(range(z*nf,(z+1)*nf));D=hm-k*k*v;rows=[]
 for a,j in enumerate(js):
  n=[0]*M;n[a]=1;r=index[tuple(n)];block=C[r*nf:(r+1)*nf,cols]
  comm=D@j-j@D;formula=-(comm+j)/(2*math.sqrt(2))
  square=(np.linalg.norm(comm,'fro')**2+np.linalg.norm(j,'fro')**2)/8
  rows.append(dict(mode=a,block_frobenius_norm=float(np.linalg.norm(block,'fro')),square_formula=float(square),matrix_identity_residual=float(np.max(abs(block-formula))),square_identity_residual=float(abs(np.linalg.norm(block,'fro')**2-square)),cross_trace_abs=float(abs(np.trace(comm.conj().T@j)))))
  assert rows[-1]['matrix_identity_residual']<1e-13 and rows[-1]['square_identity_residual']<1e-13
  assert rows[-1]['block_frobenius_norm']>.1
 # Diagnostic only: source invariance for all energy-conserving unitaries need not imply affinity.
 h=np.diag([0.,1.,1.,2.]);j=h@h;swap=np.eye(4);swap[[1,2]]=swap[[2,1]]
 assert np.linalg.norm(swap.T@j@swap-j)==0
 residual=np.linalg.norm(np.diag(j)-np.column_stack([np.diag(h),np.ones(4)])@np.linalg.lstsq(np.column_stack([np.diag(h),np.ones(4)]),np.diag(j),rcond=None)[0])
 assert residual>.9
 return dict(status='PASS_ADMISSION_SCOPE',old_source_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),low_block_exact_scope='vacuum to one quantum; all intermediate occupations included; not a full-space truncation estimate',matrix_dimension=len(H),kappa=k,blocks=rows,all_conserved_conversions_do_not_force_affine_source_example=dict(H=[0,1,1,2],J=[0,1,1,4],affine_least_squares_residual=float(residual)),new_research_round=False,new_axioms=0,goal_complete=False)
def same(a,b):
 if isinstance(a,dict):
  assert a.keys()==b.keys()
  for k in a:same(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):same(x,y)
 elif isinstance(a,float):assert abs(a-b)<1e-12,(a,b)
 else:assert a==b,(a,b)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=calculate();out=HERE/'results.json'
 if args.write:
  with out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
 else:same(r,json.loads(out.read_text(encoding='utf-8')))
 print(json.dumps(dict(status=r['status'],source_commutator_low_block_norms=[x['block_frobenius_norm'] for x in r['blocks']],new_round=False)))
