"""Round 1188: exact same-event slicing checks; no physical instrument simulation."""
from fractions import Fraction as F
from pathlib import Path
import argparse, json
import numpy as np

Z=F(0); O=F(1)
def mat(rows): return [[F(x) for x in row] for row in rows]
def eye(n=3): return [[F(i==j) for j in range(n)] for i in range(n)]
def add(a,b): return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def scale(c,a): return [[c*x for x in row] for row in a]
def mul(a,b): return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Z) for j in range(len(b[0]))] for i in range(len(a))]
def tr(a): return sum((a[i][i] for i in range(len(a))),Z)
def inv(a):
 n=len(a); w=[row[:]+eye(n)[i] for i,row in enumerate(a)]
 for j in range(n):
  k=next(k for k in range(j,n) if w[k][j]);w[j],w[k]=w[k],w[j]
  d=w[j][j];w[j]=[x/d for x in w[j]]
  for k in range(n):
   if k!=j:
    d=w[k][j];w[k]=[x-d*y for x,y in zip(w[k],w[j])]
 return [row[n:] for row in w]
def outer(v): return [[x*y for y in v] for x in v]
def invariants(hinv,k):
 b=mul(hinv,k);return tr(b),tr(mul(b,b))
def normalized_constraint(hinv,k,a,b,lam,ambient,rho,u):
 t,n=invariants(hinv,k); ll=lam/(3*lam-1)
 pm=scale(1/a,add(mul(hinv,k),scale(-ll*t,eye())))
 back=scale(a,add(pm,scale(-lam*tr(pm),eye())))
 assert back==mul(hinv,k)
 kin=a*(tr(mul(pm,pm))-lam*tr(pm)**2)
 assert kin==(n-ll*t*t)/a
 scalar_R=ambient+n-t*t
 return kin-b*scalar_R+2*u+rho

def run():
 source={'matter':F(7,5),'clock':F(13,17),'control':F(2,9),'record':F(5,11),'binding':F(-3,13)}
 rho=sum(source.values(),Z); u=F(19,23); ambient=F(11,7);s=F(2,7)
 k0=mat([[F(2,3),F(1,5),F(-1,7)],[F(1,5),F(-3,4),F(2,9)],[F(-1,7),F(2,9),F(5,6)]])
 hs=[eye(),mat([[1,F(1,2),0],[F(1,2),1,0],[0,0,2]])]
 pars=[(F(2),F(1,2),F(1,2)),(F(3,2),F(4,7),F(2,5)),(F(5,3),F(7,11),F(3,5))]
 records=[]
 for hi,h in enumerate(hs):
  hI=inv(h); v=[h[i][0] for i in range(3)];w=[h[i][1] for i in range(3)]
  rank=outer(v);free=add(rank,scale(-1,outer(w)))
  for shape,A in [('rank_one',rank),('traceless',free)]:
   ta,aa=invariants(hI,A)
   assert (ta==1 and aa==1) if shape=='rank_one' else (ta==0 and aa>0)
   for a,b,lam in pars:
    vals=[normalized_constraint(hI,add(k0,scale(z*s,A)),a,b,lam,ambient,rho,u) for z in [-1,0,1]]
    diff=vals[0]+vals[2]-2*vals[1];ll=lam/(3*lam-1)
    formula=2*s*s*((1/a-b)*aa+(b-ll/a)*ta*ta)
    assert diff==formula
    shifted=[normalized_constraint(hI,add(k0,scale(z*s,A)),a,b,lam,ambient+13,rho-19,u+23) for z in [-1,0,1]]
    assert shifted[0]+shifted[2]-2*shifted[1]==diff
    records.append({'metric':hi,'shape':shape,'a':str(a),'b':str(b),'lambda_H':str(lam),'second_difference':str(diff)})
 # Two shape equations are independent in (1/a-b, b-lambda_L/a).
 coeff=mat([[1,1],[2,0]])
 assert coeff[0][0]*coeff[1][1]-coeff[0][1]*coeff[1][0]==-2
 # Hamilton lambda=1/3 has trace(K)=a(1-3 lambda)trace(P)=0.
 assert 1-3*F(1,3)==0
 eps=F(1,1000); amp=F(1,2)
 loose=F(16,5)*eps
 lam_bad=F(51,100); c_bad=(2*lam_bad-1)/(3*lam_bad-1)
 bad_diff=2*amp*amp*c_bad
 assert abs(bad_diff)>4*eps
 d=8*eps; lower=(1+d)/(2+3*d);upper=(1-d)/(2-3*d)
 assert F(2,5)<lower<F(1,2)<upper<F(3,5)
 assert F(1,2)-lower<=loose and upper-F(1,2)<=loose
 # Exact endpoint/max denominator bound; samples supplement algebra, not prove interval theorem.
 assert max(abs(1-3*x) for x in [F(2,5),F(3,5)])==F(4,5)
 for j in range(41):
  lam=F(2,5)+F(j,200)
  assert abs((2*lam-1)/(3*lam-1))>=F(5,2)*abs(lam-F(1,2))
 # Coherent classical-label processing keeps the complete record and auxiliary.
 V=np.zeros((32,8),dtype=np.int64)
 for d0 in range(2):
  for r in range(4): V[(d0*4+r)*4+r//2,d0*4+r]=1
 assert np.array_equal(V.T@V,np.eye(8,dtype=np.int64))
 X=np.array([[0,1],[1,0]],dtype=np.int64);Zp=np.diag([1,-1])
 for Q in [X,Zp]: assert np.array_equal(np.kron(Q,np.eye(16,dtype=np.int64))@V,V@np.kron(Q,np.eye(4,dtype=np.int64)))
 # Inverse on the image restores all coherent input, so arbitrary external reference is retained.
 assert np.array_equal(V.T@V,np.eye(8,dtype=np.int64))
 post_lam=F(2,5);post_gap=2*amp*amp*(2*post_lam-1)/(3*post_lam-1)
 assert post_gap==F(-1,2)
 return {
  'round':1188,'date':'2026-10-14','status':'PASS_CONDITIONAL_SAME_HISTORY_SLICING',
  'groups':{
   '1_full_source_same_event_central_difference':{'pass':True,'exact_tests':len(records),'source_components':{k:str(v) for k,v in source.items()},'total_source':str(rho),'nonzero_ambient_gauss_projection':str(ambient),'nonzero_constant_u':str(u),'records':records,'source_constant_and_ambient_shifts_cancel':True},
   '2_two_shapes_degeneracy_and_finite_error':{'pass':True,'shape_equation_determinant':'-2','exact_selection':'lambda_H=1/2 and a*b=1; a and constant term remain free','lambda_H_one_third':'trace velocity image zero; nonzero-trace slicing family unavailable','per_slice_residual_example':str(eps),'loose_lambda_radius':str(loose),'exact_lambda_lower':str(lower),'exact_lambda_upper':str(upper),'excluded_lambda_example':str(lam_bad),'central_difference_example':str(bad_diff),'allowed_total_residual':str(4*eps),'physical_instrument_tolerance_certified':False},
   '3_postprocessing_does_not_imply_canonical_admission':{'pass':True,'input_dimension':8,'output_dimension':32,'isometry_exact':True,'unknown_data_and_reference_retained_on_full_output':True,'original_record_and_auxiliary_retained':True,'lambda_not_selected_by_label_isometry':str(post_lam),'candidate_geometric_constraint_second_difference':str(post_gap),'full_axiom_countermodel_claimed':False}
  },
  'scope_limits':['CT1188 is an unadopted candidate connecting record slices to the same canonical geometric dynamics.','Smooth effective Lorentz history and a declared two-derivative quadratic metric Hamiltonian are inputs.','Complete physical source must be the same normal projection; nonminimal mixed-momentum couplings need a different calculation.','Finite jet/source dictionary errors must enter the residual budget; no point measurement or actual quantum instrument is certified.','Gauss-Codazzi and ADM reconstruction are reused mature tools; no new generic HDA theorem.','No full joint quantum geometry, nonlinear GR construction or dark-sector prediction delivered.'],
  'new_adopted_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables_completed':0}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 result=run();path=Path(__file__).with_name('results.json')
 if args.write:
  with path.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
  print('WROTE 1188 results: 3 groups; conditional same-history slicing only')
 else:
  assert json.loads(path.read_text(encoding='utf-8-sig'))==result
  print('READ_ONLY_MATCH 1188: 3 groups; full-source central differences and conditional selection')
if __name__=='__main__':main()
