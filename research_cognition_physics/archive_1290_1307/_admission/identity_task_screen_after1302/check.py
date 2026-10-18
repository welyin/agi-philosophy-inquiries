"""Count-zero identity screen; exact polynomials, no new anomaly classification."""
import argparse,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def add(a,b):
 return [ (a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b))) ]
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def data():
 weights=[6,3,3,2,1,1]; f=[1,-1,-1,-3,3,3]
 q=[[1,1],[-4,-1],[2,-1],[-3,-3],[6,3],[0,3]]
 mixed2=[0,0,0];mixed1=[0,0]
 for w,fi,p in zip(weights,f,q):
  mixed2=add(mixed2,[w*fi*c for c in mul(p,p)])
  mixed1=add(mixed1,[w*fi*fi*c for c in p])
 values={'SU3_squared_F':2*f[0]+f[1]+f[2],
         'SU2_squared_F':3*f[0]+f[3],
         'F_cubic':sum(w*fi**3 for w,fi in zip(weights,f)),
         'F_gravity':sum(w*fi for w,fi in zip(weights,f))}
 assert all(v==0 for v in values.values())
 assert mixed2==[0,0,0] and mixed1==[0,0]
 sources=['README.md','../../_shared/global_form_corollary_after1299.md',
          '../../1302/proof.md','../../1302/acceptance.json',
          '../../../archive_1046_1062/_admission/physical_subsystems_after1055/selection.md']
 return dict(status='PASS_COUNT_ZERO_IDENTITY_SCREEN',new_science_rounds=0,new_axioms=0,
             native_identity='3(B-L)',weights=weights,F=f,q_m_coefficients=q,
             local_coefficients=values,q_squared_F_polynomial=mixed2,q_F_squared_polynomial=mixed1,
             representative_q={str(m):[a+b*m for a,b in q] for m in (0,-1)},
             no_selection_of_standard_branch_by_these_coefficients=True,
             no_full_quantum_or_O_certificate=True,
             source_bindings={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in sources},
             check_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
def main():
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');args=p.parse_args()
 r=data();out=HERE/'results.json'
 if args.write:
  with out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
  print('SAVED_COUNT_ZERO_IDENTITY_SCREEN')
 else:
  assert json.loads(out.read_text(encoding='utf-8-sig'))==r
  print('READ_ONLY_MATCH PASS_COUNT_ZERO_IDENTITY_SCREEN')
if __name__=='__main__':main()
