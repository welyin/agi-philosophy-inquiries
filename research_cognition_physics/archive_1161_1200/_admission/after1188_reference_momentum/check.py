"""Exact local chart/clock jets and metric velocity response; not a quantum model."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
Z=F(0)
def mat(rows):return [[F(x) for x in row] for row in rows]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def sc(c,a):return [[c*x for x in r] for r in a]
def mul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Z) for j in range(len(b[0]))] for i in range(len(a))]
def tr(a):return sum((a[i][i] for i in range(len(a))),Z)
def inv(a):
 n=len(a); w=[r[:]+eye(n)[i] for i,r in enumerate(a)]
 for j in range(n):
  k=next(k for k in range(j,n) if w[k][j]);w[k],w[j]=w[j],w[k]
  d=w[j][j];w[j]=[x/d for x in w[j]]
  for k in range(n):
   if k!=j:
    d=w[k][j];w[k]=[x-d*y for x,y in zip(w[k],w[j])]
 return [r[n:] for r in w]
def dot(a,b):return sum((x*y for r,s in zip(a,b) for x,y in zip(r,s)),Z)
def invs(hi,a):
 b=mul(hi,a);return tr(b),tr(mul(b,b))
def ref(hi,k,c,alpha,beta):
 t,n=invs(hi,add(k,sc(-1,c)));return alpha*n+beta*t*t
def p_ref(hi,k,c,alpha,beta):
 b=add(k,sc(-1,c));return add(sc(alpha,mul(mul(hi,b),hi)),sc(beta*invs(hi,b)[0],hi))
def jet(h,k,c0,k0,q):
 # Flat physical metric -dt^2+h_ij dx^i dx^j, t=tau+f(x).
 # At p, df=0 and d2 f=k. The same T=q[t+(c0-k0)_ij x^i x^j/2].
 g=[[Z]*4 for _ in range(4)];g[0][0]=F(-1)
 for i in range(3):
  for j in range(3):g[i+1][j+1]=h[i][j]
 gi=inv(g);dg=[[[Z]*4 for _ in range(4)] for _ in range(4)]
 for d in range(3):
  for i in range(3):dg[d+1][0][i+1]=dg[d+1][i+1][0]=-k[d][i]
 gamma=[[[sum((gi[a][d]*(dg[b][d][c]+dg[c][d][b]-dg[d][b][c]) for d in range(4)),Z)/2 for c in range(4)] for b in range(4)] for a in range(4)]
 pulled=add(c0,add(k,sc(-1,k0)))
 # T_0=q, T_i=0. Normalize u=-dT/q; spatial q derivatives drop at p.
 covh=[[q*pulled[i][j]-q*gamma[0][i+1][j+1] for j in range(3)] for i in range(3)]
 b=sc(-1/q,covh)
 assert b==add(k0,sc(-1,c0))
 assert [[gamma[0][i+1][j+1] for j in range(3)] for i in range(3)]==k
 return pulled,b

def run():
 hs=[eye(3),mat([[1,F(1,2),0],[F(1,2),1,0],[0,0,2]])]
 k0=mat([[F(2,3),F(1,5),0],[F(1,5),F(-3,7),F(2,9)],[0,F(2,9),F(5,11)]])
 c0=mat([[F(1,7),F(-1,11),0],[F(-1,11),F(2,5),F(1,13)],[0,F(1,13),F(-1,3)]])
 alpha,beta,s=F(2,5),F(-1,7),F(1,2)
 cases=[]; derivative_checks=0
 for idx,h in enumerate(hs):
  hi=inv(h)
  v=[h[i][0] for i in range(3)];w=[h[i][1] for i in range(3)]
  rank=[[x*y for y in v] for x in v];tf=add(rank,sc(-1,[[x*y for y in w] for x in w]))
  for name,a in [('rank_one',rank),('traceless',tf)]:
   ks=[add(k0,sc(z*s,a)) for z in [-1,0,1]]
   cs=[]; full=[];bare=[];mix=[];prs=[]
   for k in ks:
    c,b=jet(h,k,c0,k0,F(7,3));cs.append(c)
    f=ref(hi,k,c,alpha,beta);bareval=ref(hi,k,sc(0,c),alpha,beta)
    full.append(f);bare.append(bareval);mix.append(f-bareval)
    pr=p_ref(hi,k,c,alpha,beta);prs.append(pr)
    # Quadratic density: exact symmetric difference equals directional derivative.
    for e in [rank,tf,mat([[0,1,0],[1,0,0],[0,0,0]])]:
     d=F(1,97)
     dl=(ref(hi,add(k,sc(d,e)),c,alpha,beta)-ref(hi,add(k,sc(-d,e)),c,alpha,beta))/(2*d)
     assert dl==2*dot(pr,e);derivative_checks+=1
   delta=lambda vals:vals[0]+vals[2]-2*vals[1]
   assert full[0]==full[1]==full[2] and prs[0]==prs[1]==prs[2]
   ta,aa=invs(hi,a);expected=2*s*s*(alpha*aa+beta*ta*ta)
   assert delta(bare)==expected and delta(mix)==-expected and delta(full)==0
   if idx==0:assert expected==({'rank_one':F(9,70),'traceless':F(2,5)}[name])
   # With the EH comparison input, total momentum changes only by A-tr(A)h^-1.
   def total(k,c):
    return add(add(mul(mul(hi,k),hi),sc(-invs(hi,k)[0],hi)),p_ref(hi,k,c,alpha,beta))
   change=add(total(ks[2],cs[2]),sc(-1,total(ks[1],cs[1])))
   expected_change=sc(s,add(mul(mul(hi,a),hi),sc(-ta,hi)))
   assert change==expected_change
   cases.append({'metric':idx,'shape':name,'bare_second_difference':str(delta(bare)),'omitted_terms_second_difference':str(delta(mix)),'full_reference_second_difference':str(delta(full))})
 # Clock-adapted reading of the comparison action; deliberately not CT across cuts.
 a=1/(1+alpha);ll=(1-beta)/(1+alpha);lh=ll/(3*ll-1)
 assert (a,ll,lh)==(F(5,7),F(40,49),F(40,71))
 return {'status':'PASS_LOCAL_REFERENCE_MOMENTUM_AUDIT','numbered_round':None,'groups':{
  'chart_hessian_transport':{'pass':True,'metrics':2,'shapes':2,'cuts_per_shape':3,'clock_gradient':'7/3','actual_reference_jet_fixed':True,'cases':cases},
  'metric_velocity_dictionary':{'pass':True,'exact_directional_derivatives':derivative_checks,'reference_momentum_invariant':True,'clock_adapted_apparent_a':str(a),'clock_adapted_apparent_lambda_L':str(ll),'clock_adapted_apparent_lambda_H':str(lh),'curvature_coefficient_in_EH_comparison':'1'}},
 'limits':['Local jet and action/metric-velocity identities only; no complete Hamiltonian constraint cancellation.','EH is an explicit comparison input, not generated.','No arbitrary-reference quantum implementation, stability, empirical calibration or full-axiom countermodel.','1188 remains valid in its declared pure-metric/non-mixing class.','No new numbered round, axiom, scientific calibration or GR deliverable.']}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run();p=Path(__file__).with_name('results.json')
 if args.write:
  with p.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
 else:assert json.loads(p.read_text(encoding='utf-8-sig'))==r
 print('PASS: local chart and metric-velocity identities; 2 unnumbered checks')
