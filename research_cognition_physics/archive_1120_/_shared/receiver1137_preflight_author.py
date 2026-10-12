"""Read-only receiver qualification; no numbered round or scientific trial.
--write-results writes only the matching preflight JSON.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
# The old source digest is read and reported, not inferred from a model name.
FROZEN={
 '1133/proof.md':'476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
 '1134/proof.md':'4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666',
 '1135/proof.md':'f3577777e8983a02d02bf66ef5f869b207ac4f94df1610202c1ce6e995e17068',
 '1136/proof.md':'360c450f85f40798c6fb9524b00124364b604434d814c770503e0e97e2196462'}

def pack(d):
 return dict(fractions={k:str(v) for k,v in d.items()},decimal={k:float(v) for k,v in d.items()})

def calculate():
 sha={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
 assert sha==FROZEN
 # Analytic pointwise q''' bound using e > 8/3, followed by tensor Sobolev.
 third=sum(F(c)*n**n*F(3,8)**n for c,n in [(24,3),(60,4),(48,5),(8,6)])
 p3square=F(3*32*3000**2+18*320000*128+6*128**3)
 assert third<3000 and p3square<41000**2
 radius=F(13,4)
 j2g=radius**2*1800+2*radius*20+F(1,4)
 j2hg=radius**2*(41000+1800)+2*radius*1801+F(21,4)
 assert j2g<19200 and j2hg<464000
 q1=4*(5872+21*66); q2=4*(464000+2*5872*66+21*19200)
 n2=2*(19200+66**2); y2=4*(19200+66**2)
 lam=F(11,7); kap=F(11,14)
 c1=q1+16*lam*132+16*kap*264
 c2=q2+16*lam*48000+16*kap*96000
 assert c1<36000 and c2<10**7
 M=I=mo=mc=F(10**80); t=F(4); a=F(1,10**9); sig=F(1,10**30)
 l1=40/a; l2=2000/a**2; C1=F(36000); C2=F(10**7)
 drot=(t*l2+2*C1*t*t*l1+F(3,2)*C2*t*t+C1*C1*t**3)/(2*I)
 dorient=t*C1*a
 assert drot<F(41,10**60)
 A=F(10**34); eps=F(1,10**5); C=F(33536000000000000000)
 v=16*lam*C+F(558,7); mmin=3*M-v
 cy=F(24780800000000000000000,49)+16316
 dy=cy*(2*t*sig+t*t/(2*M*sig)); dd=28*t/(2*M*M*sig*sig)
 ds=3*t/(eps*A)*(lam*C*512+42*8/F(1,10**6)+(lam+kap)*512)
 dm=F(11,8750); dc=F(42,10**6)
 rows={}
 for key,b1,b2 in [('code',268,28000),('far_local_menu',600,400000)]:
  p1=F(b1)+94*A; p2=F(b2)+188*b1*A+12000*A*A
  a1=1/sig+p1; a2=1/sig**2+2*p1/sig+p2; bmix=p1/sig+p2
  bt=t*a2/(2*mmin)+t*(a2+2*bmix)/(2*M)+cy*a1*t*t/(2*M)
  task=dm+dc+ds+dy+dd+bt+drot+dorient
  causal=2*(dy+dd+bt)
  assert bt<F(9,10**8) and task<F(145,10**5) and causal<F(2,10**7)
  rows[key]=pack(dict(p1=p1,p2=p2,A1=a1,A2=a2,B2=bmix,
          delta_BT=bt,complete_task_error=task,finite_causal_error=causal))
 angular=pack(dict(q_third_derivative_upper=third,p3_norm_square_bound=p3square,
    J2_g=j2g,J2_hg=j2hg,Q0_first=q1,Q0_second=q2,N_second=n2,Y_second=y2,
    receiver_first=c1,receiver_second=c2,C1=C1,C2=C2,
    rotor_L_norm=l1,rotor_L2_norm=l2,rotor_kinetic_norm=l2/(2*I),
    rotor_mean_kinetic_upper=l1*l1/(2*I),delta_rotor=drot,delta_orientation=dorient,
    canonical_spin_bound=F(5*10**10),actual_rotor_impulse_norm=4*C1*t))
 common=pack(dict(I=I,mO=mo,mC=mc,M=M,m_min=mmin,V_all=v,C_Y=cy,
    delta_source=ds,delta_Y=dy,delta_DCP=dd,delta_capture=dc,delta_material=dm,
    rest_mass_without_old_D=3*M+1))
 # Independent one-dimensional quadrature of Haar radial normalisation and
 # angular moments. It does not replace the stated analytic upper bounds.
 numeric=[]
 af=float(a)
 for n in [128,256]:
  x,w=np.polynomial.legendre.leggauss(n); u=(x+1)/2; w=w/2
  z=1/(1-u*u); q=np.exp(-z)
  q1n=-2*u*z*z*q
  q2n=(4*u*u*z**4-2*z*z-8*u*u*z**3)*q
  haar=(2/math.pi)*np.sin(af*u/2)**2*af
  zn=float(np.dot(w,haar*q*q)); l1n=float(np.dot(w,haar*q1n*q1n)/zn)
  radial=q2n+af/np.tan(af*u/2)*q1n
  l2n=float(np.dot(w,haar*radial*radial)/zn)
  assert zn>af**3/6144 and math.sqrt(l1n)<40 and math.sqrt(l2n)<2000
  numeric.append(dict(nodes=n,normalization_integral=zn,
      L_norm_times_radius=math.sqrt(l1n),L2_norm_times_radius_squared=math.sqrt(l2n)))
 assert abs(numeric[0]['L2_norm_times_radius_squared']-numeric[1]['L2_norm_times_radius_squared'])<1e-5
 return dict(status='PASS',scope='Author qualification gate 1 only; no round number assigned.',
    new_axioms=0,science_experiments=0,frozen_sources=sha,angular=angular,
    common=common,task_rows=rows,Haar_quadrature=numeric,
    source_Z_selector_branch=0,old_D_dimension=288,
    old_clock_Mr_at_T4_claimed=False,wrong_orientation_and_boost_independent_gate_pending=True)

def compare(x,y):
 if isinstance(x,dict):
  assert x.keys()==y.keys()
  for k in x: compare(x[k],y[k])
 elif isinstance(x,list):
  assert len(x)==len(y)
  for a,b in zip(x,y): compare(a,b)
 elif isinstance(x,float): assert math.isclose(x,y,rel_tol=1e-7,abs_tol=1e-12),(x,y)
 else: assert x==y,(x,y)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
 out=HERE/'receiver1137_preflight_author.json'; data=calculate()
 if args.write_results:out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:compare(data,json.loads(out.read_text(encoding='utf-8')))
 print('receiver preflight author PASS: frozen identities, rational angular/task budgets, normal Haar packet. No numbered round.')
