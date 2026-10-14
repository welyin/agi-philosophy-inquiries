"""Finite original-source reference/boost witness. Default read-only.
Exact rational interval arithmetic certifies a frozen rational mode; NumPy
is used only by make_witness.py to discover that mode, never for this proof.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import isqrt, factorial
import json,argparse
HERE=Path(__file__).resolve().parent
D=10**30
def pt(x):return (F(x),F(x))
def out(a,b):return (F((a*D).__floor__(),D),F((b*D).__ceil__(),D))
def add(a,b):return out(a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 t=[x*y for x in a for y in b];return out(min(t),max(t))
def sc(a,x):return mul(a,pt(x))
def div(a,b):
 assert b[0]>0
 return mul(a,(1/b[1],1/b[0]))
def sq(a):
 if a[0]<=0<=a[1]:return out(F(0),max(a[0]*a[0],a[1]*a[1]))
 return out(min(a[0]*a[0],a[1]*a[1]),max(a[0]*a[0],a[1]*a[1]))
def sqrt(a):
 assert a[0]>=0
 lo=isqrt(a[0].numerator*D*D//a[0].denominator)
 hi=isqrt(a[1].numerator*D*D//a[1].denominator)+1
 return F(lo,D),F(hi,D)
def exp(a):
 if a==pt(0):return pt(1)
 if a[1]<=0:return div(pt(1),exp(neg(a)))
 assert a[0]>=0 and a[1]<25
 s=pt(1);t=pt(1);n=140
 for j in range(1,n+1):t=sc(mul(t,a),F(1,j));s=add(s,t)
 rem=a[1]**(n+1)/factorial(n+1)/(1-a[1]/(n+2))
 return out(s[0],s[1]+rem)
def atan_inv(d,n):
 x=F(1,d);s=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
 r=x**(2*n+1)/F(2*n+1)
 return out(min(s,s+(-1)**n*r),max(s,s+(-1)**n*r))
PI=sub(sc(atan_inv(5,60),16),sc(atan_inv(239,20),4))
def trig(x,cos=False):
 n=45;x2=sq(x);s=pt(0);t=pt(1) if cos else x
 for k in range(n):
  s=add(s,t)
  d=(2*k+1)*(2*k+2) if cos else (2*k+2)*(2*k+3)
  t=sc(mul(t,x2),-F(1,d))
 power=2*n if cos else 2*n+1
 rem=max(abs(x[0]),abs(x[1]))**power/F(factorial(power))
 return out(s[0]-rem,s[1]+rem)
def pair(a):return {'lo':str(a[0]),'hi':str(a[1]),'float_lo':float(a[0]),'float_hi':float(a[1])}
def upperabs(a):return max(abs(a[0]),abs(a[1]))
def cadd(a,b):return add(a[0],b[0]),add(a[1],b[1])
def cmul(a,b):return sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0]))
def cnorm(a):return add(sq(a[0]),sq(a[1]))
def run():
 ell=F(16,65);eps=F(3,20);nu=F(7,20);om=F(1,2);m=F(2)
 z0=F(989835554078098,10**15)
 sins=[trig(sc(PI,F(2*i,65))) for i in range(33)]
 ks=[sc(sq(s),1/ell**2) for s in sins]
 n0=pt(1+eps*z0);v0=pt(1+nu*z0)
 total=pt(0);weights=0
 for i,j,k in combinations_with_replacement(range(33),3):
  perm=1 if i==k else 3 if i==j or j==k else 6
  weight=perm*(1 if i==0 else 2)*(1 if j==0 else 2)*(1 if k==0 else 2)
  weights+=weight
  q=add(add(ks[i],ks[j]),ks[k])
  num=add(sc(n0,eps*m*m),sc(mul(v0,q),nu))
  e=sqrt(add(sc(sq(n0),m*m),mul(sq(v0),q)))
  total=add(total,sc(div(num,e),weight))
 assert weights==65**3
 b=sub(sc(total,-F(2,65**3)),pt(eps*om))
 f=div(neg(b),sqrt(add(pt(om*om),sq(b))))
 # Global contraction: |b'|<=4m(eps-nu)^2/(3sqrt(3)Nmin vmin^2);
 # outer derivative <=1/om, sqrt(3)>=3/2.
 lip=4*m*(eps-nu)**2/(3*F(3,2)*(1-eps)*(1-nu)**2*om)
 assert lip<F(2,5)
 root_error=upperabs(sub(f,pt(z0)))/(1-lip)
 assert root_error<F(1,10**12)
 zi=(z0-F(1,10**12),z0+F(1,10**12))
 mu=div(sc(add(pt(1),sc(zi,eps)),m*ell),add(pt(1),sc(zi,nu)))
 # On each finite M block, the full finite-K boost error is <=(M-.5)/K.
 # This interval widening is an operator-error proxy including v_K/v,
 # not an assertion that the physical B_K equals B(muK).
 K=10**14
 muK=(mu[0]-F(1,K),mu[1]+F(1,K))
 seed=json.loads((HERE/'mode_witness.json').read_text(encoding='utf-8'))
 M=seed['sites'];u=[F(x) for x in seed['vector']];q=sum((x*x for x in u),F(0))
 assert M==32 and len(u)==2*M and q>0
 mu0=F(seed['mu']);lam=F(seed['eigenvalue'])
 bu=[mu0*F(i//2*2+1,2)*(1 if i%2==0 else -1)*u[i] for i in range(2*M)]
 for site in range(M-1):
  for spin in (0,1):
   a=2*site+spin;bidx=2*(site+1)+(1-spin);t=F(site+1,2)
   bu[a]+=t*u[bidx];bu[bidx]+=t*u[a]
 res=[bu[i]-lam*u[i] for i in range(2*M)]
 mean=sum((u[i]*res[i] for i in range(2*M)),F(0))/q
 r2=sum((x*x for x in res),F(0))/q
 edge=F(M*M,4)*(u[-1]**2+u[-2]**2)/q
 dmu=max(abs(muK[0]-mu0),abs(muK[1]-mu0))
 shift=dmu*F(2*M-1,2)
 r_up=sqrt(pt(r2))[1]+shift
 # |F'|<=2, sup |F''|/2<=2 for F(x)=fermi(2*pi*x).
 # Infinite Jacobi B is essentially self-adjoint (Carleman).
 endpoint_bound=2*(abs(mean)+shift)+2*(r_up*r_up+edge)
 occ=pt(0)
 for n in range(65):
  angle=sc(PI,F(2*n,65))
  if n>32:angle=sub(angle,sc(PI,2))
  sn=trig(angle);co=trig(angle,True)
  # transformed real B gauge: w_i=i^i*u_i; Fourier phase exp(-ik i)*i^i
  root=(sn,co);power=(pt(1),pt(0));a=(pt(0),pt(0));bb=(pt(0),pt(0))
  for j in range(M):
   a=cadd(a,(sc(power[0],u[2*j]),sc(power[1],u[2*j])))
   bb=cadd(bb,(sc(power[0],u[2*j+1]),sc(power[1],u[2*j+1])))
   power=cmul(power,root)
  aa=cnorm(a);cc=cnorm(bb)
  cross=add(mul(a[0],bb[0]),mul(a[1],bb[1]))
  numerator=add(mul(mu,sub(aa,cc)),sc(mul(sn,cross),2))
  e=sqrt(add(sq(mu),sq(sn)))
  occ=add(occ,sc(sub(add(aa,cc),div(numerator,e)),F(1,2*65)*1/q))
 thermal=div(pt(1),add(pt(1),exp(mul(sc(PI,2),pt(lam)))))
 # Torus versus auxiliary infinite-line covariance, same mu, no change in actual object.
 eta=F(3,10);ep=exp(pt(eta));em=div(pt(1),ep)
 sinh=sc(sub(ep,em),F(1,2));cosh=sc(add(ep,em),F(1,2))
 strip=add(pt(F(1,2)),sc(div(add(mu,cosh),sqrt(sub(sq(mu),sq(sinh)))),F(1,2)))
 ap=pt(0);am=pt(0);pp=pt(1);pm=pt(1)
 for j in range(M):
  a=sqrt(pt((u[2*j]**2+u[2*j+1]**2)/q))
  ap=add(ap,mul(a,pp));am=add(am,mul(a,pm))
  pp=mul(pp,ep);pm=mul(pm,em)
 tail=exp(pt(-eta*65))
 alias=div(sc(mul(mul(mul(strip,tail),ap),am),2),sub(pt(1),tail))
 gap=occ[0]-thermal[1]-endpoint_bound-alias[1]
 assert gap>F(81,10000)
 return {'schema':'original_parent_reference_certificate_v1','status':'PASS',
 'groups':3,'pi_interval':pair(PI),'fixed_point':{'grid_weights':weights,'global_lipschitz_upper':str(lip),
 'map_at_seed':pair(f),'root_error_upper':str(root_error),'root_enclosure':pair(zi),'mu':pair(mu)},
 'witness':{'sites':M,'finite_K':str(K),'norm_squared':str(q),'rayleigh_residual_squared':str(r2),
 'far_end_residual_squared':str(edge),'finite_and_extended_boost_error_upper':float(endpoint_bound),
 'source_occupation':pair(occ),'boost_occupation':pair(thermal),
 'torus_alias_upper':float(alias[1]),'certified_difference_lower':str(gap),
 'certified_difference_float_lower':float(gap),'threshold':'.0081'},
 'scope':'Specified source-CAR preparation and instantaneous full-parent energy compression only. Not a full-parent vacuum, stable state, physical boost instrument, finite-time heat flux, area entropy, or GR theorem.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
 result=run()
 if a.write:
  with (HERE/'certificate_results.json').open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
 print(json.dumps(result,ensure_ascii=False,indent=2))
