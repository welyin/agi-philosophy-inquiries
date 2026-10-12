"""Read-only dynamics qualification; no numbered round or scientific trial.

Only --write-results writes the matching JSON. Continuous domain claims are
proved in the matching Markdown; finite diagnostics do not replace them.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FROZEN = {
 '1129/proof.md':'c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa',
 '1133/proof.md':'476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
 '1134/proof.md':'4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666',
 '1137/proof.md':'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb',
 '1137/results.json':'1fb7f63266bc240d80ab880af707414d4e41a92c2075e14b69459e53b819ceb4',
 '_shared/after1137_local_balance_preflight.md':'95915f871317e8638d23c3ee29cec779a0ff52aad5f76873a0ff896f6a2fcbe0',
}


def pack(values):
 return dict(fractions={k:str(v) for k,v in values.items()},
             decimal={k:float(v) for k,v in values.items()})


def exact_budget():
 lam,kap=F(11,7),F(11,14)
 eps,A,M,sig,T=F(1,10**5),F(10**34),F(10**80),F(1,10**30),F(4)
 eta_lap=F(36080,7)
 chi_raw=F(8*104*88,21); chi_lap=F(3500)
 assert chi_raw<chi_lap
 q0_lap=F(4*(42800+2*1801*20+21*1800))
 n_lap=F(2*(1800+20**2)); y_lap=2*n_lap
 hold=chi_lap*42+2*68*8884+q0_lap
 record=lam*(16*chi_lap+2*16*68*40+16*n_lap)
 cont=kap*(16*chi_lap+2*16*68*80+16*y_lap)
 c2rec=hold+record+cont
 c2eta=16*lam*64*eps**-3*eta_lap
 c2safe=F(9*10**21)
 assert c2rec==F(18148168,7) and c2eta+c2rec<c2safe
 cy=F(24780800000000000000000,49)+16316
 v=16*lam*33536000000000000000+F(558,7)
 mmin=3*M-v
 assert mmin>2*M
 rows={}
 for name,b1,b2 in [('original_code',268,28000),('far_menu',600,400000)]:
  p1=F(b1)+94*A;p2=F(b2)+188*b1*A+12000*A*A
  a2=sig**-2+2*p1/sig+p2
  p1t=p1+cy*T
  p2t=p2+2*cy*p1*T+cy**2*T**2+c2safe*T
  dg=T*a2/(2*mmin)
  dr=(T*p2+cy*p1*T**2+cy**2*T**3/3+c2safe*T**2/2)/(2*M)
  assert p2t<F(1201*10**69)
  assert dg+dr<F(3201,10**11)
  rows[name]=pack(dict(p1=p1,p2=p2,A2=a2,Pf_norm_at_T=p1t,Pf2_norm_at_T=p2t,
                       external_residual_integral=dg,recoil_residual_integral=dr,
                       complete_HL_difference=dg+dr,safe_HL_difference=F(3201,10**11)))
 return dict(constants=pack(dict(eta_Laplacian_L1=eta_lap,chi_raw_Laplacian_L1=chi_raw,
      chi_Laplacian_L1=chi_lap,Q0_Laplacian=q0_lap,N_Laplacian=n_lap,Yg_Laplacian=y_lap,
      receiver_hold=hold,receiver_record=record,receiver_continuation=cont,
      C2_receiver=c2rec,C2_handoff=c2eta,C2_sum=c2eta+c2rec,C2_safe=c2safe,
      C_Y=cy,V_bound=v,m_min=mmin,T=T,A=A,M=M,epsilon=eps,sigma_Y=sig)),rows=rows)


def radial_diagnostic():
 rows=[]
 for nodes in [256,512]:
  x,w=np.polynomial.legendre.leggauss(nodes);u=(x+1)/2;w=w/2
  bump=np.exp(-1/(u*(1-u)));ib=float(np.dot(w,bump))
  rp=bump/ib
  rpp=rp*(1-2*u)/(u*u*(1-u)**2)
  r=4+u
  eta_lap=float(4*math.pi*np.dot(w,r*r*np.abs(-rpp-2*rp/r)))
  z=1/(1-u*u);q=np.exp(-z)
  qp=-2*u*z*z*q
  qpp=(4*u*u*z**4-2*z*z-8*u*u*z**3)*q
  norm=float(4*math.pi*np.dot(w,u*u*q))
  chi_grad=float(4*math.pi*np.dot(w,u*u*np.abs(qp))/norm)
  chi_lap=float(4*math.pi*np.dot(w,u*u*np.abs(qpp+2*qp/u))/norm)
  total_rpp=float(np.dot(w,np.abs(rpp)))
  assert total_rpp<16 and eta_lap<36080/7
  assert norm>1/8 and chi_grad<68 and chi_lap<3500
  rows.append(dict(nodes=nodes,R_prime_max=float(rp.max()),R_second_L1=total_rpp,
       eta_Laplacian_L1=eta_lap,chi_normalization=norm,
       chi_gradient_L1=chi_grad,chi_Laplacian_L1=chi_lap))
 # Absolute values have interior cusps: check convergence at its natural rate,
 # while analytic inequalities, not quadrature, certify the constants.
 assert abs(rows[0]['eta_Laplacian_L1']-rows[1]['eta_Laplacian_L1'])<.2
 return rows


def exp_h(a,t):
 w,v=np.linalg.eigh(a)
 return (v*np.exp(-1j*t*w))@v.conj().T


def sqrt_h(a):
 w,v=np.linalg.eigh(a)
 assert w.min()>0
 return (v*np.sqrt(w))@v.conj().T


def correlated_dynamics_diagnostic():
 # A finite 9-momentum x 2-source x 3-total-momentum matrix illustrates the
 # comparison on a correlated input. It has no claim to be a normal QFT model.
 p=np.arange(-4.,5.);ks=np.array([-1.3,.2,2.1]);M=80.
 pf=np.kron(np.diag(p),np.eye(2))
 source=np.array([[.5,.2],[.2,.8]],complex)
 sx=np.array([[0,1],[1,0]],complex)
 kernel=np.exp(-(p[:,None]-p[None,:])**2/2)/8
 v=.3*np.kron(kernel,sx)
 free=3*M*np.eye(18)+np.kron(np.diag(np.sqrt(1+p*p)),np.eye(2))+np.kron(np.eye(9),source)
 l=free+v
 recoil=np.kron(np.diag(np.sqrt(M*M+p*p)-M),np.eye(2))
 mint=l+recoil
 hblocks=[sqrt_h(mint@mint+k*k*np.eye(18)) for k in ks]
 CY=float(np.linalg.norm(pf@v-v@pf,2))
 second=pf@(pf@v-v@pf)-(pf@v-v@pf)@pf
 C2=float(np.linalg.norm(second,2))
 lower=3*M-float(np.linalg.norm(v,2))
 assert np.linalg.norm(l@recoil-recoil@l)>1e-5
 logical=np.array([[1,1j],[1j,1]],complex)/math.sqrt(2)
 wave=np.exp(-p*p/5);wave/=np.linalg.norm(wave)
 blocks=[]
 for k in ks:
  block=np.einsum('p,sq->psq',wave,logical)
  phase=np.exp(1j*.7*k*p[:,None]*np.array([1.,-1.])[None,:])
  blocks.append((phase[:,:,None]*block).reshape(18,2)/math.sqrt(3))
 j=np.concatenate(blocks,axis=0)
 assert np.linalg.norm(j.conj().T@j-np.eye(2))<1e-12
 pfull=np.kron(np.eye(3),pf)
 p1=float(np.linalg.norm(pfull@j,2));p2=float(np.linalg.norm(pfull@pfull@j,2))
 k2=np.kron(np.diag(ks*ks),np.eye(18))
 a2=float(np.linalg.norm(k2@j,2))
 amp=j[:,0].reshape(3,18);reduced=amp@amp.conj().T
 purity=float(np.trace(reduced@reduced).real)
 assert purity<.95  # no replacement by an independent K packet
 rows=[]
 for t in [.1,.8,2.,4.]:
  evolved_l=np.concatenate([exp_h(l,t)@q for q in blocks],axis=0)
  evolved_h=np.concatenate([exp_h(h,t)@q for h,q in zip(hblocks,blocks)],axis=0)
  difference=float(np.linalg.norm(evolved_h-evolved_l,2))
  predicted=t*a2/(2*lower)+(t*p2+CY*p1*t*t+CY*CY*t**3/3+C2*t*t/2)/(2*M)
  p2actual=float(np.linalg.norm(pfull@pfull@evolved_l,2))
  p2bound=p2+2*CY*p1*t+CY*CY*t*t+C2*t
  assert difference<=predicted+1e-10 and p2actual<=p2bound+1e-10
  gram=float(np.linalg.norm(evolved_h.conj().T@evolved_h-np.eye(2)))
  assert gram<1e-12
  rows.append(dict(time=t,actual_complete_map_difference=difference,
       Duhamel_upper=predicted,actual_Pf2_map_norm=p2actual,Pf2_upper=p2bound,
       complete_output_gram_residual=gram))
 return dict(dimension=54,logical_input_dimension=2,source_dimension=2,
       total_momentum_reduced_purity=purity,C_Y=CY,C2=C2,rows=rows,
       scope='Finite correlated mechanism only; continuous domain proof is in Markdown.')


def calculate():
 sha={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in FROZEN}
 assert sha==FROZEN
 return dict(status='PASS',scope='Dynamics qualification only; no numbered round.',
       new_axioms=0,science_experiments=0,frozen_sources=sha,
       exact=exact_budget(),radial_quadrature=radial_diagnostic(),
       finite_correlated_diagnostic=correlated_dynamics_diagnostic(),
       selector_b=0,maximum_time=4,all_normal_sources_and_outputs_retained=True,
       actual_H_Pf2_domain_claimed=False,L_Pf2_domain_proved=True,
       local_energy_flux_or_readout_domain_certified=False)


def compare(a,b,path=()):
 if isinstance(a,dict):
  assert a.keys()==b.keys()
  for k in a:compare(a[k],b[k],path+(k,))
 elif isinstance(a,list):
  assert len(a)==len(b)
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+(str(i),))
 elif isinstance(a,float):
  absolute=1e-10 if 'finite_correlated_diagnostic' in path else 0.
  assert math.isclose(a,b,rel_tol=1e-8,abs_tol=absolute),(path,a,b)
 else:assert a==b,(path,a,b)


if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
 args=parser.parse_args();data=calculate();out=HERE/'local1138_dynamics_preflight.json'
 if args.write_results:out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 else:compare(data,json.loads(out.read_text(encoding='utf-8')))
 print('Dynamics qualification PASS: 6 frozen sources; exact C2 and H/L budgets; radial derivatives; correlated finite diagnostic. No round assigned; local E/J/readout domains remain separate.')
