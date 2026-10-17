"""1254: exact rational certificate plus finite coherent-state diagnostics.
Default is read-only replay. --write-once saves a new result exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, math
import numpy as np
P=Path(__file__).resolve().parent

def rat(x):
 return {'numerator':x.numerator,'denominator':x.denominator}

def compute():
 xmin=F(99,10); bmax=F(21,10); qmax=1+bmax/2*(1/xmin-F(1,10))
 amax=bmax/(2*xmin*xmin)*F(501,500)
 alow=F(11,10)/(2*xmin*xmin)*F(501,500)
 assert qmax==1+F(7,6600)<F(501,500)
 assert amax==F(10521,980100)<F(11,1000)
 assert 10-amax/2>xmin and amax<F(11,1000)
 low_drop=alow/2
 assert low_drop==F(167,59400)
 E=lambda B:B-B*B/20
 assert 2*E(bmax)/xmin<F(2,5)
 clock2=F(600121,360000)
 assert clock2<F(169,100)
 gap=F(1,200)*F(10,13)**2-low_drop
 assert gap==F(1477,10038600)>F(1,7000)
 agap=2*xmin*gap
 assert agap==F(1477,507000)
 elo=E(F(11,10)); ehi=E(F(2)); egap=ehi-elo
 assert elo==F(2079,2000) and ehi==F(9,5) and egap==F(1521,2000)
 kgap=48*(ehi*ehi-elo*elo)/20**6
 assert kgap==F(25913277,16000000000000)
 weights=[F(1,8),F(3,8),F(3,8),F(1,8)]
 six=[F(0) for _ in range(6)]
 for e,w in enumerate(weights):six[e]+=w/2;six[e+2]+=w/2
 assert six==[F(x,16) for x in [1,3,4,4,3,1]]
 assert sum(weights)==sum(six)==1
 leak_bounds=[]
 for b,c2 in [(1,F(19,81)),(2,F(9,16))]:
  low=c2/(2*(1+c2))
  assert low==({1:F(19,200),2:F(9,50)}[b])
  leak_bounds.append(dict(b=b,derivative_square_lower=rat(c2),uniform_probability_lower_hbar1=rat(low)))
 # Floating point diagnosis only: no claim of certified ODE or quantum evolution.
 samples=[]; plus=np.array([1.,1.])/math.sqrt(2); minus=np.array([1.,-1.])/math.sqrt(2)
 for b in [1,2]:
  for N in [320,1280,5120,81920]:
   eps=1/N**2; probs=[]; residuals=[]
   def ps(lam):
    w=lam/10; u=math.sqrt(2*w-w*w)
    return 10*(math.atanh(u)-u)
   for e,w in enumerate(weights):
    mplus=18+2*math.sqrt(19)+e;mminus=mplus+2
    dp=ps(b+mminus/N)-ps(b+mplus/N)
    r=math.exp(-dp*dp/(4*eps))
    fplus=np.array([1.,0.]);fminus=np.array([r,math.sqrt(max(0.,1-r*r))])
    state=(np.kron(plus,fplus)+np.kron(minus,fminus))/math.sqrt(2)
    pR=float(np.sum(np.abs(state[2:])**2))
    residuals.append(abs(pR-(1-r)/2));probs.append(float(w)*pR)
    assert abs(np.vdot(state,state)-1)<2e-14
   prob=sum(probs)
   assert prob>float({1:F(19,200),2:F(9,50)}[b])
   samples.append(dict(b=b,N=N,p_R_nonvac_initial=round(prob,12),swap_formula_residual_below_1e_14=max(residuals)<1e-14))
 return dict(round=1254,status='PASS',scientific_groups=3,
   exact_classical=dict(window_proper=[0,1],window_boundary=[0,1],minimum_N=320,qmax=rat(qmax),amax=rat(amax),low_drop=rat(low_drop),clock2=rat(clock2),radius_gap=rat(gap),area_gap=rat(agap)),
   exact_independent_geometry=dict(energy_gap=rat(egap),compactness_gap=rat(egap/10),scaled_K_gap=rat(kgap)),
   exact_preparation=dict(FK_weights=[rat(x) for x in weights],six_source_weights=[rat(x) for x in six],uniform_R_leak_bounds=leak_bounds),
   floating_diagnostics=samples,
   limitations=['Classical branch histories, not a quantum propagator.','Coherent preparation test, not later R transfer.','No actual collective or local geometric readout contract.','EH model input; CT not adopted; no new axiom.'])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write-once',action='store_true');args=ap.parse_args();r=compute();out=P/'results.json'
 if args.write_once:
  with out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
  print('WROTE_ONCE_1254')
 else:
  assert json.loads(out.read_text('utf-8'))==r
  print('READ_ONLY_MATCH_1254')
 print('PASS: 3 scoped groups; exact rational bounds and separate floating preparation diagnostics')
if __name__=='__main__':main()
