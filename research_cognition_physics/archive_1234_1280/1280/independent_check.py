"""1280 independent verification; no author code import.
Analytic arguments carry continuous claims. Fraction arithmetic certifies budgets;
finite tensor jets and RK/Gauss calculations are separately labelled diagnostics.
Default read-only match; --write creates results exclusively.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sumf(xs):return sum(xs,F())
def matzero():return [[F() for _ in range(4)] for _ in range(4)]
def mm(a,b):return [[sumf(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def fractions():
    amax=F(1,1000); r0=1-2*amax
    dc=F(3,4)*amax/r0**2
    de=1-r0**2
    total=dc+de*(1+dc);target=F(1,200)
    assert dc==F(375,498002) and de==F(999,250000)
    assert total==F(591628623,124500500000)<target
    assert target-total==F(30873877,124500500000)
    assert 1/r0**2==F(250000,249001)
    # Coherent stresses need not agree, despite equal KG operators on R=0.
    # Jet f=f'=1,f''=0, R_uu=2*(f')^2/36000 gives this improvement.
    improvement=(F(1,18000)-2)/6
    assert improvement!=0
    return {'r_pulse_lower':str(r0),'delta_C':str(dc),'delta_eta':str(de),
      'relative_bound':str(total),'relative_bound_decimal':float(total),
      'strict_target':str(target),'strict_margin':str(target-total),
      'clock_ratio_upper':str(1/r0**2),'clock_ratio_excess':str(1/r0**2-1),
      'B_positive_lower':str(4*r0**2),
      'coherent_improvement_nonzero_jet_example':str(improvement)}

def tensor_jet(r,p,pp):
    # Derive connection and curvature from metric and first/second metric jets.
    g=matzero();g[0][1]=g[1][0]=F(-1);g[2][2]=g[3][3]=r*r
    gi=matzero();gi[0][1]=gi[1][0]=F(-1);gi[2][2]=gi[3][3]=1/(r*r)
    assert mm(g,gi)==[[F(i==j) for j in range(4)] for i in range(4)]
    dg=[matzero() for _ in range(4)]
    dg[0][2][2]=dg[0][3][3]=2*r*p
    ddg=[[matzero() for _ in range(4)] for __ in range(4)]
    ddg[0][0][2][2]=ddg[0][0][3][3]=2*(p*p+r*pp)
    dgi=[[[ -x for x in row] for row in mm(mm(gi,dg[l]),gi)] for l in range(4)]
    C=[[[F() for _ in range(4)] for __ in range(4)] for ___ in range(4)]
    dC=[[[[F() for _ in range(4)] for __ in range(4)] for ___ in range(4)] for ____ in range(4)]
    for a in range(4):
      for b in range(4):
       for c in range(4):
        C[a][b][c]=sumf(gi[a][e]*(dg[b][e][c]+dg[c][e][b]-dg[e][b][c]) for e in range(4))/2
        for l in range(4):
         dC[l][a][b][c]=sumf(dgi[l][a][e]*(dg[b][e][c]+dg[c][e][b]-dg[e][b][c])+gi[a][e]*(ddg[l][b][e][c]+ddg[l][c][e][b]-ddg[l][e][b][c]) for e in range(4))/2
    R={}
    for a in range(4):
     for b in range(4):
      for c in range(4):
       for d in range(4):
        R[a,b,c,d]=dC[c][a][d][b]-dC[d][a][c][b]+sumf(C[a][c][e]*C[e][d][b]-C[a][d][e]*C[e][c][b] for e in range(4))
    Ric=[[sumf(R[a,b,a,d] for a in range(4)) for d in range(4)] for b in range(4)]
    expected=matzero();expected[0][0]=-2*pp/r
    assert Ric==expected
    scalar=sumf(gi[a][b]*Ric[a][b] for a in range(4) for b in range(4))
    assert scalar==0 and mm(mm(Ric,gi),Ric)==matzero()
    assert all(C[0][a][b]==0 for a in range(4) for b in range(4))
    # Weyl tensor, using R=0; this checks all 256 components at the exact jet.
    for a in range(4):
     for b in range(4):
      for c in range(4):
       for d in range(4):
        lower=sumf(g[a][e]*R[e,b,c,d] for e in range(4))
        W=lower-(g[a][c]*Ric[d][b]-g[a][d]*Ric[c][b]-g[b][c]*Ric[d][a]+g[b][d]*Ric[c][a])/2
        assert W==0
    # For Ric=c(u)l*l and parallel l=du, Box Ric=c'' g^{uu} l*l=0.
    assert gi[0][0]==0
    return {'r':str(r),'rprime':str(p),'rsecond':str(pp),'Ricci_uu':str(Ric[0][0]),
      'all_Weyl_components_zero':True,'Ricci_square_zero':True,'parallel_du':True}

def bump(u):
    y=2*u+3
    if abs(y)>=1:return 0.0,0.0
    v=1-y*y;f=math.exp(-(y*y)/v)
    return f,-4*y*f/(v*v)

def rhs(u,z):
    r,p,eta=z[:3];f,fp=bump(u);q=fp*fp/36000
    s=p*f+r*fp;J=(-u)*s*s;B=(-u)*r*r*fp*fp
    M=(-eta*r*r)*s*s;C=(r*p+(-u)*q*r*r)*f*f
    Eg=(2*u+2)*r*r*fp*fp/18000-2*(u+2)*p*p
    return np.array([p,-q*r,1/(r*r),fp*fp,B,J,M,C,Eg],dtype=float)

def rk(n):
    z=np.zeros(9);z[0]=1;z[2]=-1;h=-1.0/n
    minratio=1.0;maxratio=0.0
    for i in range(n):
      u=-1+i*h
      ratio=(-z[2]*z[0]*z[0])/(-u)
      minratio=min(minratio,ratio);maxratio=max(maxratio,ratio)
      k1=rhs(u,z);k2=rhs(u+h/2,z+h*k1/2);k3=rhs(u+h/2,z+h*k2/2);k4=rhs(u+h,z+h*k3)
      z+=h*(k1+2*k2+2*k3+k4)/6
    r,p,eta=map(float,z[:3]);E,B,J,M,C,Eg=map(float,-z[3:])
    area=1-r*r; uc=-2-r/p
    residuals={'J_minus_B_minus_C':J-B-C,'Raychaudhuri_balance':area-B/18000-Eg}
    assert max(map(abs,residuals.values()))<1e-10
    assert E>4 and E<36 and uc<-500
    assert minratio>float(F(499,500)**2)-1e-12 and maxratio<1+1e-12
    assert abs((M-B)/B)<float(F(591628623,124500500000))
    return {'steps':n,'E':E,'r_minus2':r,'rprime_minus2':p,'eta_minus2':eta,
      'past_Rosen_zero':uc,'B':B,'J':J,'M_eta':M,'C_integral':C,
      'relative_actual_difference':(M-B)/B,'area_delta':area,'Egeom':Eg,
      'sampled_weight_ratio_min':minratio,'sampled_weight_ratio_max':maxratio,
      'identity_residuals':residuals}

def gauss(n):
    a,w=np.polynomial.legendre.leggauss(n)
    y=np.array([bump(-1.5+0.5*float(t))[1]**2 for t in a])
    return {'order':n,'E':float(np.dot(w,y)/2)}

def calculate():
    guards=fractions()
    tensors=[tensor_jet(*xs) for xs in ((F(1,2),F(2,3),F(-3,5)),(F(499,500),F(1,1000),F(-1,36000)),(F(3,2),F(-4,7),F(5,9)))]
    numerical=[rk(n) for n in (4096,8192)]
    g=[gauss(n) for n in (128,256,512)]
    assert abs(numerical[-1]['E']-g[-1]['E'])<1e-9
    names=['admission.md','fixed_spec.md','proof.md','check.py','results.json','independent_check.py','../1279/proof.md','../1279/results.json']
    return {'round':1280,'status':'PASS_WITH_SCOPE','independent_groups':1,
      'strict_fraction_certificate':guards,
      'curvature_exact_jet_diagnostics':{'continuous_identity_proved_in_review_not_from_samples':True,'jets':tensors},
      'numerical_diagnostics':{'rigorous_floating_error':False,'author_implementation_imported':False,'RK4':numerical,'Gauss':g},
      'reviewed_scope':{'natural_eta_reference_new':True,'transverse_topology':'R^2, response densities or finite transverse null subalgebra',
        'bulk_reference_Hadamard':True,'absolute_minimal_reference_stress_zero_in_fixed_scheme':True,
        'local_normal_coherent_state':True,'global_finite_energy_Fock_state':False,
        'actual_observer_KMS_certified':False,'full_finite_resource_readout_certified':False,'whole_goal_completed':False},
      'bindings':{name:sha(HERE/name) for name in names}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate();p=HERE/'independent_results.json'
    if args.write:
      with p.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
      mode='CREATED_ONCE'
    else:
      assert json.loads(p.read_text('utf-8'))==result,'saved result differs'
      mode='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'mode':mode,'groups':1,
      'bound':result['strict_fraction_certificate']['relative_bound'],
      'diagnostic':result['numerical_diagnostics']['RK4'][-1]},ensure_ascii=False))
