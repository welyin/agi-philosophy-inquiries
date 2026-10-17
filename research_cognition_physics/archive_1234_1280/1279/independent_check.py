"""1279 independent audit: analytic Fraction guards + independent RK/Gauss diagnostics.
No import of the author implementation. Floats are diagnostics, never certified bounds.
Default read-only match; --write creates the independent result exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def frac(x):return str(x)

def analytic():
    U=F(500); E_low=F(4); E_high=F(36); Aq=E_high/F(36000)
    # Independent bump energy proof: y=2x-1, t=y/sqrt(1-y^2).
    # E=8 int_R t^2(1+t^2)^(3/2)exp(-2t^2)dt.
    # Bound (1+t^2)^(3/2)<=(1+t^2)^2; Gaussian moment ratios:
    moments={0:F(1),2:F(1,4),4:F(3,16),6:F(15,64)}
    polynomial=moments[2]+2*moments[4]+moments[6]
    E_upper_alternative=8*polynomial*F(3,2) # sqrt(pi/2)<3/2.
    assert polynomial==F(55,64) and E_upper_alternative<E_high
    r_global=1-U*Aq; r_pulse=1-2*Aq
    assert r_global==F(1,2) and r_pulse==F(499,500)
    # Endpoint Cauchy: sup f^2<=E/4 and int_0^1 f^2<=E/4.
    C_over_AE=F(1,4)+F(2,4)
    assert C_over_AE==F(3,4)
    delta_C=C_over_AE*Aq/(r_pulse*r_pulse)
    delta_w=F(4)/(U*U)
    total=delta_C+delta_w*(1+delta_C)
    target=F(1,1250)
    assert delta_C==F(375,498002)
    assert delta_w==F(1,62500)
    assert total==F(23935877,31125125000)<target
    margin=target-total
    assert margin==F(964223,31125125000)
    assert 1/(1-delta_w)==F(62500,62499)
    # s*(8piG)=2pi, represented as the coefficient of pi.
    beta_g=F(2,36000); s_over_pi=F(36000)
    assert beta_g*s_over_pi==2
    return {'method':'Continuous bounds are proved analytically; exact rational constants checked here.',
      'E_lower':frac(E_low),'E_upper_working':frac(E_high),
      'independent_Gaussian_moment_polynomial':frac(polynomial),
      'independent_Gaussian_energy_upper':frac(E_upper_alternative),
      'Aq_upper':frac(Aq),'r_whole_domain_lower':frac(r_global),
      'r_pulse_lower':frac(r_pulse),'C_upper_over_Aq_E':frac(C_over_AE),
      'delta_C':frac(delta_C),'delta_w':frac(delta_w),
      'relative_bound':frac(total),'relative_bound_decimal':float(total),
      'strict_target':frac(target),'strict_margin':frac(margin),
      'clock_rate_upper':frac(1/(1-delta_w)),
      'B_positive_lower':frac(r_pulse*r_pulse*E_low),
      'eight_pi_G':frac(beta_g),'s_over_pi':frac(s_over_pi),
      'new_independent_groups':1}

def bump(u):
    y=2*(u+2)-1
    if abs(y)>=1:return 0.0,0.0
    v=1-y*y
    f=math.exp(-(y*y)/v)
    return f,-4*y*f/(v*v)

def rhs(u,state):
    r,p=state[:2];f,fp=bump(u);q=fp*fp/36000.0
    w=500.0*math.atanh(u/500.0)
    xi=-w*(1-(u/500.0)**2)
    psi_p=p*f+r*fp
    energy=fp*fp
    B=(-u)*r*r*energy
    J=(-u)*psi_p*psi_p
    M=xi*psi_p*psi_p
    C=(r*p+(-u)*q*r*r)*f*f
    Eg=(2*u+2)*r*r*energy/18000.0-2*(u+2)*p*p
    return np.array([p,-q*r,energy,B,J,M,C,Eg],dtype=float)

def rk(n):
    h=-1.0/n;state=np.zeros(8);state[0]=1.0
    for i in range(n):
        u=-1.0+i*h
        k1=rhs(u,state);k2=rhs(u+h/2,state+h*k1/2)
        k3=rhs(u+h/2,state+h*k2/2);k4=rhs(u+h,state+h*k3)
        state += h*(k1+2*k2+2*k3+k4)/6
    r,p=map(float,state[:2]);E,B,J,M,C,Eg=map(float,-state[2:])
    area=1-r*r
    residual=area-B/18000.0-Eg
    identities={'J_minus_B_minus_C':J-B-C,'finite_Raychaudhuri_balance':residual}
    assert max(abs(v) for v in identities.values())<1e-10
    assert 0<C and 0<B and abs(M-B)/B<float(F(23935877,31125125000))
    assert r>=float(F(499,500)) and r+(-500+2)*p>0.5
    return {'steps':n,'E':E,'r_minus2':r,'rprime_minus2':p,
      'r_at_minus500_by_exact_linear_extension':r+(-500+2)*p,
      'B':B,'J':J,'M':M,'C_independent_integral':C,'relative_actual_difference':(M-B)/B,
      'area_delta':area,'Egeom_independent_integral':Eg,'identity_residuals':identities}

def gauss(n):
    xs,weights=np.polynomial.legendre.leggauss(n)
    fp=np.array([bump(-1.5+0.5*float(t))[1] for t in xs])
    return {'order':n,'E':float(np.dot(weights,fp*fp)/2)}

def calculate():
    a=analytic(); g=[gauss(n) for n in (128,256,512)]
    d=[rk(n) for n in (4096,8192)]
    for row in g:assert 4<row['E']<36
    assert max(row['E'] for row in g)-min(row['E'] for row in g)<1e-9
    assert abs(d[-1]['E']-g[-1]['E'])<1e-9
    changes={key:abs(d[-1][key]-d[0][key]) for key in ('E','r_minus2','rprime_minus2','B','J','M','C_independent_integral','area_delta','Egeom_independent_integral')}
    # Bind evidence, do not use author results to construct this calculation.
    names=['admission.md','fixed_spec.md','proof.md','check.py','results.json','independent_check.py']
    bindings={n:sha(HERE/n) for n in names}
    return {'round':1279,'status':'PASS_WITH_SCOPE','independent_groups':1,
      'analytic_fraction_certificate':a,
      'diagnostic_only':{'floating_error_certified':False,'method':'Independent direct backward RK4 and Gauss-Legendre; no author import.',
        'gauss':g,'rk':d,'coarse_fine_absolute_changes':changes},
      'source_bindings':bindings,
      'scope':{'same_fixed_background_reference_subtraction':True,
        'normal_coherent_null_algebra':True,'null_zero_modes_excluded':True,
        'absolute_reference_source_certified':False,'bulk_Hadamard_state_certified':False,
        'actual_clock_instrument_certified':False,'complete_goal_certified':False}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate();path=HERE/'independent_results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text('utf-8'))==result,'saved independent result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'mode':mode,'independent_groups':1,
       'exact_relative_bound':result['analytic_fraction_certificate']['relative_bound'],
       'fine_diagnostic':result['diagnostic_only']['rk'][-1]},ensure_ascii=False))
