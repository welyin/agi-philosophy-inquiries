"""1279: rational finite bounds and explicitly non-certifying numerical diagnostics.
Default replays saved results without mutation; --write creates them exactly once.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def bump(u):
    x=u+2.0
    if x<=0.0 or x>=1.0: return 0.0,0.0
    h=4*x*(1-x)
    exponent=1-1/h
    if exponent < -740: return 0.0,0.0
    f=math.exp(exponent)
    return f,f*(1-2*x)/(4*x*x*(1-x)*(1-x))

def derivative(u,y):
    f,fp=bump(u)
    r,rp=y[:2]; q=fp*fp/36000; z=-u
    xi=500*math.atanh(z/500)*(1-z*z/250000)
    p=rp*f+r*fp
    return np.array([rp,-q*r,fp*fp,z*r*r*fp*fp,z*p*p,xi*p*p,
                     (r*rp+z*q*r*r)*f*f,
                     (2*u+2)*r*r*(2*q)-(u+2)*2*rp*rp],dtype=float)

def integrate(n):
    y=np.zeros(8);y[0]=1
    step=-1/n
    for i in range(n):
        u=-1+i*step
        k1=derivative(u,y)
        k2=derivative(u+step/2,y+step*k1/2)
        k3=derivative(u+step/2,y+step*k2/2)
        k4=derivative(u+step,y+step*k3)
        y+=step*(k1+2*k2+2*k3+k4)/6
    r,rp=y[:2]
    E,B,J,M,C,Eg=-y[2:]
    area=1-r*r
    return dict(steps=n,E=float(E),B=float(B),J=float(J),M=float(M),C=float(C),
                r_pulse_past=float(r),rprime_pulse_past=float(rp),
                r_at_minus_U=float(r-498*rp),area_change=float(area),geometry_correction=float(Eg),
                rescaling_identity_residual=float(J-B-C),
                raychaudhuri_identity_residual=float(area-B/18000-Eg),
                relative_modular_source_difference=float((M-B)/B))

def compute():
    A=F(1,1000);U=F(500)
    rmin=1-U*A;rpulse=1-2*A
    dc=3*A/(4*rpulse*rpulse);dw=4/(U*U)
    total=dc+dw*(1+dc);target=F(1,1250)
    assert rmin==F(1,2) and rpulse==F(499,500)
    assert dc==F(375,498002) and dw==F(1,62500)
    assert total==F(23935877,31125125000)<target
    assert target-total==F(964223,31125125000)
    assert F(8,3)**2>F(64,9)-F(1,100) # lower bound on e follows from its positive series
    e_series=sum((F(1,math.factorial(i)) for i in range(5)),F())
    assert e_series>F(8,3) and 16/e_series<6
    gauss=[]
    for n in (128,256,512):
        z,w=np.polynomial.legendre.leggauss(n)
        val=sum(float(wi)*bump(float(zi/2-1.5))[1]**2/2 for zi,wi in zip(z,w))
        gauss.append(dict(nodes=n,E=val))
    rows=[integrate(n) for n in (4096,8192)]
    for row in rows:
        assert 4<row['E']<36
        assert row['r_at_minus_U']>float(rmin)
        assert abs(row['relative_modular_source_difference'])<float(total)
        assert abs(row['rescaling_identity_residual'])<1e-9
        assert abs(row['raychaudhuri_identity_residual'])<1e-10
    bindings=['admission.md','fixed_spec.md','proof.md','check.py',
              '../_shared/ct_hg_adoption_after1278.md','../_shared/jacobson_reconstruction_after1278.md']
    return dict(schema='round1279_null_reference_matching_v1',status='PASS_WITH_SCOPE',round=1279,
                author_check_groups=1,new_adopted_axioms=0,
                strict_certificate=dict(method='Analytic inequalities plus exact Fraction arithmetic; not inferred from samples',
                    E_lower='4',E_strict_upper='36',A_upper=str(A),U=str(U),
                    r_global_lower=str(rmin),r_pulse_lower=str(rpulse),
                    relative_rescaling_bound=str(dc),relative_window_bound=str(dw),
                    relative_total_bound=str(total),target=str(target),strict_margin=str(target-total),
                    clock_rate_excess_bound=str(dw/(1-dw))),
                numerical_diagnostics=dict(scope='Uncertified floating RK4/Gauss diagnostics only; differences are not rigorous numerical error bars',
                    gauss=gauss,rk4=rows,
                    RK_B_difference=abs(rows[0]['B']-rows[1]['B'])),
                physical_scope=dict(reference='Specified w-null vacuum, coherent difference on one fixed geometry',
                    bulk_Hadamard_proved=False,absolute_reference_source_proved=False,
                    actual_clock_instrument_proved=False,full_goal_completed=False),
                source_url='https://arxiv.org/html/1502.04106',
                bindings={p:sha(HERE/p) for p in bindings})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=compute();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        saved=json.loads(path.read_text('utf-8'))
        assert saved==result,'Saved result differs'
        mode='READ_ONLY_MATCH'
    row=result['numerical_diagnostics']['rk4'][-1]
    print(json.dumps(dict(status=result['status'],mode=mode,
        bound=result['strict_certificate']['relative_total_bound'],diagnostic_E=row['E'],
        diagnostic_relative_difference=row['relative_modular_source_difference'],
        full_goal_completed=False),ensure_ascii=False))
