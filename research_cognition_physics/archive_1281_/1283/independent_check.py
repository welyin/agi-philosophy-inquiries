"""1283 independent audit: radial chi-square moments and exact group jets.
No author/old calculation imports. All diagnostics use Fraction arithmetic.
Default is read-only replay; --write creates independent_results.json once.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]

def poch(a, n):
    out = F(1)
    for j in range(n):
        out *= a + j
    return out

def chi_moment(variance, dimension, n):
    return (2 * variance) ** n * poch(F(dimension, 2), n)

def radial_moments(s, t):
    # A=|y|^2, B=z_parallel^2, C=|z_perp|^2 are independent
    # scaled chi-square variables of dimensions 4,1,3.
    A = [chi_moment(s, 4, n) for n in range(5)]
    B = [chi_moment(t, 1, n) for n in range(5)]
    C = [chi_moment(t, 3, n) for n in range(5)]
    zeta = (A[1] - B[1] - C[1]) / 2
    zeta2 = (A[2] + B[2] + C[2] - 2*A[1]*B[1]
             - 2*A[1]*C[1] + 2*B[1]*C[1]) / 4
    zeta_D = (A[2]*C[1] - A[1]*B[1]*C[1] - A[1]*C[2]) / 2
    D = [A[n]*C[n] for n in range(5)]
    assert D[1] == 12*s*t and D[2] == 360*s*s*t*t
    assert zeta2 == 6*(s*s+t*t)-8*s*t
    assert zeta_D == 36*s*t*(s-t)
    for n in range(1, 5):
        assert D[n] == (4*s*t)**n*poch(F(2),n)*poch(F(3,2),n)
    return dict(zeta=zeta,zeta2=zeta2,zeta_D=zeta_D,D=D)

# Exact complex pairs; no floating-point differentiation.
def cmul(a,b):
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cscale(c,a):
    return (c*a[0],c*a[1])
def cadd(a,b):
    return (a[0]+b[0],a[1]+b[1])
def cabs2(a):
    return a[0]*a[0]+a[1]*a[1]
def creprod(a,b):
    return a[0]*b[0]+a[1]*b[1]

def group_energy(s,t,kappa,theta,m):
    # Along exp(h J_a), zeta(h)=cos(h/2) zeta+sin(h/2) d_a.
    # Group-dependent wavefunction is exp((q-i theta*kappa)*zeta).
    # Summing its exact second jets gives L psi_theta / psi_theta.
    z = ((s-t)/(4*s*t),-theta*kappa)
    p = cscale(F(3,4),z)
    d = cscale(F(-1,4),cmul(z,z))
    mean = cadd(cscale(m['zeta'],p),cscale(m['D'][1],d))
    assert mean[1] == 0
    squared = (cabs2(p)*m['zeta2']
               +2*creprod(p,d)*m['zeta_D']
               +cabs2(d)*m['D'][2])
    return mean[0],squared

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def qm(a,b):
    # Quaternion multiplication, using the oriented unit table.
    cyclic={(1,2):3,(2,3):1,(3,1):2}
    out=[F(0)]*4
    for i in range(4):
        for j in range(4):
            if i == 0: k,sign=j,1
            elif j == 0: k,sign=i,1
            elif i == j: k,sign=0,-1
            elif (i,j) in cyclic: k,sign=cyclic[(i,j)],1
            else: k,sign=cyclic[(j,i)],-1
            out[k] += sign*a[i]*b[j]
    return out

def exact_group_jet():
    x=[F(2,3),F(-1,5),F(3,4),F(5,7)]
    w=[F(-2,7),F(4,9),F(1,6),F(-3,5)]
    zeta=dot(x,w)
    determinant=dot(x,x)*dot(w,w)-zeta*zeta
    first=[];second=[]
    for a in range(1,4):
        e=[F(0)]*4;e[a]=F(1)
        jw=[v/2 for v in qm(e,w)]
        jjw=[v/2 for v in qm(e,jw)]
        first.append(dot(x,jw));second.append(dot(x,jjw))
        assert jjw == [-v/4 for v in w]
    assert sum(v*v for v in first)==determinant/4
    assert sum(second)==-3*zeta/4
    return dict(zeta=str(zeta),determinant=str(determinant),
                first=[str(v) for v in first],
                laplacian=str(sum(second)),residual='0')

def response_case(s,t,kappa,akin,b,c,theta):
    m=radial_moments(s,t)
    L0,L20=group_energy(s,t,kappa,F(0),m)
    assert L0==3*(s-t)**2/(16*s*t)
    assert L20==(9*(s-t)**2/(64*s*t)
                 +45*(s-t)**4/(512*s*s*t*t))
    # Matter gradient energy: F=kappa|z|^2, so |grad_z F|^2=4 kappa^2 |z|^2.
    matter_second=akin*4*kappa*kappa*chi_moment(t,4,1)
    base_other=akin*(1/s+1/t+F(1,2))+kappa*chi_moment(t,4,1)
    energies=[];group_rows=[]
    for th in (F(0),theta,2*theta):
        L,L2=group_energy(s,t,kappa,th,m)
        e=base_other+matter_second*th*th+b*L+c*L2
        energies.append(e);group_rows.append([str(L),str(L2)])
    quartic=(kappa*kappa/4)**2*m['D'][2]
    A2=(16*akin*kappa*kappa*t+3*b*kappa*kappa*s*t
         +9*c*kappa*kappa*(5*(s*s+t*t)-6*s*t)/16)
    for th,e in zip((F(0),theta,2*theta),energies):
        assert e==energies[0]+A2*th*th+c*quartic*th**4
    residual=energies[2]-4*energies[1]+3*energies[0]
    assert residual==12*c*quartic*theta**4
    assert A2>0
    # Reality gives exact evenness, checked at an independent negative phase.
    assert group_energy(s,t,kappa,-theta,m)==group_energy(s,t,kappa,theta,m)
    return dict(parameters={k:str(v) for k,v in dict(s=s,t=t,kappa=kappa,akin=akin,b=b,c=c,theta=theta).items()},
        moments={k:str(m[k]) for k in ('zeta','zeta2','zeta_D')},
        D_moments=[str(v) for v in m['D'][1:]],
        L0=str(L0),L20=str(L20),group_L_and_L2=group_rows,
        energies_minus_same_W=[str(v) for v in energies],
        A2=str(A2),C4=str(quartic),doubling_residual=str(residual))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def calculate():
    fixed=(F(2),F(1),F(4,5),F(1,2),F(3,10))
    rows=[response_case(*fixed,c,F(1,2)) for c in (F(0),F(1,10))]
    rows.append(response_case(F(3,2),F(2,5),F(7,6),F(2,3),F(4,7),F(1,9),F(2,3)))
    assert F(rows[0]['A2'])==F(784,125)
    assert F(rows[1]['A2'])==F(337,50)
    assert F(rows[1]['doubling_residual'])==F(1728,625)
    denominator=12*F(rows[1]['C4'])*F(1,2)**4
    eta=F(1,10000);eta_R=eta+4*eta+3*eta
    assert denominator==F(3456,125)
    assert eta_R/denominator==F(1,34560)
    # Sharp triangle-budget signs for the three separate reports.
    errors=[eta,-eta,eta]
    assert errors[2]-4*errors[1]+3*errors[0]==eta_R
    # Deletion of finite-amplitude exactness: positive c remains O(theta^2).
    c=F(1,10);theta=F(1,2);C4=F(rows[1]['C4'])
    weak_remainders=[c*C4*(theta/2**j)**2 for j in range(3)]
    assert weak_remainders[1]==weak_remainders[0]/4
    # The exact response does not select b; these are analytic source curves.
    s,t,kappa,akin,_=fixed
    A2_other_b=16*akin*kappa*kappa*t+3*F(7,10)*kappa*kappa*s*t
    assert A2_other_b!=F(rows[0]['A2'])
    # Highest-symbol diagnostics through radial moments, not author Wick polynomials.
    moments=radial_moments(s,t)
    high=[]
    for n in (2,3,4):
        coefficient=(kappa*kappa/4)**n*moments['D'][n]
        assert coefficient==(kappa*kappa*s*t)**n*poch(F(2),n)*poch(F(3,2),n)>0
        high.append(dict(n=n,unit_leading_coefficient=str(coefficient)))
    paths=[HERE/'proof.md',HERE/'check.py',HERE/'results.json',HERE/'NEXT.md',
           HERE.parent/'research_note_1283.md',Path(__file__).resolve(),
           RESEARCH/'archive_554_584/research_note_557.md',HERE.parent/'1282/proof.md']
    return dict(round=1283,status='PASS_WITH_SCOPE',independent_groups=1,new_axioms=0,
        method='exact conditional chi-square moments and quaternion geodesic jets; no author calculation import',
        exact_group_jet=exact_group_jet(),source_response_cases=rows,
        finite_error=dict(each_energy_error=str(eta),weighted_error=str(eta_R),
            denominator=str(denominator),c_bound_when_reported_residual_zero=str(eta_R/denominator)),
        deletion_checks=dict(positive_c_has_leading_quadratic_Taylor=True,
            quartic_correction_div_theta_squared=[str(v) for v in weak_remainders],
            distinct_b_also_passes_exact_quadratic=True,other_b_A2=str(A2_other_b),
            degenerate_kappa_zero_would_remove_detector=True),
        polynomial_highest_moments=high,
        scope=dict(same_normal_Gauss_preparation=True,original_W_is_common_finite_term=True,
            c_L2_baseline_retained=True,effective_phase_permission_explicit=True,
            exact_open_interval_is_stronger_than_finite_sample_reports=True,
            finite_error_bound_only_degree_two_family=True,
            numerical_response_is_all_unknown_inputs=False,
            autonomous_controller_or_energy_meter_constructed=False,
            phase_source_work_is_total_control_cost=False,
            gauge_group_or_all_dynamics_selected=False),
        source_sha256={str(p.relative_to(RESEARCH)).replace('\\','/'):sha(p) for p in paths})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true')
    args=ap.parse_args();data=calculate();dest=HERE/'independent_results.json'
    if args.write:
        with dest.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(data,out,ensure_ascii=False,indent=2);out.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(dest.read_text('utf-8'))==data,'Independent result/source binding mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=data['status'],mode=mode,round=1283,
        exact_cases=len(data['source_response_cases']),C4=data['source_response_cases'][1]['C4'],
        doubling_residual=data['source_response_cases'][1]['doubling_residual'],
        conditional_c_bound=data['finite_error']['c_bound_when_reported_residual_zero']),
        ensure_ascii=False))

if __name__=='__main__':
    main()
