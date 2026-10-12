"""Nonformal finite actual-BT defect qualification; default is read-only.
The exact budgets audit the stated analytic estimates, not their domain proof.
The finite matrix and scalar diagnostics do not simulate a continuum source.
Only --write-results writes this qualification's JSON; no round is reserved.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
RESULT=HERE/'after1141_defect_results.json'
FROZEN={'1128/proof.md': 'fc9963b8f77d0a60e2233acbcc68e02e583af345a5653d1d4431ea619be0fcd2', '1129/proof.md': 'c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa', '1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d', '1134/proof.md': '4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666', '1137/proof.md': 'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb', '1138/proof.md': '5668334517745145330f0504905c11bdaf70b96a374080180177c0ddf5272110', '1141/proof.md': 'd71543b186660711ac1602a15495c5700ba2430577366d3057f70249fb650c96', '1141/results.json': '62d67da39ec6e6c016f5d6150760ed21b84729a96dcbc0ee5b81f9891919376a', '1141/acceptance.json': 'fb78daab0c9d9d75102b74d2770f9919494d0cae230f82a72e004aa9a12a1df7', '_shared/after1137_local_balance_preflight.md': '95915f871317e8638d23c3ee29cec779a0ff52aad5f76873a0ff896f6a2fcbe0', '_shared/after1141_BT_power_qualification.md': '3269ee80aa0ce3ad55db4bf621f5d312f7628a7c88a671a0b2e187ce0a2d13f6', '_shared/after1141_BT_power_check.py': '5430100453dca062fa4df3ca40aeb0e8ce5858aa92ab4bca474733f23caf30ac', '_shared/after1141_BT_power_results.json': '9768b0529ba3bd49a81db67868ecd15c933eddaa5a3c4c81ee73da459790e223', '_shared/after1141_defect_direct_qualification.md': '740df64ce95b6358a2e46593569856ce1c57fdee66cc2a518664f4b2ae874cb7'}

def exact(x):
    x=F(x)
    with localcontext() as ctx:
        ctx.prec=45
        val=str(Decimal(x.numerator)/Decimal(x.denominator))
    return {'numerator':str(x.numerator),'denominator':str(x.denominator),'decimal':val}

def sqrt_upper(x, digits=30):
    x=F(x); scale=10**digits
    n=isqrt(x.numerator*scale*scale//x.denominator)+1
    bound=F(n,scale)
    assert bound*bound>x
    return bound

def rational_budgets():
    A,M,eps,T,sigma=F(10**34),F(10**80),F(1,10**5),F(5,4),F(1,10**30)
    Cs,lam=F(33536000000000000000),F(11,7)
    v=16*lam*Cs+F(558,7); mstar=3*M-v
    CY=F(24780800000000000000000,49)+16316
    C2x=F(9*10**21); C1=4*CY; C2=2*C2x+20*CY*CY/mstar
    assert v<10**21 and mstar>2*M
    assert C1<3*10**21 and C2<2*10**22
    p1=268+94*A; p2=28000+188*268*A+12000*A*A
    A1=1/sigma+p1; A2=1/sigma**2+2*p1/sigma+p2
    P1=F(10**36); P2=F(1201,1000)*10**72
    p_actual=p1+C1*T
    p2_actual=p2+2*C1*p1*T+C1*C1*T*T+C2*T
    assert p_actual<P1 and p2_actual<P2 and A1<P1 and A2<P2

    # epsilon^(-7/2) <= epsilon^(-3)*317: certified by squaring.
    assert F(317)**2>1/eps
    eps_72=eps**-3*317
    Vh=F(210000)
    assert 400*524<Vh
    c_kernel=12*Vh*(66*eps**-3+30*eps_72)
    d_kernel=12*Vh*(30*eps_72+142*eps**-4)
    c,d=F(10**27),F(10**30)
    assert c_kernel<c and d_kernel<d

    D=F(12,10)*10**21; q=D/A; nu=100*eps**-4
    assert 64*2146<F(300)**2*2
    assert F(128)+F(16,7)<131
    assert 131**2*2146<F(5000)**2*2
    p_comparison=p1+CY*T
    r=A/M*(300*(p_comparison+A1)+9000)
    assert r<F(1,10**7)
    source_derivative=q*sqrt_upper(2+600*nu)
    bootstrap_derivative=source_derivative+2*r
    assert bootstrap_derivative<F(295,1000)
    bootstrap_e=T*bootstrap_derivative
    assert bootstrap_e<F(4,10)

    Bstar=F(1074,1000)*10**71
    assert Bstar<(33*A)**2
    Bc_expression=Bstar+2*(33*A)**2+1073*A*A+94*A
    Bc=F(44,10)*10**71
    assert Bc_expression<Bc
    delta_exact=T*v*A2/(2*mstar*mstar)
    delta_safe=T*v*A2/(8*M*M)
    assert delta_exact<delta_safe<F(1581,1000)*F(1,10**68)
    e_actual=F(4,10)+2*delta_safe*Bc
    assert e_actual<15000

    assert 600*nu<(4*10**12)**2
    local_column_base=2*sqrt_upper(15000)+4*10**12
    assert local_column_base<5*10**12
    # The proof's interaction contribution is 4 sqrt(v Bc s).
    assert 16*v*Bc<(10**47)**2
    # A common conservative mixed-column constant for Pf and K.
    mixed_column=2*(sqrt_upper(P2*Bstar)+33*A*P1)+4*10**12*P1
    assert mixed_column<F(15,10)*10**72
    Lc=F(16,10)*10**73
    assert 8*(F(15,10)*10**72)+400*A*P1==Lc
    g0=F(1,10**8)
    assert A2/(2*mstar)<g0

    B2src=lam*32768*Cs+42*128*10**12+(lam+F(11,14))*32768
    ccr=4*eps**-3*B2src/(A*A)
    NP=D/A*(47*A+230*eps**-2)+ccr
    Nc=F(2*10**98)
    assert c*Bstar+d+NP<Nc
    # Each coefficient of the polynomial domination of ||Pf g Phi_s||.
    pg_constant=sqrt_upper(A2*P2)+5*CY
    pg_sqrt_squared=A2*(2*C1*P1+C2)
    pg_linear=sqrt_upper(A2)*C1
    assert pg_constant<2*10**72
    assert pg_sqrt_squared<(10**65)**2
    assert pg_linear<4*10**57

    # Exact integrals of 1 and s, with Gamma(3/2)<1 for sqrt(s).
    # M^(3/2) is rational here, since M=10^80.
    M32=F(10**120); assert M32*M32==M**3
    recoil=6*Lc*(F(5*10**12)/(2*M)+600*A/(4*M*M))
    outer=6*Lc*(F(5*10**12)/(4*M)+1800*A/(16*M*M)+F(10**47)/(8*M32))
    gCo=2*g0*Nc/(4*M)
    gCR=6*Nc/M*(F(2*10**72)/(4*M)+F(10**65)/(8*M32)+F(4*10**57)/(16*M*M))+3*Nc*g0/(4*M*M)
    assert recoil<3*10**6 and outer<2*10**6
    assert gCo<2*10**10 and gCR<7*10**10
    point=recoil+outer+gCo+gCR
    integrated=T*point
    assert point<10**11 and integrated<2*10**11

    values=locals()
    names='A M eps T v mstar CY C2x C1 C2 p1 p2 A1 A2 P1 P2 p_actual p2_actual Vh c_kernel d_kernel c d D q nu p_comparison r source_derivative bootstrap_derivative bootstrap_e Bstar Bc_expression Bc delta_exact delta_safe e_actual local_column_base mixed_column Lc g0 NP Nc pg_constant pg_sqrt_squared pg_linear recoil outer gCo gCR point integrated'.split()
    return {'status':'PASS','arithmetic':'Fraction plus squared rational upper bounds',
            'analytic_graph_and_domain_proof_verified_by_code':False,
            'qualification_intermediate_constants_require_formal_clarification':True,
            'reviewed_graph_intermediate_c_formula':'12 Vh (66 epsilon^-3 + 30 epsilon^-7/2)',
            'reviewed_graph_intermediate_d_formula':'12 Vh (30 epsilon^-7/2 + 142 epsilon^-4)',
            'values':{name:exact(values[name]) for name in names}}

def add(a,b):
    return [[a[i][j]+b[i][j] for j in range(len(a))] for i in range(len(a))]

def scale(a,s):
    return [[F(s)*x for x in row] for row in a]

def mul(a,b):
    n=len(a)
    return [[sum((a[i][k]*b[k][j] for k in range(n)),F(0)) for j in range(n)] for i in range(n)]

def comm(a,b): return add(mul(a,b),scale(mul(b,a),-1))
def anti(a,b): return add(mul(a,b),mul(b,a))
def diag(xs): return [[F(x) if i==j else F(0) for j,x in enumerate(xs)] for i in range(len(xs))]
def zero(a): return all(x==0 for row in a for x in row)
def sylvester(diagonal,b): return [[b[i][j]/(diagonal[i]+diagonal[j]) for j in range(len(b))] for i in range(len(b))]
def encode_matrix(a): return [[str(x) for x in row] for row in a]

def exact_matrix_identity():
    # C quantities below omit one common i.  All arrays are exact rational.
    # Pf=K means SY=m0 I and [K-Pf,E]=0; all three axes contribute.
    E=[[F(1),F(2)],[F(2),F(-1)]]
    K=[diag([0,F(8,3)]),diag([0,F(4,3)]),diag([0,F(8,3)])]
    Pf=K; mass=diag([8,F(15,2)]); S=diag([3,5]); SY=diag([3,3])
    H=diag([8,F(17,2)]); g=add(H,scale(mass,-1))
    Ho=add(add(mass,scale(S,-1)),SY)
    I=diag([1,1]); z=diag([0,0])
    Ksq=z
    for k in K: Ksq=add(Ksq,mul(k,k))
    assert mul(S,S)==add(scale(I,9),Ksq)
    assert mul(H,H)==add(mul(mass,mass),Ksq)
    assert zero(comm(SY,E)) and all(zero(comm(k,mass)) for k in K)
    assert all(zero(comm(add(k,scale(p,-1)),E)) for k,p in zip(K,Pf))
    C=[comm(k,E) for k in K]; CR=comm(S,E); Co=comm(Ho,E)
    assert not zero(CR) and all(not zero(cj) for cj in C)
    internal_rhs=z; external_rhs=z
    for k,p,cj in zip(K,Pf,C):
        internal_rhs=add(internal_rhs,anti(p,cj))
        external_rhs=add(external_rhs,anti(k,cj))
    assert anti(S,CR)==internal_rhs
    assert sylvester([3,5],internal_rhs)==CR
    direct=comm(add(H,scale(Ho,-1)),E)
    old_rhs=add(add(anti(mass,CR),external_rhs),scale(anti(g,Co),-1))
    assert anti(H,direct)==old_rhs
    reordered=add(CR,sylvester([8,F(17,2)],add(add(external_rhs,scale(anti(g,Co),-1)),scale(anti(g,CR),-1))))
    assert direct==reordered
    # Retain a falsification witness: dropping internal recoil changes the result.
    no_recoil=sylvester([8,F(17,2)],add(external_rhs,scale(anti(g,Co),-1)))
    assert no_recoil!=direct
    return {'status':'PASS','arithmetic':'exact Gaussian-imaginary commutators represented by rational coefficients',
            'dimension':2,'nonzero_K_commutator_axes':3,'all_K_commute_with_mass':True,
            'SY_commutes_with_E':True,'K_minus_Pf_commutes_with_E':True,
            'source_mass_positive_gap':True,
            'internal_recoil':encode_matrix(CR),'direct_defect':encode_matrix(direct),
            'reordered_defect':encode_matrix(reordered),'dropping_recoil_fails':True,
            'continuum_domain_proof':False}

def scalar_diagnostics():
    rng=np.random.default_rng(1142)
    max_r_ratio=max_s_ratio=max_h_ratio=max_s_semigroup=max_h_semigroup=0.0
    negative_h=0; cases=0
    # Dimensionless moderate parameters. No huge-source numerical simulation.
    for j in range(512):
        M=float([0.5,1.0,3.0,10.0][j%4])
        m0=M*(1+0.25*rng.random())
        p=rng.normal(size=3)*M*3
        K=rng.normal(size=3)*M*9
        if j%16==0: p=np.zeros(3); K=np.array([-60*M,0.,0.])
        k=rng.normal(size=3)*M*(0.01+0.9*rng.random())
        if j%16==0: k=np.array([0.1*M,0.,0.])
        omega=float(np.linalg.norm(k)); ef=2*M+np.linalg.norm(p)+rng.random()*M
        S=float(np.hypot(m0,np.linalg.norm(p)))
        Sp=float(np.hypot(m0,np.linalg.norm(p+k)))
        dS=(2*float(p@k)+omega*omega)/(Sp+S)
        mass=ef+S; dmass=omega+dS
        H=float(np.hypot(mass,np.linalg.norm(K)))
        Hp=float(np.hypot(mass+dmass,np.linalg.norm(K+k)))
        dH=(2*mass*dmass+dmass*dmass+2*float(K@k)+omega*omega)/(Hp+H)
        dR=dH-omega
        bound=omega/M*(np.linalg.norm(p)+(2/3)*np.linalg.norm(K)+(5/3)*omega)
        assert abs(dR)<=bound*(1+1e-13)
        assert abs(dS)<=omega*(1+1e-13) and abs(dH)<=3*omega*(1+1e-13)
        assert S>=M and Sp>=M and H>=3*M and Hp>=3*M
        max_r_ratio=max(max_r_ratio,abs(dR)/bound)
        max_s_ratio=max(max_s_ratio,abs(dS)/omega)
        max_h_ratio=max(max_h_ratio,abs(dH)/(3*omega))
        negative_h+=int(dH<0)
        for u in [1e-4,0.01,0.1,1.,3.]:
            t=u/M
            # Stable exp differences: exp(-min*t)*[-expm1(-abs(delta)*t)].
            ds_exp=np.exp(-min(S,Sp)*t)*(-np.expm1(-abs(dS)*t))
            dh_exp=np.exp(-min(H,Hp)*t)*(-np.expm1(-abs(dH)*t))
            bs=t*omega*np.exp(-M*t)
            bh=3*t*omega*np.exp(-3*M*t)
            assert ds_exp<=bs*(1+1e-13) and dh_exp<=bh*(1+1e-13)
            max_s_semigroup=max(max_s_semigroup,float(ds_exp/bs))
            max_h_semigroup=max(max_h_semigroup,float(dh_exp/bh))
        cases+=1
    assert negative_h>0
    # Rounded diagnostics avoid platform-specific BLAS final-bit result changes.
    rd=lambda x:round(float(x),11)
    return {'status':'PASS','seed':1142,'free_charge_cases':cases,
            'positive_semigroup_samples':5*cases,'negative_free_H_increments_retained':negative_h,
            'max_finite_difference_R_to_bound':rd(max_r_ratio),
            'max_abs_S_increment_to_omega':rd(max_s_ratio),
            'max_abs_H_increment_to_3omega':rd(max_h_ratio),
            'max_S_semigroup_difference_to_bound':rd(max_s_semigroup),
            'max_H_semigroup_difference_to_bound':rd(max_h_semigroup),
            'analytic_continuum_or_all_spectrum_proof':False,
            'huge_source_state_simulated':False}

def run():
    observed={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in FROZEN}
    assert observed==FROZEN,'Qualification/source version changed; explicit review required'
    acc=json.loads((BASE/'1141/acceptance.json').read_text(encoding='utf-8'))
    assert acc['status']=='PASS' and acc['round']==1141
    return {'qualification':'after1141_actual_BT_defect','status':'PASS',
            'formal_round_reserved':False,'new_cognitive_axioms':0,'scientific_calibration_increment':0,
            'groups':{
                'frozen_source_identity':{'status':'PASS','sha256':observed,'round1141_acceptance_verified':True},
                'exact_fraction_budgets':rational_budgets(),
                'exact_constrained_matrix_sylvester':exact_matrix_identity(),
                'finite_scalar_diagnostics':scalar_diagnostics()},
            'scope':{'candidate_proof_requires_independent_analytic_domain_review':True,
                     'code_alone_establishes_actual_BT_balance':False,
                     'fixed_original_region_and_short_window_only':True,
                     'actual_relative_energy_and_flux_bounds_migrated':False,
                     'full_material_local_balance_or_stress_proved':False,
                     'finite_source_subaccount_candidate_only':True}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); data=run()
    if args.write_results:
        RESULT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert RESULT.exists(),'Saved result missing; generate explicitly with --write-results'
        assert data==json.loads(RESULT.read_text(encoding='utf-8')),'Saved result differs from recomputation'
    vals=data['groups']['exact_fraction_budgets']['values']
    print(json.dumps({'status':'PASS','mode':'write-results' if args.write_results else 'read-only',
                      'groups':len(data['groups']),'frozen_sources':len(FROZEN),
                      'candidate_formal_acceptance_claimed':False,
                      'bootstrap_e':vals['bootstrap_e']['decimal'],
                      'actual_e':vals['e_actual']['decimal'],
                      'defect_integrated_budget':vals['integrated']['decimal']},ensure_ascii=False,indent=2))
if __name__=='__main__': main()
