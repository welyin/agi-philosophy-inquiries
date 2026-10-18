"""1303: compact continuum flavor phase interface.
Exact rational inequalities certify the complete time window.
Quadrature and Decimal matrix arithmetic are diagnostics, not proofs of integration.
Default/--check read only; --write exclusively creates results.json.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import math
import numpy as np
HERE = Path(__file__).resolve().parent
OUT = HERE/'results.json'


def rational_certificate():
    a,L,mu,S = F(1,32),F(1,4),F(1,4096),F(55,256)
    assert S < F(1,4) and 3*mu < 1
    det = (1**3-3**3)*(2**5-3**5)-(1**5-3**5)*(2**3-3**3)
    assert det == 888
    K = 888*L*a*a/F(163840)
    assert K == F(111,83886080) and 2*mu*mu < K/2
    assert F(2188,12288) < F(1,5) and F(2315,12288) < F(1,5)
    assert F(26,8)+F(242,256*9) < 4
    assert F(19,8)+F(211,256*9) < 3
    assert F(7,5)+F(1,25*81) < 2
    assert F(1,16)-F(1,256)-F(1,12288) > F(1,32)
    assert 9*(512/(a*a)+1+9*mu*mu) < 2173**2
    # e^(m^2 S/4) <= 1/(1-1/16) < 2 supplies the analytic tail bound.
    assert 1/(1-F(1,16)) < 2
    c = K*a**6*mu**8/288
    singular = a**3*mu**3/384
    # eta=1/20, log(12/eta)=log(240)<6, so this integer budget is sufficient.
    need = 1728/(c*c)
    repetitions = (need.numerator+need.denominator-1)//need.denominator
    return {'a':str(a),'L':str(L),'mu':str(mu),'s_squared_max':str(S),
            'mass_determinant':det,'K':str(K),
            'tail_margin_K_over_2_minus_2_mu_squared':str(K/2-2*mu*mu),
            'cycle_strict_lower_bound':str(c),'cycle_lower_float':float(c),
            'amplitude_singular_strict_lower_bound':str(singular),
            'relative_field_energy_upper':2173,
            'per_probability_calibration_budget':str(c/24),
            'per_probability_statistical_budget':str(c/24),
            'repetitions_per_known_input_probability_eta_0_05':str(repetitions),
            'all_time_window_proven_analytically':True,
            'practical_precision_or_autonomous_device_claimed':False}


def gauss(n,lo,hi):
    x,w=np.polynomial.legendre.leggauss(n)
    return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2


def bump(r,a):
    r=np.asarray(r)
    z=np.maximum(1-(r/a)**2,1e-300)
    return np.where(r<a,np.exp(-1/z),0.)


def continuum_moments(order,t):
    a,L,p=1/32,1/4,1.
    radial,wr=gauss(order,0,a)
    b=bump(radial,a)
    norm=math.sqrt(4*math.pi*np.sum(wr*radial**2*b*b))
    g=b/norm
    B0=4*math.pi*np.sum(wr*radial**2*g)
    B2=4*math.pi/3*np.sum(wr*radial**4*g)
    Z=4*math.pi*np.sum(wr*radial**2*g*np.sinc(p*radial/math.pi))
    # Stable convergent series for d/dp sinc(pr), |pr| <= 1/32.
    dp=np.zeros_like(radial)
    for j in range(1,8):
        dp += (-1)**j*(2*j)*p**(2*j-1)*radial**(2*j)/math.factorial(2*j+1)
    Zprime=4*math.pi*np.sum(wr*radial**2*g*dp)
    G,Gprime=Z*Z,2*Z*Zprime
    # The radial autocorrelation is an actual continuum convolution, not a lattice.
    rr,ww=gauss(order,0,2*a)
    nodes,weights=np.polynomial.legendre.leggauss(order)
    rho=[]
    for r in rr:
        lo=np.abs(radial-r)
        hi=np.minimum(a,radial+r)
        width=np.maximum(hi-lo,0)
        vv=lo[:,None]+width[:,None]*(nodes[None,:]+1)/2
        inner=np.sum(vv*bump(vv,a)/norm*weights[None,:],axis=1)*width/2
        rho.append(2*math.pi/r*np.sum(wr*radial*g*inner))
    rho=np.array(rho)
    u,wu=gauss(order,-1,1)
    ss=t*t-(L*L+rr[:,None]**2+2*L*rr[:,None]*u[None,:])
    phase=np.exp(-1j*p*rr[:,None]*u[None,:])
    measure=2*math.pi*(ww*rr**2*rho)[:,None]*wu[None,:]
    raw=[complex(np.sum(measure*phase*ss**ell)) for ell in range(5)]
    coeff=[(-1)**ell*raw[ell]/(2**(2*ell+3)*math.pi*math.factorial(ell)*math.factorial(ell+1)) for ell in range(5)]
    lhs=(coeff[0]*coeff[1].conjugate()).imag
    rhs=-L*G*Gprime/(256*math.pi**2)
    assert lhs > 0 and rhs > 0
    return {'order':order,'t':t,'B0':float(B0),'B2':float(B2),'W':float(B0*B0),
            'G':float(G),'Gprime':float(Gprime),'raw':raw,'coeff':coeff,
            'phase_area':float(lhs),'phase_area_fourier':float(rhs),
            'moment_zero_relative_error':float(abs(raw[0]-G)/G),
            'phase_area_relative_error':float(abs(lhs-rhs)/rhs)}


def dc(z):
    if isinstance(z,tuple): return z
    if isinstance(z,complex): return Decimal(str(z.real)),Decimal(str(z.imag))
    return Decimal(str(z)),Decimal(0)

def add(a,b): return a[0]+b[0],a[1]+b[1]
def neg(a): return -a[0],-a[1]
def sub(a,b): return add(a,neg(b))
def mul(a,b): return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def conj(a): return a[0],-a[1]
def scale(a,b): return a[0]*b,a[1]*b
def norm2(a): return a[0]*a[0]+a[1]*a[1]


def matrix_amplitude(u,fs):
    n=len(u)
    return [[sumc([mul(mul(u[b][j],fs[j]),conj(u[a][j])) for j in range(n)]) for a in range(n)] for b in range(n)]

def sumc(xs):
    out=dc(0)
    for x in xs: out=add(out,x)
    return out


def cycle(a):
    probs=[[norm2(z) for z in row] for row in a]
    return probs[1][0]+probs[2][1]+probs[0][2]-probs[1][2]-probs[2][0]-probs[0][1]


def flavor_diagnostics(mom):
    with localcontext() as ctx:
        ctx.prec=90
        mu=Decimal(1)/Decimal(4096)
        root3=Decimal(3).sqrt()
        w=(Decimal(-1)/2,root3/2)
        wp=[dc(1),w,conj(w)]
        u=[[scale(wp[(a*j)%3],1/root3) for j in range(3)] for a in range(3)]
        cs=[dc(z) for z in mom['coeff']]
        fs=[sumc([scale(c,(mu*j)**(2*ell+3)) for ell,c in enumerate(cs)]) for j in [1,2,3]]
        amp=matrix_amplitude(u,fs)
        C=cycle(amp)
        J=mul(mul(mul(u[0][0],u[1][1]),conj(u[0][1])),conj(u[1][0]))[1]
        area=mul(sub(fs[0],fs[2]),conj(sub(fs[1],fs[2])))[1]
        assert C>0 and area>0 and abs(C-12*J*area)<abs(C)*Decimal('1e-55')
        # Model-independent reference baselines cancel without a physical empty code.
        baselines=[Decimal(1)/4,Decimal(1)/3,Decimal(2)/5]
        pp=[[baselines[b]+norm2(amp[b][a]) for a in range(3)] for b in range(3)]
        six=(pp[1][0]-pp[1][2])+(pp[2][1]-pp[2][0])+(pp[0][2]-pp[0][1])
        assert abs(six-C)<Decimal('1e-86')
        conjC=cycle(matrix_amplitude([[conj(z) for z in row] for row in u],fs))
        assert abs(conjC+C)<abs(C)*Decimal('1e-55')
        o=[[dc(Decimal(int(i==j))-Decimal([1,2,3][i]*[1,2,3][j])/7) for j in range(3)] for i in range(3)]
        realC=cycle(matrix_amplitude(o,fs))
        no_carrier=cycle(matrix_amplitude(u,[(f[0],Decimal(0)) for f in fs]))
        degenerate=cycle(matrix_amplitude(u,[fs[0],fs[0],fs[2]]))
        for z in [realC,no_carrier,degenerate]: assert abs(z)<Decimal('1e-82')
        phase=Decimal('0.7')
        # A genuinely complex U(2) diagnostic; no plane-wave propagator is substituted.
        phasec=dc(complex(math.cos(float(phase)),math.sin(float(phase))))
        s=Decimal(1)/Decimal(2).sqrt()
        u2=[[dc(s),scale(phasec,s)],[scale(conj(phasec),-s),dc(s)]]
        a2=matrix_amplitude(u2,fs[:2])
        assert abs(norm2(a2[0][1])-norm2(a2[1][0]))<Decimal('1e-82')
        # Fixed mass/electric-role diagnostics only: continuum propagation is above.
        uf=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
        hf=uf@np.diag([1.,4.,9.])@uf.conj().T
        of=np.eye(3)-np.outer([1.,2.,3.],[1.,2.,3.])/7
        hr=of@np.diag([1.,4.,9.])@of.T
        assert min(abs(hf[i,j]) for i in range(3) for j in range(i))>0.1
        assert min(abs(hr[i,j]) for i in range(3) for j in range(i))>0.1
        lower=F(rational_certificate()['cycle_strict_lower_bound'])
        assert C>Decimal(lower.numerator)/Decimal(lower.denominator)
        # The diagnostic coefficient series uses five terms; proof tail bounds all terms.
        return {'C':str(C),'J':str(J),'spectral_triangle':str(area),
                'six_probability_baseline_cancellation':True,
                'synthetic_baselines_are_not_computed_vacuum_expectations':True,
                'conjugate_C':str(conjC),'real_mixing_zero':str(realC),
                'common_phase_zero':str(no_carrier),'two_equal_masses_zero':str(degenerate),
                'two_flavor_normal_reciprocity':True,
                'minimum_amplitude_eigenvalue_abs':str(min(norm2(z).sqrt() for z in fs)),
                'complex_and_real_common_flavor_irreducibility':True,
                'synthetic_baseline_probabilities_decimal':[[str(z) for z in row] for row in pp]}


def packed(mom):
    return {k:([[z.real,z.imag] for z in v] if k in ['raw','coeff'] else v) for k,v in mom.items()}


def payload():
    exact=rational_certificate()
    # Three fixed times diagnose the proven uniform window; two resolutions check quadrature.
    times=[7/16,15/32,1/2]
    coarse=[continuum_moments(40,t) for t in times]
    fine=[continuum_moments(64,t) for t in times]
    for low,high in zip(coarse,fine):
        assert high['moment_zero_relative_error']<2e-6
        assert high['phase_area_relative_error']<2e-5
        assert abs(high['phase_area']-low['phase_area'])/high['phase_area']<3e-5
    diagnostics=[flavor_diagnostics(m) for m in fine]
    sources=['proof.md','NEXT.md','../research_note_1303.md','../_shared/round1303_admission.md',
             '../1293/proof.md','../1295/proof.md','../1300/proof.md',
             '../_shared/nondegenerate_joint_corollary_after1300.md',
             '../../archive_1009_1043/research_note_1011.md']
    return {'status':'PASS_AUTHOR_COMPACT_FLAVOR_PHASE_INTERFACE',
            'exact_uniform_certificate':exact,'continuum_quadrature_coarse':[packed(m) for m in coarse],
            'continuum_quadrature_fine':[packed(m) for m in fine],
            'finite_flavor_diagnostics':diagnostics,
            'matrix_series_terms':5,'quadrature_is_not_uniform_proof':True,
            'new_axioms':0,'exactly_three_generations_claimed':False,
            'physical_CP_or_T_instrument_claimed':False,'full_quantum_SM_claimed':False,
            'old_1295_or_1300_precision_inherited':False,
            'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
            'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('CREATED',result['status'])
    else:
        expected=json.loads(OUT.read_text(encoding='utf-8'))
        if result!=expected: raise SystemExit('MISMATCH: diagnostics or bindings changed')
        print('READ_ONLY_MATCH',result['status'])
    print('cycle_lower',result['exact_uniform_certificate']['cycle_lower_float'])
    print('cycles',[x['C'] for x in result['finite_flavor_diagnostics']])

if __name__=='__main__': main()
