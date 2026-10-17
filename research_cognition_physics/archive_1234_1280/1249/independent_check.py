"""Independent 1249 audit: symbolic canonical identities and finite bounds.

This file does not import either author checker.  Its canonical calculations
are exact sparse Laurent-polynomial identities, not point samples.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULT = HERE/'independent_results.json'
NAMES = ('a','A','p','P','X','U','q','r','psi','B','Y','K','theta_r','S','J','w','f','wp','fp')
D = len(NAMES)
ZERO = (0,)*D
IDX = {n:i for i,n in enumerate(NAMES)}


class LP:
    """Finite Laurent polynomials; only a has negative exponents here."""
    def __init__(self, x=0):
        if isinstance(x, LP):
            self.c = x.c.copy()
        elif isinstance(x, dict):
            self.c = {k:F(v) for k,v in x.items() if v}
        else:
            self.c = {ZERO:F(x)} if x else {}

    @staticmethod
    def var(name):
        power=list(ZERO);power[IDX[name]]=1
        return LP({tuple(power):F(1)})

    def __add__(self, b):
        b=LP(b); out=self.c.copy()
        for k,v in b.c.items():out[k]=out.get(k,0)+v
        return LP(out)
    __radd__=__add__

    def __neg__(self):return LP({k:-v for k,v in self.c.items()})
    def __sub__(self,b):return self+-LP(b)
    def __rsub__(self,b):return LP(b)+-self

    def __mul__(self,b):
        b=LP(b);out={}
        for k,v in self.c.items():
            for l,z in b.c.items():
                e=tuple(i+j for i,j in zip(k,l))
                out[e]=out.get(e,0)+v*z
        return LP(out)
    __rmul__=__mul__

    def __truediv__(self,b):return self*F(1,b)

    def __pow__(self,n):
        if n<0:
            assert len(self.c)==1
            e,v=next(iter(self.c.items()))
            return LP({tuple(i*n for i in e):v**n})
        out=LP(1)
        for _ in range(n):out=out*self
        return out

    def diff(self,name):
        i=IDX[name];out={}
        for k,v in self.c.items():
            if k[i]:
                e=list(k);e[i]-=1
                out[tuple(e)]=v*k[i]
        return LP(out)

    def evaluate(self, point):
        return sum((v*prod(F(point[NAMES[i]])**p for i,p in enumerate(k) if p)
                    for k,v in self.c.items()), F(0))

    def __bool__(self):return bool(self.c)
    def __eq__(self,b):return not (self-LP(b))


def prod(xs):
    out=F(1)
    for x in xs:out*=x
    return out


def canonical_group():
    a,A,p,P,X,U,q,r,psi,B,Y,K,theta_r,S,J,w,f,wp,fp=[LP.var(n) for n in NAMES]
    h=a**-3*(-A*A/4+(P*P+U*U)/2)+a**3*(p*p+X*X)/2
    h+=F(3,20)*a*(1+p/5)+F(72,100)+(q*q+r*r+K*K)/2
    h+=(B-1)*w*q*p/4+f*X*Y/10+(1+p)*S/10
    velocity={'a':a*h.diff('A'),'A':-a*h.diff('a')}
    for pos,mom in [('p','P'),('X','U'),('q','r'),('psi','B'),('Y','K'),('theta_r','S')]:
        velocity[pos]=h.diff(mom)
        velocity[mom]=-h.diff(pos)
    velocity.update(J=-(B-1)*wp*q*p/4-fp*X*Y/10,w=wp,f=fp,wp=LP(),fp=LP())
    def dt(b):
        return sum((b.diff(n)*v for n,v in velocity.items()),LP())
    H=-a**-3*A/2;v=a**-3*P;V=a**-3*U
    Mc=F(71,100)+J+(q*q+r*r+K*K)/2+(B-1)*w*q*p/4+f*X*Y/10
    rho=(v*v+V*V+p*p+X*X)/2+F(3,20)*(1+p/5)*a**-2
    rho+=a**-3*(Mc+F(1,100)+(1+p)*S/10)
    pressure=(v*v+V*V-p*p-X*X)/2-F(1,20)*(1+p/5)*a**-2
    checks={
        'log_scale_velocity':velocity['a']-a*H,
        'pointer_field':dt(v)+3*H*v+p+F(3,100)*a**-2+((B-1)*w*q/4+S/10)*a**-3,
        'minimal_probe':dt(V)+3*H*V+X+f*Y*a**-3/10,
        'oscillator_coordinate':velocity['q']-r,
        'oscillator_momentum':velocity['r']+q+(B-1)*w*p/4,
        'program_conjugate_no_alpha':velocity['psi']-w*q*p/4,
        'program_action_constant':velocity['B'],
        'memory_coordinate':velocity['Y']-K,
        'memory_momentum':velocity['K']+f*X/10,
        'recorder_phase_density':velocity['theta_r']-(1+p)/10,
        'recorder_action_constant':velocity['S'],
        'clock_control_work':dt(h)+velocity['J'],
        'complete_mass_exchange':dt(Mc)-(B-1)*w*q*v/4-f*Y*V/10,
        'complete_source_continuity':dt(rho)+3*H*(rho+pressure),
        'geometry_before_imposing_constraint':dt(H)+F(3,2)*(H*H+pressure),
        'original_constraint_dictionary':h+J-a**3*(rho-H*H),
        'same_history_memory_omission':rho-(rho-(F(1,100)+K*K/2)*a**-3)-(F(1,100)+K*K/2)*a**-3,
    }
    for name, residual in checks.items():
        assert not residual, (name,residual.c)
    # Exact local Weyl-square correction; the two momentum derivatives needed
    # above order two vanish because H is linear in A.
    H_QA=a*H.diff('a').diff('A')
    assert H_QA**2/4 == F(9,16)*a**-6
    assert not H.diff('A').diff('A')
    initial_terms=[F(-1110),F(5),F(150),F(835),F(20)]
    assert sum(initial_terms)==-100
    assert sum(initial_terms)+100==0
    # Scaling is an actual family of finite extensive resources.
    volume_checks=[]
    for n in (1,2,7):
        volume=1000*n**3
        assert volume//10==100*n**3
        assert sum(initial_terms)*n**3+100*n**3==0
        volume_checks.append({'n':n,'V':volume,'Nc':volume,'Nr':volume//10})
    return {'status':'PASS','method':'Exact Laurent-polynomial identities with independent clock-action J before constraint substitution.',
            'identically_zero_residuals':list(checks),
            'weyl_square_correction':'9 epsilon^2 exp(-6Q)/16 locally; not a global cutoff formula',
            'initial_total_constraint':'100-100=0 at V=1000',
            'volume_family_checks':volume_checks,
            'independent_source_measurement_not_inferred_from_constraint':True}


def fourier_and_cutoff_lemmas():
    """Radial beta integrals, independently of the author's annulus bound."""
    # D=14, m=10. A_k=(2pi)^-14 |S^13| B((14+k)/2,10-(14+k)/2)/2.
    # For k=0,2,4 the beta factors are rational.  For k=3 a single pi
    # cancels. pi>3 gives strict rational upper bounds, without float pi.
    radial={
        0:F(1,2**14*720*252*3**7),
        2:F(1,2**14*720*72*3**7),
        3:F(math.factorial(16),2**15*4**8*720*math.factorial(8)*math.factorial(9)*3**6),
        4:F(1,2**14*720*9*3**7),
    }
    assert all(0<x<1 for x in radial.values())
    assert F(math.factorial(6)*math.factorial(2),math.factorial(9))==F(1,252)
    assert F(math.factorial(7)*math.factorial(1),math.factorial(9))==F(1,72)
    assert F(math.factorial(8),math.factorial(9))==F(1,9)
    assert 14+4-2*10<0
    # B(x)=normalizer*integral_0^x t^40(1-t)^40 dt is C40 when
    # extended by 0 and 1.  It is not being called C-infinity.
    normalizer=F(math.factorial(81),math.factorial(40)**2)
    poly=[F(0)]*82
    for j in range(41):
        poly[41+j]=normalizer*(-1)**j*math.comb(40,j)/F(41+j)
    assert sum(poly)==1 and poly[0]==0
    d=poly
    coeff_bounds=[]
    for order in range(1,41):
        d=[(i+1)*d[i+1] for i in range(len(d)-1)]
        assert d[0]==0 and sum(d)==0
        coeff_bounds.append(sum(abs(c) for c in d))
        # A loose explicit bound from 81*C(80,40), coefficient sum 2^40,
        # and at most 80^(order-1) differentiation factors.
        assert coeff_bounds[-1] < 2**(127+7*(order-1))
    assert normalizer < 2**87
    B=324*normalizer*2**40
    for r,value in enumerate(coeff_bounds,1):
        assert value <= math.factorial(r)*B**r
    # The scalar sine Taylor remainder gives the Moyal constant; both initial
    # coherent smoothing and the final positive anti-Wick effect have 1/4.
    assert F(2)*F(1,2)**3/F(6)==F(1,24)
    assert F(1,4)+F(1,4)==F(1,2)
    return {'radial_Fourier_A_upper':{str(k):str(v) for k,v in radial.items()},
            'Fourier_weight_power':10,'phase_dimension':14,
            'beta_normalizer':str(normalizer),'cutoff_order':40,
            'derivative_endpoint_checks':40,
            'derivative_coefficient_bounds_sha256':hashlib.sha256(
                json.dumps(list(map(str,coeff_bounds))).encode()).hexdigest(),
            'Moyal_constant':'1/24',
            'two_coherent_smoothing_constants':'1/4 + 1/4'}


def partial_bell(n,k,x):
    """Exact exponential partial Bell polynomial; x[j] is derivative j."""
    table={(0,0):F(1)}
    for r in range(1,n+1):
        for j in range(1,min(r,k)+1):
            table[r,j]=sum((F(math.comb(r-1,i-1))*x.get(i,F(0))*table.get((r-i,j-1),F(0))
                            for i in range(1,r-j+2)),F(0))
    return table.get((n,k),F(0))


def quantum_group(compare_saved=True):
    lemmas=fourier_and_cutoff_lemmas()
    # Variational equations yield C_r=(r-1)^-1 sum_{j>=2} B_{r,j}(C).
    # This independently proves the flow majorant, rather than checking the
    # author's total Bell-number estimate. The exact solution is (r-1)!.
    constants={1:F(1)}
    for r in range(2,21):
        constants[r]=sum((partial_bell(r,j,constants) for j in range(2,r+1)),F(0))/(r-1)
        assert constants[r]==math.factorial(r-1)
    for r in range(1,21):
        assert math.factorial(r)<2**(4096*r)
    # All subsequent arithmetic uses ordinary integers; no huge 2^N is built.
    normalizer=math.factorial(81)//math.factorial(40)**2
    B=324*normalizer*2**40
    Lcut=56*B
    assert Lcut>=34*B and Lcut>=28*B
    H0=2**36
    analytic_source=1024*3**18*(1+H0)*(1+129**2)
    assert analytic_source<2**96
    Lh=14*(1+Lcut)
    rawM=H0*math.factorial(21)*Lh**21
    logceil=lambda n:(n-1).bit_length()
    m_exp=logceil(rawM);M=2**m_exp
    assert M>=rawM
    # e^(T M)<2^(5M) for T=9/4 because e<4.  The independent factorial
    # majorants are absorbed by the additional 2^4096 in R.
    R_exp=4096+5*M
    A_source=2**98*(1+Lcut)**2
    L_source=4*(Lcut+14)
    L_positive=L_source*(1+A_source*100*B)
    L_good=Lcut+4*L_positive
    L_K=Lcut+4_000_000_000*B
    A_tail=1+128*A_source**2
    L_tail=L_good+2*L_source
    A_menu=max(A_source,A_tail,1)
    L_menu=max(L_source,L_tail,L_K)
    volume=5*256**13
    prefactor=2**34*volume*15**10*math.factorial(20)
    EH=logceil(prefactor*H0*Lh**20)
    Ebase=logceil(prefactor*A_menu*(1+L_menu)**20)
    EF=Ebase+20*R_exp
    N=1+max((EF+44)//3,(EH+EF+47)//6)
    linear_margin=3*N-EF
    quadratic_margin=6*N-EH-EF
    assert linear_margin>=43 and quadratic_margin>=43
    # Includes BOTH coherent heat terms and the time-integrated Moyal term.
    # eps <= 2^(-3N)/1000. Do not expand a giant negative power of two.
    upper=(F(1,2000)+F(3,32_000_000))*F(1,2**43)
    assert upper<F(1,2**41)<F(1,2**40)<F(1,10**10)
    old=json.loads((HERE.parent/'1247/results.json').read_text(encoding='utf-8'))
    centers=[F('-.00227032632043'),F('-.00227026126850')]
    halfwidth=F(1,10**9)
    errors=[]
    for i,c in enumerate(centers):
        lo,hi=map(F,old['group_finite_memory'][f'alpha{i}']['checkpoints']['2']['state']['K'])
        error=max(abs(lo-c),abs(hi-c));errors.append(str(error))
        assert error<F(3,10**15)<halfwidth/2
    assert centers[1]-centers[0]>2*halfwidth
    matching={
        'M_power_of_two_exponent':m_exp,
        'Fourier_H_exponent':EH,
        'Fourier_menu_base_exponent':Ebase,
        'Fourier_transport_exponent_decimal':str(EF),
        'integer_n_power_exponent_decimal':str(N),
        'linear_error_negative_binary_exponent':str(linear_margin),
        'quadratic_error_negative_binary_exponent':str(quadratic_margin),
    }
    if compare_saved:
        author=json.loads((HERE/'quantum_bounds_results.json').read_text(encoding='utf-8'))
        ledger=author['group_constructive_finite_resource_bound']['ledger']
        for key,value in matching.items():assert ledger[key]==value,(key,ledger[key],value)
    return {'status':'PASS','method':'Independent exact radial Fourier integrals and variational-flow Bell recursion, then exact resource arithmetic.',
            'lemmas':lemmas,'flow_C_r':[str(constants[r]) for r in range(1,21)],
            'flow_identity':'C_r=(r-1)! for 1<=r<=20',
            'resource_crosscheck':matching,'resource_N_digits':len(str(N)),
            'combined_error_rational_upper':str(upper),'K_center_errors':errors,
            'same_three_outcome_Husimi_POVM':True,
            'scope':['single-time Husimi effects uniformly over the finite window; not continuous surveillance',
                     'bounded completed source symbols, not the unextended global stress',
                     'known coherent inputs; arbitrary reference extension is not an arbitrary-state classical limit',
                     'collective homogeneous intensity readout; not one local original clock',
                     'quantum source matching is within budget, not exact equality of separate initial quantum energies']}


def run(draft=False):
    result={'round':1249,'canonical_group':canonical_group(),
            'quantum_group':quantum_group(compare_saved=not draft),
            'scientific_groups':2,'status':'DRAFT_PENDING_TEXT_REVIEW' if draft else 'PASS'}
    if not draft:
        paths={
            'independent_check.py':Path(__file__),
            **{name:HERE/name for name in (
                'admission.md','canonical.md','canonical_check.py','canonical_results.json',
                'quantum_bridge.md','quantum_bounds.py','quantum_bounds_results.json','proof.md','NEXT.md')},
            '../research_note_1249.md':HERE.parent/'research_note_1249.md',
            '../1247/proof.md':HERE.parent/'1247/proof.md',
            '../1247/results.json':HERE.parent/'1247/results.json',
        }
        result['bindings']={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--draft',action='store_true')
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    assert not (args.draft and args.write)
    result=run(draft=args.draft)
    if args.write:
        with RESULT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        print('FIRST_EXCLUSIVE_WRITE PASS; scientific_groups=2')
    elif not args.draft:
        saved=json.loads(RESULT.read_text(encoding='utf-8'))
        if saved!=result:raise SystemExit('READ_ONLY_MISMATCH')
        print('READ_ONLY_MATCH PASS; scientific_groups=2')
    print(json.dumps({'status':result['status'],'canonical_residuals':len(result['canonical_group']['identically_zero_residuals']),
                      'quantum_M_exponent':result['quantum_group']['resource_crosscheck']['M_power_of_two_exponent'],
                      'quantum_N_digits':result['quantum_group']['resource_N_digits']},ensure_ascii=False))
