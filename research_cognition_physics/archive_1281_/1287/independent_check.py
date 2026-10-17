"""1287 independent exact checks; no import of author code.
Combined Bessel series, integer-complex Clifford algebra, and independent
first-order dual arithmetic for the stress divergence. --write is exclusive;
default recomputes and compares without modifying any file.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
SIG = (-1, 1, 1, 1)
M2 = ((Q(5,2), Q(-3,2)), (Q(-3,2), Q(5,2)))


def rat(x):
    x = Q(x)
    return f'{x.numerator}/{x.denominator}'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def mul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def matadd(a,b):
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def scale(a,z):
    return [[z*x for x in row] for row in a]


def dagger(a):
    return [[complex(a[j][i]).conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def dirac_check():
    # All matrix entries and unnormalized spinors below are Gaussian integers.
    # Their arithmetic is exact in binary complex arithmetic at these sizes.
    pauli = ([[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]])
    gamma0 = [[int(i==j)*(1 if i<2 else -1) for j in range(4)] for i in range(4)]
    gammas = [gamma0]
    for s in pauli:
        gammas.append([[s[i][j-2] if i<2 and j>=2 else
                         -s[i-2][j] if i>=2 and j<2 else 0 for j in range(4)] for i in range(4)])
    identity = [[int(i==j) for j in range(4)] for i in range(4)]
    for a in range(4):
        for b in range(4):
            ant = matadd(mul(gammas[a],gammas[b]),mul(gammas[b],gammas[a]))
            assert ant == scale(identity,-2*SIG[a] if a==b else 0)
    alphas = [mul(gamma0,g) for g in gammas[1:]]
    assert all(dagger(g)==g for g in [gamma0]+alphas)
    hd = matadd(scale(alphas[2],4),scale(gamma0,3))
    u = [[2],[0],[1],[0]]
    assert mul(hd,u)==scale(u,5)
    norm = mul(dagger(u),u)[0][0]
    assert norm == 5
    bar = mul(dagger(u),gamma0)
    current = [mul(mul(bar,scale(gammas[a],SIG[a])),u)[0][0] for a in range(4)]
    derivatives = [-5j,0,0,4j]
    raw = [[(-1j/2)*(current[a]*derivatives[b]+current[b]*derivatives[a])
            for b in range(4)] for a in range(4)]
    expected = [[25,0,0,-20],[0,0,0,0],[0,0,0,0],[-20,0,0,16]]
    assert raw == expected
    p_lower = [-5,0,0,4]
    assert raw == [[p_lower[a]*p_lower[b] for b in range(4)] for a in range(4)]

    return {'clifford_pairs_exact':16, 'mass':'3/1', 'momentum_z':'4/1',
            'positive_energy':'5/1', 'T00':'5/1', 'T03':'-4/1', 'T33':'16/5',
            'plane_wave_stress_equals_momentum_outer_product_over_energy':True, 'floating_quadrature_used':False}


class Dual:
    """Exact first-order jet in four coordinates."""
    def __init__(self,v,d=None):
        self.v = Q(v)
        self.d = tuple(Q(x) for x in (d if d is not None else [0]*4))
    def __add__(self,o):
        if not isinstance(o,Dual): o=Dual(o)
        return Dual(self.v+o.v,[a+b for a,b in zip(self.d,o.d)])
    __radd__=__add__
    def __mul__(self,o):
        if not isinstance(o,Dual): o=Dual(o)
        return Dual(self.v*o.v,[self.v*b+o.v*a for a,b in zip(self.d,o.d)])
    __rmul__=__mul__
    def __neg__(self): return self*(-1)
    def __sub__(self,o): return self+(-o if isinstance(o,Dual) else -Q(o))


def tensor_jets(sign,mass):
    values = [Q(1)*sign,Q(2)*sign]
    first = [[1,2,0,-1],[-1,1,3,2]]
    first = [[sign*Q(v) for v in row] for row in first]
    second = [[[Q(0) for b in range(4)] for a in range(4)] for i in range(2)]
    # A local polynomial solution jet, with spatial Hessian zero, for full M2.
    for i in range(2):
        second[i][0][0] = -sum(M2[i][j]*values[j] for j in range(2))
    field = [Dual(values[i],first[i]) for i in range(2)]
    gradient = [[Dual(first[i][a],second[i][a]) for a in range(4)] for i in range(2)]
    lag_bracket = sum(SIG[a]*gradient[i][a]*gradient[i][a] for i in range(2) for a in range(4))
    lag_bracket += sum(mass[i][j]*field[i]*field[j] for i in range(2) for j in range(2))
    return [[sum(gradient[i][a]*gradient[i][b] for i in range(2))
             -(SIG[a]*Q(1,2)*lag_bracket if a==b else Dual(0))
             for b in range(4)] for a in range(4)]


def stress_check():
    full = tensor_jets(1,M2)
    neg = tensor_jets(-1,M2)
    divergence = [sum(SIG[a]*full[a][b].d[a] for a in range(4)) for b in range(4)]
    assert divergence == [0]*4
    assert all((full[a][b].v,full[a][b].d)==(neg[a][b].v,neg[a][b].d)
               for a in range(4) for b in range(4))
    omitted = tensor_jets(1,((Q(5,2),0),(0,Q(5,2))))
    defect = [sum(SIG[a]*omitted[a][b].d[a] for a in range(4)) for b in range(4)]
    assert defect == [Q(-3,2),Q(-15,2),Q(-9,2),Q(0)]
    assert all(abs(full[a][b].v)<=full[0][0].v for a in range(4) for b in range(4))
    return {'method':'first-order dual arithmetic of Hilbert tensor on one exact local KG polynomial jet',
            'full_divergence':[rat(x) for x in divergence],
            'omit_mass_mixing_defect':[rat(x) for x in defect],
            'opposite_coherent_sign_tensor_and_derivative_identical':True,
            'continuous_Ward_identity_requires_analytic_proof':True}


def series_sum(x,n):
    # Combined kernel, not the author's separate mass inequalities.
    return sum(Q((-1)**j*(4**(j+1)-1),2**(2*j+1)*factorial(j)*factorial(j+1))*x**j
               for j in range(n+1))


def calculate():
    eigenvectors = ([1,1],[1,-1])
    for vector,eigenvalue in zip(eigenvectors,(1,4)):
        assert [sum(M2[i][j]*vector[j] for j in range(2)) for i in range(2)] == [eigenvalue*x for x in vector]
    comparison = [[M2[i][j]-int(i==j) for j in range(2)] for i in range(2)]
    assert [sum(row) for row in comparison] == [0,0]
    a,d,tmin,tmax,zeta = Q(1,32),Q(1,4),Q(7,16),Q(1,2),Q(1,256)
    rlo,rhi = d-2*a,d+2*a
    xlo,xhi = tmin**2-rhi**2,tmax**2-rlo**2
    assert (rlo,rhi,xlo,xhi)==(Q(3,16),Q(5,16),Q(3,32),Q(55,256))
    assert rhi<tmin and xhi<Q(1,4)
    # coefficient ratio = [(4^(n+2)-1)/(4*(4^(n+1)-1))]/((n+1)(n+2));
    # 4*(4^(n+2)-1) <= 5*4*(4^(n+1)-1) for n>=0.
    # Thus ratio* x <= (5/4)*x/2 <= 5/32 on x<=1/4.
    ratio_bound=Q(5,32)
    lower=Q(3,2)-Q(15,16)*xhi
    assert lower==Q(5319,4096)>Q(5,4)
    brackets=[]
    for x in (xlo,Q(1,8),xhi):
        lo,hi=series_sum(x,3),series_sum(x,4)
        assert lower<lo<hi<Q(3,2)
        brackets.append({'s_squared':rat(x),'lower':rat(lo),'upper':rat(hi)})
    # Alternative bump route: ||b/Z||_2^2 <= sup(b)/Z
    # <= 6 exp(1/3)/(pi*a^3) < 4/a^3 because exp(1/3)<2, pi>3.
    exp_third_upper=1/(1-Q(1,3))
    assert exp_third_upper<Q(2)
    l2=4/a**3
    variance=Q(3,8)*l2
    exponent=zeta*zeta*variance/2
    assert (l2,variance,exponent)==(131072,49152,Q(3,8))
    damping_lower=1-exponent
    assert damping_lower==Q(5,8)>Q(1,2)
    phase_upper=zeta/(4*3)
    assert phase_upper<1 and 1-phase_upper**2/6>Q(1,2)
    gap=Q(1,2)*Q(1,2)*zeta*Q(5,128)
    assert gap==Q(5,131072)
    error=Q(1,10**6)
    robust=gap-2*error
    assert robust==Q(74029,2048000000)>0
    energy_kg=l2/2
    energy_tag=Q(3)+Q(3,2)/6
    assert energy_kg+energy_tag==Q(262157,4)
    assert Q(21,2)<16
    author=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    verified=[]
    for relative,digest in author['source_bindings'].items():
        assert sha(ROOT/relative)==digest,relative
        verified.append(relative)
    assert Q(author['record_certificate']['whole_window_strict_lower'])==gap
    assert author['scope']['initial_all_fields_in_B_blind'] is False
    assert author['scope']['actual_terminal_instrument_constructed'] is False
    files=[HERE/'proof.md',HERE/'check.py',HERE/'results.json',HERE/'NEXT.md',
           HERE.parent/'research_note_1287.md',
           HERE.parent/'_shared/source_scope_and_round1287_admission.md',
           ROOT/'research_cognition_physics/archive_935_955/research_note_942.md',
           Path(__file__).resolve()]
    return {'round':1287,'status':'PASS_WITH_SCOPE','independent_groups':1,'new_axioms':0,
            'method':'Independent combined Bessel series, exact finite Clifford/Dirac plane wave, Fraction dual Hilbert-tensor differentiation; no author import',
            'source_bindings':{p.relative_to(ROOT).as_posix():sha(p) for p in files},
            'verified_author_source_bindings':verified,
            'mass_and_comparison':{'eigenvalues':['1/1','4/1'],'common_comparison_kernel':True},
            'causal_bessel':{'distance':[rat(rlo),rat(rhi)],'time':[rat(tmin),rat(tmax)],
              's_squared':[rat(xlo),rat(xhi)],'global_term_ratio_upper':rat(ratio_bound),
              'combined_numerator_uniform_lower':rat(lower),'three_brackets':brackets},
            'normalization_noise_resources':{'bump_L2_squared_strict_upper':rat(l2),
              'variance_strict_upper':rat(variance),'exponent_strict_upper':rat(exponent),
              'damping_lower':rat(damping_lower),'tag_energy_upper':rat(energy_tag),
              'total_energy_strict_upper':rat(energy_kg+energy_tag),'tag_derivative_norm_squared':'21/2',
              'source_smearing_constant_strict_upper':'65540/1'},
            'signal':{'whole_window_strict_gap':rat(gap),'gap_decimal_diagnostic':float(gap),
              'optional_each_probability_error':rat(error),'optional_gap':rat(robust)},
            'dirac_signature':dirac_check(),'local_source':stress_check(),
            'analytic_scope':{'unknown_Q_R_isometry':True,'initial_second_KG_and_R_product':True,
              'all_fields_in_B_initially_blind':False,'source_compression_scalar_for_unread_mean_only':True,
              'terminal_effect_is_actual_instrument':False,'self_gravitating_flat_solution':False,
              'deletion_is_same_state_quench':False,'finite_samples_prove_continuous_field_claims':False}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    payload=calculate()
    target=HERE/'independent_results.json'
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(payload,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        print('PASS_WITH_SCOPE: EXCLUSIVE_FIRST_WRITE')
    else:
        assert payload==json.loads(target.read_text(encoding='utf-8')),'Saved independent payload differs'
        print('PASS_WITH_SCOPE: READ_ONLY_MATCH')
    print('whole_window_gap > '+payload['signal']['whole_window_strict_gap'])
    print('combined_Bessel_numerator >= '+payload['causal_bessel']['combined_numerator_uniform_lower'])
    print('results_sha256 = '+sha(target))


if __name__=='__main__':
    main()
