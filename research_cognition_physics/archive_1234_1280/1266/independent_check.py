"""Independent 1266 arithmetic and finite algebra; no author code imports.
Default is read-only. --write exclusively creates independent_results.json.
The continuous Gaussian propagation and source claims require the reviewed proof.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]

def tr(a):
    return [list(x) for x in zip(*a)]

def diagonal(v):
    return [[F(v[i]) if i==j else F(0) for j in range(len(v))] for i in range(len(v))]

def calculate():
    # Reconstruct the old ledger from its primitive factorial/Fourier bounds,
    # then compare with its saved certificate. Do not import either author script.
    ceil_log2=lambda x:(x-1).bit_length()
    BS=324*(factorial(81)//factorial(40)**2)*(1<<40)
    LC=56*BS
    LH=14*(1+LC)
    H0=1<<36
    Mexp=ceil_log2(H0*factorial(21)*LH**21)
    M=1<<Mexp
    AS=(1<<98)*(1+LC)**2
    LS=4*(LC+14)
    LP=LS*(1+100*AS*BS)
    LXI=LC+4*LP
    LK=LC+4*10**9*BS
    AE=1+128*AS**2
    LE=LXI+2*LS
    AST=max(AS,AE,1)
    LST=max(LS,LE,LK)
    CF=(1<<34)*5*256**13*15**10*factorial(20)
    EH=ceil_log2(CF*H0*LH**20)
    E0=ceil_log2(CF*AST*(1+LST)**20)
    EF=E0+20*(4096+5*M)
    N=1+max((EF+44)//3,(EH+EF+47)//6)
    ledger=json.loads((STAGE/'1249/quantum_bounds_results.json').read_text(encoding='utf-8-sig'))['group_constructive_finite_resource_bound']['ledger']
    assert (Mexp,EH,E0)==(3059,3093,14008)
    assert EF-100*M==95928
    assert N==int(ledger['integer_n_power_exponent_decimal'])
    assert EF==int(ledger['Fourier_transport_exponent_decimal'])
    assert EH==int(ledger['Fourier_H_exponent'])
    assert E0==int(ledger['Fourier_menu_base_exponent'])

    # Gaussian even moments by the independent integration-by-parts recurrence.
    moments={0:F(1)}
    for j in range(2,7,2): moments[j]=(j-1)*moments[j-2]
    assert [moments[j] for j in (2,4,6)]==[1,3,15]
    cG2=moments[6]/(factorial(3)**2*2**3)
    assert cG2==F(5,96)<F(1,16)
    # Time length 9/4 and two nonorthogonal coefficients have sum <=sqrt(2).
    assert cG2*2*F(9,4)**2<1
    # A diagnostic illustrating why the two comparison vectors need not be
    # orthogonal or normalized in the middle of the gate (not a simulation).
    c,s,overlap=F(3,5),F(4,5),F(3,5)
    middle_norm_sq=c*c+s*s+2*c*s*overlap
    assert middle_norm_sq==F(197,125)>1
    assert middle_norm_sq<=(abs(c)+abs(s))**2<2

    eta=F(1,10**9)
    delta_center=F(3,10**15)
    assert 2*eta-delta_center>eta
    c0=F('-0.00227032632043'); c1=F('-0.00227026126850')
    assert c1-c0>6*eta
    interval=[(F(15,16),F(17,16)),(F(31,16),F(33,16))]
    six=sorted((a+d,b+d) for d in (0,4,8) for a,b in interval)
    assert all(six[i][1]<six[i+1][0] for i in range(5))
    assert max(b for a,b in six)==F(161,16)<61<64
    beta_coefficient=(F(16**2)+eta**-2)/2
    assert beta_coefficient<2**60
    assert 128*(1<<8)==1<<15
    # From beta<epsilon 2^(10M+68), pi sqrt(2 beta)
    # has squared constant 2*pi^2*2^68 <2^74.
    assert 2*F(22,7)**2<2**6
    hres=EH+15*M+12
    gres=5*M+37
    d_factor=EH+15*M+16
    post_factor=EH+15*M+18
    source_factor=E0+EH+15*M+22
    assert gres<hres and hres+1<d_factor
    assert 5*M+9<d_factor and d_factor+1<post_factor
    # Vector-source comparison before reading includes |b_c|, so it does not
    # use a normalized intermediate comparison vector or cancel cross terms.
    gate_source_1=E0+2+d_factor
    gate_source_2=E0+5*M+4
    assert gate_source_2<gate_source_1
    assert gate_source_1+1<source_factor
    # After the cq instrument: 2||b||D, one Gaussian fluctuation, and AW heat.
    post_source_1=E0+2+post_factor
    post_source_2=E0+5*M+4
    heat_source=E0-2
    assert max(post_source_2,heat_source)<post_source_1
    assert post_source_1+2==source_factor

    exponents={
        'state':d_factor-2*N,
        'cq_poststate':post_factor-2*N,
        'bounded_source':source_factor-2*N,
        'blank_bad_rate':10*M+15-4*N,
        'blank_K_bad_rate':10*M+68-4*N,
        'interaction_mass':9-4*N,
    }
    required={'state':-100,'cq_poststate':-100,'bounded_source':-80,
              'blank_bad_rate':-200,'blank_K_bad_rate':-200,'interaction_mass':-80}
    assert all(exponents[k]<v for k,v in required.items())
    assert 2*F(1,2**100)+F(1,2**200)<F(1,2**98)
    assert F(98,1000)-F(1,100)>F(2,25)
    assert (F(1,2**79)/F(2,25))**2==F(625,2**160)<F(1,2**150)

    # Exact finite flag identities. Blank=0, flags=1,2, outside=3.
    Pblank=diagonal([1,0,0,0]); eye=diagonal([1,1,1,1])
    flags=[]
    for j in (1,2):
        U=diagonal([1,1,1,1])
        U[0][0]=U[j][j]=F(0)
        U[j][0]=F(1);U[0][j]=F(-1)
        Q=diagonal([int(k==j) for k in range(4)])
        assert mm(tr(U),U)==eye
        assert mm(mm(tr(U),Q),U)==Pblank
        # Removing blank preparation: incoming same flag returns to blank.
        assert U[j][j]==0 and U[0][j]==-1
        flags.append({'flag':j,'pullback_equals_blank':True,'nonblank_deletion':True})
    # AW single-momentum convolution has variance epsilon/2. Spectral losses
    # use a 1/32 mask margin for B and eta for K, not AW=spectral equality.
    Btail=F(1,2*F(1,32)**2)
    Ktail=F(1,2*eta**2)
    assert Btail==512 and Ktail==5*10**17<2**59
    assert LC+128*BS<LK
    assert 4*N>101
    eq=F(1,2**41)
    assert eq+F(1,2**42)==F(3,2**42)<F(1,2**40)
    assert 2*F(3,2**42)==F(3,2**41)<F(1,2**39)
    assert 4*N-EF>42 and 8*N-EH-EF>42

    ta=F(17,8);tau=F(1,128);tb=ta+tau;tf=tb+tau
    assert (tb,tf)==(F(273,128),F(137,64))
    assert F(2)<ta<tb<tf<F(9,4)
    # Integrate the cubic beta pulse exactly, without evaluating at a grid.
    pulse_integral=140*(F(1,4)-F(3,5)+F(3,6)-F(1,7))
    peak=F(140,64)/tau
    assert pulse_integral==1 and peak==280
    assert peak*F(11,7)==440<512
    assert tf-ta==F(1,64)
    # The numerical upper bound is intensive; physical energy scales V0*s.
    deps=[Path(__file__).resolve(),STAGE/'1249/quantum_bounds_results.json',
          STAGE/'1249/quantum_bridge.md',STAGE/'1249/canonical.md',
          STAGE/'1259/proof.md',
          STAGE/'_admission/finite_semiclassical_budget_gate_after1256.md',
          STAGE/'_admission/causal_block_after1265/proof.md']
    return {
        'round':1266,'status':'PASS_INDEPENDENT_1266_WITH_SCOPE','independent_groups':1,
        'method':'Independent primitive-constant reconstruction, exact Gaussian-moment recurrence, rational matrix identities and exponent inequalities; no author-code imports',
        'old_ledger_reconstructed':{'M_exponent':Mexp,'EH':EH,'E0':E0,'EF_offset_after_100M':EF-100*M,
            'N_digits':len(str(N)),'N_decimal_sha256':hashlib.sha256(str(N).encode()).hexdigest()},
        'continuous_gate_constants':{'sixth_moment':str(moments[6]),'cG_squared':str(cG2),
            'nonorthogonal_diagnostic_norm_squared':str(middle_norm_sq),'beta_coefficient':str(beta_coefficient),
            'B_largest_interval_endpoint':str(six[-1][1]),'pulse_integral':str(pulse_integral),
            'pulse_peak':str(peak),'gate_begin':str(ta),'gate_end':str(tb),'flag_end':str(tf)},
        'negative_exponents':{k:str(v) for k,v in exponents.items()},
        'bounded_source_proof_checks':{'norm_comparison_retains_abs_center':True,
            'nonorthogonal_cross_terms_not_discarded':True,'cq_trace_distance_factor_two_retained':True,
            'AW_Weyl_heat_difference_retained':True,'full_source_identity_requires_common_linear_map':True},
        'flag_spectral_checks':{'matrix_checks':flags,'B_convolution_tail_coefficient':str(Btail),
            'K_convolution_tail_coefficient':str(Ktail),'no_actual_blank_projection':True},
        'scope':['Only two specified coherent input preparations have the Gaussian-trajectory bound.',
                 'CPTP reference extension does not enlarge this preparation theorem.',
                 'Source coverage is the 48 symbols declared in 1249 section 3.2 and their declared masks/tails; bare B, psi, Y are not included.',
                 'The exact instrument acts on normal states; it need not preserve unbounded conjugate moments.',
                 'Completed bounded sources and intensive units only; no physical whole-block implementation, unbounded stress tail, or complete GR certification.',
                 'Analytic Gaussian/Duhamel estimates are reviewed separately; this script is not a wavefunction simulation.'],
        'input_sha256':{str(p.relative_to(STAGE)).replace('\\','/'):sha(p) for p in deps},
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    out=calculate()
    path=HERE/'independent_results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(out,f,indent=2,ensure_ascii=False,sort_keys=True);f.write('\n')
        receipt='FIRST_WRITE_EXCLUSIVE'
    else:
        assert json.loads(path.read_text(encoding='utf-8'))==out,'READ_ONLY_MISMATCH'
        receipt='READ_ONLY_MATCH'
    print(json.dumps({'round':1266,'status':out['status'],'receipt':receipt,
        'independent_groups':1,'result_sha256':sha(path)},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
