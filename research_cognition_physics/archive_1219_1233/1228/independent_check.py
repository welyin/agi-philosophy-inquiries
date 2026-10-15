"""1228 independent checks: exact arithmetic guards and complete finite CAR fixture.
No author code or author results are imported. The fixture is not the Dirac vacuum.
The continuum identities, domains and effective-region permission require the review.
Default only replays saved results; --write-results refuses an existing result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
OUT=HERE/"independent_results.json"
SOURCES=(
    "archive_1120_1158/1126/proof.md",
    "archive_1219_/1227/proof.md",
    "archive_1219_/1227/results.json",
    "archive_1219_/_admission/after1227_reference_flow_priority.md",
)
PI_LO=F(157,50)
PI_HI=F(22,7)
K=10**6
C=16*K
T=2
DELTA=F(1,10**24)
SIGMA=F(1,10**120)
M=10**260
L=16/DELTA**2
A0=2**28

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def show(x): return {"exact":str(x),"display":float(x)}
def poly_der(p): return [F(i)*p[i] for i in range(1,len(p))]
def poly_mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j]+=a*b
    return out
def poly_int(p):
    return sum((2*a/F(i+1) for i,a in enumerate(p) if i%2==0),F(0))
def packet():
    b=[F(0)]*17
    for j in range(9): b[2*j]=F((-1)**j*math.comb(8,j))
    b2=poly_mul(b,b); norm=poly_int(b2)
    m2=poly_int([F(0)]*2+b2)/norm
    m4=poly_int([F(0)]*4+b2)/norm
    d1=poly_int(poly_mul(poly_der(b),poly_der(b)))/norm
    d4=b
    for _ in range(4):d4=poly_der(d4)
    rho4=sum((2*abs(a)/F(i+1) for i,a in enumerate(d4)),F(0))/poly_int(b)
    assert norm==F(4294967296,9917826435)
    assert m2==F(1,35) and m4==F(3,35*37) and d1==F(44,5)
    assert rho4==F(6891885,16)<A0
    rlo=1-(2*PI_HI)**2*m2/2
    rhi=1-(2*PI_LO)**2*m2/2+(2*PI_HI)**4*m4/24
    assert F(43,100)<rlo<rhi<F(59,100)
    assert 1+3*d1<F(53,10)**2
    return dict(norm=show(norm),moment2=str(m2),moment4=str(m4),
                derivative_energy=str(d1),rho_derivative4_L1_upper=show(rho4),
                r_lower=show(rlo),r_upper=show(rhi),
                actual_packet_regular_class="C7 and H8, not C-infinity",
                quadrature_used=False)

D=16
def zero():return [[F(0) for _ in range(D)] for _ in range(D)]
def eye():
    a=zero()
    for i in range(D):a[i][i]=F(1)
    return a
def transpose(a):return [list(x) for x in zip(*a)]
def madd(a,b,sign=1):return [[x+sign*y for x,y in zip(r,s)] for r,s in zip(a,b)]
def mm(a,b):
    out=zero()
    for i in range(D):
        for k in range(D):
            if a[i][k]:
                for j in range(D):
                    if b[k][j]:out[i][j]+=a[i][k]*b[k][j]
    return out
def annihilator(k):
    out=zero()
    for s in range(D):
        if (s>>k)&1:
            out[s^(1<<k)][s]=F((-1)**((s&((1<<k)-1)).bit_count()))
    return out
def number(coeff,norm,ann):
    c=zero()
    for k,w in enumerate(coeff):
        for i in range(D):
            for j in range(D):c[i][j]+=w*ann[k][i][j]
    a=mm(transpose(c),c)
    return [[x/norm for x in row] for row in a]
def finite_car():
    ann=[annihilator(k) for k in range(4)]
    I=eye();Z=zero()
    for i in range(4):
        for j in range(4):
            anti=madd(mm(ann[i],transpose(ann[j])),mm(transpose(ann[j]),ann[i]))
            assert anti==(I if i==j else Z)
    nf=number([1,0,1,0],F(2),ann)
    ng=number([3,4,3,4],F(50),ann)
    assert mm(nf,nf)==nf and mm(ng,ng)==ng
    # Orthogonal PE words are 0000,0011,1100,1111; no environment is erased.
    kraus=[mm(p2,p1) for p1 in [madd(I,nf,-1),nf]
           for p2 in [madd(I,ng,-1),ng]]
    gram=zero()
    for k in kraus:gram=madd(gram,mm(transpose(k),k))
    assert gram==I  # all 16 input columns, therefore every passive reference extension
    vacuum=12  # fill modes 2,3; S=diag(0,0,1,1) in this finite fixture
    p=[sum((k[i][vacuum]**2 for i in range(D)),F(0)) for k in kraus]
    assert p==[F(17,50),F(4,25),F(4,25),F(17,50)]
    alpha=F(3,5);cfg=F(3,10)
    formula=F(1,4)-cfg*cfg+alpha*alpha/2
    assert formula==p[3]
    repeat=[mm(p2,p1) for p1 in [madd(I,nf,-1),nf]
            for p2 in [madd(I,nf,-1),nf]]
    rp=[sum((k[i][vacuum]**2 for i in range(D)),F(0)) for k in repeat]
    assert rp==[F(1,2),F(0),F(0),F(1,2)]
    # Tracing E creates orthogonal record blocks; it changes global purity.
    dephased_purity=sum(x*x for x in p)
    assert dephased_purity==F(353,1250)<1
    return dict(fixture_one_particle_dimension=4,fixture_Fock_dimension=16,
                full_PE_dimension=16,all_input_isometry_columns=16,
                joint_order=["00","01","10","11"],
                joint_probabilities=[str(x) for x in p],
                repeated_mode_probabilities=[str(x) for x in rp],
                p11_formula=str(formula),
                complete_output_purity="1",after_record_dephasing_purity=str(dephased_purity),
                Dirac_vacuum_replaced_by_fixture=False)

def flow():
    kappa=2*PI_HI
    # Independent broader norm constants than the author: B<5.3 and lab<7.4 tau.
    # Boost residual integrated as 17 + 3 K s, rather than its interval supremum.
    ew=(17*kappa+F(3,2)*kappa*kappa)/K
    et=F(74,10)*2*kappa/K
    spec=F(53,10)/K
    target=(1-F(59,100)**2)/4
    lower=target-spec-2*(ew+et)
    assert lower>F(1624,10000)
    assert F(3)*F(53,10)+F(1,2)<17
    assert F(93,5)<F(44,10)**2 and F(44,5)<3**2
    assert kappa/K<F(1,100000)
    # e^eta <= 1/(1-eta), and ||h f||^2=K^2+137/5.
    energy_square=(F(K*K)+F(137,5))/(1-kappa/K)**2
    assert energy_square<(2*K)**2
    assert F(2)*PI_HI*(2*K+14)<C
    return dict(K=K,boost_parameter_upper=str(kappa/K),
                spectral_error_upper=show(spec),
                boost_vector_error_upper=show(ew),lab_vector_error_upper=show(et),
                comparison_gap_lower=show(target),ideal_event_gap_lower=show(lower),
                both_exact_mode_occupancies="1/2",
                physical_same_second_time=True,
                boost_equals_actual_accelerated_detector=False)

def source_guards():
    # Each 4x4 spectral projection has rank 2. Two flavors times rank 2
    # cancels the 1/4 of ((E-E')/2)^2, leaving D4^2 (E+E')^-6.
    assert F(2*2,4)==1
    # (2*pi)^-3 times integral(1+p^2)^-3 = pi^2/4.
    hs_ratio=F(8)/(32*PI_LO)
    author_broader_ratio=F(8)/(16*PI_LO)
    assert hs_ratio<author_broader_ratio<1
    assert 16*16*3<28*28
    assert 2*PI_HI*(2*K+14)<C
    return dict(two_flavors_times_spin_rank_over_energy_factor="1",
                spatial_Fourier_convention="unitary; convolution coefficient (2 pi)^(-3/2)",
                momentum_integral="pi^2/4",
                independent_pair_HS_square_over_D4square_upper=show(hs_ratio),
                author_wider_pair_ratio_upper=show(author_broader_ratio),
                TQ_vacuum_vector_norm_upper=A0,
                commutator_constant=C,
                CAR_commutator_sign="[T(Q),c(u)] = -c(A_Q u)",
                full_Fock_TQ_bounded_claim=False,
                future_source_meaning="terminal state followed by a specified free reference continuation",
                actual_autonomous_future_source_history_certified=False)

def transport():
    ep=2*C*DELTA
    ec=L*T*(SIGMA+F(T,2)/M/SIGMA)
    RH=F(2*C);RQ=F(A0)+RH
    RHA=F(2*T*C)/DELTA;RQA=F(A0)+RHA
    qerr=(RQA+RQ)*ec+2*RQ*ep
    herr=(RHA+RH)*ec+2*RH*ep
    per_menu_event=2*(ep+ec)
    both_menus=2*per_menu_event
    base=F(flow()["ideal_event_gap_lower"]["exact"])
    actual_gap=base-both_menus
    assert qerr<F(1,10**6) and herr<F(1,10**6)
    assert actual_gap>F(16,100)
    assert 4*qerr<F(1,10**6)
    return dict(pulse_vector_error=show(ep),quantum_clock_vector_error=show(ec),
                actual_TQ_vector_norm_bound=str(RQA),ideal_TQ_vector_norm_bound=str(RQ),
                TQ_branch_absolute_error=show(qerr),
                HF_total_or_branch_absolute_error=show(herr),
                both_menu_probability_error=show(both_menus),
                actual_joint_gap_lower=show(actual_gap),
                record_branches=4,weighted_transport_not_trace_distance=True,
                comparison_order=["actual two clocks","finite scalar pulses","complete instantaneous instruments"],
                scalar_area_replaces_quantum_operator_norm_integral=False)

def clocks():
    assert F(T,2)/M/SIGMA<SIGMA
    residual=L*(3*SIGMA+L*T*T/F(2*M))
    nbar=(SIGMA*L*T+L*T*T/F(4*M)/SIGMA)**2
    assert residual<F(1,10**60) and nbar<F(1,10**140)
    # g(n)<=n(2-log n), g increasing, log(10)<3.
    entropy_upper=2*F(1,10**140)*(2+140*3)
    assert entropy_upper<F(1,10**136)
    initial=F(M)+F(1,4*M)/SIGMA**2+2
    assert initial>10**260 and initial<10**260+3
    assert 2*DELTA<F(8)*PI_LO/K
    return dict(clock_count=2,M=str(M),sigma=str(SIGMA),pulse_width=str(DELTA),
                final_time=T,total_force_bound=str(L),
                total_initial_clock_energy_exact=str(initial),
                pointer_baseline_energy=4,
                endpoint_interaction_sum_upper=show(residual),
                each_back_free_mode_occupation_upper=show(nbar),
                total_clock_entropy_upper=show(entropy_upper),
                total_clock_energy_reaction="|Delta E_C1+C2 + Delta E_F| <= endpoint bound",
                each_clock_loss_sign_claim=False,
                terminal_clock_Gaussian_or_independent_assumed=False,
                finite_region_effective_gate_permission_is_input=True,
                pointwise_Lorentz_implementation_certified=False,
                actual_Unruh_thermometer_certified=False,
                same_Dirac_independent_area_entropy_bridge_certified=False)

def calculate():
    groups={
        "packet_and_real_cosine_moments":packet(),
        "full_finite_CAR_record_fixture":finite_car(),
        "independent_high_carrier_flow_gap":flow(),
        "normal_Wick_source_domain_guards":source_guards(),
        "three_process_weighted_transport":transport(),
        "positive_clock_resource_bounds":clocks(),
    }
    return dict(schema="1228-independent-v1",status="PASS",all_checks_passed=True,
                check_groups=len(groups),groups=groups,
                independent_code_sha256=sha(Path(__file__)),
                source_sha256={p:sha(RESEARCH/p) for p in SOURCES},
                scope="Independent rational guards and an exact finite CAR fixture; continuum domains and effective gate permissions are reviewed separately.")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--write-results",action="store_true")
    args=p.parse_args()
    data=calculate()
    if args.write_results:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write("\n")
    saved=json.loads(OUT.read_text(encoding="utf-8"))
    if saved!=data:raise AssertionError("saved independent results differ")
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,
                          check_groups=data["check_groups"]),ensure_ascii=False))

if __name__=="__main__":main()
