"""1229 independent audit: rational domain/budget guards and finite algebra fixtures.
No author code/results imported. Infinite-algebra statements are analytically
reviewed; finite matrices do not approximate the Dirac vacuum or its entropy.
Default read-only replay. --write-results exclusively creates a missing result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
OUT=HERE/"independent_results.json"
SOURCES=(
 "archive_1120_1158/1126/proof.md",
 "archive_1219_/1228/proof.md",
 "archive_1219_/1228/results.json",
 "archive_1219_/1228/NEXT.md",
 "archive_1219_/_admission/after1228_geometric_information_admission.md",
)
K=10**6; B2=F(137,5); C=16*K; CB=128*K
PI=F(22,7); TP=2; TF=4
DELTA=F(1,10**24); SIGMA=F(1,10**120); M=10**260
L=16/DELTA**2
QMIN=F(1,5*10**12)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(x):return {"exact":str(x),"display":float(x)}
def clock_spread(t):return SIGMA+F(t,2)/M/SIGMA
def clock_error(t):return L*t*clock_spread(t)

def wedge_and_boost():
    # Exact interval endpoints; |t-1| is convex, so endpoint minima suffice.
    margins=[F(z+8)-abs(F(t-1)) for z in [-1,5] for t in [0,4]]
    assert min(margins)==4
    support_first_min=1-max(abs(t-1) for t in [0,4])
    support_second_min=1-max(abs(t-3) for t in [0,4])
    assert support_first_min==support_second_min==-2
    assert support_first_min+8==6
    # h^2=p^2+1 and [h,z]=-i alpha_z => [h,k]=-i p_z.
    # Thus Heisenberg k(t)=k+t p, and to compare the same slice charge,
    # the Schrodinger coefficient of p is 1-T.
    coeff={str(t):1-t for t in [0,1,2,4]}
    assert coeff=={"0":1,"1":0,"2":-1,"4":-3}
    boost_mode=14*(K+F(53,10))+F(1,2)
    assert B2<F(53,10)**2 and boost_mode<15*K
    assert 2*PI*15*K<CB
    return dict(full_window=[0,4],preparation_cut=2,gate_centers=[1,3],
                physical_second_support_z=[-1,5],physical_wedge_margin=4,
                common_Cauchy_support_min_z=-2,Cauchy_wedge_margin=6,
                Schrodinger_p_coefficients=coeff,
                commutator="[h,k_z]=-i p_z",
                common_charge="k_z+8h at t=1",
                boost_mode_norm_upper=val(boost_mode),
                per_gate_boost_commutator_upper=CB,
                Gaussian_time_tails_included=True,
                infinite_future_same_wedge_claim=False)

def zmat(d):return [[F(0) for _ in range(d)]for _ in range(d)]
def eye(d):
    a=zmat(d)
    for i in range(d):a[i][i]=F(1)
    return a
def mt(a):return [list(x)for x in zip(*a)]
def plus(a,b,s=1):return [[u+s*v for u,v in zip(r,t)]for r,t in zip(a,b)]
def mm(a,b):
    return [[sum((a[i][k]*b[k][j]for k in range(len(b))),F(0))
             for j in range(len(b[0]))]for i in range(len(a))]
def ann(k,n):
    d=2**n;a=zmat(d)
    for s in range(d):
        if (s>>k)&1:a[s^(1<<k)][s]=F((-1)**((s&((1<<k)-1)).bit_count()))
    return a
def finite_car():
    a0=ann(0,2);a1=ann(1,2);I=eye(4)
    nu=mm(mt(a0),a0);nv=mm(mt(a1),a1)
    c=plus(a0,a1)
    nf=[[x/2 for x in row]for row in mm(mt(c),c)]
    nf0=plus(I,nf,-1)
    rhs=[[x/2 for x in row]for row in plus(nu,nv)]
    lhs=plus(mm(mm(nf,nu),nf),mm(mm(nf0,nu),nf0))
    assert lhs==rhs and mm(nf,nf)==nf
    assert plus(mm(nf,nf),mm(nf0,nf0))==I
    examples=[]
    for s in [F(1,16),F(1,4),F(1,2),F(15,16)]:
        weights=[s*(1-s),s*s,(1-s)*(1-s),s*(1-s)]
        rho=zmat(4)
        for i,p in enumerate(weights):rho[i][i]=p
        out=plus(mm(mm(nf,rho),nf),mm(mm(nf0,rho),nf0))
        p=sum(mm(out,nu)[i][i]for i in range(4))
        before=sum(mm(rho,nu)[i][i]for i in range(4))
        assert before==s and p==F(1,2)
        examples.append(dict(vacuum_comparator_population=str(s),
                             after_nonselective_record=str(p)))
    return dict(exact_operator_dephasing_identity=True,all_input_columns=4,
                preserved_full_PE_isometry=True,examples=examples,
                actual_continuum_state_replaced=False,
                old_PE_record_used_as_field_channel=False)

def record_kl():
    smin=1/(4*(K*K+B2));smax=B2/(K*K)
    alpha=2*(clock_error(TP)+C*DELTA)
    gamma=2*(clock_error(TF)+2*C*DELTA)
    eta=alpha+gamma
    zeta=2*(clock_error(TF)+C*DELTA)
    assert QMIN<smin and smax<1-QMIN
    assert eta<F(1,10**15) and zeta<F(1,10**16) and zeta<QMIN
    # e<11/4 from Taylor through 5, with ratio <=1/7 on the remaining tail.
    eupper=sum((F(1,math.factorial(j))for j in range(6)),F(0))
    eupper+=F(1,math.factorial(6))/(1-F(1,7))
    assert eupper<F(11,4)
    assert smax+zeta<(F(4,11))**24
    # exp(0.7)>2 proves log 2<0.7, without floating logarithms.
    assert sum((F(7,10)**j/F(math.factorial(j))for j in range(6)),F(0))>2
    assert F(1,2)-max(alpha,gamma)>F(49,100)
    information_lower=F(49,100)*24-F(7,10)
    assert information_lower==F(553,50)>11
    # eta log(1/eta) monotone here, log10<3. log(1/qmin)<40 since e>8/3.
    assert (F(8,3))**40>1/QMIN
    continuity=F(1,10**15)*(2+15*3)+2*F(1,10**15)*40+zeta/(QMIN-zeta)
    assert continuity<F(1,1000)
    return dict(vacuum_population_lower=val(smin),vacuum_population_upper=val(smax),
                common_reference_probability_floor=str(QMIN),
                preparation_probability_error=val(alpha),terminal_record_error=val(gamma),
                actual_record_to_fixed_field_channel_error=val(eta),
                vacuum_readout_error=val(zeta),
                finite_binary_KL_transport_upper=val(continuity),
                exact_field_and_actual_record_KL_lower=val(information_lower),
                second_clock_initial_independence_at_T2_assumed=False,
                infinite_algebra_relative_entropy_continuity_used=False)

def auxiliary_moments():
    d4=clock_spread(TF)
    assert d4<2*SIGMA
    # [N,g(t+nu a+nu* a*)] = g'(-nu a+nu* a*)-|nu|^2 g''.
    c=16*PI*SIGMA/DELTA**2+64*PI*SIGMA**2/DELTA**3
    x=c*TF
    assert x<F(1,10**69)
    # Differential Gronwall ||(N+1)psi||<=exp(cT)<=1/(1-cT).
    second_moment_upper=1/(1-x)**2
    assert second_moment_upper<4
    nbar=(SIGMA*L*TF+L*TF*TF/F(4*M)/SIGMA)**2
    assert nbar<F(1,10**140)
    clock_entropy_upper=2*F(1,10**140)*(2+140*3)
    assert clock_entropy_upper<F(1,10**136)
    assert F(28,10)+clock_entropy_upper<F(29,10)
    # Sharper log2<.7 is strict, so 4 log2+tiny<2.8 can be certified using log2<.694.
    assert sum((F(347,500)**j/F(math.factorial(j))for j in range(9)),F(0))>2
    total_entropy=4*F(347,500)+clock_entropy_upper
    assert total_entropy<F(28,10)
    # A simple deletion shows norm convergence alone does not transport N.
    deletion=[dict(n=n,trace_distance=str(F(1,n)),mean_N=1,mean_N2=n)
              for n in [10,100,1000]]
    return dict(independent_commutator_graph_constant=val(c),
                uniform_Nplus1_second_moment_upper=val(second_moment_upper),
                retained_coarse_uniform_guard=4,
                tail_expectation_bound="E[N 1_(N>R)] <= 4/R",
                each_actual_clock_number_upper=val(nbar),
                full_auxiliary_entropy_upper=val(total_entropy),
                faithful_reference="exp(-N)/Z times I_PE/16; Z=(1-exp(-1))^-2",
                reference_is_actual_bath=False,
                fixed_reference_embedding_correction="log(Z/Z_L)",
                moment_deletion_examples=deletion,
                entropy_continuity_used_for_lsc_limit=False)

def boost_geometry():
    ep=C*DELTA;ec=clock_error(TP)
    RH=F(2*C);RB=F(2*CB)
    RHA=F(4*C)/DELTA;RBA=F(4*CB)/DELTA
    errH=(RHA+RH)*ec+2*RH*ep
    errB=(RBA+RB)*ec+2*RB*ep
    errK=2*PI*errB
    assert errH<F(1,10**6) and errB<F(1,10**6) and errK<F(1,10**6)
    # The proof uses the still wider certified 1e-6 norm/expectation envelope.
    EH=RH+F(1,10**6);Babs=RB+F(1,10**6)
    eps=F(1,10**12)
    geom=eps*(2*PI*Babs)+(1+eps)*2*PI*(2*eps)*EH
    assert geom<F(1,100)
    endpoint=L*(3*SIGMA+L*TF*TF/F(2*M))
    initial_tail=L*SIGMA
    assert endpoint<F(1,10**60) and 0<initial_tail<endpoint
    positive_boost=F(11)-F(28,10)
    assert positive_boost==F(41,5)
    return dict(HF_weighted_error=val(errH),boost_weighted_error=val(errB),
                modular_energy_weighted_error=val(errK),
                uniform_geometry_error_upper=val(geom),
                strict_modular_energy_positive_margin=str(positive_boost),
                initial_interaction_upper=val(initial_tail),
                initial_plus_terminal_interaction_upper=val(endpoint),
                free_initial_energy="M + 1/(4 M sigma^2) + 6",
                total_initial_energy_adds_V0=True,
                geometry_errors="uniform origin/time shifts and uniform area multiplier only",
                pointwise_T00_operator_bound_claim=False,
                total_apparatus_boost_conservation_claim=False)

def entropy(rho):
    e=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert float(e.min())>-1e-12
    return float(-sum(x*math.log(x)for x in e if x>1e-14))
def rel(rho,sigma):
    e,v=np.linalg.eigh(sigma)
    assert float(e.min())>0
    logs=(v*np.log(e))@v.conj().T
    return -entropy(rho)-float(np.trace(rho@logs).real)
def marginals(rho):
    a=rho.reshape(2,2,2,2)
    return np.einsum("iaja->ij",a),np.einsum("iaib->ab",a)
def finite_entropy_chain():
    omega=np.diag([.9,.1]).astype(complex)
    tau=np.diag([1,math.exp(-1)]).astype(complex)
    tau/=np.trace(tau)
    blank=np.diag([1,0]).astype(complex)
    rho0=np.kron(omega,blank);reference=np.kron(omega,tau)
    KF=-np.diag(np.log([.9,.1]));N=np.diag([0,1])
    swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    cnot=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],complex)
    cases={"swap":swap,"partial_swap":(np.eye(4)-1j*swap)/math.sqrt(2),"record_only":cnot}
    initialD=rel(rho0,reference);out={}
    assert initialD>0  # pure auxiliary is not the faithful reference
    for name,U in cases.items():
        assert np.linalg.norm(U.conj().T@U-np.eye(4))<1e-13
        rho=U@rho0@U.conj().T
        rF,rE=marginals(rho)
        DF=rel(rF,omega);DE=rel(rE,tau);DJ=rel(rho,reference)
        I=entropy(rF)+entropy(rE)-entropy(rho)
        deltaK=float(np.trace((rF-omega)@KF).real)
        deltaN=float(np.trace((rE-blank)@N).real)
        joint_balance=DJ-initialD-deltaK-deltaN
        chain=DJ-DF-DE-I
        reduced=DF+I-deltaK-entropy(rE)
        assert max(abs(joint_balance),abs(chain),abs(reduced))<1e-12
        out[name]=dict(Dfield=DF,mutual_information=I,aux_entropy=entropy(rE),
                       delta_K=deltaK,
                       joint_balance_residual=joint_balance,chain_residual=chain,
                       reduced_balance_residual=reduced)
    assert out["swap"]["delta_K"]<0<out["swap"]["Dfield"]
    assert abs(out["record_only"]["Dfield"])<1e-13
    assert out["record_only"]["mutual_information"]>.3
    return dict(initial_relative_entropy_to_faithful_reference=initialD,
                cases=out,omitting_aux_entropy_fails=True,
                fixture_not_infinite_wedge_entropy_computation=True,
                pure_aux_not_substituted_for_faithful_tau=True)

def calculate():
    groups={
      "wedge_time_and_boost_sign":wedge_and_boost(),
      "exact_CAR_field_channel":finite_car(),
      "vacuum_probability_and_finite_KL":record_kl(),
      "uniform_auxiliary_domain_guards":auxiliary_moments(),
      "same_state_boost_geometry_energy":boost_geometry(),
      "general_initial_state_entropy_chain_fixture":finite_entropy_chain(),
    }
    return dict(schema="1229-independent-v1",status="PASS",all_checks_passed=True,
                check_groups=len(groups),groups=groups,
                independent_code_sha256=sha(Path(__file__)),
                source_sha256={p:sha(RESEARCH/p)for p in SOURCES},
                scope="Exact finite guards and distinct algebra fixtures; domain/lsc arguments require analytic review, not a vacuum truncation.")

def main():
    p=argparse.ArgumentParser();p.add_argument("--write-results",action="store_true")
    args=p.parse_args();result=calculate()
    if args.write_results:
        with OUT.open("x",encoding="utf-8",newline="\n")as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write("\n")
    if json.loads(OUT.read_text(encoding="utf-8"))!=result:
        raise AssertionError("Saved independent result mismatch; no files changed.")
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,
                          check_groups=result["check_groups"]),ensure_ascii=False))
if __name__=="__main__":main()
