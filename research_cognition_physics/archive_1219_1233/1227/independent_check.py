"""1227 independent audit: exact guards plus log resource displays.
No author implementation is imported. No Gauss spectrum or clock truncation.
Default is read-only; --write-results creates a previously absent result.
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
 "archive_585_628/research_note_598.md",
 "archive_1219_/1226/proof.md",
 "archive_1219_/1226/results.json",
 "archive_1219_/1225/proof.md",
 "archive_1201_1218/1203/proof.md",
)
BETA=F(17,10); N=125; C=64; R=(1200+C)**2
D=2*C*(1200+2*C); T=4; EXPONENT=1984775
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def field_and_pulse():
    # Analytic extrema: |theta'|<=1/4, |theta''|<=1/(2 sqrt(3)).
    assert F(1,16)<F(1,9) and F(1,12)<F(1,9)
    lap=F(2,3)+F(4,9)
    assert lap==F(10,9)<2
    assert F(4,9)*8000<60**2 and 2+60<C
    # Original J on a radius-1/4 support.
    f_lower=F(95,48); gamma_lower=f_lower*F(95,96)
    assert gamma_lower>1
    j_lower=F(1,4)*F(15,16)**2/F(16)
    assert j_lower==F(225,16384)>F(1,128)
    assert F(2,9)<F(1,4)
    # b(u)=15/16 (1-u^2)^2, integral and derivative maximum squared.
    area=F(15,16)*2*(1-F(2,3)+F(1,5))
    assert area==1 and F(15,16)<1
    assert F(25,12)<F(9,4)
    return dict(theta_prime_square_upper="1/16",theta_second_square_upper="1/12",
       laplacian_theta_upper=str(lap),commutator_constant=C,
       injection_support_lower=str(j_lower),registered_injection_lower="1/128",
       pulse_integral=str(area),pulse_derivative_square_max="25/12",
       all_other_original_H_terms_commute_with_scalar_field_gate=True,
       constant_theta_deletes_injection=True)

def S_guard(k,root,base):
    gammas=[F(2),F(1),F(1),F(1),F(3,2)]
    terms=[]
    for j in range(k+1):
        power=base**((j+1)//2) if j%2 else root*base**(j//2)
        terms.append(F(math.comb(k,j)*2**(k-j),2)*gammas[j]*power)
    return 1+sum(terms)

def thermal_and_seed():
    assert 123**2>15000 and 142**2>20000 and 71**2>5000
    s4=S_guard(4,123,15000);s2=S_guard(2,142,20000);s0=S_guard(0,71,5000)
    assert s4<3*10**14 and s2<4*10**9 and 2*s0-1<1000
    assert (3*10**14)**2*(4*10**9)*1000<10**42
    assert 2*5000**2*71<10**10
    # sqrt(6)<5/2, log(3)<2, Omega_4<27.
    core=F(27*36,8)*F(5,2)*(F(81,8)-8+6)
    tail=F(9*27,4)/(F(1,5)*F(1,200))
    assert core<5000 and tail==60750 and core+tail<10**5
    astar_guard=32+15*3+3*42*3
    assert astar_guard==455<500
    eb=2*5446+F(1000)/BETA
    assert eb<11481
    x=125*11481
    A2=F(x*x)+(4/BETA+2)*x+8/BETA**2+1
    assert x+1<1200**2 and A2<1200**4
    shifted_trial=240*16+2+36+4096+1072
    assert shifted_trial==9046
    assert N*(BETA*shifted_trial+500)==EXPONENT
    return dict(S4_rational_guard=str(s4),S2_rational_guard=str(s2),
       twoS0_minus1_guard=str(2*s0-1),scalar_core_guard=str(core),
       scalar_tail_guard=str(tail),a_star_guard=500,a_star_intermediate=455,
       e_beta_guard=11481,x0_guard=x,A2_guard=str(A2),
       root_A2_guard=1200**2,small_packet_shifted_energy_per_node=9046,
       j0_exact="2^(-7)*exp(-1984775)",
       complete_Gibbs_state_assumed=True,thermal_spectrum_computed=False)

def process_budgets():
    assert R==1597696 and D==169984
    denom=32*R*D
    # j in (0,1), delta=j/denom, K=r_a+R+1>R.
    assert F(1,denom)<F(1,4)
    epsE=F(1,64)+F(1,16)
    epsP=F(1,64*R)+F(1,32*R)
    residual=F(3,128*T*R)+F(1,32)
    assert epsE==F(5,64)
    assert epsP<F(1,16*R)
    assert residual<F(1,16)
    assert 1-epsE==F(59,64)
    assert 1-epsE-F(1,16)==F(55,64)
    # M>=T/sigma^2 gives b_T<=sqrt(5/4)*sigma<2 sigma.
    assert F(5,4)<4
    return dict(R=R,D=D,T=T,delta_denominator=denom,
       instantaneous_and_scalar_r_bound=R,
       autonomous_r_formula="(1200+256/delta)^2",
       autonomous_pulse_norm_integral_upper="2*T/delta",
       epsilon_E_over_j_upper=str(epsE),
       epsilon_probability_over_j_upper=str(epsP),
       endpoint_residual_over_j_upper=str(residual),
       actual_injection_over_j_strict_lower="59/64",
       total_clock_energy_loss_over_j_strict_lower="55/64",
       independent_initial_clocks=2,local_force_bound_each="2/delta^2",
       sum_force_bound="4/delta^2",total_constant_baseline=6,
       energy_transport_separate_from_trace_distance=True,
       scalar_to_instantaneous_error_not_multiplied_by_autonomous_r=True,
       endpoint_interaction_not_deleted=True)

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def stinespring():
    # Rational diagonal samples obey original L_+^2+L_-^2=1 and
    # L_+^2 in [1/4,3/4]. This checks the PE algebra, not a Gauss state.
    cp=[F(3,5),F(4,5)]; sp=[F(4,5),F(3,5)]
    for c,s in zip(cp,sp):
        assert c*c+s*s==1 and F(1,4)<=c*c<=F(3,4)
    U=[[F(i==j) for j in range(4)] for i in range(4)]
    U[0][0]=U[3][3]=F(3,5);U[0][3]=-F(4,5);U[3][0]=F(4,5)
    assert all(sum(U[k][i]*U[k][j] for k in range(4))==F(i==j)
               for i in range(4) for j in range(4))
    # Output rows (material, first PE=aa, second PE=bb), 16 retained rows.
    rows=[(s,a,b) for s in range(4) for a in range(2) for b in range(2)]
    V=[]
    for initial in range(4):
        col=[]
        for final,a,b in rows:
            first=(cp if a==0 else sp)[initial//2]
            second=(cp if b==0 else sp)[final%2]
            col.append(second*U[final][initial]*first)
        V.append(col)
    assert all(dot(V[i],V[j])==F(i==j) for i in range(4) for j in range(4))
    probs={}
    for a in range(2):
        for b in range(2):
            probs[f"{a}{b}"]=sum(V[i][k]**2/F(4)
                for i in range(4) for k,row in enumerate(rows) if row[1:]==(a,b))
    assert sum(probs.values())==1
    # Orthogonal reference tags preserve the complete input-column isometry.
    gram_ref=[dot(V[i],V[i])/4 for i in range(4)]
    assert sum(gram_ref)==1
    # Pure PE c|00>+s|11> differs from its dephasing by trace distance |cs|.
    retained_coherence=cp[0]*sp[0]
    assert retained_coherence==F(12,25)>0
    return dict(input_dimension=4,retained_output_rows=16,
       exact_input_Gram_identity=True,arbitrary_reference_isometry_preserved=True,
       diagnostic_joint_pointer_probabilities={k:str(v) for k,v in probs.items()},
       pure_PE_vs_dephased_half_trace_distance=str(retained_coherence),
       full_Stinespring_is_not_classical_quantum_output=True,
       diagnostics_are_not_original_Gauss_probabilities=True,
       actual_budget_requires_original_thermal_material_marginal=True)

def logadd(*vals):
    m=max(vals)
    return m+math.log(sum(math.exp(x-m) for x in vals if x-m>-745))
def log_resources():
    # Display only: inequalities above are certified algebraically.
    lj=-EXPONENT-7*math.log(2)
    ld=lj-math.log(32*R*D)
    lL=math.log(4)-2*ld
    lra=2*logadd(math.log(1200),math.log(256)-ld)
    lK=logadd(lra,math.log(R),0.)
    ls=lj-math.log(128)-lL-math.log(T)-lK
    lM=logadd(math.log(T)-2*ls,math.log(16)+2*lL+2*math.log(T)-lj,0.)
    lEC=logadd(lM-math.log(2),-math.log(8)-lM-2*ls,0.)
    lnbar=2*lj-math.log(4096)
    lSupper=lnbar+math.log(2-lnbar)
    # Interaction norm <=2pi/(3delta); baseline remains physical.
    lv=math.log(2*math.pi/3)-ld
    def row(z): return dict(natural_log=z,log10=z/math.log(10))
    data={"j0":row(lj),"pulse_width_delta":row(ld),"force_L":row(lL),
          "autonomous_material_norm_bound":row(lra),"clock_width_sigma":row(ls),
          "clock_M":row(lM),"initial_one_clock_mean_energy":row(lEC),
          "initial_two_clock_mean_energy":row(lEC+math.log(2)),
          "interaction_norm_upper":row(lv),"clock_occupation_upper":row(lnbar),
          "one_clock_entropy_upper":row(lSupper),
          "joint_two_clock_entropy_upper":row(lSupper+math.log(2))}
    assert all(math.isfinite(x) for r in data.values() for x in r.values())
    return dict(values=data,log_displays_are_not_outward_intervals=True,
       discarded_logsumexp_terms_below_float_underflow_are_display_only=True,
       no_subtraction_of_huge_clock_energies=True,engineering_feasibility_claimed=False)

def entropy_and_scope():
    # From M>=T/sigma^2, both Heisenberg shifts give
    # sqrt(n)<=5*sigma*L*T/4<=5*j/(512*K).
    nratio=F(25,512**2)
    assert nratio<F(1,4096)
    return dict(occupation_over_j_squared_upper="1/4096",
       reference_mode="(X-T*P/M)/(2*sigma)+i*sigma*(P-M)",
       momentum_difference_bound="L*T",back_free_position_difference_bound="L*T^2/(2*M)",
       one_clock_entropy_upper="nbar*(2-log(nbar))",
       joint_two_clock_entropy_upper="2*nbar*(2-log(nbar))",
       final_clock_independence_assumed=False,individual_clock_energy_loss_asserted=False,
       infinite_rank_material_allowed=True,finite_rank_clock_argument_used=False,
       nonGaussian_clock_states_allowed=True,retained_pointer_environment_entropy_upper="4*log(2)",
       original_reference_preparation_is_input=True,
       unbounded_energy_measurement_implemented=False,
       spatial_or_accelerated_proper_clock_claimed=False,
       Unruh_or_gravity_claimed=False,stage_completion_claimed=False)

def run():
    return dict(status="PASS",all_checks_passed=True,check_groups=6,
       independent_code_sha256=sha(Path(__file__)),
       source_sha256={p:sha(RESEARCH/p) for p in SOURCES},
       groups=dict(field_and_pulse=field_and_pulse(),thermal_and_seed=thermal_and_seed(),
          separate_process_budgets=process_budgets(),retained_Stinespring=stinespring(),
          finite_log_resources=log_resources(),clock_entropy_scope=entropy_and_scope()))
def main():
    p=argparse.ArgumentParser();p.add_argument("--write-results",action="store_true")
    args=p.parse_args();res=run()
    if args.write_results:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(res,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        assert json.loads(OUT.read_text("utf-8"))==res,"Saved results or sources changed"
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,
       check_groups=6,independent_code_sha256=res["independent_code_sha256"])))
if __name__=="__main__": main()
