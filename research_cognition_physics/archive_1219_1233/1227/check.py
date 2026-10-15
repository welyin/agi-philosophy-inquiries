"""1227: same untruncated Gauss material, two records, clock energy and entropy.
Four groups. Exact rational inequalities certify the coarse constants.
Floating values are log-scale displays / pointer identities, not full dynamics.
Default is read-only; --write-results creates a new receipt only.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TARGET=HERE/"results.json"
BETA=Q(17,10)
N=125
C=64
T=4
R=(1200+C)**2
D=2*C*(1200+2*C)
BEXP=N*(BETA*9046+500)
SOURCES=[
    "archive_554_584/574/joint_curved_quantum_source.py",
    "archive_585_628/research_note_598.md",
    "archive_585_628/research_note_603.md",
    "archive_702_741/research_note_723.md",
    "archive_702_741/research_note_725.md",
    "archive_1120_1158/1128/proof.md",
    "archive_1201_1218/1203/proof.md",
    "archive_1219_/1225/proof.md",
    "archive_1219_/1226/proof.md",
    "archive_1219_/1226/check.py",
    "archive_1219_/1226/results.json",
    "archive_1219_/_admission/after1226_stage_priority_review.md",
    "archive_1219_/1227/proof.md",
    "archive_1219_/1227/check.py",
]


def exact_material_guards():
    # Inherited trace formulas: exact broad inequalities replace numerical rounding.
    assert Q(17,20)/4000 > Q(1,5000)
    assert Q(17,20)/4 > Q(1,5)
    assert 4*5000**5 < 10**20  # 2*5000^(5/2) < 10^10.
    core=Q(27*9*81,4)
    tail=Q(9*27,4)*5*200
    assert core<5000 and tail==60750 and core+tail<10**5
    s4=82*15000**3
    s2=10*20000**2
    assert s4<3*10**14 and s2<=4*10**9
    assert 1+2*71<1000 and 71**2>5000
    assert (3*10**14)**2*(4*10**9)*1000<10**42
    assert 32+45+3*126==455<500
    ebar=2*5446+Q(1000,1)/BETA
    assert ebar<11481
    x=N*11481
    m2=Q(x*x)+(4/BETA+2)*x+8/(BETA*BETA)+1
    assert x==1435125 and x+1<1200**2 and m2<1200**4
    assert Q(10,9)<2 and (C-2)**2>Q(4*8000,9)
    assert Q(1,12)<Q(1,9)  # derivative envelopes <= 1/3.
    assert Q(95,48)*Q(95,96)>1
    assert Q(225,16384)>Q(1,128)
    small_trial=3840+2+36+4096+1072
    assert small_trial==9046 and BEXP==1984775
    previous=json.loads((ROOT/"archive_1219_/1226/results.json").read_text("utf8"))
    thermal=previous["groups"]["uniform_thermal_comparison"]["rows"]
    assert thermal[1]["e_beta_per_node"]<11481
    assert max(thermal[1]["a_star_at_half_beta"],thermal[2]["a_star_at_half_beta"])<500
    return dict(status="PASS",beta=str(BETA),N=N,
                trace_a_upper=500,trace_log_budget_rational=455,
                scalar_core_guard=str(core),scalar_tail_guard=str(tail),
                e_beta_rational_upper=str(ebar),chosen_e_upper=11481,
                A_mean_upper=x+1,A_second_moment_upper_rational=str(m2),
                A_second_norm_upper=1200**2,
                material_commutator_constant=C,
                positive_injection_formula="j0 = 2^(-7) exp(-1984775)",
                narrow_trial_shifted_energy_per_node=small_trial,
                no_full_Gibbs_spectrum_or_moment_evaluation=True)


def theta(s):
    return .5*math.acos(math.sin(s)/2)


def pointer_and_field():
    # Actual P/E logical subspace, all four physical pointer states retained.
    Y=np.zeros((4,4),dtype=complex)
    Y[0,3]=-1j;Y[3,0]=1j
    blank=np.array([1.,0.,0.,0.],dtype=complex)
    rows=[];max_error=0.
    for s in (-1.2,-.7,0.,.4,1.2):
        th=theta(s)
        U=np.eye(4,dtype=complex)
        for j in (0,3): U[j,j]=math.cos(th)
        U[0,3]=-math.sin(th);U[3,0]=math.sin(th)
        out=U@blank
        expected=np.array([math.sqrt(.5+math.sin(s)/4),0,0,
                           math.sqrt(.5-math.sin(s)/4)],dtype=complex)
        err=float(np.linalg.norm(out-expected))
        max_error=max(max_error,err)
        assert err<5e-16 and np.linalg.norm(U.conj().T@U-np.eye(4))<1e-15
        # Real original target ray phi=(0,0,0,0,s), with w=1.
        F=2-s*s/6
        gamma=F*(1-s*s/12)
        den=1-math.sin(s)**2/4
        dth=-math.cos(s)/(4*math.sqrt(den))
        ddth=3*math.sin(s)/(16*den**1.5)
        delta_s=-F*s/6
        delta_th=ddth*gamma+dth*delta_s
        J=gamma*dth*dth/2
        assert abs(dth)<=1/3 and abs(ddth)<=1/3 and abs(delta_th)<2
        assert 0<J<.25
        rows.append(dict(s=s,theta=th,plus_probability=float(abs(out[0])**2),
                         minus_probability=float(abs(out[3])**2),
                         Gamma_s=gamma,Delta_s=delta_s,Delta_theta=delta_th,
                         original_energy_injection_density_w1=J))
    assert abs(np.vdot(blank,Y@blank))==0
    # Constant angle deletes actual field sensitivity and the original positive cost.
    assert abs(rows[2]["plus_probability"]-rows[3]["plus_probability"])>.09
    return dict(status="PASS",pointer_dimensions="P1 E1 P2 E2 = 2^4",
                pointer_bare_positive_energy=4,rows=rows,maximum_Stinespring_error=max_error,
                constant_angle_deleted_field_injection=0,
                constant_angle_deleted_field_readout_dependence=True,
                full_state="coherent Stinespring with E retained; CQ records only after matching query/trace",
                not_a_finite_matrix_replacement_of_material_or_clock=True)


def logadd(*xs):
    m=max(xs)
    return m+math.log(sum(math.exp(x-m) for x in xs))


def parameter_logs():
    assert R==1597696 and D==169984
    j=-float(BEXP)-7*math.log(2)
    delta=j-math.log(32*R*D)
    L=math.log(4)-2*delta
    ra=2*logadd(math.log(1200),math.log(256)-delta)
    S=logadd(ra,math.log(R+1))
    sigma=j-math.log(128)-L-math.log(T)-S
    M=logadd(math.log(T)-2*sigma, math.log(16)+2*L+2*math.log(T)-j, 0.)
    single_clock_energy=logadd(M-math.log(2),-math.log(8)-M-2*sigma,0.)
    initial_clock_energy=math.log(2)+single_clock_energy
    return dict(log_j0=j,log_delta=delta,log_force_bound=L,
                log_actual_material_A_norm_bound=ra,
                log_r_a_plus_R_plus_one=S,log_sigma=sigma,log_M=M,
                log_single_clock_energy=single_clock_energy,
                log_initial_two_clock_energy=initial_clock_energy)


def error_certificate():
    # Exact symbolic implications hold for every j0 in (0,1), hence for fixed (14).
    assert Q(1,32*R*D)<Q(1,4)
    assert Q(1,16)+Q(1,64)==Q(5,64)
    endpoint=Q(3,128*T)+Q(1,32)
    assert endpoint<Q(1,16)
    probability_ratio=Q(1,32*R)+Q(1,64*R)
    assert probability_ratio<Q(1,16*R)
    assert 1-Q(5,64)==Q(59,64)
    assert 1-Q(5,64)-Q(1,16)==Q(55,64)
    logs=parameter_logs()
    assert logs["log_delta"]<math.log(.25)
    assert logs["log_M"]>=math.log(T)-2*logs["log_sigma"]-1e-8
    assert logs["log_M"]>=math.log(16)+2*logs["log_force_bound"]+2*math.log(T)-logs["log_j0"]-1e-8
    return dict(status="PASS",R=R,D=D,T=T,
                fixed_parameters=[
                    "delta=j0/(32 R D)","L=4/delta^2","r_a=(1200+256/delta)^2",
                    "sigma=j0/(128 L T (r_a+R+1))",
                    "M=T/sigma^2+16 L^2 T^2/j0+1"],
                classical_pulse_total_area=2,
                actual_operator_time_bound="2 T/delta (not area 2)",
                pulse_material_energy_error_over_j0="<= 1/16",
                clock_tracking_material_energy_error_over_j0="< 1/64",
                total_material_energy_error_over_j0="< 5/64",
                endpoint_interaction_bound_over_j0=str(endpoint)+" < 1/16",
                total_record_halftrace_over_j0="< 1/(16 R)",
                actual_material_gain_lower_over_j0="59/64",
                actual_clock_loss_lower_over_j0="55/64",
                number_of_independently_prepared_local_clocks=2,
                force_bounds="L1=L2=2/delta^2; L1+L2=L; sum Lj^2 <= L^2",
                energy_identity="Delta EH + Delta E_C1C2 + Delta E_PE + Delta V = 0",
                parameter_log_displays=logs,
                displays_are_not_interval_certification=True)


def entropy_resources():
    logs=parameter_logs()
    # M>=T/sigma^2 makes the second annihilator term <= one quarter the first.
    assert Q(25,512**2)<Q(1,4096)
    logn=2*logs["log_j0"]-math.log(4096)
    assert logn<0
    logS=math.log(2)+logn+math.log(2-logn)
    log_actual_n_envelope=2*(math.log(Q(5,4))+logs["log_sigma"]+
                            logs["log_force_bound"]+math.log(T))
    assert log_actual_n_envelope<logn
    ln10=math.log(10)
    records={
        "j0":logs["log_j0"]/ln10,
        "pulse_width_delta":logs["log_delta"]/ln10,
        "clock_spatial_width_sigma":logs["log_sigma"]/ln10,
        "clock_mass_parameter":logs["log_M"]/ln10,
        "initial_two_clock_energy":logs["log_initial_two_clock_energy"]/ln10,
        "each_clock_number_upper":logn/ln10,
        "joint_two_clock_entropy_upper":logS/ln10,
    }
    return dict(status="PASS",
                each_clock_number_bound="n_Cj <= (sigma L T + L T^2/(4 M sigma))^2 <= j0^2/4096",
                joint_two_clock_entropy_bound="S(C1 C2) <= 2 nbar (2-log nbar)",
                no_final_clock_independence_assumption=True,
                no_finite_rank_assumption=True,
                actual_clock_may_be_non_Gaussian=True,
                initial_two_clock_energy="2 (M/2 + 1/(8 M sigma^2) + 1)",
                pointer_environment_entropy_upper_display=4*math.log(2),
                material_entropy_finite_from="finite actual A mean and inherited positive-temperature trace",
                log10_displays=records,
                no_empirical_precision_or_practical_feasibility_claim=True,
                all_four_pointer_environment_systems_retained=True)


def run():
    groups=dict(material_commutator_and_thermal_moments=exact_material_guards(),
                original_field_Stinespring_identity=pointer_and_field(),
                two_comparisons_and_full_energy_certificate=error_certificate(),
                infinite_rank_clock_entropy_and_resources=entropy_resources())
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    # No hidden CR/backspace etc in the new proof; Markdown display math preserved.
    proof=(HERE/"proof.md").read_bytes()
    assert not [x for x in proof if x<32 and x not in (9,10)]
    return dict(round=1227,status="PASS",check_groups=4,tests_run=4,failures=0,errors=0,
                source_sha256=hashes,groups=groups,
                scope="one finite sufficient same-untruncated-material two-record clock certificate",
                claims_not_made=["full dynamics numerical simulation","practical clock specification",
                    "three-dimensional clock position","accelerated proper time","Unruh response",
                    "full Lorentz-causal certification of the original Gauss graph","gravitational backreaction",
                    "autonomous preparation of all hardware","stage completion"])


def main():
    p=argparse.ArgumentParser();p.add_argument("--write-results",action="store_true")
    args=p.parse_args();d=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            json.dump(d,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        assert json.loads(TARGET.read_text("utf8"))==d,"Saved receipt differs; no files changed"
    print(json.dumps(dict(round=1227,status="PASS",check_groups=4,
        receipt_mode="created" if args.write_results else "read_only_verified",
        log10_j0=d["groups"]["infinite_rank_clock_entropy_and_resources"]["log10_displays"]["j0"],
        no_actual_Gauss_or_clock_truncation=True),ensure_ascii=False))


if __name__=="__main__":
    main()

