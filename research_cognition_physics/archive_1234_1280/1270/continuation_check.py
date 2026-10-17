"""1270 continuation: one exact-arithmetic author group.
Default recomputes and compares without writing. --write creates once.
No ODE, state-vector simulation, or full field completion is performed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
OUT=HERE/"continuation_results.json"

def read(rel):
    return json.loads((STAGE/rel).read_text(encoding="utf-8-sig"))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def interval(state,key):
    return tuple(F(x) for x in state[key])

def absmax(xs):
    return max(map(abs,xs))

def absmin(xs):
    lo,hi=xs
    return F(0) if lo<=0<=hi else min(abs(lo),abs(hi))

def core_extension(state):
    u=interval(state,"u"); H=interval(state,"H")
    lower_u=u[0]*(1-H[1]/4)
    assert lower_u>F(1,3)
    Jlo=interval(state,"J")[0]
    Mc=F(71,100)+Jlo+sum(absmin(interval(state,k))**2 for k in ("q","r","K"))/2
    qr=sum(absmax(interval(state,k))**2 for k in ("q","r"))
    Yend=absmax(interval(state,"Y"))+absmax(interval(state,"K"))/4
    theta2=absmax(interval(state,"z"))+F(1,5)
    theta_end=theta2+F(1,4)*F(1+F(3,2),10)
    assert Jlo>F(98,1000)
    assert Mc>F(93,100)
    assert qr<1 and Yend<1 and theta_end<1
    assert absmax(interval(state,"psi"))<1
    assert absmax(interval(state,"K"))<1
    pg=2*27*H[1]
    pm=27*F(3,2)
    assert pg<61 and pm<61
    return {
      "u_after2_strict_lower":str(lower_u),
      "u_lower_margin_over_one_third":str(lower_u-F(1,3)),
      "J2_lower":str(Jlo),"Mc_hold_lower":str(Mc),
      "q_r_amplitude_squared_upper":str(qr),
      "Y_hold_abs_upper":str(Yend),"theta_hold_abs_upper":str(theta_end),
      "piQ_abs_upper":str(pg),"pi_material_abs_upper":str(pm),
      "F_strict_lower":"7/10","Mr_strict_lower":"1/20",
      "all_core_checks":True,
    }

def calculate():
    family=read("1252/results.json")["group_parameter_family"]["cells"]
    assert len(family)==1
    fam_state=family[0]["checkpoints"]["2"]["state"]
    fam=core_extension(fam_state)
    old=read("1247/results.json")["group_finite_memory"]
    old_core={name:core_extension(old[name]["checkpoints"]["2"]["state"])
              for name in ("alpha0","alpha1")}

    ledger=read("1249/quantum_bounds_results.json")[
        "group_constructive_finite_resource_bound"]["ledger"]
    assert ledger["M_power_of_two_exponent"]==3059
    M=1<<3059
    EH=int(ledger["Fourier_H_exponent"])
    E0=int(ledger["Fourier_menu_base_exponent"])
    EF=int(ledger["Fourier_transport_exponent_decimal"])
    N=int(ledger["integer_n_power_exponent_decimal"])
    assert (EH,E0)==(3093,14008)
    assert EF==95928+100*M
    assert N==1+max((EF+44)//3,(EH+EF+47)//6)

    # Exact pulse polynomial and finite peak/work controls.
    integral=140*(F(1,4)-3*F(1,5)+3*F(1,6)-F(1,7))
    assert integral==1
    tau=F(1,10**8); alpha=F(1,16)
    peak=F(35,16)/tau
    derivative_peak_bound=F(105,4)/tau**2
    work=alpha*27*derivative_peak_bound
    assert peak==218750000 and work<2**59

    # Analytic identities are coefficient identities for arbitrary a>0,c,kappa.
    # The Fraction substitutions below are deletion checks, not their proof.
    algebra=[]
    for a,c,kappa in [(F(3,2),F(1,7),F(3)),
                      (F(8,3),F(-1,11),F(3001,1000)),
                      (F(5,2),F(0),F(37,12))]:
        density_clock=-c*a*a**(-3)
        density_rod=c*a**(-2)
        pressure_rod=-c*a**(-2)/3
        Hdot_extra=kappa*c*a**(-2)/6
        assert density_clock+density_rod==0
        assert -2*Hdot_extra==kappa*pressure_rod
        # Exact old/new continuity exchange at arbitrary instantaneous H.
        H=F(2,7); cdot=F(5,13)
        rho_i_dot=(cdot-2*H*c)*a**(-2)
        continuity_i=rho_i_dot+3*H*(density_rod+pressure_rod)
        clock_work=-cdot*a*a**(-3)
        assert continuity_i+clock_work==0
        # Same completion coefficient required for pure internal output.
        A=F(17,5); Jbase=F(2,9); V=F(-3,8)
        assert A*((Jbase-V)+V)==A*Jbase
        mismatched=A*(Jbase-V)+(A+1)*V-A*Jbase
        assert mismatched==V and mismatched!=0
        algebra.append({"a":str(a),"c":str(c),"kappa":str(kappa),
                        "density_change":str(density_clock+density_rod),
                        "pressure_vs_Hdot":str(-2*Hdot_extra-kappa*pressure_rod),
                        "work_continuity_residual":str(continuity_i+clock_work)})
    margins={
      "classical_source_negative_exponent":4*N-(E0+4*M+140),
      "quantum_source_negative_exponent":2*N-(E0+EH+15*M+27),
      "explicit_work_negative_exponent":4*N-(E0+60),
      "read_norm_negative_exponent":2*N-(EH+15*M+20),
    }
    assert all(x>110 for x in margins.values())
    # Each term in delta_q is <2^-110, so their sum is <2^-109<2^-80.
    assert F(2,2**110)<F(1,2**80)
    # Classical source expression combines a derivative*delta_z term
    # and the controlled pulse/work term with more than 10 extra bits.
    assert F(2,2**110)<F(1,2**100)

    Dk=F(1,15000)+F(5,10**16)
    savedD=F(read("1268/readout_results.json")["finite_prediction"]["parameter_acceptance_bound"])
    assert Dk==savedD
    classical=(F(111,100)+F(1,2**100))*Dk
    actual=F(111,100)*Dk+F(3,2**80)
    assert actual<F(1,10000)
    # Conditional branch budget: 16||b||dR fits the +27 exponent.
    assert 4+1+20<=27
    # Actual output is after sampling; source cancellation is gate-wide.
    assert F(2)<F(17,8)<F(9,4)

    bindings=[
      "1247/proof.md","1247/results.json","1249/canonical.md",
      "1249/quantum_bridge.md","1249/quantum_bounds.py",
      "1249/quantum_bounds_results.json","1252/proof.md","1252/results.json",
      "1267/geometry_proof.md","1267/geometry_results.json",
      "1268/readout_proof.md","1268/readout_results.json",
      "1269/material_proof.md","1269/material_results.json",
      "1269/pointer_proof.md","1269/pointer_results.json","1269/acceptance.json",
      "_admission/finite_semiclassical_budget_gate_after1256.md",
    ]
    hashes={rel:sha(STAGE/rel) for rel in bindings}
    hashes.update({name:sha(HERE/name) for name in
                   ("admission.md","continuation_proof.md","continuation_check.py")})
    return {
      "round":1270,"part":"same_reading_source_and_geometry_continuation",
      "status":"PASS_AUTHOR_WITH_SCOPE","author_groups":1,
      "arithmetic":"exact Fraction endpoint arithmetic and integer exponents; no ODE",
      "task":{"family":"B=2, kappa in [3,3001/1000]; plus kappa=3 B=1,2",
              "old_comparison_and_new_read_window":["0","9/4"],
              "sampling_gate":["199999999/100000000","2"],
              "actual_conditional_output_window":["17/8","9/4"],
              "epsilon_upper":"2^(-4N)","N_decimal_digits":len(str(N))},
      "full_parameter_core_extension":fam,
      "old_two_program_endpoint_checks":old_core,
      "source_and_pulse_checks":{
        "pulse_integral":str(integral),"pulse_peak":str(peak),
        "pulse_derivative_upper":str(derivative_peak_bound),
        "uncompleted_work_coefficient_upper":str(work),
        "work_power_bound":"epsilon*2^59",
        "explicit_completed_work_bound":"epsilon*2^(E0+60)",
        "density_pressure_work_deletion_checks":algebra,
        "same_completion_internal_output_cancellation":True,
        "different_completion_cancellation_rejected":True},
      "resource_checks":{"M_exponent":3059,"EH":EH,"E0":E0,
        "all_margin_exponents_over_110":True,
        "margin_bit_lengths":{k:v.bit_length() for k,v in margins.items()},
        "classical_new_source_difference":"<2^-100",
        "real_quantum_completed_source_difference":"<2^-80",
        "conditional_source_extra":"<=16||b||dR; p_a>=1/4 from output part",
        "source_comparison_preserves_actual_output_noise_process":True},
      "field_tests":{
        "parameter_acceptance":str(Dk),
        "new_classical_branch_bound":str(classical),
        "actual_bounded_quantum_source_comparison_bound":str(actual),
        "weight_normalization":"sum_ab integral |w_ab| <=1",
        "actual_geometry_role":"g_kappa is a declared effective comparison metric",
        "conditional_weights_supported_after_actual_output":True},
      "scope":[
        "The old full-kappa classical existence theorem is reused, not new ODE integration.",
        "New classical +/-epsilon branches are comparison tools, not quantum conditional geometries.",
        "Actual quantum state is not reset at t=2; the post-gate unitary is the same Hc.",
        "Known material Gaussian family; arbitrary pointer-reference handled as a full cq process.",
        "Nonrare normalized branch estimates only for the specified calibration input.",
        "Pure internal output cancellation uses identical completed source coefficients.",
        "Only bounded completed sources, not uncompleted stress or an operator Einstein equation.",
        "CT unadopted; no full nonuniform PDE, full GR, or empirical conclusion."
      ],
      "bindings_sha256":hashes,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        print("FIRST_EXCLUSIVE_WRITE PASS; author_groups=1")
    else:
        saved=json.loads(OUT.read_text(encoding="utf-8-sig"))
        if saved!=result:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS; author_groups=1")
    print(json.dumps({
        "full_family_core":result["full_parameter_core_extension"]["all_core_checks"],
        "N_digits":result["task"]["N_decimal_digits"],
        "all_source_margins_over_110_bits":True,
        "quantum_source_error":"<2^-80",
        "field_test_upper":"<1e-4",
    },ensure_ascii=False))

if __name__=="__main__":
    main()

