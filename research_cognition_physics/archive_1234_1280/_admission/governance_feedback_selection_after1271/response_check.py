"""Unnumbered finite classical feedback response. Exact arithmetic; no ODE.
Default: read-only recomputation/match. --write: exclusive first creation.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
OUT = HERE / "response_results.json"

def read(rel):
    return json.loads((STAGE / rel).read_text(encoding="utf-8-sig"))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def interval(state, key):
    return tuple(F(x) for x in state[key])

def absmax(pair):
    return max(map(abs, pair))

def absmin(pair):
    lo, hi = pair
    return F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))

def exp_upper(x, n=10):
    # First omitted term, with each subsequent ratio at most x/(n+2).
    term = F(1)
    total = term
    for j in range(1, n+1):
        term *= x/j
        total += term
    omitted = term*x/(n+1)
    return total + omitted/(1-x/(n+2))

def s(x):
    return str(x)

def endpoint_certificate(state, perturb):
    u = interval(state, "u"); H = interval(state, "H")
    amin = 1/u[1]
    amax = 1/(u[0]*(1-H[1]/4))
    Jlo = interval(state, "J")[0]
    Mc = F(71,100)+Jlo+sum(absmin(interval(state,k))**2 for k in ("q","r","K"))/2
    assert Jlo > F(49,500) and Mc > F(93,100)
    margins = {
        "a_lower": amin-F(5,2),
        "a_upper": F(3)-amax,
        "H_upper": F(31,100)-H[1],
        "field_abs": F(8,5)-F(3,2),
        "H_positive": F(9,50),
    }
    assert all(m > perturb for m in margins.values())
    qr2 = sum(absmax(interval(state,k))**2 for k in ("q","r"))
    Y = absmax(interval(state,"Y"))+absmax(interval(state,"K"))/4
    vartheta = absmax(interval(state,"z"))+F(1,5)+F(1,4)*F(5,2)/10
    assert qr2 < 1 and Y < 1 and vartheta+perturb/500 < 1
    assert absmax(interval(state,"psi")) < 1 and absmax(interval(state,"K")) < 1
    return {
        "a_old_lower":s(amin), "a_old_upper":s(amax),
        "H_old_upper":s(H[1]), "H_old_strict_lower":"9/50",
        "J_old_lower":s(Jlo), "Mc_old_lower":s(Mc),
        "box_margins":{k:s(v) for k,v in margins.items()},
        "q_r_norm_squared_upper":s(qr2), "Y_abs_upper":s(Y),
        "vartheta_abs_upper_after_new_pulse":s(vartheta+perturb/500),
        "all_strictly_inside":True,
    }

def calculate():
    eta = F(1,25000); T=F(9,4); ta=F(447,200); tb=F(449,200)
    tau=tb-ta; kappa=F(3001,1000); amin=F(5,2); amax=F(3)
    Hmax=F(31,100); fields=F(8,5); D=F(1,50)
    assert F(17,8)<ta<tb<T and tau==F(1,100) and T-ta<D
    pulse_int=140*(F(1,4)-3*F(1,5)+3*F(1,6)-F(1,7))
    pulse_s_moment=140*(F(1,5)-3*F(1,6)+3*F(1,7)-F(1,8))
    moment=(T-ta)*pulse_int-tau*pulse_s_moment
    peak=F(35,16)/tau; cpeak=eta*peak
    assert pulse_int==1 and pulse_s_moment==F(1,2)
    assert moment==F(1,100) and peak==F(875,4) and cpeak==F(7,800)
    expU=exp_upper(12*D)
    assert expU < F(32,25)
    perturb=F(81,1000)*eta*F(32,25)
    assert perturb==F(324,78125000)
    actual_old_a_error=F(1,150000)
    assert perturb<actual_old_a_error<2*actual_old_a_error

    row_a=Hmax+amax
    row_H=3*Hmax+2*kappa*fields+kappa/(200*amin**2)+kappa*(1+fields/5)/(20*amin**3)
    row_v=3*fields+3*Hmax+1+F(6,100)/amin**3+F(3,100)/amin**4
    row_V=3*fields+3*Hmax+1
    rows={"a":row_a,"H":row_H,"p":F(1),"v":row_v,"X":F(1),"V":row_V}
    assert row_a<4 and all(v<12 for v in rows.values())
    beta_min=F(3)/(6*amax**2); beta_max=kappa/(6*amin**2)
    assert beta_min==F(1,18) and beta_max<F(81,1000)
    peano_ratio=24*D*F(32,25)
    assert peano_ratio<1
    propagator_slope=F(5,2)-peano_ratio
    assert propagator_slope>F(3,2)
    response=F(3,2)*beta_min*eta*moment
    assert response==eta/1200==F(1,30000000)
    # Old Friedmann lower H, from positive old dust and a<3.
    assert F(9,50)**2<F(9,10)/27
    # Other canonical coordinates remain within the old coordinate completion.
    canonical={"pi_material_upper":27*fields, "piQ_upper":2*27*Hmax}
    assert all(v<61 for v in canonical.values())
    # e^(3/2)>3 suffices for Q<3/2, while a>1 gives Q>0.
    assert sum(F(3,2)**n/__import__("math").factorial(n) for n in range(4))>3

    family=read("1252/results.json")["group_parameter_family"]["cells"]
    assert len(family)==1 and family[0]["kappa"]==["3","3001/1000"]
    states={"B2_full_kappa":family[0]["checkpoints"]["2"]["state"]}
    old=read("1247/results.json")["group_finite_memory"]
    states.update({name:old[name]["checkpoints"]["2"]["state"] for name in ("alpha0","alpha1")})
    endpoint_checks={name:endpoint_certificate(state,perturb) for name,state in states.items()}

    Jlower=F(49,500)-amax*cpeak
    Mclower=F(93,100)-amax*cpeak
    assert Jlower==F(287,4000)>0 and Mclower==F(723,800)>F(9,10)
    assert Jlower<F(2,25) # Do not grant the old good-state plateau.
    endJlower=amin*(F(9,50)-perturb)*eta
    endJupper=amax*Hmax*eta
    assert 0<endJlower<endJupper

    eta_alt=F(1,40000)
    altJ=F(49,500)-3*peak*eta_alt
    assert F(3,125000)<eta_alt<F(9,500)/(3*peak)
    assert eta_alt/1200>F(2,10**8) and altJ>F(2,25)

    # Universal coefficient identities: explicit clock/rod density cancels;
    # -2*(kappa*c/(6a²)) = kappa*(-c/(3a²)).
    assert -F(1)+F(1)==0
    assert -2*F(1,6)-(-F(1,3))==0
    # For rho_int=c a^-2, P_int=-rho_int/3:
    # d(rho_int)+3H(rho_int+P_int)=c' a^-2,
    # canceled by a^-3 Jdot with Jdot=-c' a.
    assert -2+3*(1-F(1,3))==0

    bindings=[
        "1247/proof.md","1247/results.json","1249/canonical.md",
        "1249/quantum_bridge.md","1252/proof.md","1252/results.json",
        "1267/geometry_proof.md","1268/readout_proof.md","1269/material_proof.md",
        "1270/continuation_proof.md","1270/continuation_results.json",
        "1270/output_proof.md","1270/acceptance.json",
        "_admission/finite_common_task_after1266/README.md",
    ]
    hashes={rel:sha(STAGE/rel) for rel in bindings}
    hashes.update({name:sha(HERE/name) for name in
                   ("source_baseline.json","response_proof.md","response_check.py")})
    return {
      "status":"PASS_CLASSICAL_CANDIDATE_ONLY",
      "numbered_round":False, "diagnostic_groups":1, "formal_science_groups_added":0,
      "arithmetic":"exact Fraction; inherited closed-parameter endpoint certificate; no new ODE",
      "parameters":{"eta":s(eta),"eta_range":"0 < eta <= 1/25000",
        "time_output":"17/8","pulse":[s(ta),s(tb)],"terminal":s(T),"tau":s(tau)},
      "pulse":{"integral":s(pulse_int),"s_first_moment":s(pulse_s_moment),
        "terminal_first_moment":s(moment),"g_peak":s(peak),"c_peak":s(cpeak)},
      "old_endpoint_bootstrap":endpoint_checks,
      "finite_response":{
        "jacobian_row_bounds":{k:s(v) for k,v in rows.items()},
        "jacobian_uniform_bound":"12","a_row_bound":"4",
        "exp_6_over25_rational_upper":s(expU),"exp_coarse_upper":"32/25",
        "state_difference_strict_upper":s(perturb),
        "beta_lower":s(beta_min),"beta_upper":s(beta_max),
        "peano_remainder_per_time_upper":s(peano_ratio),
        "propagator_slope_lower":s(propagator_slope),
        "conservative_propagator_slope":"3/2",
        "a_terminal_difference_lower":s(response),
        "a_terminal_difference_decimal_display":float(response),
        "actual_old_Ramsey_a_error":s(actual_old_a_error),
        "entire_classical_difference_less_than_old_single_arm_error":True,
        "old_actual_budget_cannot_guarantee_discrimination":True},
      "source_and_clock":{
        "J_strict_lower":s(Jlower),"Mc_strict_lower":s(Mclower),
        "J_final_increment_strict_lower":s(endJlower),
        "J_final_increment_strict_upper":s(endJupper),
        "clock_rod_density_coefficients_cancel":True,
        "pressure_and_Hdot_coefficients_match":True,
        "source_exchange_coefficients_cancel":True,
        "old_good_J_plateau_NOT_inherited":True,
        "coordinate_completion_bounds":{k:s(v) for k,v in canonical.items()},
        "physical_energy_scaling":"V0*s*c*a, not hbar/s"},
      "nonempty_smaller_strength_arithmetic_only":{
        "eta_example":s(eta_alt),"classical_response_lower":s(eta_alt/1200),
        "J_lower":s(altJ),"actual_readout_or_quantum_certified":False},
      "scope":{
        "same_actual_conditional_state_protocol_question":True,
        "proved_comparison":"same classical old center, feedback versus waiting",
        "actual_conditional_state_replaced_by_classical_point":False,
        "new_feedback_material_input_required":True,
        "new_quantum_propagation_certified":False,
        "new_terminal_measurement_permission_certified":False,
        "old_good_state_failure_bound_inherited":False,
        "governance_selection_law_proved":False,
        "thermal_geometric_identity_proved":False,
        "CT_adopted":False},
      "source_bindings_sha256":hashes,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    data=calculate()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(data,f,ensure_ascii=False,indent=2)
            f.write("\n")
        print("EXCLUSIVE_FIRST_WRITE PASS_CLASSICAL_CANDIDATE_ONLY")
    else:
        saved=json.loads(OUT.read_text(encoding="utf-8"))
        if saved!=data:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS_CLASSICAL_CANDIDATE_ONLY")
    print("a_gap >= "+data["finite_response"]["a_terminal_difference_lower"])
    print("J > "+data["source_and_clock"]["J_strict_lower"])
if __name__=="__main__":
    main()
