"""1266: finite same-H Gaussian gate, instrument, and bounded-source budget.
Default recomputes and matches results.json without writing.
--write-results creates results.json exclusively.
No wave-function simulation or new classical integration is performed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def calculate():
    old = json.loads((STAGE/"1249/quantum_bounds_results.json").read_text(encoding="utf-8-sig"))
    ledger = old["group_constructive_finite_resource_bound"]["ledger"]
    M = 1 << int(ledger["M_power_of_two_exponent"])
    EH = int(ledger["Fourier_H_exponent"])
    E0 = int(ledger["Fourier_menu_base_exponent"])
    EF = int(ledger["Fourier_transport_exponent_decimal"])
    N = int(ledger["integer_n_power_exponent_decimal"])
    assert M == 1 << 3059
    assert (EH,E0) == (3093,14008)
    assert EF == 95928+100*M
    ceildiv = lambda a,b: (a+b-1)//b
    assert N == 1+max(ceildiv(EF+42,3),ceildiv(EH+EF+42,6))
    assert 3*N-EF > 42 and 6*N-EH-EF > 42

    eta = F(1,10**9)
    center_error = F(3,10**15)
    assert 2*eta-center_error > eta
    intervals = [(F(15,16),F(17,16)),(F(31,16),F(33,16))]
    translated = sorted((lo+d,hi+d) for d in (0,4,8) for lo,hi in intervals)
    assert all(translated[i][1] < translated[i+1][0] for i in range(5))
    assert max(hi for _,hi in translated) == F(161,16)
    assert F(161,16)<61<64
    for center,(lo,hi) in zip((F(1),F(2)),intervals):
        assert center-lo == hi-center == F(1,16)

    ta,tau = F(17,8),F(1,128)
    tb,tf = ta+tau,ta+2*tau
    assert (tb,tf)==(F(273,128),F(137,64))
    assert F(2)<ta<tb<tf<F(9,4)
    assert F(21,400)<ta-2
    gpeak=F(35,16)/tau
    assert gpeak==280
    assert gpeak*F(11,7)==440

    # c_G^2 = 15/(144*2) = 5/96.  This is the exact Gaussian sixth-moment constant.
    cG2=F(5,96)
    assert cG2<F(1,16)
    assert cG2*2*F(9,4)**2<1
    assert 14<16
    # log(2)>1/2: exp((9/4)M)<2^(5M); no floating exponentiation.
    assert F(9,4)<F(5,2)
    beta_coefficient=(16**2+eta**-2)/2
    assert beta_coefficient < 2**60
    assert 2*F(22,7)**2 < 25
    hres_exp=EH+15*M+12
    gres_exp=5*M+37
    assert hres_exp>gres_exp
    d_factor_exp=EH+15*M+16
    assert hres_exp+1 < d_factor_exp
    assert 5*M+9 < d_factor_exp
    post_factor_exp=EH+15*M+18

    # epsilon <= 2^(-4N), sqrt(epsilon) <= 2^(-2N).
    exponents={
        "state_error_upper_power2":d_factor_exp-2*N,
        "post_instrument_half_trace_upper_power2":post_factor_exp-2*N,
        "bounded_source_error_upper_power2":E0+EH+15*M+22-2*N,
        "gaussian_blank_failure_upper_power2":10*M+15-4*N,
        "gaussian_blank_and_K_failure_upper_power2":10*M+68-4*N,
        "gate_mass_upper_power2":9-4*N,
    }
    thresholds={
        "state_error_upper_power2":-100,
        "post_instrument_half_trace_upper_power2":-100,
        "bounded_source_error_upper_power2":-80,
        "gaussian_blank_failure_upper_power2":-200,
        "gaussian_blank_and_K_failure_upper_power2":-200,
        "gate_mass_upper_power2":-80,
    }
    for key,value in exponents.items(): assert value<thresholds[key],key
    # Source sum: 2 F0 dpost + first Gaussian moment + AW ordering error.
    leading_source_exp=E0+post_factor_exp+2
    assert leading_source_exp>=E0+5*M+4
    assert leading_source_exp>=E0
    assert leading_source_exp+2==E0+EH+15*M+22
    # 2d + beta_B < 2^-98; clock norm bound gives 625*2^-160 < 2^-150.
    assert 2*F(1,2**100)+F(1,2**200)<F(1,2**98)
    assert (F(1,2**79)/F(2,25))**2 == F(625,2**160)
    assert 625<2**10

    deps=[
        HERE/"admission.md", HERE/"proof.md", Path(__file__).resolve(),
        STAGE/"1249/canonical.md", STAGE/"1249/quantum_bridge.md",
        STAGE/"1249/quantum_bounds_results.json",
        STAGE/"1259/proof.md",
        STAGE/"_admission/finite_semiclassical_budget_gate_after1256.md",
        STAGE/"_admission/causal_block_after1265/proof.md",
        STAGE/"_admission/causal_block_after1265/results.json",
    ]
    return {
        "round":1266,"status":"PASS_AUTHOR_WITH_SCOPE","author_groups":1,
        "arithmetic":"exact integers and Fraction; no new ODE or quantum wave-function simulation",
        "source_bindings":{str(p.relative_to(STAGE)).replace("\\","/"):sha(p) for p in deps},
        "inherited_constants":{"M_definition":"2^3059","EH":EH,"E0":E0,
            "EF_definition":"95928+100*2^3059","N_decimal":str(N)},
        "gate_geometry":{"ta":str(ta),"duration":str(tau),"tb":str(tb),"flag_end":str(tf),
            "normalized_pulse_peak":str(gpeak),"B_max":str(F(161,16)),
            "K_safe_margin":str(eta),"same_chi_core":True},
        "gaussian_remainder":{"cG_squared":str(cG2),"S_upper":"2^(5*M+4)",
            "two_component_coefficient_sum":"at most sqrt(2); no orthogonality assumption",
            "beta_coefficient":str(beta_coefficient),
            "actual_generator":"Op_epsilon(chi*h) + epsilon*g*G; no holding completion replacement"},
        "authoritative_negative_binary_exponents":{k:str(v) for k,v in exponents.items()},
        "reported_strict_bounds":{"state_error":"2^-100","post_instrument_half_trace":"2^-100",
            "bounded_source_and_tail_expectation_error":"2^-80",
            "incorrect_flag_probability":"2^-98",
            "gate_time_clock_spectral_failure_J_le_0.01":"2^-150"},
        "scope":[
            "two original specified coherent preparations, not arbitrary unknown initial states",
            "full finite continuous gate and all terminal instrument outcomes retained",
            "bounded completed source menus; no unbounded stress tail or psi^2 domain claim",
            "J_read+V=-Hc is a common reduced constraint ledger, not a physical energy identification",
            "clock positivity is a specified-time state-dependent probability bound",
            "actual block boundary, clock trigger and complete new participant sources remain explicit inputs",
            "no new recognition axiom, no general relativity completion claim",
        ],
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    current=calculate()
    destination=HERE/"results.json"
    if args.write_results:
        with destination.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(current,f,ensure_ascii=False,indent=2,sort_keys=True);f.write("\n")
        receipt="FIRST_WRITE_EXCLUSIVE"
    else:
        saved=json.loads(destination.read_text(encoding="utf-8"))
        assert saved==current,"Saved results mismatch"
        receipt="READ_ONLY_MATCH"
    print(json.dumps({"round":1266,"status":current["status"],"receipt":receipt,
        "author_groups":1,"result_sha256":sha(destination)},ensure_ascii=False,indent=2))
if __name__=="__main__": main()

