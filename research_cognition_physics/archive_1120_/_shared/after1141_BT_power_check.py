"""Nonformal after1141 qualification: exact BT power migration budgets.
Default is read-only. Only --write-results writes this qualification's JSON.
No new round, scientific calibration, continuum simulation, or BT balance proof.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RESULT = HERE / "after1141_BT_power_results.json"
FROZEN = {'1128/proof.md': 'fc9963b8f77d0a60e2233acbcc68e02e583af345a5653d1d4431ea619be0fcd2', '1129/proof.md': 'c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa', '1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d', '1134/proof.md': '4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666', '1137/proof.md': 'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb', '1138/proof.md': '5668334517745145330f0504905c11bdaf70b96a374080180177c0ddf5272110', '1141/proof.md': 'd71543b186660711ac1602a15495c5700ba2430577366d3057f70249fb650c96', '1141/results.json': '62d67da39ec6e6c016f5d6150760ed21b84729a96dcbc0ee5b81f9891919376a', '1141/acceptance.json': 'fb78daab0c9d9d75102b74d2770f9919494d0cae230f82a72e004aa9a12a1df7', '_shared/after1137_local_balance_preflight.md': '95915f871317e8638d23c3ee29cec779a0ff52aad5f76873a0ff896f6a2fcbe0', '_shared/after1141_BT_power_qualification.md': '3269ee80aa0ce3ad55db4bf621f5d312f7628a7c88a671a0b2e187ce0a2d13f6'}


def exact(x):
    x = F(x)
    with localcontext() as ctx:
        ctx.prec = 50
        decimal = str(Decimal(x.numerator) / Decimal(x.denominator))
    return {"numerator": str(x.numerator), "denominator": str(x.denominator),
            "decimal": decimal}


def from_saved(x):
    return F(int(x["numerator"]), int(x["denominator"]))


def run():
    observed = {p: hashlib.sha256((BASE / p).read_bytes()).hexdigest()
                for p in FROZEN}
    assert observed == FROZEN, "Frozen qualification/source file changed"
    acceptance = json.loads((BASE / "1141/acceptance.json").read_text(encoding="utf-8"))
    assert acceptance["status"] == "PASS" and acceptance["round"] == 1141
    assert acceptance["core_sha256"]["proof.md"] == FROZEN["1141/proof.md"]
    assert acceptance["core_sha256"]["results.json"] == FROZEN["1141/results.json"]

    A, M, eps, sigma, dc, T = F(10**34), F(10**80), F(1,10**5), F(1,10**30), F(1,10**6), F(5,4)
    lam, kap, Cs, D = F(11,7), F(11,14), F(33536000000000000000), F(12,10)*10**21
    v = 16*lam*Cs+F(558,7)
    mstar = 3*M-v
    assert mstar > 2*M
    CY = F(24780800000000000000000,49)+16316
    p1 = 268+94*A
    p2 = 28000+188*268*A+12000*A*A
    A1 = 1/sigma+p1
    A2 = 1/sigma**2+2*p1/sigma+p2
    B2_mom = p1/sigma+p2
    dBT = T*A2/(2*mstar)+T*(A2+2*B2_mom)/(2*M)+CY*A1*T*T/(2*M)
    assert dBT < F(26,10**9)

    F_second, S_second = F(32768), F(128)/dc**2
    B2_source = lam*F_second*Cs+42*S_second+(lam+kap)*F_second
    q = 8*eps**-3
    ccr_half = q*B2_source/(2*A*A)
    assert B2_source < F(1727,1000)*10**24
    assert ccr_half < F(6908,1000)*F(1,10**29)

    Bstar = F(1074,1000)*10**71
    Ho_mean_upper = 1073*A*A+300+1/(M*sigma**2)+2*v
    assert Ho_mean_upper < Bstar
    Vw, nu_upper = F(524), 100*eps**-4
    C = 47*A+230*eps**-2
    assert 2*Bstar < (47*A)**2
    assert Vw*nu_upper < (230*eps**-2)**2
    assert 2*Bstar+Vw*nu_upper < C*C
    NP = D/A*C+ccr_half
    point_error = 2*dBT*NP
    integral_error = T*point_error
    assert NP < F(5641,1000)*10**22
    assert point_error < 3*10**15
    assert integral_error < 4*10**15

    r1141 = json.loads((BASE / "1141/results.json").read_text(encoding="utf-8"))
    old_budget = r1141["groups"]["exact_rational_budget"]
    original_D = from_saved(old_budget["force_column_D"])
    assert original_D < D
    original_C = from_saved(old_budget["free_column_Cw"])
    assert original_C == C
    ho_integral = C*original_D/A*T+(original_D/A*T)**2/2
    assert ho_integral == from_saved(old_budget["absolute_work_integral_bound"])
    actual_bare_integral = ho_integral+integral_error
    assert actual_bare_integral < F(72,10)*10**22 < 10**23

    external_g_exact = T*v*A2/(2*mstar*mstar)
    external_g_safe = T*v*A2/(2*(2*M)**2)
    assert external_g_exact < external_g_safe < F(1581,1000)*F(1,10**68)

    vals = dict(A=A,M=M,epsilon=eps,sigma_Y=sigma,delta_c=dc,T=T,
                v_all=v,m_star=mstar,D_safe=D,D_from_round1141=original_D,
                F_second_upper=F_second,S_second_upper=S_second,
                B2_source=B2_source,CCR_half=ccr_half,
                H_actual_source_mean_bound=Bstar,Ho_source_mean_bound=Ho_mean_upper,
                common_source_mean_bound=Bstar,common_pi_column_C=C,
                common_power_graph_NP=NP,p1=p1,p2=p2,A1=A1,A2=A2,
                B2_momentum=B2_mom,dBT_short=dBT,
                pointwise_power_migration=point_error,
                integrated_power_migration=integral_error,
                Ho_absolute_power_integral=ho_integral,
                actual_H_same_bare_power_integral_bound=actual_bare_integral,
                external_g_interaction_state_exact=external_g_exact,
                external_g_interaction_state_safe=external_g_safe)
    return {
        "qualification": "after1141_BT_power",
        "status": "PASS",
        "formal_round_reserved": False,
        "scientific_calibration_increment": 0,
        "new_cognitive_axioms": 0,
        "groups": {
            "frozen_source_identity": {"status":"PASS","sha256":observed,
                                       "round1141_acceptance_verified":True},
            "exact_fraction_budgets": {"status":"PASS",
                                       "values":{k:exact(x) for k,x in vals.items()}},
            "scope": {"status":"PASS",
                      "original_normal_inputs_and_quantum_tails_retained":True,
                      "all_QR_uniform":True,
                      "fixed_lab_region_common_LLP_transport":True,
                      "actual_H_same_bare_power_expectation_migration":True,
                      "actual_BT_regional_balance_proved":False,
                      "actual_Ho_regional_energy_bound_migrated":False,
                      "new_energy_instrument_permission":False,
                      "continuum_physics_simulation":False,
                      "external_g_bound_replaces_free_reference_analysis":False}
        }
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    data = run()
    if args.write_results:
        RESULT.write_text(json.dumps(data, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        assert RESULT.exists(), "Saved result missing; generate only with --write-results"
        saved = json.loads(RESULT.read_text(encoding="utf-8"))
        assert data == saved, "Saved result does not match exact recomputation"
    vals=data["groups"]["exact_fraction_budgets"]["values"]
    print(json.dumps({"status":"PASS","mode":"write-results" if args.write_results else "read-only",
                      "qualification":"after1141_BT_power","groups":len(data["groups"]),
                      "frozen_sources":len(FROZEN),"formal_round_reserved":False,
                      "dBT_short":vals["dBT_short"]["decimal"],
                      "power_migration":vals["pointwise_power_migration"]["decimal"],
                      "integrated_migration":vals["integrated_power_migration"]["decimal"],
                      "actual_BT_regional_balance_proved":False}, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
