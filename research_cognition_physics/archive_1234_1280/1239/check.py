"""1239: complete homogeneous conversion source -> old finite record budget.

Fraction inequalities certify the analytic bounds. RK4/Hermite values are
diagnostics, not a substitute for the pressure-work/ODE comparison proof.
Default replay is read-only; --write creates results.json exactly once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
FUEL = 1e-5
END = .82
EMISSIONS = np.array([0., .331, .728])
DISTANCE = .03


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def rat(value):
    return {"fraction": str(value), "decimal": float(value)}


def det(e, r):
    return (e[1]-e[0])*(r[2]-r[0])-(e[2]-e[0])*(r[1]-r[0])


def exact_certificates():
    fuel = F(1, 100000)
    d = F(3, 100)
    max_e = F(91, 125)
    # All factors are independent domain endpoints, not fitted diagnostics.
    terms = [4*d*max_e/3, 4*d*d*F(6, 5)/3,
             4*d*d*max_e/9, 16*d**3/27]
    derivative = sum(terms)
    need(derivative == F(24115, 781250), "four rational derivative bounds")
    for E0 in (F(1), F(4, 5)):
        need(E0-fuel > F(4, 5)**2, "both lower h exceed .8")
        need(E0 <= 1 and E0 > fuel, "positive finite complete energy window")
    h_shift = 5*fuel/8
    delta_g = derivative*h_shift
    need(delta_g == F(4823, 25000000000), "uniform geometric return budget")
    delta_D = 2*max_e*delta_g
    old_path = HERE.parent/"1236/results.json"
    old = json.loads(old_path.read_text(encoding="utf-8"))
    old_gap = F(old["group_3"]["analytic_Peano_gap_lower"]["fraction"])
    old_D = F(old["group_2"]["D"]["fraction"])
    coefficient = F(old["group_2"]["linear_error_coefficient"]["fraction"])
    need(old_gap == F(40383498710737, 1194393600000000000), "frozen source gap")
    need(old_D == -F(5529, 12500000), "frozen six-record determinant")
    need(2*delta_D < old_gap, "all conversion histories retain ideal source separation")
    eta = F(1, 1250000)
    total = eta + delta_g
    need(total < F(1, 1000000), "remaining readout contract plus derived geometry")
    flat_error = coefficient*(2*total) + 8*(2*total)**2
    source_error = 10*total + 16*total**2
    need(flat_error < -old_D, "arbitrary flat ideal two-box exclusion")
    need(source_error < old_gap, "old source gap survives both complete budgets")
    e = [F(0), F(331, 1000), F(728, 1000)]
    r = [(q+F(1, 50))**3-1 for q in (F(1), F(11, 10), F(6, 5))]
    need(det(e, r) == old_D and max(r) < F(82, 100), "same finite menu and proper-time window")
    # Exact polynomial identity for constraint preservation. The erroneous
    # source contribution is independent of H, U, and scale-factor samples.
    for a, H, E, U, q in ((F(7, 5), F(3, 7), F(9, 10), F(1, 1000), F(1, 2000)),
                           (F(6, 5), F(2, 5), F(4, 5), F(1, 2000), F(1, 4000))):
        Hdot = -F(2, 3)*(E+U/3)/a**3
        correct = 2*H*Hdot - F(4, 9)*(-H*U-3*H*E)/a**3
        wrong = 2*H*Hdot - F(4, 9)*(q-H*U-3*H*E)/a**3
        need(correct == 0 and wrong == -4*q/(9*a**3), "exact complete/missing-fuel source identity")
    return dict(status="PASS",certificate_type="analytic inequalities with exact rational endpoints",
                fuel_per_cell=rat(fuel),derivative_terms=[rat(x) for x in terms],
                return_derivative_bound=rat(derivative),h_change_bound=rat(h_shift),
                per_return_geometry_budget=rat(delta_g),per_branch_D_geometry_budget=rat(delta_D),
                ideal_source_gap_after_conversion=rat(old_gap-2*delta_D),
                remaining_readout_contract=rat(eta),total_six_record_budget=rat(total),
                flat_two_box_error=rat(flat_error),flat_certified_margin=rat(-old_D-flat_error),
                source_two_branch_error=rat(source_error),source_certified_margin=rat(old_gap-source_error),
                all_conversions_with_nonnegative_fuel_covered=True,
                actual_remaining_readout_accuracy_obtained=False)


def rhs(state, E0, wrong=False):
    # Simultaneous Raychaudhuri + matter integration; H is never reset from C=0.
    a, H, fuel, U, work, conformal, bad_integral = state
    E = E0-FUEL+fuel+U
    transfer = fuel  # gamma=1, no external time-dependent program
    return np.array([a*H, -(2/3)*(E+U/3)/a**3,
                     0. if wrong else -transfer, transfer-H*U,
                     H*U, 1/a, -4*transfer/(9*a**3) if wrong else 0.])


def integrate(E0, steps, wrong=False):
    ts = np.linspace(0., END, steps+1)
    values = np.empty((steps+1, 7))
    values[0] = [1., (2/3)*np.sqrt(E0), FUEL, 0., 0., 0., 0.]
    dt = END/steps
    for j in range(steps):
        v = values[j]
        k1 = rhs(v, E0, wrong)
        k2 = rhs(v+dt*k1/2, E0, wrong)
        k3 = rhs(v+dt*k2/2, E0, wrong)
        k4 = rhs(v+dt*k3, E0, wrong)
        values[j+1] = v+dt*(k1+2*k2+2*k3+k4)/6
    return ts, values


def conformal_at(t, ts, values):
    # Cubic Hermite interpolation with exact ODE endpoint derivative 1/a.
    i = min(int(np.searchsorted(ts, t, side="right")-1), len(ts)-2)
    i = max(0, i)
    dt = ts[i+1]-ts[i]
    u = (t-ts[i])/dt
    f0, f1 = values[i, 5], values[i+1, 5]
    d0, d1 = 1/values[i, 0], 1/values[i+1, 0]
    return ((2*u**3-3*u*u+1)*f0 + (u**3-2*u*u+u)*dt*d0
            +(-2*u**3+3*u*u)*f1 + (u**3-u*u)*dt*d1)


def returns(ts, values):
    result = []
    for e in EMISSIONS:
        target = conformal_at(float(e), ts, values)+2*DISTANCE
        lo, hi = float(e), END
        need(conformal_at(hi, ts, values)>target, "finite return exists inside domain")
        for _ in range(54):
            mid = (lo+hi)/2
            if conformal_at(mid, ts, values)<target:
                lo = mid
            else:
                hi = mid
        result.append((lo+hi)/2)
    return np.array(result)


def dust_returns(h):
    z = 1+h*EMISSIONS
    d = DISTANCE
    return EMISSIONS+2*d*z**(2/3)+(4/3)*h*d*d*z**(1/3)+(8/27)*h*h*d**3


def nonzero_conversion():
    rows = []
    for E0 in (1., .8):
        ts, v = integrate(E0, 4096)
        ts0, v0 = integrate(E0, 2048)
        a, H, fuel, U, work, conformal, _ = v.T
        E = E0-FUEL+fuel+U
        C = H*H-(4/9)*E/a**3
        source = np.max(np.abs(E+work-E0))
        split = np.max(np.abs(fuel+U+work-FUEL))
        fuel_exact = np.max(np.abs(fuel-FUEL*np.exp(-ts)))
        lower = (1+np.sqrt(E0-FUEL)*ts)**(2/3)
        upper = (1+np.sqrt(E0)*ts)**(2/3)
        violations = max(0., float(np.max(lower-a)), float(np.max(a-upper)))
        rr = returns(ts, v)
        coarse = returns(ts0, v0)
        ref = dust_returns(np.sqrt(E0))
        low = dust_returns(np.sqrt(E0-FUEL))
        need(min(H)>0 and min(fuel)>0 and U[-1]>0 and work[-1]>0,
             "nonzero positive conversion and pressure work")
        need(min(E)>=E0-FUEL-1e-14 and max(E)<=E0+1e-14, "finite complete energy interval")
        need(source<3e-14 and split<3e-18 and fuel_exact<3e-18, "complete source accounting")
        need(np.max(np.abs(C))<5e-14, "independently integrated constraint preservation")
        need(violations<5e-14, "analytic scale-factor envelopes")
        need(np.max(np.abs(rr-coarse))<5e-12, "two-grid return diagnostic stability")
        need(min(ref-rr)>0 and min(rr-low)>0, "strict nonzero return shift inside analytic bracket")
        need(max(ref-rr)<1.9292e-7, "analytic budget contains autonomous example")
        rows.append(dict(E0=E0,gamma=1.0,steps=[2048,4096],
                         final_fuel=float(fuel[-1]),final_comoving_radiation=float(U[-1]),
                         final_pressure_work=float(work[-1]),final_comoving_energy=float(E[-1]),
                         source_account_max_error=float(source),fuel_split_max_error=float(split),
                         exact_decay_max_error=float(fuel_exact),constraint_max_error=float(np.max(np.abs(C))),
                         envelope_violation=float(violations),returns=rr.tolist(),
                         original_dust_returns=ref.tolist(),lower_envelope_returns=low.tolist(),
                         derived_geometry_shifts=(ref-rr).tolist(),
                         return_grid_difference=float(np.max(np.abs(rr-coarse))),
                         converted_D=float(det(EMISSIONS,rr)),unconverted_D=float(det(EMISSIONS,ref))))
    need(rows[1]["converted_D"]-rows[0]["converted_D"]>3.3e-5,
         "actual diagnostic source contrast retains conservative old gap")
    return dict(status="PASS",interpretation="diagnostic only; exact bounds proved independently",cases=rows)


def missing_fuel_deletion():
    ts, v = integrate(1., 4096, wrong=True)
    a, H, fuel, U, work, _, integrated_bad = v.T
    E = 1-FUEL+fuel+U
    C = H*H-(4/9)*E/a**3
    need(np.max(np.abs(C-integrated_bad))<5e-14,
         "wrong-source constraint drift equals complete analytic residual")
    need(C[-1]<-1e-6 and integrated_bad[-1]<-1e-6,
         "finite wrong-source residual is nonzero")
    # Global comoving energy + pressure work now gains the unchanged fuel rate.
    added = E+work-1
    need(np.max(np.abs(added-FUEL*ts))<3e-14, "explicit double-counted fuel energy")
    need(np.max(np.abs(fuel-FUEL))==0., "the deletion really leaves fuel untouched")
    return dict(status="PASS",steps=4096,
                final_constraint_residual=float(C[-1]),
                independently_integrated_residual=float(integrated_bad[-1]),
                residual_identity_max_error=float(np.max(np.abs(C-integrated_bad))),
                spurious_comoving_energy_gain=float(added[-1]),
                source_double_counting_max_error=float(np.max(np.abs(added-FUEL*ts))),
                refitting_constraint_would_move_failure_to_continuity=True)


def calculate():
    result = dict(status="PASS",science_groups=3,
                  group_1=exact_certificates(),group_2=nonzero_conversion(),
                  group_3=missing_fuel_deletion(),
                  scope=dict(CT="candidate hypothesis, not adopted",
                             effective_geometric_branch_and_initial_constraint_are_inputs=True,
                             autonomous_covariant_transfer_is_new_effective_input=True,
                             only_homogeneous_geometry_component_certified=True,
                             actual_directed_signal_and_recoil=False,
                             actual_trigger_record_and_quantum_code=False,
                             remaining_readout_eta_is_conditional=True,
                             independent_horizon_entropy_matching=False,
                             full_GR_proved=False))
    files = [HERE/name for name in ("admission.md","proof.md","check.py","NEXT.md")]
    files += [HERE.parent/"1236"/name for name in ("proof.md","results.json")]
    files += [RESEARCH/"archive_956_989/research_note_964.md",
              RESEARCH/"archive_1161_1200/1188/proof.md",HERE.parent/"1235/proof.md"]
    result["source_sha256"]={str(p.relative_to(RESEARCH)).replace("\\","/"):
                             hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=HERE/"results.json"
    if args.write:
        if path.exists():
            raise FileExistsError("results.json exists; default is read-only replay")
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        mode="CREATED"
    else:
        need(result==json.loads(path.read_text(encoding="utf-8")),"read-only results/source-hash match")
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,science_groups=3,
                          geometry_return_bound=result["group_1"]["per_return_geometry_budget"]["decimal"],
                          source_margin=result["group_1"]["source_certified_margin"]["decimal"])))


if __name__=="__main__":
    main()
