"""1230: fixed finite Dirac geometry and actual readable mutual information.
Exact Fraction inequalities certify the margins; NumPy checks only finite CAR.
Default/--verify is read-only. --write-results creates, never overwrites, receipt.
No full-field simulation, sampled experiment, cutoff, area law, or Unruh device.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = HERE / "results.json"
GAMMA = F(1, 10**8)
B = F(1, 100)
D = F(1)
DELTA = F(1, 10**24)
SIGMA = F(1, 10**120)
M = F(10**260)
T = F(2)
C = F(3000)
A0 = F(2**28)
SOURCES = [
    "03_认知操作五公理与三维关系空间_阶段论文.md",
    "04_认知时空操作公理与洛伦兹运动学_阶段论文.md",
    "05_认知动力学公理与狭义相对论基础动力学_阶段论文.md",
    "archive_1120_1158/1126/proof.md",
    "archive_1219_/1227/proof.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1228/check.py",
    "archive_1219_/1228/results.json",
    "archive_1219_/_admission/after1229_geometry_roles_and_area_priority.md",
    "archive_1219_/1230/proof.md",
    "archive_1219_/1230/check.py",
]


def poly_integral(n, extra_power=0):
    return sum((F(2 * (-1)**j * math.comb(n, j), 2*j+extra_power+1)
                for j in range(n+1)), F(0))


def atan_interval(q, n):
    s = sum((F((-1)**j, 2*j+1)*q**(2*j+1) for j in range(n)), F(0))
    nxt = F((-1)**n, 2*n+1)*q**(2*n+1)
    return min(s, s+nxt), max(s, s+nxt)


def kernel_and_geometry():
    j = poly_integral(8)
    i0 = poly_integral(16)
    q = j*j/i0
    d1 = 256*poly_integral(14, 2)/i0
    assert j > 0 and i0 > 0 and d1 == F(44, 5)
    l1, u1 = atan_interval(F(1, 5), 15)
    l2, u2 = atan_interval(F(1, 239), 5)
    pi_lo, pi_hi = F(314159, 100000), F(22, 7)
    assert pi_lo < 16*l1-4*u2 < 16*u1-4*l2 < pi_hi
    # Heat kernel: (4 pi t)^(-3/2)/sqrt(pi t) = t^-2/(8 pi^2).
    # P_- contributes 1/2 and radial derivative of exp(-r^2/(4t)) contributes d/2.
    heat_prefactor_without_pi = F(1, 8)
    dirac_prefactor_without_pi = heat_prefactor_without_pi/F(4)
    assert dirac_prefactor_without_pi == F(1, 32)
    assert 8*dirac_prefactor_without_pi == F(1, 4)  # K2/r^2 conversion.
    intervals = []
    rows = []
    for a in (F(1, 10), F(1, 5)):
        am, ap, bm, bp = a-GAMMA, a+GAMMA, B-GAMMA, B+GAMMA
        dm, dp = D-GAMMA-2*bp, D+GAMMA+2*bp
        r2 = dp*dp + 8*ap*ap
        assert 0 < dm and r2 < 4
        lo = am*am*bm*q**3*dm*(2/r2**2-F(1, 2)/r2)/(4*pi_hi**2)
        hi = ap*ap*bp*q**3*dp*2/dm**4/(4*pi_lo**2)
        h2max = 1+d1*(2/am**2+1/bm**2)
        assert 0 < lo < hi < F(1, 2)
        assert h2max < 400**2
        intervals.append((lo, hi))
        rows.append(dict(a_nominal=float(a),area_nominal=float(4*a*a),
                         area_interval=[float(4*am*am),float(4*ap*ap)],
                         actual_mode_geometry_intervals=dict(
                             ax=[float(am),float(ap)],ay=[float(am),float(ap)],
                             b=[float(bm),float(bp)],D=[float(D-GAMMA),float(D+GAMMA)]),
                         longitudinal_gap_lower=float(dm),distance_squared_upper=float(r2),
                         c_lower=float(lo),c_upper=float(hi),
                         field_mode_h_squared_upper=float(h2max),
                         exact_c_lower=str(lo),exact_c_upper=str(hi)))
    assert 2*pi_hi*(400+14) < C
    # Optional Schwinger quadrature is only a diagnostic of the analytic bounds.
    x, w = np.polynomial.laguerre.laggauss(80)
    diagnostics = []
    for r in (.98, 1.1, 1.2):
        alpha = r*r/4
        val = float(np.dot(w, x*np.exp(-alpha/x))/(8*alpha**2))
        lower = 2/r**4-.5/r**2
        upper = 2/r**4
        assert lower < val < upper
        diagnostics.append(dict(r=r,K2_over_r2_diagnostic=val,
                                analytic_lower=lower,analytic_upper=upper))
    group = dict(status="PASS",mass=1,flavors=2,used_flavor=1,
                 envelope_order=8,regularity="C7 and H8, not C-infinity",
                 J_exact=str(j),I0_exact=str(i0),Q_exact=str(q),
                 derivative_norm_ratio_exact=str(d1),
                 Fourier_kernel_G="K1(r)/(2*pi^2*r)",
                 covariance_kernel="d_z*K2(r)/(4*pi^2*r^2)",
                 Schwinger_covariance_prefactor_without_pi=str(dirac_prefactor_without_pi),
                 K2_over_r2_bounds="2/r^4-1/(2*r^2) <= K2(r)/r^2 <= 2/r^4",
                 pi_certified_interval=[str(pi_lo),str(pi_hi)],
                 geometry_tolerance_exact=str(GAMMA),
                 interval_applies_to_mode_normalization_and_support=True,
                 menus=rows,schwinger_diagnostic_only=diagnostics,
                 certificate="analytic continuum kernel plus exact rational geometry bounds")
    return group, intervals


def annihilator(j):
    a = np.zeros((4,4), complex)
    for s in range(4):
        if (s >> j) & 1:
            a[s ^ (1 << j),s] = (-1)**((s & ((1 << j)-1)).bit_count())
    return a


def full_car_records():
    a, b = annihilator(0), annihilator(1)
    nf, ng = a.conj().T@a, b.conj().T@b
    eye = np.eye(4)
    assert np.linalg.norm(a@b.conj().T+b.conj().T@a) == 0
    assert np.linalg.norm(nf@ng-ng@nf) == 0
    v = np.zeros((64,4),complex)
    for x in (0,1):
        for y in (0,1):
            op = (nf if x else eye-nf)@(ng if y else eye-ng)
            for i in range(4):
                for j in range(4):
                    v[(i*4+3*x)*4+3*y,j] = op[i,j]
    assert np.linalg.norm(v.conj().T@v-eye) == 0
    residual = 0.
    cases = []
    for c in (3e-6,.125,.3):
        p = .25-c*c
        rho = np.diag([p,.25+c*c,.25+c*c,p]).astype(complex)
        rho[1,2], rho[2,1] = -1j*c,1j*c
        assert np.min(np.linalg.eigvalsh(rho)) > 0
        cfg = np.trace(rho@a.conj().T@b)
        joint = v@rho@v.conj().T
        diag = np.real(np.diag(joint)).reshape(4,4,4)
        probs = [float(diag[:,3*x,3*y].sum()) for x,y in ((0,0),(0,1),(1,0),(1,1))]
        expected = [p,.25+c*c,.25+c*c,p]
        err = max(abs(cfg-1j*c),max(abs(x-y) for x,y in zip(probs,expected)),
                  abs(np.trace(joint)-1))
        residual = max(residual,float(err))
        assert err < 1e-14
        cases.append(dict(c=c,probabilities=probs))
    return dict(status="PASS",field_restriction_dimension=4,PE_dimension=16,
                full_joint_dimension=64,mode_algebras_commute_evenly=True,
                cases=cases,max_float_residual=residual,
                ideal_probabilities="(.25-c^2,.25+c^2,.25+c^2,.25-c^2)",
                full_continuum_not_replaced_by_finite_matrix=True,
                both_environment_registers_retained=True)


def clock_source_budget():
    force = 16/DELTA**2
    bT = SIGMA+T/(2*M*SIGMA)
    epsC, epsp = force*T*bT, 2*C*DELTA
    e = 2*(epsC+epsp)
    rh, rq = 2*C, A0+2*C
    rha, rqa = 2*T*C/DELTA, A0+2*T*C/DELTA
    errH = (rha+rh)*epsC+2*rh*epsp
    errQ = (rqa+rq)*epsC+2*rq*epsp
    endpoint = force*(3*SIGMA+force*T*T/(2*M))
    n = (SIGMA*force*T+force*T*T/(4*M*SIGMA))**2
    energy = 2*(M/2+1/(8*M*SIGMA**2)+1)
    assert e < F(121,10**22)
    assert errH < F(721,10**19) and errQ < F(323,10**14)
    assert endpoint < F(1,10**60) and n < F(1,10**140)
    assert 2*F(1,10**140)*(2+140*3) < F(1,10**136)
    group = dict(status="PASS",same_center_times=[1,1],T=2,clock_count=2,
                 M="1e260",sigma="1e-120",delta="1e-24",
                 per_gate_H_and_Q_commutator_upper=int(C),
                 field_mode_H_norm_upper=400,force_upper=float(force),
                 clock_to_classical_vector_error=float(epsC),
                 simultaneous_pulse_to_joint_instrument_vector_error=float(epsp),
                 any_event_probability_error=float(e),
                 branch_stress_error_upper=float(errQ),field_energy_error_upper=float(errH),
                 classical_and_ideal_field_H_norm_upper=float(rh),
                 actual_field_H_norm_upper=float(rha),
                 initial_two_clock_energy_display=float(energy),pointer_baseline_energy=4,
                 initial_and_terminal_interaction_sum_upper=float(endpoint),
                 field_energy_gain_upper=float(rh+errH),
                 reaction="|Delta E_clocks + Delta E_field| <= endpoint bound",
                 each_clock_ref_free_occupation_upper=float(n),
                 both_clock_entropy_upper="1e-136",PE_entropy_upper="4 log 2",
                 stress_Q="same 1228 compact spacetime T00 smearing",
                 source_is_same_terminal_state_free_continuation=True,
                 not_actual_driven_future_source_history=True,
                 all_apparatus_boost_or_proper_time=False)
    return group,e,energy+4


def information_response(intervals,e,per_trial_energy):
    def mi_bounds(bounds, event_error):
        lo,hi=bounds
        xi=2*event_error+event_error**2
        assert lo**2 > xi and event_error < F(1,2)
        low=8*(lo**2-xi)**2
        high=(hi**2+xi)**2/(F(1,4)-event_error**2)**2
        return low,high,xi
    actual=[mi_bounds(row,e) for row in intervals]
    gap=actual[1][0]-actual[0][1]
    assert gap > F(4,10**22)
    # Optional empirical table interval, not a performed experiment.
    s=F(1,10**14)
    ntrials=4*10**28
    alpha=F(1,100)
    assert 2*ntrials*s*s == 8
    exp8lower=sum((F(8)**j/math.factorial(j) for j in range(20)),F(0))
    assert exp8lower > 16/alpha
    empirical=[mi_bounds(row,e+2*s) for row in intervals]
    empirical_gap=empirical[1][0]-empirical[0][1]
    assert empirical_gap > F(3,10**22)
    assert actual[0][0]>0 and actual[1][0]>0
    return dict(status="PASS",logarithm="natural",quantity="actual binary pointer mutual information",
                actual_MI_intervals=[[float(x),float(y)] for x,y,_ in actual],
                actual_covariance_error_upper=float(actual[0][2]),
                actual_large_minus_small_MI_lower=float(gap),
                exact_actual_gap_lower=str(gap),certified_gap_greater_than="4e-22 nats",
                ideal_region_MI_lower_bounds=[float(8*x**4) for x,_ in intervals],
                actual_MI_not_automatically_bounded_by_original_region_MI=True,
                region_bound_from_actual="I_region >= 8 max(|d_actual|-xi,0)^2",
                optional_sampling=dict(per_menu_trials=ntrials,independent_preparations_input=True,
                    per_cell_frequency_error=str(s),predeclared_failure_probability=float(alpha),
                    union_failure_bound="16 exp(-8) < .01",
                    empirical_MI_intervals=[[float(x),float(y)] for x,y,_ in empirical],
                    empirical_gap_lower=float(empirical_gap),actual_trials_performed=False,
                    all_PE_registers_retained=8*ntrials,all_clock_pairs_retained=2*ntrials,
                    total_preparation_baseline_energy_display=float(2*ntrials*per_trial_energy),
                    no_free_erasure_or_reused_blank_records=True))


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-results",action="store_true")
    mode.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    kernel,intervals=kernel_and_geometry()
    car=full_car_records()
    budget,e,energy=clock_source_budget()
    info=information_response(intervals,e,energy)
    result=dict(round=1230,status="PASS",check_groups=4,
                groups=dict(continuum_kernel_and_independent_geometry=kernel,
                            full_CAR_and_entity_records=car,
                            actual_finite_information_response=info,
                            same_source_and_complete_resources=budget),
                scope=dict(normal_same_massive_Dirac_vacuum=True,
                    continuum_positive_and_negative_spectra_retained=True,
                    geometry_precedes_entropy=True,finite_region_effective_gate_input=True,
                    same_reference_units_slice_and_geometry_role_contract_input=True,
                    microscopic_pointwise_device=False,accelerated_thermometer=False,
                    full_region_entropy_area_law=False,universal_area_coefficient=False,
                    dynamical_area_backreaction=False,all_apparatus_boost_conservation=False,
                    empirical_experiment=False,stage_complete=False),
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES})
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
    else:
        saved=json.loads(TARGET.read_text(encoding="utf8"))
        if result != saved:
            raise SystemExit("Saved receipt differs; no files changed.")
    print(json.dumps(dict(round=1230,status="PASS",check_groups=4,
        actual_MI_gap_lower=info["actual_large_minus_small_MI_lower"],
        field_energy_error_upper=budget["field_energy_error_upper"],
        branch_stress_error_upper=budget["branch_stress_error_upper"],
        read_only=not args.write_results),ensure_ascii=False))


if __name__ == "__main__":
    main()
