"""1232 independent finite-measure and full-probability audit.

No author code/results imported. All accepted parameter margins use Fraction;
finite thermal matrices and atomic spectra are regression fixtures, not Dirac
spectral measurements or an autonomous accelerated detector.
Default/--verify is read-only; --write-results exclusively creates a receipt.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = HERE / "independent_results.json"
SOURCES = (
    "archive_1120_1158/1126/proof.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1229/proof.md",
    "archive_1219_/_admission/stage_scope_audit_after1230.md",
    "archive_1219_/1232/independent_check.py",
)
T = F(128)
R = F(1, 4)
CORE = F(1, 64)
LAMBDA = F(1, 2**16)
EPS_OP = F(1, 10**10)
RHO = F(1, 100)
C0 = F(23, 24)**4


def view(x):
    return {"exact": str(x), "display": float(x)}


def cos_interval(x):
    lower = sum(((-1)**j*x**(2*j)/math.factorial(2*j) for j in range(8)), F(0))
    upper = lower+x**16/math.factorial(16)
    return lower, upper


def sin_interval(x):
    assert 0 <= x <= 1
    lower = sum(((-1)**j*x**(2*j+1)/math.factorial(2*j+1)
                 for j in range(8)), F(0))
    upper = lower+x**17/math.factorial(17)
    return lower, upper


def sign_and_KMS():
    # Environment H=diag(0,1), rho_beta=diag(2/3,1/3), beta=ln 2.
    x = np.array([[0., 1.], [1., 0.]])
    rho = np.diag([2/3, 1/3])
    diagnostics = []
    for t in (.125, .75, 1.25):
        u = np.diag([1., np.exp(1j*t)])
        w = np.trace(rho@(u@x@u.conj().T)@x)
        spectral = (2*np.exp(-1j*t)+np.exp(1j*t))/3
        assert abs(w-spectral) < 1e-14
        diagnostics.append(dict(t=t, residual=float(abs(w-spectral))))
    plus_weight, minus_weight = F(1, 3), F(2, 3)
    assert minus_weight == 2*plus_weight
    return dict(status="PASS", correlation="omega(B_t B_0)",
                Fourier_sign="W(t)=integral exp(+i nu t) mu(dnu)",
                detailed_balance="mu(-dnu)=exp(beta nu) mu(dnu)",
                excitation_filter="abs(chi_hat(E-nu))^2",
                finite_thermal_fixture_not_Dirac=True,
                positive_frequency_weight=str(plus_weight),
                negative_frequency_weight=str(minus_weight),
                diagnostics=diagnostics)


def tent_fourier_sq(frequency, duration):
    x = frequency*duration/2
    if x == 0:
        return duration*duration
    return duration*duration*(math.sin(x)/x)**4


def band_guard_and_atomic_regression():
    assert 1-F(1, 2)**2/6 == F(23, 24)
    assert 1/T <= R < 1
    area = C0*T*T*CORE
    tail = 16/(R**4*T*T)
    u = tail/area
    assert area > 0 and tail == F(1, 4)
    # Different exact atomic KMS fixture from a Dirac continuum calculation.
    atoms = {1: F(1, 32), -1: F(1, 16),
             2: F(1, 48), -2: F(1, 12), 0: F(5, 96)}
    assert sum(atoms.values()) == F(1, 4)
    for nu in (1, 2):
        assert atoms[-nu] == 2**nu*atoms[nu]
        assert atoms[nu] >= CORE
    data = []
    for e in (1, 2):
        fp = sum(float(weight)*tent_fourier_sq(e-nu, float(T))
                 for nu, weight in atoms.items())
        fm = sum(float(weight)*tent_fourier_sq(-e-nu, float(T))
                 for nu, weight in atoms.items())
        reweighted = sum(float(weight)*2**nu*tent_fourier_sq(e-nu, float(T))
                         for nu, weight in atoms.items())
        assert abs(fm-reweighted) < 1e-10
        assert fp >= float(area) and fm >= float(area)
        beta = math.log(2)
        ratio_low = math.exp(beta*(e-float(R)))/(1+float(u))
        ratio_high = math.exp(beta*(e+float(R)))+float(u)
        assert ratio_low <= fm/fp <= ratio_high
        bias = abs(math.log(fm/fp)-beta*e)
        assert bias <= beta*float(R)+math.log1p(float(u))
        data.append(dict(gap=e, response_plus=fp, response_minus=fm,
                         log_bias_diagnostic= bias))
    return dict(status="PASS", exact_core_filter_factor=str(C0),
                leading_response_lower=view(area), each_tail_upper=view(tail),
                reflected_measure_total_mass_bound_used=True,
                tail_over_core=view(u),
                atomic_fixture={str(k): str(v) for k, v in atoms.items()},
                atomic_fixture_not_Dirac_data=True,
                float_regressions_not_continuum_certificates=data), area, u


def full_probability_contract(area, u):
    x = LAMBDA*T
    r3_upper = x**3/(6*(1-x*x/20))
    dyson = 2*x*r3_upper+r3_upper**2
    qmin = LAMBDA**2*area
    all_error = dyson+EPS_OP
    assert 0 < x < 1
    assert all_error < RHO*qmin
    assert LAMBDA**2 <= min(1/T**2, RHO*area/T**4)
    assert EPS_OP <= RHO*qmin/2
    # ln 2 = 2 atanh(1/3), geometric remainder proves ln 2 < .7.
    z = F(1, 3)
    ln2_lower = 2*sum((z**(2*j+1)/(2*j+1) for j in range(8)), F(0))
    ln2_upper = ln2_lower+2*z**17/(17*(1-z*z))
    assert ln2_upper < F(7, 10)
    log_total = ln2_upper*R+u+2*RHO/(1-RHO)
    assert log_total < F(1, 5)
    assert qmin-all_error > 0
    # Certificate is positive in finite parameters, not a populated apparatus class.
    return dict(status="PASS", x=view(x), odd_Dyson_remainder_upper=view(r3_upper),
                complete_probability_error_upper=view(all_error),
                leading_probability_lower=view(qmin),
                actual_probability_lower_if_contract_holds=view(qmin-all_error),
                operation_error_input=view(EPS_OP),
                relative_error_allowance=str(RHO),
                joint_log_ratio_error_upper=view(log_total),
                ln2_interval=[view(ln2_lower), view(ln2_upper)],
                fixed_gap_menu=[1, 2], fixed_core_mass=str(CORE),
                specified_Dirac_core_mass_certified=False,
                actual_motion_clock_source_contract_constructed=False)


def all_order_odd_flip():
    # Commuting pure-flip complete unitary exp(-ix X_D tensor X_F).
    xmat = np.array([[0., 1.], [1., 0.]])
    gen = np.kron(xmat, xmat)
    ident = np.eye(4)
    target = np.diag([0., 0., 1., 1.])
    init = np.array([1., 0., 0., 0.])
    cases = []
    for x in (F(1, 8), F(1, 3), F(1)):
        sl, su = sin_interval(x)
        cl, cu = cos_interval(x)
        r3 = x**3/(6*(1-x*x/20))
        dyson = 2*x*r3+r3*r3
        assert dyson <= x**4/2
        assert 0 < sl <= su <= x
        actual_probability_difference_upper = x*x-sl*sl
        assert actual_probability_difference_upper <= dyson
        u = math.cos(float(x))*ident-1j*math.sin(float(x))*gen
        output = u@init
        probability = float(np.vdot(output, target@output).real)
        assert np.linalg.norm(u.conj().T@u-ident) < 1e-14
        assert float(sl*sl)-1e-15 <= probability <= float(su*su)+1e-15
        cases.append(dict(x=str(x), full_probability=probability,
                          rational_complete_error_guard=view(dyson)))
    return dict(status="PASS", all_order_time_ordered_bound="sum odd n>=3 x^n/n!",
                generic_x_le_1_bound="d_D <= x^4/2",
                pure_flip_structure_required=True,
                full_environment_kept=True,
                finite_complete_unitary_fixture_only=True, cases=cases)


def deletion_checks():
    # A zero-frequency atom is KMS for every beta, but finite even filtering
    # gives F+=F- at every nonzero gap (when the filter is nonzero).
    sl, _ = sin_interval(F(1, 2))
    zero_atom_response_lower = F(1, 4)*(2*sl)**4
    assert zero_atom_response_lower > 0
    zero_atom_ratio = F(1)
    assert zero_atom_ratio != 2
    # A fixed absolute device error can erase an arbitrarily small hot/cold ratio.
    qp, qm, eps = F(1, 10**14), F(2, 10**14), F(1, 10**10)
    pa = ma = eps
    assert abs(pa-qp) <= eps and abs(ma-qm) <= eps and ma/pa == 1
    # Diagonal first, flip second: U=exp(-il X)exp(-il Z).
    # Flip amplitude has even part -sin(l)^2, so the odd-only bound is inapplicable.
    xmat = np.array([[0., 1.], [1., 0.]])
    zmat = np.diag([1., -1.])
    init = np.array([1., 0.])
    l = .125
    def amplitude(v):
        ux = math.cos(v)*np.eye(2)-1j*math.sin(v)*xmat
        uz = math.cos(v)*np.eye(2)-1j*math.sin(v)*zmat
        return (ux@uz@init)[1]
    even_sum = amplitude(l)+amplitude(-l)
    assert abs(even_sum+2*math.sin(l)**2) < 1e-14 and abs(even_sum) > .01
    # Rescaling the geometric rate while keeping detector E fixed rescales beta.
    scale = F(2)
    assert 2**1 == 2**(F(1, 2)*scale)
    return dict(status="PASS", zero_frequency_KMS_without_core_mass=True,
                zero_atom_response_lower=view(zero_atom_response_lower),
                zero_atom_ratio=1,
                fixed_absolute_error_can_destroy_ratio=True,
                diagonal_interaction_creates_second_order_flip=True,
                even_amplitude_sum_diagnostic=[float(even_sum.real), float(even_sum.imag)],
                wrong_clock_rate_effect="beta -> beta/c",
                no_deletion_is_a_counterexample_to_the_entire_cognitive_program=True)


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write-results", action="store_true")
    modes.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    g1 = sign_and_KMS()
    g2, area, u = band_guard_and_atomic_regression()
    g3 = full_probability_contract(area, u)
    g4 = all_order_odd_flip()
    g5 = deletion_checks()
    result = dict(round=1232, status="PASS", check_groups=5,
                  groups=dict(KMS_sign_and_actual_excitation_convention=g1,
                              finite_switch_band_bound_and_regression=g2,
                              simultaneous_finite_conditional_budget=g3,
                              complete_odd_flip_Dyson_bound=g4,
                              necessary_condition_deletions=g5),
                  analytic_review_only=dict(
                      massive_pair_boost_absolute_continuity=True,
                      positive_Lebesgue_point_existence=True,
                      no_assumed_pointwise_continuous_density=True,
                      specified_Dirac_numerical_spectrum=False,
                      finite_comparison_field_energy_from_bounded_commutator=True),
                  scope=dict(conditional_complete_probability_Unruh=True,
                             scalar_spectrum_substituted=False,
                             actual_autonomous_accelerated_apparatus=False,
                             actual_source_internalization=False,
                             TG2_or_stage_completion=False),
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                 for p in SOURCES})
    if args.write_results:
        with TARGET.open("x", encoding="utf8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
    elif json.loads(TARGET.read_text(encoding="utf8")) != result:
        raise SystemExit("Saved independent receipt differs; no files changed.")
    print(json.dumps(dict(round=1232, status="PASS", check_groups=5,
                          probability_error_upper=g3["complete_probability_error_upper"]["display"],
                          joint_log_error_upper=g3["joint_log_ratio_error_upper"]["display"],
                          read_only=not args.write_results), ensure_ascii=False))


if __name__ == "__main__":
    main()
