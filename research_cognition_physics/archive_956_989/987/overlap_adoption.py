"""Recompute overlap diagnostics from frozen 978/985/986 evidence.

This is an adoption audit, not a new simulation of a joined physical model.
No cutoff, runtime installation, or writes occur without --write.
"""
from pathlib import Path
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / "overlap_adoption_results.json"


def read(path):
    return json.loads(path.read_text("utf-8-sig"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(value):
    return {"exact": str(value), "decimal": float(value)}


def run():
    r978 = read(STAGE / "978/moving_complete_source_results.json")
    r985 = read(STAGE / "985/two_active_sources_results.json")
    r986 = read(STAGE / "986/metric_dipole_bridge_results.json")
    for value in (r978, r985, r986):
        assert value["all_scientific_checks_passed"]
    sigma = r978["parameters"]["position_sigma"]
    k = r986["geometry_input"]["TT_frequency"]
    z = k * sigma
    e985 = Fraction(r985["uniform_budget"]["rounded_state_upper"]["exact"])
    e986 = Fraction(r986["infinite_domain_bound"]["reported_state_upper"]["exact"])
    overlap = e985 + e986
    assert overlap == Fraction(860567, 500000000)
    assert 2 * e986 == Fraction(67, 250000000)

    # Preserve exp(-800) as a nonzero decimal; float64 would underflow.
    with localcontext() as ctx:
        ctx.prec = 70
        zd = Decimal(str(z))
        mean = (-zd * zd / 2).exp()
        oscillatory_second = (-2 * zd * zd).exp()
        mean_string, second_correction_string = str(mean), str(oscillatory_second)
    second = .5 * (1 + math.exp(-2 * z * z))
    difference = second - 2 * float(mean_string) + 1
    assert z == 20 and abs(second - .5) < 1e-15
    assert abs(difference - 1.5) < 1e-15

    # Independent Gaussian quadrature in units of the COM standard deviation.
    # At z=20, the numerical mean only tests an absolute error, not exp(-200).
    x = np.linspace(-12., 12., 131073)
    density = np.exp(-x*x/2) / math.sqrt(2*math.pi)
    dx = x[1]-x[0]
    quadrature = []
    for frequency in (.25, 1., 3., z):
        f = np.cos(frequency*x)
        integ = lambda a: dx * (np.sum(a) - .5*(a[0]+a[-1]))
        numeric = [integ(density*f), integ(density*f*f), integ(density*(f-1)**2)]
        a = math.exp(-frequency*frequency/2)
        b = .5*(1+math.exp(-2*frequency*frequency))
        expected = [a, b, b-2*a+1]
        error = max(abs(i-j) for i,j in zip(numeric, expected))
        assert error < 2e-14
        quadrature.append(dict(k_sigma=frequency, values=[float(i) for i in numeric],
                               max_absolute_error=float(error)))

    g = r986["inherited"]["g"]
    d = (1-1/math.sqrt(1.16))/2
    c_second = np.array([1+3*g*g*d, 1+g*g*d])
    lambda2 = 4*math.pi*r985["parameters"]["G"]
    source_checks = []
    for phi in r978["parameters"]["potentials"]:
        lapse = 1+phi
        constant = lambda2*lapse*lapse*c_second
        profiled = second*constant
        source_checks.append(dict(lapse=lapse,
            linear_vertex_force_second_moment_constant=constant.tolist(),
            linear_vertex_force_second_moment_profiled=profiled.tolist(),
            absolute_difference=(constant-profiled).tolist(),
            relative_difference=(1-profiled/constant).tolist()))

    paths = [Path(__file__), HERE/"drafts/overlap_decision.md"] + [STAGE/p for p in (
        "978/moving_complete_source_results.json", "985/two_active_sources_results.json",
        "986/metric_dipole_bridge_results.json", "research_note_953.md",
        "research_note_985.md", "research_note_986.md", "987/drafts/STATUS.md",
        "981/drafts/common_parent_contract_v1.md")]
    return dict(date="2026-10-07", kind="overlap_adoption_audit",
        formal_reports=986, cumulative_numbered_test_groups=3771,
        new_numbered_scientific_groups=0, all_audit_calculations_passed=True,
        common_record=dict(state_trace_distance_upper=exact(overlap),
            bounded_effect_probability_difference_upper=exact(overlap),
            two_input_contrast_difference_upper=exact(2*overlap),
            inherited_isometry_bounds_not_added_as_full_physical_errors=True,
            common_preparation_and_read_map_required=True),
        conditional_profile=dict(input="u(X_z)=cos(k X_z); centered COM Gaussian",
            profile_is_additional_diagnostic_input=True,
            full_profile_was_not_specified_in_986=True,
            sigma_COM=sigma, k=k, k_sigma=z,
            mean_exact="exp(-200)", mean_decimal=mean_string,
            second_exact="(1+exp(-800))/2", exp_minus_800=second_correction_string,
            second_float=second, squared_difference_from_one_float=difference,
            intrinsic_body_size_not_identified_with_COM_spread=True,
            source_mean_zero_hides_second_moment_difference=True,
            linear_source_second_moment_by_lapse=source_checks,
            independent_quadrature=quadrature,
            tiny_mean_certified_by_formula_not_floating_quadrature=True),
        spectator_comparison=dict(position_kinetic_energy_included=False,
            uniform_in_bounded_real_profile_abs_le_one=True,
            state_to_original_upper=exact(e986),
            constant_vs_profiled_state_upper=exact(2*e986),
            proves_join_with_985_motion=False),
        established_EFT=dict(reference="https://arxiv.org/html/hep-th/0409156v2",
            adopted="Eq.(51) binding energy belongs in long-wave lapse coupling",
            full_binary_k_r=[k*r for r in r978["parameters"]["radii"]],
            binary_long_wave_formula_directly_applicable_to_986_mode=False,
            internal_noncommuting_material_matching_automatically_done=False),
        decision=dict(old_record_overlap_certified=True,
            constant_spatial_source_transfer_certified=False,
            actual_finite_record_failure_demonstrated=False,
            all_goals_counterexample=False, full_U1_completed=False,
            full_goal_completed=False, UV_completion_required=False,
            optimize_this_candidate_next=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual: compare(actual[key], saved[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a,b in zip(actual,saved): compare(a,b)
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, rel_tol=1e-12, abs_tol=1e-14)
    else:
        assert actual == saved, (actual, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open("x", encoding="utf-8") as dest:
            json.dump(result,dest,ensure_ascii=False,indent=2); dest.write("\n")
    else:
        compare(result,read(TARGET))
    print(json.dumps({key:result[key] for key in (
        "all_audit_calculations_passed", "common_record", "spectator_comparison", "decision")},
        ensure_ascii=False,indent=2))
