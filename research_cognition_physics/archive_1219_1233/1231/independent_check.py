"""1231 independent audit of a finite effective scattering boundary.

No author code or results imported. Fraction calculations certify inequalities;
finite matrices below only verify the stated conditional logical identities.
Neither an actual phase reader nor a local finite-time collision is constructed.
Default/--verify is read-only; --write-results creates only a missing receipt.
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
    "archive_1106_1119/1118/proof.md",
    "archive_1106_1119/1119/proof.md",
    "archive_1219_/1219/candidate_proof.md",
    "archive_1219_/_admission/after1230_finite_motion_thermal_admission.md",
    "archive_1219_/1231/independent_check.py",
)
SOURCE_MASS = 10**14
WIDTH = F(1, 100)
SOURCE_WIDTH = F(1, 10**7)
BAD_WIDTH = F(1, 10**16)
# Independent rational guards tighter than the author's A<3, sqrt35>5.
A_UPPER = F(297, 100)
SQRT35_LOWER = F(591, 100)
RMAX = F(8, SOURCE_MASS)
RDER = F(16, SOURCE_MASS)


def show(value):
    return {"exact": str(value), "display": float(value)}


def even_integral(n, extra=0):
    return sum((F(2*(-1)**j*math.comb(n, j), 2*j+extra+1)
                for j in range(n+1)), F(0))


def exp_interval(x):
    assert 0 <= x <= 1
    n = 20
    lower = sum((x**j/math.factorial(j) for j in range(n+1)), F(0))
    next_term = x**(n+1)/math.factorial(n+1)
    upper = lower+next_term/(1-x/(n+2))
    return lower, upper


def sinh_cosh_bounds(x):
    lo, hi = exp_interval(x)
    return (lo-1/lo)/2, (hi-1/hi)/2, (hi+1/lo)/2


def normalized_packets():
    i0 = even_integral(16)
    dnorm = 256*even_integral(14, 2)/i0
    second = even_integral(16, 2)/i0
    assert dnorm == F(44, 5)
    assert second == F(1, 35)
    assert A_UPPER**2 > dnorm and SQRT35_LOWER**2 < 35
    support = F(1, 2)+WIDTH+SOURCE_WIDTH
    assert support < F(52, 100)
    _, sinh_upper, cosh_upper = sinh_cosh_bounds(F(52, 100))
    assert sinh_upper < 1 and cosh_upper < 2
    return dict(status="PASS", I0=show(i0), derivative_norm_squared=str(dnorm),
                variance_scaled=str(second), A_upper=str(A_UPPER),
                sqrt35_lower=str(SQRT35_LOWER), relative_support=show(support),
                sinh_support_upper=show(sinh_upper),
                cosh_support_upper=show(cosh_upper),
                normal_compact_packets=True, same_source_for_unknown_input=True)


def scatter_lightcone(mu, source_mass, q, v):
    plus = mu*q+source_mass*v
    minus = mu/q+source_mass/v
    boost_square = plus/minus
    return boost_square/q, boost_square/v, plus, minus


def exact_scattering():
    fixtures = []
    for mu in (F(3), F(4)):
        for ell in (F(17), F(SOURCE_MASS)):
            for q, v in ((F(1), F(5, 3)), (F(7, 5), F(6, 5)),
                         (F(2, 3), F(11, 7))):
                qp, vp, plus, minus = scatter_lightcone(mu, ell, q, v)
                qq, vv, pp, mm = scatter_lightcone(mu, ell, qp, vp)
                assert mu*qp+ell*vp == plus and mu/qp+ell/vp == minus
                assert qq == q and vv == v and pp == plus and mm == minus
                assert (qp/vp) == v/q
                # e^R agrees with the source's exact rapidity change.
                x = q/v
                er = (1+mu*x/ell)/(1+mu/(ell*x))
                assert vp/v == er
                rprime = mu*x/(ell+mu*x)+mu/(ell*x+mu)
                assert rprime >= 0
                # In (d,xi): [[-1,0],[R'(d),1]], exactly determinant -1.
                determinant = F(-1)*1-F(0)*rprime
                assert determinant == -1
                energy = mu*(qp+1/qp)/2
                momentum = mu*(qp-1/qp)/2
                assert energy**2-momentum**2 == mu**2
                fixtures.append(dict(mu=str(mu), source_mass=str(ell),
                                     input_exp_rapidities=[str(q), str(v)],
                                     exact_fourP_and_involution=True))
    return dict(status="PASS", exact_rational_fixtures=fixtures,
                general_measure="d eta d xi = d d d xi",
                Jacobian_determinant="-1", coherent_pullback_unitary=True,
                transverse_momenta_zero_in_declared_1plus1_model=True,
                scattering_phase_zero_is_input=True,
                local_continuous_collision_constructed=False)


def coherent_transport():
    target_width_term = 2*SOURCE_WIDTH*A_UPPER/(SQRT35_LOWER*WIDTH)
    recoil_target_term = RMAX*A_UPPER/WIDTH
    recoil_source_term = RMAX*A_UPPER/SOURCE_WIDTH
    epsilon = target_width_term+recoil_target_term+recoil_source_term
    assert epsilon < F(15, 10**6)
    difference = F(2, SOURCE_MASS)*A_UPPER*(1/WIDTH+1/SOURCE_WIDTH)
    loss = difference**2/2
    assert loss < F(181, 10**15)
    source_hi = F(1, 2)+SOURCE_WIDTH
    _, _, source_cosh = sinh_cosh_bounds(source_hi)
    _, _, target_cosh = sinh_cosh_bounds(WIDTH)
    energy_hi = SOURCE_MASS*source_cosh+4*target_cosh
    assert energy_hi < 2*SOURCE_MASS+5
    return dict(status="PASS",
                complete_vector_error_terms=[show(x) for x in
                    (target_width_term, recoil_target_term, recoil_source_term)],
                complete_isometry_error=show(epsilon),
                half_diamond_bound=show(epsilon),
                all_internal_inputs_and_passive_references=True,
                complete_motional_branch_vector_difference=show(difference),
                joint_motion_overlap_loss=show(loss),
                overlap_not_claimed_between_pure_source_marginals=True,
                initial_active_energy_upper=show(energy_hi),
                finite_energy_moments_from_compact_support=True), epsilon, loss


def narrow_source_failure():
    _, sinh_hi, _ = sinh_cosh_bounds(F(511, 1000))
    exp_hi = exp_interval(F(511, 1000))[1]
    sinh_lo = sinh_cosh_bounds(F(489, 1000))[0]
    assert sinh_hi < F(535, 1000)
    assert exp_hi < 2 and sinh_lo > F(489, 1000)
    # No numerical subtraction of nearly equal rapidities at order 1/2.
    gap = (8*F(489, 1000)/(1+F(8, SOURCE_MASS))
           -6*F(535, 1000))/SOURCE_MASS-2*BAD_WIDTH
    assert gap > F(6, 10**15)
    assert F(1, 2)-WIDTH-BAD_WIDTH > F(489, 1000)
    assert F(1, 2)+WIDTH+BAD_WIDTH < F(511, 1000)
    return dict(status="PASS", bad_source_width=str(BAD_WIDTH),
                disjoint_source_support_gap=show(gap),
                exact_joint_motion_overlap=0,
                reduced_internal_coherence_erased=True,
                complete_global_output_kept=True,
                same_exact_fourP_law=True,
                failure_scope="changed source preparation in this boundary model")


def free_phase_and_position(epsilon):
    width_out = WIDTH+2*SOURCE_WIDTH+RMAX
    assert width_out < F(10001, 10**6)
    phase_error = 2*width_out**2/(1-width_out**2)
    complete_error = epsilon+phase_error
    assert complete_error < F(215, 10**6)
    xa = ((1+RDER)*A_UPPER/WIDTH+RDER*A_UPPER/SOURCE_WIDTH+F(1, 2))/3
    xb = ((2+RDER)*A_UPPER/WIDTH
          +(1+RDER)*A_UPPER/SOURCE_WIDTH+F(1, 2))/SOURCE_MASS
    assert xa < 101 and xb < F(31, 10**8)
    # Alternating cosine series: lower terminates with a negative term.
    cos_lower = sum((F((-1)**k, math.factorial(2*k)) for k in range(8)), F(0))
    conditional_gap = cos_lower-2*complete_error
    assert conditional_gap > F(539, 1000)
    return dict(status="PASS", output_rapidity_half_width=show(width_out),
                proper_parameter=1, center_rapidity=1,
                free_generator="u_star dot P_A = mu cosh(eta-1)",
                free_phase_error=show(phase_error),
                complete_error_after_free_phase=show(complete_error),
                target_position_norm_upper=show(xa),
                source_position_norm_upper=show(xb),
                conditional_X_and_rapidity_event_gap=show(conditional_gap),
                X_effect_requires_additional_permission=True,
                actual_288dim_Mr_clock_accelerated=False,
                contact_period_proper_time_constructed=False,
                asymptotic_position_not_identified_with_actual_CO_meter=True)


def coherence_and_phase_freedom(loss):
    # The matrix only checks logical consequences of the analytically derived
    # overlap, not an assigned source Gram matrix or a simulated source.
    plus = np.array([1., 1.])/math.sqrt(2)
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    pop = np.abs(plus)**2
    assert np.array_equal(np.abs(z@plus)**2, pop)
    assert abs(np.vdot(plus, x@plus)-1) < 1e-14
    assert abs(np.vdot(z@plus, x@(z@plus))+1) < 1e-14
    # Global environment orthogonality proved above gives exact dephasing.
    # Population is identical, but the original |+> effect changes by 1/2.
    rho = np.array([[F(1, 2), F(1, 2)], [F(1, 2), F(1, 2)]], dtype=object)
    bad = np.array([[F(1, 2), F(0)], [F(0), F(1, 2)]], dtype=object)
    effect = rho
    p_good = np.trace(rho@effect)
    p_bad = np.trace(bad@effect)
    assert p_good-p_bad == F(1, 2)
    return dict(status="PASS", phase_zero_and_pi_keep_populations=True,
                exact_bad_source_conditional_effect_difference="1/2",
                good_reduced_dephasing_half_diamond_upper=show(loss/2),
                arbitrary_phase_not_determined_by_conservation=True,
                actual_phase_reader_not_constructed=True,
                source_not_discarded_in_full_model=True)


def main():
    parser = argparse.ArgumentParser()
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write-results", action="store_true")
    options.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    g1 = normalized_packets()
    g2 = exact_scattering()
    g3, epsilon, loss = coherent_transport()
    g4 = narrow_source_failure()
    g5 = free_phase_and_position(epsilon)
    g6 = coherence_and_phase_freedom(loss)
    result = dict(round=1231, status="PASS", check_groups=6,
                  groups=dict(normal_wavepacket_and_hyperbolic_bounds=g1,
                              exact_lightcone_shear_and_measure=g2,
                              same_source_complete_coherent_transport=g3,
                              source_support_decoherence_countercase=g4,
                              free_phase_position_and_conditional_prediction=g5,
                              phase_freedom_and_logical_scope=g6),
                  scope=dict(effective_asymptotic_boundary=True,
                             normal_finite_source=True,
                             coherent_mass_preparation_additional_input=True,
                             same_1219_internal_mass_spectrum=True,
                             actual_1219_X_reader_permission=False,
                             finite_contact_worldtube=False,
                             actual_accelerated_thermometer=False,
                             same_Dirac_field_coupling=False,
                             Unruh_or_TG2_completion=False),
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                 for p in SOURCES})
    if args.write_results:
        with TARGET.open("x", encoding="utf8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
    elif json.loads(TARGET.read_text(encoding="utf8")) != result:
        raise SystemExit("Saved independent receipt differs; no files changed.")
    print(json.dumps(dict(round=1231, status="PASS", check_groups=6,
                          complete_vector_error_upper=float(epsilon),
                          joint_motion_overlap_loss_upper=float(loss),
                          narrow_source_gap=float(g4["disjoint_source_support_gap"]["display"]),
                          read_only=not args.write_results), ensure_ascii=False))


if __name__ == "__main__":
    main()
