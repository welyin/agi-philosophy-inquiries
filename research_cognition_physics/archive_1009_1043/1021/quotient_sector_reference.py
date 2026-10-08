"""1021: quotient intertwining is not a choice of the full thermal sector.

Finite Z12 -> Z2 calibration on a three-edge oriented ring. No SM spectrum,
continuum reconstruction or physical measuring apparatus is claimed. The
Gauss basis is independently enumerated before the 12-dimensional reduction.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "quotient_sector_reference_results.json"


def norm(a):
    return float(np.linalg.norm(a))


def hermitian_function(a, function):
    values, vectors = np.linalg.eigh(a)
    return (vectors * function(values)) @ vectors.conj().T


def trace_real(a):
    value = np.trace(a)
    assert abs(value.imag) < 2e-11
    return float(value.real)


def gauss_calibration():
    ring = np.array([[-1, 0, 1], [1, -1, 0], [0, 1, -1]], dtype=int)
    tree = np.array([[-1, 0], [1, -1], [0, 1]], dtype=int)

    def invariant_flows(incidence, order):
        return [list(q) for q in itertools.product(range(order), repeat=incidence.shape[1])
                if np.all((incidence @ np.array(q)) % order == 0)]

    cover = invariant_flows(ring, 12)
    quotient = invariant_flows(ring, 2)
    kernel_ring = invariant_flows(ring, 6)
    kernel_tree = invariant_flows(tree, 6)
    assert cover == [[q, q, q] for q in range(12)]
    assert quotient == [[0, 0, 0], [1, 1, 1]]
    assert len(kernel_ring) == 6 and kernel_tree == [[0, 0]]
    assert [row for row in cover if all(q % 6 == 0 for q in row)] == [[0, 0, 0], [6, 6, 6]]
    return dict(oriented_ring_incidence=ring.tolist(), tree_incidence=tree.tolist(),
                cover_enumerated_electric_triples=12**3, cover_Gauss_basis=cover,
                quotient_Gauss_basis=quotient, kernel_ring_cycle_flows=kernel_ring,
                kernel_tree_cycle_flows=kernel_tree, cover_dimension=12,
                quotient_dimension=2, ring_cycle_sector_count=6,
                tree_cycle_sector_count=1, full_configuration_matrix_constructed=False)


def model(lam):
    q = np.arange(12)
    lap = np.diag(2 - 2*np.cos(2*np.pi*q/12))
    # Use exact special values at the quotient-supported electric modes.
    lap[0, 0] = 0.; lap[6, 6] = 4.
    shift = np.roll(np.eye(12), 6, axis=0)
    kappa = np.array([.7, 1.1, .9])*np.exp(-2*lam)
    total = float(kappa.sum()); v = .4
    H = total*lap + v*(np.eye(12)-shift)
    source = -2*total*lap; second = 4*total*lap
    W = np.eye(12)[:, [0, 6]]
    X = np.array([[0., 1.], [1., 0.]])
    lap_q = np.diag([0., 4.])
    Hq = total*lap_q + v*(np.eye(2)-X)
    source_q = -2*total*lap_q; second_q = 4*total*lap_q
    deck = np.diag(np.exp(2j*np.pi*q/6))
    return dict(lambda_value=float(lam), kappa=kappa, K=total, v=v,
                L=lap, X6=shift, H=H, S=source, H2=second, W=W,
                X=X, Lq=lap_q, Hq=Hq, Sq=source_q, H2q=second_q, T=deck)


def unitary(H, duration):
    return hermitian_function(H, lambda values: np.exp(-1j*duration*values))


def process_calibration(m):
    H, S, H2, W, Hq = (m[k] for k in ("H", "S", "H2", "W", "Hq"))
    identity = np.eye(12); iq = np.eye(2)
    errors = dict(isometry=norm(W.T@W-iq), Hamiltonian=norm(H@W-W@Hq),
                  source=norm(S@W-W@m["Sq"]), second_source=norm(H2@W-W@m["H2q"]),
                  deck_commutation=norm(H@m["T"]-m["T"]@H),
                  quotient_support_deck_identity=norm(m["T"]@W-W))
    kraus = [hermitian_function(identity/2 + sign*m["X6"]/4, np.sqrt) for sign in (1, -1)]
    kraus_q = [hermitian_function(iq/2 + sign*m["X"]/4, np.sqrt) for sign in (1, -1)]
    errors["full_instrument_completeness"] = norm(sum(k.conj().T@k for k in kraus)-identity)
    errors["quotient_instrument_completeness"] = norm(sum(k.T@k for k in kraus_q)-iq)
    errors["individual_Kraus_intertwining"] = max(norm(k@W-W@kq) for k, kq in zip(kraus, kraus_q))
    times = [.23, .71]
    U = [unitary(H, t) for t in times]; Uq = [unitary(Hq, t) for t in times]
    errors["waiting_intertwining"] = max(norm(u@W-W@uq) for u, uq in zip(U, Uq))
    lift = np.kron(W, iq)
    histories = []; total_probability = 0.; sum_effect = np.zeros((2, 2), dtype=complex)
    for first, second in itertools.product(range(2), repeat=2):
        K = kraus[second]@U[1]@kraus[first]@U[0]@W
        Kq = kraus_q[second]@Uq[1]@kraus_q[first]@Uq[0]
        # Normalized Choi vectors retain an untouched two-dimensional reference.
        c = K.reshape(-1)/math.sqrt(2); cq = Kq.reshape(-1)/math.sqrt(2)
        choi = np.outer(c, c.conj()); choi_q = np.outer(cq, cq.conj())
        p = float(np.vdot(c, c).real); total_probability += p
        sum_effect += Kq.conj().T@Kq
        histories.append(dict(first_outcome=first, second_outcome=second,
            full_Kraus_intertwining_error=norm(K-W@Kq),
            reference_Choi_error=norm(choi-lift@choi_q@lift.conj().T),
            branch_probability_on_maximally_entangled_input=p))
    errors["history_completeness"] = norm(sum_effect-iq)
    assert max(errors.values()) < 3e-13
    assert max(h["reference_Choi_error"] for h in histories) < 3e-13
    assert max(h["full_Kraus_intertwining_error"] for h in histories) < 3e-13
    assert abs(total_probability-1) < 3e-13
    assert norm(H@S-S@H) > 1
    return dict(lambda_value=m["lambda_value"], edge_kappa=m["kappa"].tolist(),
                total_K=m["K"], v=m["v"], errors=errors, histories=histories,
                history_total_probability=total_probability,
                Hamiltonian_source_commutator_norm=norm(H@S-S@H),
                unknown_inputs_and_untouched_reference_retained=True,
                claim_for_arbitrary_cover_sector_input=False,
                analytic_all_time_all_finite_descending_histories=True)


def thermal(H, beta, S=None, H2=None):
    values, vectors = np.linalg.eigh(H)
    raw = np.exp(-beta*values); Z = float(raw.sum()); p = raw/Z
    rho = (vectors*p)@vectors.conj().T
    logp = -beta*values-math.log(Z)
    log_rho = (vectors*logp)@vectors.conj().T
    out = dict(rho=rho, log_rho=log_rho, Z=Z, F=-math.log(Z)/beta,
               energies=values, populations=p)
    if S is not None:
        mean = trace_real(rho@S)
        Se = vectors.conj().T@S@vectors
        lm = np.empty((len(p), len(p)))
        for i in range(len(p)):
            for j in range(len(p)):
                difference = float(logp[j]-logp[i])
                lm[i, j] = p[i] if abs(difference) < 1e-13 else p[i]*math.expm1(difference)/difference
        km = float(np.sum(lm*np.abs(Se)**2).real)-mean**2
        ordinary = trace_real(rho@S@S)-mean**2
        second_mean = trace_real(rho@H2)
        out.update(mean_source=mean, KM_covariance=km, ordinary_variance=ordinary,
                   mean_second_source=second_mean, F2=second_mean-beta*km,
                   ordinary_variance_wrong_F2=second_mean-beta*ordinary)
        assert km > -1e-11 and ordinary+1e-11 >= km
    return out


def twirl(rho, T):
    total = np.zeros_like(rho, dtype=complex)
    for power in range(6):
        action = np.linalg.matrix_power(T, power)
        total += action@rho@action.conj().T/6
    return total


def coherence_calibration(m):
    # A genuinely coherent state containing every deck sector.
    psi = np.arange(1, 13, dtype=float) + 1j*np.array([1, -2, 1, 3, -1, 2, -3, 1, 2, -1, 3, -2])
    psi /= np.linalg.norm(psi)
    pure = np.outer(psi, psi.conj()); averaged = twirl(pure, m["T"])
    masked = pure.copy()
    for i, j in itertools.product(range(12), repeat=2):
        if i % 6 != j % 6:
            masked[i, j] = 0
    weights_before = [float(np.trace(pure[np.ix_([r, r+6], [r, r+6])]).real) for r in range(6)]
    weights_after = [float(np.trace(averaged[np.ix_([r, r+6], [r, r+6])]).real) for r in range(6)]
    P0 = m["W"]@m["W"].T
    conditioned = P0@pure@P0/weights_before[0]
    trace_distance = float(np.sum(np.abs(np.linalg.eigvalsh(averaged-conditioned)))/2)
    kraus = [hermitian_function(np.eye(12)/2+sign*m["X6"]/4, np.sqrt) for sign in (1, -1)]
    branches = [k@averaged@k.conj().T for k in kraus]
    nonselective = sum(branches)
    p_plus = trace_real(branches[0]); selected = branches[0]/p_plus
    nonselective_weights = [float(np.trace(nonselective[np.ix_([r, r+6], [r, r+6])]).real) for r in range(6)]
    selected_weights = [float(np.trace(selected[np.ix_([r, r+6], [r, r+6])]).real) for r in range(6)]
    assert max(abs(a-b) for a, b in zip(weights_before, nonselective_weights)) < 3e-14
    assert max(abs(a-b) for a, b in zip(weights_before, selected_weights)) > 1e-3
    error = norm(averaged-masked)
    assert error < 3e-14
    assert max(abs(a-b) for a, b in zip(weights_before, weights_after)) < 3e-14
    assert abs(trace_distance-(1-weights_before[0])) < 3e-14
    return dict(phase_average_mask_error=error, initial_sector_weights=weights_before,
                final_sector_weights=weights_after, initial_purity=trace_real(pure@pure),
                twirled_purity=trace_real(averaged@averaged),
                distance_to_sector0_conditioned_state=trace_distance,
                expected_distance=1-weights_before[0],
                twirl_is_sector_postselection=False,
                coherences_within_each_sector_retained=True,
                nonselective_instrument_weights=nonselective_weights,
                selected_plus_probability=p_plus, selected_plus_weights=selected_weights,
                normalized_selective_branch_can_reweight_sectors=True,
                sector_transfer_in_any_branch=False)


def thermal_case(lam, beta):
    m = model(lam); H, S, H2, W = (m[k] for k in ("H", "S", "H2", "W"))
    full = thermal(H, beta, S, H2)
    quotient = thermal(m["Hq"], beta, m["Sq"], m["H2q"])
    sigma0 = W@quotient["rho"]@W.T
    sectors = []
    for r in range(6):
        idx = np.ix_([r, r+6], [r, r+6])
        th = thermal(H[idx], beta, S[idx], H2[idx])
        c2 = math.cos(math.pi*r/6)**2
        radius = math.sqrt(4*m["K"]**2*c2+m["v"]**2)
        closed_Z = 2*math.exp(-beta*(2*m["K"]+m["v"]))*math.cosh(beta*radius)
        closed_effect = .5-m["K"]*c2/radius*math.tanh(beta*radius)
        matrix_effect = trace_real(th["rho"]@m["L"][idx]/4)
        assert abs(th["Z"]-closed_Z) < 2e-13
        assert abs(matrix_effect-closed_effect) < 2e-13
        sectors.append(dict(sector=r, Z=th["Z"], F=th["F"], weight=th["Z"]/full["Z"],
                            source=th["mean_source"], F2=th["F2"],
                            KM_covariance=th["KM_covariance"],
                            energy_effect=matrix_effect, closed_form_energy_effect=closed_effect,
                            closed_form_Z=closed_Z, spectral_radius=radius))
    weights = np.array([r["weight"] for r in sectors])
    means = np.array([r["source"] for r in sectors])
    hessians = np.array([r["F2"] for r in sectors])
    between_variance = float(weights@(means**2)-(weights@means)**2)
    mixture = float(weights@hessians-beta*between_variance)
    P0 = W@W.T; w0 = sectors[0]["weight"]
    conditioned_error = norm(P0@full["rho"]@P0/w0-sigma0)
    # Compute D on the true rank-two support. Never add epsilon eigenvalues.
    relative_entropy = trace_real(quotient["rho"]@(quotient["log_rho"]-W.T@full["log_rho"]@W))
    distance = float(np.sum(np.abs(np.linalg.eigvalsh(full["rho"]-sigma0)))/2)
    full_read = trace_real(full["rho"]@m["L"]/4)
    quotient_read = trace_real(sigma0@m["L"]/4)
    step = 2e-4
    fplus = thermal(model(lam+step)["H"], beta)["F"]
    fminus = thermal(model(lam-step)["H"], beta)["F"]
    first_fd = (fplus-fminus)/(2*step)
    second_fd = (fplus-2*full["F"]+fminus)/step**2
    errors = dict(partition_sector_sum=abs(sum(r["Z"] for r in sectors)-full["Z"]),
                  source_sector_mixture=abs(float(weights@means)-full["mean_source"]),
                  Hessian_sector_mixture=abs(mixture-full["F2"]),
                  conditioned_state=conditioned_error,
                  thermal_twirl_invariance=norm(twirl(full["rho"], m["T"])-full["rho"]),
                  relative_entropy_identity=abs(relative_entropy+math.log(w0)),
                  trace_distance_identity=abs(distance-(1-w0)),
                  finite_difference_first=abs(first_fd-full["mean_source"]),
                  finite_difference_second=abs(second_fd-full["F2"]))
    assert abs(weights.sum()-1) < 2e-13 and np.all(weights > 0)
    assert max(v for k, v in errors.items() if not k.startswith("finite_difference")) < 2e-11
    assert errors["finite_difference_first"] < 2e-6 and errors["finite_difference_second"] < 3e-6
    assert abs(full_read-quotient_read) > 1e-7
    radius0 = math.sqrt(4*m["K"]**2+m["v"]**2)
    analytic_lower = m["K"]*math.tanh(beta*radius0)/(6*radius0*math.cosh(beta*radius0))
    assert full_read-quotient_read > analytic_lower > 0
    assert full["ordinary_variance"]-full["KM_covariance"] > 1e-4
    return dict(lambda_value=float(lam), beta=float(beta), total_K=m["K"],
                full_dimension=12, quotient_dimension=2, sectors=sectors,
                full_Z=full["Z"], quotient_Z=quotient["Z"], sector0_weight=w0,
                full_free_energy=full["F"], quotient_free_energy=quotient["F"],
                relative_entropy_supported=relative_entropy,
                expected_relative_entropy=-math.log(w0), trace_distance=distance,
                full_reference_purity=trace_real(full["rho"]@full["rho"]),
                quotient_reference_purity=trace_real(quotient["rho"]@quotient["rho"]),
                full_energy_effect_probability=full_read,
                quotient_energy_effect_probability=quotient_read,
                probability_gap=full_read-quotient_read,
                analytic_single_sector_gap_lower_bound=analytic_lower,
                full_mean_source=full["mean_source"], quotient_mean_source=quotient["mean_source"],
                full_F_second=full["F2"], quotient_F_second=quotient["F2"],
                full_KM_covariance=full["KM_covariance"],
                full_ordinary_variance=full["ordinary_variance"],
                wrong_ordinary_variance_F_second=full["ordinary_variance_wrong_F2"],
                between_sector_source_variance=between_variance,
                sector_mixture_F_second=mixture, finite_difference_step=step,
                errors=errors, support_regularization_used=False)


def finite_probability_certificate():
    """Strict rational certificate; numerical quadrature is not used here."""
    def exp_partial(x, degree):
        return sum((x**j/F(math.factorial(j)) for j in range(degree+1)), F(0))

    K, v = F(27, 10), F(2, 5)
    radius_squared = 4*K*K+v*v
    assert F(27, 5)**2 < radius_squared < F(11, 2)**2
    lower_exp = exp_partial(F(108, 25), 12)
    assert lower_exp > F(197, 3)
    # Terms after degree 8 start at term 9, and each later ratio is at most x/10.
    x = F(11, 5)
    upper_exp = exp_partial(x, 8)+(x**9/F(math.factorial(9)))/(1-x/10)
    assert upper_exp < F(19, 2)
    assert exp_partial(x, 1) > 2
    lower = F(27, 55)*F(97, 100)/30
    assert lower == F(2619, 165000)
    return dict(lambda_value="0", beta="2/5", K=str(K), v=str(v),
        R0_squared=str(radius_squared), R0_lower="27/5", R0_upper="11/2",
        exp_108_over25_partial_degree=12, exp_108_over25_lower=str(lower_exp),
        exp_108_over25_threshold="197/3",
        exp_11_over5_partial_degree=8, exp_11_over5_geometric_tail_upper=str(upper_exp),
        exp_11_over5_upper_threshold="19/2", exp_11_over5_lower_threshold="2",
        tanh_54_over25_strict_lower="97/100", cosh_11_over5_strict_upper="5",
        probability_gap_strict_lower_exact=str(lower), probability_gap_strict_lower=float(lower),
        all_bounds_exact_rational=True, finite_nonzero_signal_certified=True,
        physical_thermal_preparation_and_instrument_implementation_proved=False,
        analytic_monotonicity_obligation="x/sqrt(4*K*K*x+v*v)*tanh(beta*sqrt(4*K*K*x+v*v)) is strictly increasing on [0,1] for K,v,beta>0")


def compare(fresh, saved, path="root"):
    if isinstance(fresh, float):
        assert isinstance(saved, (int, float)) and math.isclose(fresh, saved, rel_tol=5e-10, abs_tol=3e-12), path
    elif isinstance(fresh, dict):
        assert isinstance(saved, dict) and fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert isinstance(saved, list) and len(fresh) == len(saved), path
        for i, (a, b) in enumerate(zip(fresh, saved)):
            compare(a, b, path+f"[{i}]")
    else:
        assert fresh == saved, (path, fresh, saved)


def run():
    sources = [BASE/"archive_1009_/research_note_1020.md",
        BASE/"archive_1009_/1020/passivity_reference_selection_results.json",
        BASE/"archive_531_553/research_note_531.md",
        BASE/"archive_585_628/research_note_603.md",
        BASE/"archive_585_628/research_note_617.md",
        BASE/"archive_629_652/research_note_637.md",
        BASE/"archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
        BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return dict(round=1021, new_calibration_groups=1, cumulative_test_groups=3799,
        new_cognitive_axioms=0, all_scientific_calibrations_passed=True,
        scope="finite Z12 cover / Z2 quotient, fixed three-edge ring, Gauss constraints, noncommuting magnetic term, finite beta and source derivatives",
        Gauss_and_cycle_count=gauss_calibration(),
        process_intertwining=[process_calibration(model(lam)) for lam in (0., .2)],
        coherence_twirl=coherence_calibration(model(0.)),
        finite_probability_certificate=finite_probability_certificate(),
        thermal_cases=[thermal_case(lam, beta) for lam in (0., .2) for beta in (.4, 1.2, 3.)],
        analytic_obligations=[
            "invariant quotient Hilbert space is the trivial character sector, not the full cover Hilbert space",
            "Gauss constraints allow nontrivial kernel-character cycle flows on a ring but not on a tree",
            "Hamiltonian, geometry source, waits and every descending Kraus branch intertwine on the quotient sector",
            "operator intertwining preserves arbitrary quantum inputs entangled with untouched references",
            "finite-temperature full trace includes every nonempty sector with positive weight",
            "deck twirling dephases sectors and does not condition the trivial sector",
            "relative entropy and trace distance follow from disjoint supports and Gibbs sector normalization",
            "free-energy second derivative uses Kubo-Mori covariance, with an additional between-sector term"],
        quotient_group_uniquely_selected_by_matter_kernel=False,
        quotient_reference_equal_to_full_cover_Gibbs=False,
        same_sector_representation_equivalence_preserved=True,
        all_SM_full_H_spectrum_computed=False,
        continuum_or_nonperturbative_field_theory_proved=False,
        general_global_topology_classification_proved=False,
        actual_thermalization_or_measurement_device_built=False,
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\", "/"):
                                  hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); start = time.perf_counter(); result = run()
    if args.write:
        serialized = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
        with OUT.open("x", encoding="utf8") as stream:
            stream.write(serialized)
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1021, passed=result["all_scientific_calibrations_passed"],
        thermal_cases=len(result["thermal_cases"]),
        sector0_weights=[r["sector0_weight"] for r in result["thermal_cases"]],
        readout_gaps=[r["probability_gap"] for r in result["thermal_cases"]],
        maximum_F_second_crosscheck_error=max(r["errors"]["finite_difference_second"] for r in result["thermal_cases"]),
        elapsed_seconds=time.perf_counter()-start), indent=2))
