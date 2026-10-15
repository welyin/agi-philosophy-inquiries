"""Independent 1233 diagnostics, built without the author's implementation.

Finite CAR matrices check algebra only; they are not a cutoff of the Dirac sea
or of the Gaussian program clocks. The continuous domain argument is reviewed
separately in independent_review.md. --write-results is create-only.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "independent_results.json"
SOURCES = ("archive_1219_/1230/proof.md", "archive_1219_/1228/proof.md")


def ensure(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packq(x):
    x = F(x)
    return {"exact": str(x), "display": float(x)}


def display(x):
    return round(float(x), 11)


def norm(x):
    return float(np.linalg.norm(x))


def comm(a, b):
    return a @ b - b @ a


def kron(*xs):
    return reduce(np.kron, xs)


def calculate():
    # Four modes and two buffer sites keep A_chi f_j = h f_j or zero exactly.
    nm, df = 4, 16
    eye_f, eye_c = np.eye(df), np.eye(2)
    c = []
    for j in range(nm):
        op = np.zeros((df, df), complex)
        for bits in range(df):
            if (bits >> j) & 1:
                sign = (-1)**((bits & ((1 << j)-1)).bit_count())
                op[bits ^ (1 << j), bits] = sign
        c.append(op)
    cd = [v.conj().T for v in c]
    n = [cd[j] @ c[j] for j in range(nm)]
    h = np.diag([-1.3, -0.2, 0.8, 1.7])
    h += np.diag([-0.5, 0.35, -0.4], 1)
    h += np.diag([-0.5, 0.35, -0.4], -1)
    eps, eigenvectors = np.linalg.eigh(h)
    occupied = eigenvectors[:, eps < 0]
    sea = occupied @ occupied.conj().T

    def dgamma(b):
        ans = np.zeros((df, df), complex)
        for j in range(nm):
            for k in range(nm):
                ans += b[j, k] * cd[j] @ c[k]
        return ans

    def normal(b):
        return dgamma(b) - np.trace(sea @ b)*eye_f

    vacuum = np.zeros(df, complex)
    vacuum[0] = 1
    for u in occupied.T:
        creator = sum((u[j]*cd[j] for j in range(nm)),
                      np.zeros((df, df), complex))
        vacuum = creator @ vacuum
    hf = normal(h)
    max_car = 0.0
    for j in range(nm):
        for k in range(nm):
            max_car = max(max_car, norm(c[j]@cd[k]+cd[k]@c[j]
                           - (eye_f if j == k else 0*eye_f)),
                          norm(c[j]@c[k]+c[k]@c[j]))
    ensure(max_car < 1e-12, "CAR")
    ensure(abs(np.vdot(vacuum, vacuum)-1) < 1e-12, "sea normalization")
    ensure(norm(hf@vacuum) < 1e-12, "normal-ordered vacuum energy")
    ensure(np.linalg.eigvalsh(hf).min() > -1e-12, "positive excitation energy")

    trial_b = np.array([[1., .2j, .3, 0], [-.2j, .6, .1j, .2],
                        [.3, -.1j, -.4, .2j], [0, .2, -.2j, .9]])
    u = np.array([1, 1j, -2, .5j], complex)
    u /= np.linalg.norm(u)
    cu = sum((np.conj(u[j])*c[j] for j in range(nm)),
             np.zeros((df, df), complex))
    cbu = sum((np.conj((trial_b@u)[j])*c[j] for j in range(nm)),
              np.zeros((df, df), complex))
    lhs = comm(normal(trial_b), cu.conj().T@cu)
    rhs = cbu.conj().T@cu - cu.conj().T@cbu
    ensure(norm(lhs-rhs) < 1e-12, "CAR bounded commutator sign")
    ensure(np.linalg.norm(lhs, 2) <= 2*np.linalg.norm(trial_b@u)+1e-12,
           "commutator norm guard")
    g1 = {
        "passed": True, "modes": nm, "Fock_dimension": df,
        "one_particle_eigenvalues": [display(v) for v in eps],
        "filled_negative_modes": int((eps < 0).sum()),
        "CAR_and_normal_order_residuals": "<1e-12",
        "bounded_commutator_sign_checked": True,
        "continuous_trace_class_limit_proven_by_matrices": False,
    }

    # Exact logical reduction of complete four-dimensional P/E pairs.
    # K is zero on |01>,|10>; these sectors remain in the full device space.
    embed = np.zeros((4, 2))
    embed[0, 0] = embed[3, 1] = 1
    X = np.array([[0., 1.], [1., 0.]])
    K = np.pi/2 * (eye_c-X)
    Kfull = embed@K@embed.T
    inactive = np.eye(4)-embed@embed.T
    ensure(norm(Kfull@inactive) == 0, "full P/E inactive sectors")
    ensure(norm(embed.T@Kfull@embed-K) < 1e-14, "P/E reduction")
    hc1, hc2 = np.diag([1., 2.2]), np.diag([1.4, 2.3])
    gc1 = np.array([[1.1, .3], [.3, .6]])
    gc2 = np.array([[.7, -.2], [-.2, 1.2]])
    for mat in (hc1, hc2, gc1, gc2, K):
        ensure(np.linalg.eigvalsh(mat).min() > -1e-12, "positive control/interaction")
    HF = kron(hf, eye_c, eye_c, eye_c, eye_c)
    HC1 = kron(eye_f, hc1, eye_c, eye_c, eye_c)
    HC2 = kron(eye_f, eye_c, hc2, eye_c, eye_c)
    V1 = kron(n[0], gc1, eye_c, K, eye_c)
    V2 = kron(n[3], eye_c, gc2, eye_c, K)
    V = V1+V2
    dim = HF.shape[0]
    H = HF+HC1+HC2+V+4*np.eye(dim)

    chi1, chi2 = np.diag([1, 1, 0, 0]), np.diag([0, 0, 1, 1])
    aa = [(z@h+h@z)/2 for z in (chi1, chi2)]
    jj = [1j*comm(h, a) for a in aa]
    targets = [np.eye(4)[:, 0], np.eye(4)[:, 3]]
    for j in range(2):
        ensure(norm(aa[j]@targets[j]-h@targets[j]) < 1e-14, "own platform")
        ensure(norm(aa[j]@targets[1-j]) < 1e-14, "other platform")
        ensure(norm(jj[j]-.5j*comm(h@h, (chi1,chi2)[j])) < 1e-14,
               "J=i[h^2,chi]/2")
    lift = lambda b: kron(normal(b), eye_c, eye_c, eye_c, eye_c)
    QA = [lift(a) for a in aa]
    QJ = [lift(j) for j in jj]
    src = [1j*comm(v, HF) for v in (V1, V2)]
    controls = [HC1+V1, HC2+V2]
    ensure(norm(comm(V1, V2)) < 1e-12, "same-time distinct source gates commute")
    ensure(norm(QA[0]+QA[1]-HF) < 1e-12, "regional sum")
    for j in range(2):
        ensure(norm(1j*comm(H, QA[j])-QJ[j]-src[j]) < 1e-10,
               "local source and flux identity")
        ensure(norm(1j*comm(H, controls[j])+src[j]) < 1e-10,
               "control plus interaction reaction")
        ensure(norm(1j*comm(H, QA[j]+controls[j])-QJ[j]) < 1e-10,
               "region with control closes")
    g2 = {
        "passed": True, "finite_region_1_sites": [0, 1],
        "finite_region_2_sites": [2, 3], "source_modes": [0, 3],
        "full_device_dimension": 16*2*2*4*4,
        "exact_invariant_initial_block_dimension": dim,
        "four_dimensional_PE_pairs_kept": True,
        "reduction_reason": "K annihilates both unused P/E sectors; blank lies in invariant logical block",
        "all_state_operator_identities_residual": "<1e-10",
        "same_time_gate_commutation_only": True,
        "claims_different_time_commutation": False,
    }

    initial = kron(vacuum, np.array([1, 1j])/np.sqrt(2),
                   np.array([np.sqrt(.3), np.sqrt(.7)]),
                   np.array([1., 0.]), np.array([1., 0.]))
    evals, evecs = np.linalg.eigh(H)
    initial_e = evecs.conj().T@initial
    ensure(evals.min() > 6.39, "full finite H positive")
    ensure(norm(H@evecs-evecs*evals) < 1e-10, "eigensolver residual")

    def state(t):
        return evecs@(np.exp(-1j*t*evals)*initial_e)

    def expect(op, psi):
        value = np.vdot(psi, op@psi)
        ensure(abs(value.imag) < 2e-11, "Hermitian expectation")
        return float(value.real)

    base_H = expect(H, initial)
    base_V = expect(V, initial)
    base_HC = expect(HC1+HC2, initial)
    history = []
    max_source = max_missing_V = 0.0
    for time in (0., .15, .55, 1.1, 1.7):
        psi = state(time)
        delta_hf = expect(HF, psi)-expect(HF, initial)
        delta_hc = expect(HC1+HC2, psi)-base_HC
        delta_v = expect(V, psi)-base_V
        sj = [expect(v, psi) for v in src]
        flow = [expect(v, psi) for v in QJ]
        ensure(abs(expect(H, psi)-base_H) < 1e-10, "complete energy conservation")
        ensure(abs(delta_hf+delta_hc+delta_v) < 1e-10, "full energy balance")
        max_source = max(max_source, abs(sj[0]), abs(sj[1]))
        max_missing_V = max(max_missing_V, abs(delta_hf+delta_hc))
        history.append({
            "time": time, "field_energy_change": display(delta_hf),
            "controls_HC_energy_change": display(delta_hc),
            "interaction_energy_change": display(delta_v),
            "source_region_1": display(sj[0]), "source_region_2": display(sj[1]),
            "inward_flux_1": display(flow[0]), "inward_flux_2": display(flow[1]),
        })
    g3 = {
        "passed": True, "full_energy_conservation_residual": "<1e-10",
        "initial_interaction_energy_retained": display(base_V),
        "minimum_complete_H_eigenvalue": display(evals.min()),
        "history": history,
        "source_factorization_used": False,
        "reference_free_future_substituted": False,
    }

    # Finite spectral regulator analogue, with a genuinely evolving full state.
    # This is NOT the continuum momentum cutoff, but tests its CAR/normal-order algebra.
    tfinal = 1.7
    nodes, weights = np.polynomial.legendre.leggauss(32)
    times = (nodes+1)*tfinal/2
    evolved = [state(t) for t in times]
    final = state(tfinal)
    regulators = []
    for threshold in (.5, 1., 1.75, 3.):
        keep = abs(eps) <= threshold
        proj = eigenvectors[:,keep]@eigenvectors[:,keep].conj().T
        ar = proj@aa[0]@proj
        jr = proj@jj[0]@proj
        ensure(norm(1j*comm(h,ar)-jr) < 1e-12, "regulated J")
        ensure(abs(np.trace(sea@jr)) < 1e-12, "normal-order derivative constant")
        qar, qjr = lift(ar), lift(jr)
        sr = 1j*comm(V, qar)
        integrand = qjr+sr
        integral = tfinal/2*sum(float(w)*expect(integrand,p)
                               for w,p in zip(weights,evolved))
        change = expect(qar,final)-expect(qar,initial)
        ensure(abs(integral-change) < 1e-10, "finite weak/Duhamel balance")
        regulators.append({
            "threshold": threshold, "rank": int(keep.sum()),
            "energy_change": display(change), "integrated_flux_plus_source": display(integral),
            "residual": "<1e-10",
        })
    g4 = {
        "passed": True, "finite_regulator_analogue": "spectral projection of finite h",
        "regulators": regulators,
        "quadrature": "32-node Gauss-Legendre, diagnostic only",
        "proves_continuum_trace_class_or_cutoff_rate": False,
    }

    ensure(max_source > .07, "nonzero source deletion witness")
    ensure(max_missing_V > .22, "interaction-energy deletion witness")
    psi = state(.55)
    no_source_gap = expect(1j*comm(H, QA[0])-QJ[0], psi)
    no_interaction_balance = (expect(HF,psi)-expect(HF,initial)
                              +expect(HC1+HC2,psi)-base_HC)
    ensure(abs(no_source_gap) > .07, "fixed-time omitted-source failure")
    ensure(abs(no_interaction_balance) > .22, "fixed-time omitted-V failure")
    g5 = {
        "passed": True, "fixed_witness_time": .55,
        "omitting_source_local_derivative_residual": display(no_source_gap),
        "omitting_interaction_complete_energy_residual": display(no_interaction_balance),
        "full_equations_hold_at_same_time": True,
        "deletion_is_finite_algebra_not_QFT_counterexample": True,
    }

    # Exact geometry and C^2 guards for the ORIGINAL continuous packets.
    gamma = F(1, 10**8)
    half_xy = F(1,5)+gamma
    half_z = F(1,100)+gamma
    Dlo, Dhi = 1-gamma, 1+gamma
    small_xy_margin = F(3,10)-half_xy
    small_z_margin = F(1,10)-(half_z+gamma)
    ensure(small_xy_margin > 0 and small_z_margin > 0, "platform contains full packets")
    ensure(F(11,5) < 3-half_z-gamma, "first support disjoint from second packet")
    # theta'=-30u^2(1-u)^2; |theta'|<=15/8.
    derivative1 = F(15,8)
    # theta''=60u(1-u)(2u-1), max^2=100/3<36.
    second_derivative_sq = F(100,3)
    ensure(second_derivative_sq < 36, "theta second derivative guard")
    sum_inv_width_sq = F(25)+F(25)+F(100)
    c1_sq = derivative1**2*sum_inv_width_sq
    c2 = 6*sum_inv_width_sq
    ensure(c1_sq < 24**2 and c2 == 900, "C2 platform constants")
    momentum_bound_sq = 1+F(44,5)*(2/(F(1,10)-gamma)**2
                                  +1/(F(1,100)-gamma)**2)
    ensure(momentum_bound_sq < 400**2, "original packet h-domain norm")
    propagation = F(2)
    halo_xy_margin = 3-(half_xy+propagation)
    halo_z_low = 2-half_z-propagation
    halo_z_high = 2+Dhi+half_z+propagation
    ensure(halo_xy_margin > 0 and halo_z_low > F(-1,2)
           and halo_z_high < F(11,2), "all free propagation halos fit large platform")
    bA, bJ = F(400)+F(24,2), F(24)*400+F(900,2)
    # Endpoint values and first two derivatives of theta match the constants.
    theta_coeff = [F(1), F(0), F(0), F(-10), F(15), F(-6)]
    deriv1_coeff = [F(k)*theta_coeff[k] for k in range(1,6)]
    deriv2_coeff = [F(k)*deriv1_coeff[k] for k in range(1,5)]
    ensure(sum(theta_coeff)==0 and sum(deriv1_coeff)==0 and sum(deriv2_coeff)==0,
           "C2 right endpoint")
    ensure(deriv1_coeff[0]==0 and deriv2_coeff[0]==0, "C2 left endpoint")
    g6 = {
        "passed": True, "geometry_error": packq(gamma),
        "small_platform_xy_margin": packq(small_xy_margin),
        "second_platform_z_margin_at_fixed_center_3": packq(small_z_margin),
        "gradient_squared_upper": packq(c1_sq), "gradient_guard": 24,
        "theta_second_derivative_squared_max": packq(second_derivative_sq),
        "laplacian_guard": packq(c2),
        "packet_h_norm_squared_upper": packq(momentum_bound_sq),
        "packet_h_norm_guard": 400,
        "halo_xy_margin": packq(halo_xy_margin),
        "halo_z_interval": [packq(halo_z_low),packq(halo_z_high)],
        "large_platform_z_interval": ["-1/2","11/2"],
        "b_A": packq(bA), "b_J": packq(bJ),
        "C2_implies_quantified_cutoff_rate": False,
    }

    groups = {
        "CAR_normal_order": g1,
        "regional_operator_reaction": g2,
        "complete_autonomous_history": g3,
        "finite_regulator_weak_identity": g4,
        "deletion_witnesses": g5,
        "continuous_platform_and_halo_guards": g6,
    }
    ensure(all(v["passed"] for v in groups.values()), "all groups")
    return {
        "round":1233, "status":"PASS", "all_checks_passed":True,
        "check_groups":len(groups), "groups":groups,
        "scope":{
            "continuous_Dirac_proof_is_separate":True,
            "finite_diagnostic_is_QFT_cutoff":False,
            "finite_controls_are_Gaussian_clock_approximation":False,
            "actual_full_H_history_used":True,
            "new_region_energy_instrument_constructed":False,
            "whole_apparatus_four_momentum_or_boost_signed":False,
            "area_dynamics_or_Einstein_signed":False,
        },
        "independent_code_sha256":sha(Path(__file__)),
        "source_sha256":{name:sha(ROOT/name) for name in SOURCES},
    }


def main():
    args=sys.argv[1:]
    if args not in ([],["--write-results"]):
        raise SystemExit("usage: independent_check.py [--write-results]")
    current=calculate()
    if args:
        with RECEIPT.open("x",encoding="utf-8",newline="\n") as handle:
            json.dump(current,handle,ensure_ascii=False,indent=2,allow_nan=False)
            handle.write("\n")
        readonly=False
    else:
        ensure(json.loads(RECEIPT.read_text(encoding="utf-8"))==current,
               "read-only receipt mismatch")
        readonly=True
    print(json.dumps({"round":1233,"status":"PASS","all_checks_passed":True,
                      "results_match":True,"check_groups":6,"read_only":readonly}))


if __name__=="__main__":
    main()
