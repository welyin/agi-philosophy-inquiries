"""1153: finite calibrations of analytic tools; not a simulation of the field parent."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent

def opnorm(x):
    return float(np.linalg.norm(x, 2))

def exp_h(h, t=1.0):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j*t*e)) @ v.conj().T

def outer(v):
    return np.outer(v, v.conj())

def tensor(*matrices):
    out = np.ones((1, 1), dtype=complex)
    for matrix in matrices:
        out = np.kron(out, matrix)
    return out

def rational(x):
    return {"numerator": x.numerator, "denominator": x.denominator,
            "float": float(x)}

def exact_budgets():
    z = F(99, 56)
    lower = 1 - z*z/6 + z**4/120 - z**6/5040
    upper = F(99, 224)
    # P'(x)/x <= -1/3 + z^2/30 < 0 on [0,z]; omitted -x^4/840 is negative.
    derivative_majorant = -F(1, 3) + z*z/30
    assert derivative_majorant < 0 and lower > upper
    x2max = F(1, 10**8)
    p = F(1, 39999)
    d = 4*p
    rT = 2*F(22, 7)/(10**6*(1-x2max))
    record = d+rT
    clock = record+2*F(22, 7)/10**6
    history = 2*d+rT
    assert record < F(107, 10**6)
    assert clock < F(114, 10**6)
    assert history < F(207, 10**6)
    assert 1-2*record > F(999786, 10**6)
    assert 1-F(4, 10**6) > 0
    samples = []
    for g, C, nu in [(1., .5, .03), (.2, .001, .0001),
                     (.01, 1., .2), (.7, 1e-8, .001)]:
        D = 10/math.sqrt(g)+5/g
        R0 = 1/(4*D)
        lam = min(R0/1e4, math.sqrt(C*R0**4/(1e6*g)),
                  R0*math.sqrt(nu/(1e6*g)))
        T = math.pi/(C*lam*lam)
        q = lam*D
        p_actual = q/(1-q)
        r = 2*g*(lam/R0)**4/(1-(lam/R0)**2)
        actual_record = 4*p_actual+T*r
        omega_relative = 4*g*lam*lam/(R0*R0*nu)
        clock_to_T = C*lam*lam/(2*nu)
        assert R0 <= math.sqrt(g)/40
        assert 0 < lam <= R0/1e4 and math.isfinite(T)
        assert p_actual <= float(p)*(1+1e-12)
        assert T*r <= float(rT)*(1+1e-12)
        assert actual_record < .000107
        assert omega_relative <= 4e-6*(1+1e-12)
        assert clock_to_T <= 1/(3200*1e6)*(1+1e-12)
        samples.append(dict(g=g, C=C, nu=nu, R0=R0, coupling=lam,
                            duration=T, record_bound=actual_record,
                            omega_relative_bound=omega_relative,
                            clock_duration_fraction=clock_to_T))
    return dict(name="exact_positivity_and_finite_parameter_budgets", passed=True,
                bessel_polynomial_lower=rational(lower),
                bessel_lower_component_upper=rational(upper),
                positivity_margin=rational(lower-upper),
                derivative_majorant=rational(derivative_majorant),
                dressing_bound=rational(d), long_time_remainder_bound=rational(rT),
                record_bound=rational(record), actual_clock_bound=rational(clock),
                one_insertion_history_bound=rational(history),
                ideal_read_probability_gap_lower=rational(1-2*record),
                sample_parameters_not_physical_spectral_values=samples,
                actual_profile_C_not_numerically_evaluated=True)

def cluster_calibration():
    # Low rank-four cluster, odd auxiliary sector, and another even sector.
    # This is a finite diagnostic of the projection/Kato estimates, not the KG Fock space.
    I4 = np.eye(4)
    h0 = np.diag(np.repeat([0., 1., 2.], 4))
    parity = np.diag(np.repeat([1., -1., 1.], 4))
    P0 = np.diag(np.repeat([1., 0., 0.], 4))
    a = .7
    A = a*I4
    X = np.array([[0., 1.], [1., 0.]])
    Z = np.diag([1., -1.])
    # Non-product spin mixing makes the fourth-order matrix non-scalar.
    rotation = exp_h(.31*np.kron(X, X)+.23*np.kron(Z, X)+.19*np.kron(X, Z))
    B = rotation @ np.diag([.2, .4, .6, .8])
    V = np.zeros((12, 12), dtype=complex)
    V[:4, 4:8] = A
    V[4:8, :4] = A.conj().T
    V[4:8, 8:] = B
    V[8:, 4:8] = B.conj().T
    assert opnorm(parity@V@parity+V) < 1e-14
    assert opnorm(V) < 2.5  # Sufficient for the deliberately loose proof relative bound.
    second = -A@A.conj().T
    fourth = a**4*I4-(a*a/2)*(B@B.conj().T)
    fourth_offdiag = opnorm(fourth-np.diag(np.diag(fourth)))
    assert fourth_offdiag > .01
    R0 = 1/60
    def data(lam):
        H = h0+lam*V
        e, v = np.linalg.eigh(H)
        P = v[:, :4]@v[:, :4].conj().T
        delta = P-P0
        square = np.eye(12)-delta@delta
        ev, basis = np.linalg.eigh(square)
        invroot = (basis*(1/np.sqrt(ev)))@basis.conj().T
        S = P@P0+(np.eye(12)-P)@(np.eye(12)-P0)
        U = S@invroot
        K = (U.conj().T@H@U)[:4, :4]
        return H, P, U, K
    samples = []
    for fraction in [.25, .35, .5]:
        lam = fraction*R0
        H, P, U, K = data(lam)
        _, Pminus, Uminus, Kminus = data(-lam)
        p_actual = opnorm(P-P0)
        p_bound = lam*15/(1-lam*15)
        r_bound = 2*fraction**4/(1-fraction*fraction)
        remainder = K-lam*lam*second
        scaled = remainder/lam**4
        errors = dict(isometry=opnorm(U.conj().T@U-np.eye(12)),
                      intertwining=opnorm(U@P0@U.conj().T-P),
                      spectral_intertwining=opnorm(H@U[:, :4]-U[:, :4]@K),
                      projection_parity=opnorm(Pminus-parity@P@parity),
                      rotation_parity=opnorm(Uminus-parity@U@parity),
                      effective_evenness=opnorm(K-Kminus))
        assert max(errors.values()) < 2e-13
        assert p_actual <= p_bound and opnorm(remainder) <= r_bound
        assert opnorm((U-np.eye(12))[:, :4]) <= 2*p_actual*(1+1e-10)
        assert opnorm(scaled-fourth) < .0001
        assert opnorm(scaled-np.diag(np.diag(scaled))) > .01
        # Bare full isometry, retaining the auxiliary components.
        t = 1/(lam*lam)
        W0 = np.eye(12, 4)
        bare_error = opnorm(exp_h(H, t)@W0-W0@exp_h(lam*lam*second, t))
        tight_spectral_bound = 4*p_actual+t*opnorm(remainder)
        assert bare_error <= tight_spectral_bound+1e-12
        samples.append(dict(lambda_fraction=fraction, coupling=lam,
                            projection_distance=p_actual, projection_bound=p_bound,
                            errors=errors, remainder_norm=opnorm(remainder),
                            cauchy_remainder_bound=r_bound,
                            scaled_fourth_order_error=opnorm(scaled-fourth),
                            higher_order_spin_mixing=opnorm(scaled-np.diag(np.diag(scaled))),
                            bare_full_output_error=bare_error,
                            spectral_full_output_bound=tight_spectral_bound))
    return dict(name="rank_four_even_cluster_and_operator_remainder", passed=True,
                hilbert_dimension=12, cluster_rank=4,
                interaction_norm=opnorm(V), second_order_scalar=-a*a,
                fourth_order_offdiagonal_norm=fourth_offdiag, samples=samples,
                includes_high_order_spin_mixing=True,
                field_model_simulation=False)

def record_and_history():
    I = np.eye(2)
    n = np.diag([0., 1.])
    X = np.array([[0., 1.], [1., 0.]])
    pa = tensor(n, I, I, I)
    pb = tensor(I, I, n, I)
    xa = tensor(X, I, I, I)
    # Dimensionless time T=1 is solely a calibration of the code algebra.
    phase = .371
    common = .137
    K = common*np.eye(16)+phase*(pa+pb)-math.pi*pa@pb
    W = np.zeros((16, 8), complex)
    for a in range(2):
        for sa in range(2):
            for sb in range(2):
                col = (a*2+sa)*2+sb
                for b in range(2):
                    row = ((a*2+sa)*2+b)*2+sb
                    W[row, col] = 1/math.sqrt(2)
    assert opnorm(W.conj().T@W-np.eye(8)) < 1e-14
    U = exp_h(K)
    CZ = np.eye(16)-2*pa@pb
    local_phase = exp_h(common*np.eye(16)+phase*(pa+pb))
    assert opnorm(U-local_phase@CZ) < 1e-14
    beta = np.array([1., np.exp(-1j*phase)])/math.sqrt(2)
    EB = tensor(I, I, outer(beta), I)
    pulled_effect = W.conj().T@U.conj().T@EB@U@W
    read_A0 = tensor(np.diag([1., 0.]), I, I)
    read_error = opnorm(pulled_effect-read_A0)
    assert read_error < 1e-14
    # The isometry preserves off-diagonal unknown-input coherence, not only probabilities.
    ideal = local_phase@CZ@W
    assert opnorm(ideal.conj().T@ideal-np.eye(8)) < 1e-14
    insertion = exp_h(xa, math.pi/4)
    half = exp_h(K, .5)
    history = half@insertion@half@W
    assert opnorm(history.conj().T@history-np.eye(8)) < 1e-14
    assert opnorm(insertion@pa-pa@insertion) > .7
    history_change = opnorm(history-U@W)
    assert history_change > .7
    # All-spin identity follows from a complete operator check over the 16-dimensional code.
    spin_x_A = tensor(I, X, I, I)
    spin_x_B = tensor(I, I, I, X)
    spin_z_A = tensor(I, np.diag([1., -1.]), I, I)
    spin_z_B = tensor(I, I, I, np.diag([1., -1.]))
    spin_residual = max(opnorm(U@s-s@U) for s in
                        [spin_x_A, spin_x_B, spin_z_A, spin_z_B])
    assert spin_residual < 1e-14
    # A 32-dimensional full-output test of tails and one non-reset continuation.
    # Its exact parent is deliberately finite; it supplies no QFT spectral estimates.
    theta = .003
    Q = np.block([[math.cos(theta)*np.eye(16), -math.sin(theta)*np.eye(16)],
                  [math.sin(theta)*np.eye(16), math.cos(theta)*np.eye(16)]])
    high = np.diag(np.linspace(5., 7., 16))
    perturb = 1e-5*tensor(X, np.diag([1., -1.]), X, X)
    H = Q@np.block([[K+perturb, np.zeros((16, 16))],
                   [np.zeros((16, 16)), high]])@Q.T
    W0 = np.eye(32, 16)
    full_insert = np.block([[insertion, np.zeros((16, 16))],
                           [np.zeros((16, 16)), np.eye(16)]])
    actual = exp_h(H)@W0
    target = W0@U
    actual_history = exp_h(H, .5)@full_insert@exp_h(H, .5)@W0
    target_history = W0@half@insertion@half
    p = math.sin(theta)
    r = opnorm(perturb)
    bare_error = opnorm(actual-target)
    hist_error = opnorm(actual_history-target_history)
    assert bare_error <= 4*p+r and hist_error <= 8*p+r
    # Direct operator norm is the certificate for arbitrary passive reference dimensions.
    reference_error = opnorm(np.kron(actual-target, np.eye(3)))
    assert abs(reference_error-bare_error) < 2e-14
    tail_norm = opnorm(actual[16:, :])
    assert tail_norm > 1e-4
    return dict(name="complete_record_isometry_and_noncommuting_history", passed=True,
                code_dimension=16, input_unknown_flavor_and_spins_dimension=8,
                bare_B_pulled_effect_error=read_error,
                code_spin_operator_residual=spin_residual,
                insertion_noncommutator_norm=opnorm(insertion@pa-pa@insertion),
                inserted_history_differs_from_CZ=history_change,
                finite_full_parent_dimension=32,
                nonzero_retained_tail=tail_norm,
                bare_full_output_error=bare_error,
                bare_full_output_bound=4*p+r,
                history_full_output_error=hist_error,
                history_full_output_bound=8*p+r,
                passive_reference_norm_difference=abs(reference_error-bare_error),
                field_model_simulation=False, field_reset_performed=False,
                arbitrary_unknown_input_supported_by_operator_identity=True)

def compare(a, b, trail="root"):
    if isinstance(a, float):
        assert isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=2e-8, abs_tol=1e-11), (trail, a, b)
    elif isinstance(a, dict):
        assert a.keys() == b.keys(), trail
        for key in a:
            compare(a[key], b[key], trail+"."+key)
    elif isinstance(a, list):
        assert len(a) == len(b), trail
        for j, (x, y) in enumerate(zip(a, b)):
            compare(x, y, trail+"."+str(j))
    else:
        assert a == b, (trail, a, b)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    source_path = HERE/"source_hashes.json"
    if source_path.exists():
        for name, digest in json.loads(source_path.read_text(encoding="utf-8")).items():
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, name
    result = dict(round=1153, all_passed=True, numerical_check_groups=3,
                  groups=[exact_budgets(), cluster_calibration(), record_and_history()],
                  interpretation="Analytic-tool calibration; the full normal field result is proved in proof.md.",
                  physical_new_experiments=0, new_adopted_axioms=0,
                  full_goal_complete=False)
    dest = HERE/"results.json"
    if args.write_results and not dest.exists():
        dest.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        assert dest.exists(), "Run --write-results once to create the calibrated results."
        compare(result, json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="PASS", round=1153, groups=3,
                          source_hashes_checked=source_path.exists(),
                          record_bound=result["groups"][0]["record_bound"]["float"],
                          clock_bound=result["groups"][0]["actual_clock_bound"]["float"],
                          history_bound=result["groups"][0]["one_insertion_history_bound"]["float"])))
if __name__ == "__main__":
    main()
