"""1293: one bounded author group, exact constants + finite CAR/source diagnostics.
Default and --check are read-only. --write exclusively creates results.json.
No old checker import, continuum discretization, propagation scan or apparatus.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results.json"
TOL = 3e-12


def frac(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def encode(x):
    if isinstance(x, F):
        return frac(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    return x


def err(a, b=0):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))


def assert_close(a, b=0):
    e = err(a, b)
    assert e < TOL, e
    return e


def upper_exp(x, n):
    p = sum((x**j/F(math.factorial(j)) for j in range(n+1)), F(0))
    tail = x**(n+1)/F(math.factorial(n+1))
    return p+tail/(1-x/F(n+2))


def j1_over_x(q, n):
    # q=s^2; alternating series, odd n lower and even n upper.
    return sum(((-1)**j*q**j /
                F(2**(2*j+1)*math.factorial(j)*math.factorial(j+1))
                for j in range(n+1)), F(0))


def exact_window():
    a, L = F(1, 32), F(1, 4)
    t0, t1 = F(7, 16), F(1, 2)
    r0, r1 = L-2*a, L+2*a
    q0, q1 = t0*t0-r1*r1, t1*t1-r0*r0
    assert 0 < r0 <= r1 < t0
    assert (q0, q1) == (F(3, 32), F(55, 256))
    assert upper_exp(F(4, 3), 12) < 4   # b>1/4 on inner half-ball
    assert F(1)+1+F(1, 2) > 2           # e>2, hence sup b<1/2
    # Successive positive Bessel term ratio is q/[4(j+1)(j+2)].
    ratio_bound = q1/8
    assert ratio_bound < 1
    lower = j1_over_x(q1, 1)
    assert lower == F(1993, 4096) > F(31, 64)
    brackets = []
    for q in (q0, (q0+q1)/2, q1):
        lo, hi = j1_over_x(q, 5), j1_over_x(q, 6)
        assert lo > F(31, 64) and lo < hi
        brackets.append({"s_squared": q, "lower_odd5": lo, "upper_even6": hi})
    # Z>pi*a^3/24 and ||b||_2^2<Z/2 eliminate pi directly.
    stronger_amp = F(31)*a**3/3072
    public_amp = F(31, 2**27)
    public_gap = public_amp**2
    assert stronger_amp == F(31, 3*2**25) > public_amp
    assert public_gap == F(961, 2**54)
    grad_squared_bound = F(512)/a**2
    assert grad_squared_bound == 524288
    assert 1+grad_squared_bound < 725**2
    return {
        "m": 1, "a": a, "L": L, "t_window": [t0, t1],
        "distance_range": [r0, r1], "s_squared_range": [q0, q1],
        "exp_4_over_3_upper": upper_exp(F(4, 3), 12),
        "Bessel_term_ratio_upper": ratio_bound,
        "Bessel_uniform_lower": lower,
        "Bessel_public_lower": F(31, 64),
        "Bessel_rational_brackets": brackets,
        "amplitude_strict_lower_stronger": stronger_amp,
        "amplitude_strict_lower_public": public_amp,
        "probability_strict_lower_public": public_gap,
        "gradient_squared_strict_upper": grad_squared_bound,
        "relative_energy_strict_upper": 1450,
        "all_window_proof_not_sampling": True,
        "full_SM_remainder_budget_assigned": False,
    }


def annihilation(mode, n):
    size = 1 << n
    a = np.zeros((size, size), complex)
    bit = 1 << (n-1-mode)
    for col in range(size):
        if col & bit:
            preceding = (col >> (n-mode)).bit_count()
            a[col ^ bit, col] = (-1)**preceding
    return a


def trace_out(rho, dims, keep):
    a = rho.reshape(tuple(dims)+tuple(dims))
    current = list(dims)
    for i in reversed(range(len(dims))):
        if i not in keep:
            a = np.trace(a, axis1=i, axis2=i+len(current))
            current.pop(i)
    size = math.prod(current)
    return a.reshape(size, size)


def covariance(rho, annihilators):
    n = len(annihilators)
    return np.array([[np.trace(rho @ annihilators[j].conj().T @ annihilators[i])
                      for j in range(n)] for i in range(n)])


def car_diagnostics():
    # A correlated reference Slater state: the F modes are NOT initially
    # independent of the remote modes. All quantities are finite diagnostics.
    aa = [annihilation(i, 4) for i in range(4)]
    I16 = np.eye(16)
    residual = 0.0
    for i in range(4):
        for j in range(4):
            residual = max(residual,
                           assert_close(aa[i]@aa[j].conj().T+
                                        aa[j].conj().T@aa[i],
                                        I16 if i == j else np.zeros((16, 16))),
                           assert_close(aa[i]@aa[j]+aa[j]@aa[i]))
    u = np.array([3/5, 0, 4/5, 0], complex)
    v = np.array([0, 5/13, 0, 12/13], complex)
    cre = lambda x: sum(x[i]*aa[i].conj().T for i in range(4))
    vacuum = np.zeros(16, complex)
    vacuum[0] = 1
    psi = cre(u)@cre(v)@vacuum
    assert_close(np.vdot(psi, psi), 1)
    reference = np.outer(psi, psi.conj())
    P = np.outer(u, u.conj())+np.outer(v, v.conj())
    assert_close(covariance(reference, aa), P)
    sigmaE = trace_out(reference, [4, 4], [1])
    rhoF = trace_out(reference, [4, 4], [0])
    factorization_defect = err(reference, np.kron(rhoF, sigmaE))
    assert factorization_defect > .1

    # Fixed one-occupation, odd-parity code in the original two CAR modes.
    V = np.zeros((4, 2), complex)
    V[2, 0] = 1
    V[1, 1] = 1
    bell = np.array([1, 0, 0, 1], complex)/math.sqrt(2)
    rhoQR = np.outer(bell, bell.conj())
    VR = np.kron(V, np.eye(2))
    rhoFR = VR@rhoQR@VR.conj().T
    # F,R,E -> F,E,R, with all density matrices even.
    full = np.kron(rhoFR, sigmaE).reshape(4, 2, 4, 4, 2, 4)
    full = full.transpose(0, 2, 1, 3, 5, 4).reshape(32, 32)
    rhoR = trace_out(rhoQR, [2, 2], [1])
    remoteR = trace_out(full, [4, 4, 2], [1, 2])
    residual = max(residual, assert_close(remoteR, np.kron(sigmaE, rhoR)))
    recoveredFR = trace_out(full, [4, 4, 2], [0, 2])
    residual = max(residual, assert_close(VR.conj().T@recoveredFR@VR, rhoQR))
    rho_system = trace_out(full, [4, 4, 2], [0, 1])
    rhoQ = trace_out(rhoQR, [2, 2], [0])
    projection = np.diag([1, 1, 0, 0])
    Q = np.eye(4)-projection
    insertion = np.zeros((4, 4), complex)
    insertion[:2, :2] = rhoQ
    residual = max(residual, assert_close(covariance(rho_system, aa),
                                         Q@P@Q+insertion))
    parity = np.diag([(-1)**i.bit_count() for i in range(16)])
    residual = max(residual, assert_close(parity@rho_system, rho_system@parity))
    n0, n1 = aa[0].conj().T@aa[0], aa[1].conj().T@aa[1]
    double = float(np.trace(rho_system@n0@n1).real)
    wick_false = float(np.linalg.det(rhoQ).real)
    assert abs(double) < TOL and wick_false > .24

    # Unknown pure inputs and complex coherences retain the same remote state.
    tests = [np.array([1, 0], complex), np.array([0, 1], complex),
             np.array([1, 1], complex)/math.sqrt(2),
             np.array([1, 1j], complex)/math.sqrt(2)]
    for q in tests:
        r = np.outer(q, q.conj())
        out = np.kron(V@r@V.conj().T, sigmaE)
        residual = max(residual, assert_close(trace_out(out, [4, 4], [1]), sigmaE))
        predicted = Q@P@Q
        predicted[:2, :2] += r
        residual = max(residual, assert_close(covariance(out, aa), predicted))
    # A vacuum/one-fermion coherent qubit would violate parity; our code does not.
    bad = np.array([1, 0, 1, 0], complex)/math.sqrt(2)
    parF = np.diag([1, -1, -1, 1])
    badrho = np.outer(bad, bad.conj())
    parity_negative = err(parF@badrho, badrho@parF)
    assert parity_negative > .9
    return {
        "finite_modes": 4, "system_dimension": 16,
        "original_reference_factorization_defect": factorization_defect,
        "max_identity_residual": residual,
        "full_remote_R_product_verified": True,
        "logical_R_state_recovered_exactly_within_tolerance": True,
        "actual_double_occupation": double,
        "wrong_quasifree_double_occupation": wick_false,
        "vacuum_one_particle_code_parity_defect": parity_negative,
        "diagnostic_is_not_continuum_approximation": True,
    }


sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.diag([1, -1]).astype(complex)
pauli = [sx, sy, sz]


def mass_matrix(H, yn=1.0, ye=2.0):
    M = np.array([[yn*H[1].conjugate(), ye*H[0]],
                  [-yn*H[0].conjugate(), ye*H[1]]], complex)
    return np.block([[np.zeros((2, 2)), M],
                     [M.conj().T, np.zeros((2, 2))]])


def higgs_derivative(dH, yn=1.0, ye=2.0):
    # The exact parent mass map is real-linear in the four Higgs components.
    return mass_matrix(dH, yn, ye)


def blocks(mats):
    d = sum(len(x) for x in mats)
    out = np.zeros((d, d), complex)
    start = 0
    for x in mats:
        n = len(x)
        out[start:start+n, start:start+n] = x
        start += n
    return out


def background_sources():
    H0 = np.array([0, 1], complex)
    directions = [np.array(x, complex) for x in
                  ([1, 0], [1j, 0], [0, 1], [0, 1j])]
    weak = []
    for t in pauli:
        weak.append(blocks([t/2, np.zeros((2, 2))]))
    hyper = np.diag([-.5, -.5, 0, -1]).astype(complex)
    generators = weak+[hyper]
    residual = 0.0
    omitted_higgs_defects = []
    # Vary the original Higgs-dependent mass matrix, then take the background.
    for H in (H0, np.array([.2+.1j, .8-.15j])):
        HM = mass_matrix(H)
        for j, T in enumerate(generators):
            dH = 1j*(pauli[j]/2)@H if j < 3 else .5j*H
            dHM = higgs_derivative(dH)
            jdot = 1j*(T@HM-HM@T)  # current sign from proof (14)
            force = -dHM
            residual = max(residual, assert_close(dHM, jdot),
                           assert_close(jdot+force))
            omitted_higgs_defects.append(err(jdot))
    assert max(omitted_higgs_defects) > .9

    compressions = []
    for n in (3, 4):
        # Different spectator masses allowed; only the first two form the code.
        yn = [1, 1]+[j+1 for j in range(2, n)]
        ye = [2, 2]+[j+2 for j in range(2, n)]
        for momentum in ([0, 0, 0], [.2, -.3, .1], [-.4, .1, .25]):
            sigp = sum(x*s for x, s in zip(momentum, pauli))
            kinetic = np.kron(np.diag([-1, -1, 1, 1]), sigp)
            parents = [kinetic+np.kron(mass_matrix(H0, yn[j], ye[j]), np.eye(2))
                       for j in range(n)]
            h = blocks(parents)
            eig, basis = np.linalg.eigh(h)
            U = (basis*np.exp(-.37j*eig))@basis.conj().T
            J = np.zeros((8*n, 2), complex)
            J[4, 0] = 1                 # nu_R, fixed spin in family 1
            J[8+4, 1] = 1               # identical family 2
            evolving = U@J
            sources = {"metric_lapse": h,
                       "metric_spatial_kinetic": blocks([kinetic]*n)}
            for i, T in enumerate(generators):
                # Time and spatial background gauge-current vertices.
                sources[f"gauge_{i}_time"] = blocks([np.kron(T, np.eye(2))]*n)
                for a, sp in enumerate(pauli):
                    alpha = np.kron(np.diag([-1, -1, 1, 1]), sp)
                    source = alpha@np.kron(T, np.eye(2))
                    sources[f"gauge_{i}_space_{a}"] = blocks([source]*n)
            for i, dH in enumerate(directions):
                sources[f"Higgs_real_component_{i}"] = blocks([
                    -np.kron(higgs_derivative(dH, yn[j], ye[j]), np.eye(2))
                    for j in range(n)])
            for name, source in sources.items():
                code = evolving.conj().T@source@evolving
                residual = max(residual, assert_close(code, np.trace(code)/2*np.eye(2)))
                if name.startswith("Higgs_real_component_0") or name.startswith("Higgs_real_component_1"):
                    residual = max(residual, assert_close(code))
            compressions.append({"n": n, "momentum": momentum,
                                 "source_vertices": len(sources)})
    # Arbitrary flavor spurions are NOT protected by the common U(2).
    flavor = np.array([[0, 1], [1, 0]], complex)
    y_spurion = np.kron(flavor, higgs_derivative(np.array([0, 1], complex)))
    flavor_generator = np.kron(sz, np.eye(4))
    spurion_defect = err(y_spurion@flavor_generator-flavor_generator@y_spurion)
    assert spurion_defect > 1
    # Field identity including inherited discrete r and conjugated right fields.
    rL, rH, rN = 1j, -1, 1j
    rchiL = rH.conjugate()*rL
    rchiR = rN.conjugate()
    assert rchiL == rchiR == -1j
    return {"max_source_or_Ward_residual": residual,
            "source_tests": compressions,
            "drop_Higgs_exchange_max_Ward_defect": max(omitted_higgs_defects),
            "arbitrary_Yukawa_spurion_U2_defect": spurion_defect,
            "r_chiL_equals_r_chiR": "-i",
            "source_menu": ["metric", "all gauge-background components",
                            "four real Higgs components"],
            "couplings_varied_as_sources": False}


def dirac_kernel_symbol():
    residual = 0.0
    for m in (0.0, 1.0, 2.0):
        for p in ([.2, -.3, .4], [.1, .2, -.25]):
            sigp = sum(x*s for x, s in zip(p, pauli))
            h = np.block([[-sigp, m*np.eye(2)], [m*np.eye(2), sigp]])
            energy = math.sqrt(sum(x*x for x in p)+m*m)
            residual = max(residual, assert_close(h@h, energy**2*np.eye(4)))
            for t in (7/16, 15/32, 1/2):
                eig, basis = np.linalg.eigh(h)
                U = (basis*np.exp(-1j*t*eig))@basis.conj().T
                target = -1j*m*math.sin(t*energy)/energy*np.eye(2)
                residual = max(residual, assert_close(U[:2, 2:], target))
                if m == 0:
                    residual = max(residual, assert_close(U[:2, 2:]))
    return {"max_residual": residual, "m0_cross_chirality_exactly_zero": True,
            "symbol_diagnostic_not_spatial_propagation": True}


def source_bindings():
    paths = ["proof.md", "check.py", "../_shared/round1293_admission.md",
             "../../archive_702_741/research_note_730.md",
             "../../archive_797_818/research_note_798.md",
             "../../archive_797_818/research_note_813.md",
             "../../archive_797_818/research_note_818.md",
             "../../archive_1281_1287/1287/proof.md",
             "../../archive_1046_1062/1060/proof.md"]
    return {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}


def calculate():
    return encode({
        "schema": "round1293-dirac-field-record-source-v1",
        "round": 1293, "author_check_groups": 1, "verdict": "PASS_AUTHOR_SCOPED",
        "source_bindings": source_bindings(),
        "exact_finite_window": exact_window(),
        "complete_CAR_preparation": car_diagnostics(),
        "same_parent_background_sources": background_sources(),
        "Dirac_mass_transfer": dirac_kernel_symbol(),
        "scope": {
            "original_SM_field_leading_dictionary": True,
            "exact_quadratic_effective_model": True,
            "normal_finite_resource_initial_contract": True,
            "initial_remote_state_and_R_product_preserved": True,
            "arbitrary_unknown_QR_preserved_globally": True,
            "no_positive_energy_projection_claimed_compact": True,
            "encoded_mixed_state_quasifree": False,
            "all_fixed_g_A_H_mean_sources_code_scalar": True,
            "all_fluctuations_or_conditioned_readout_sources_equal": False,
            "autonomous_preparation_or_readout_device_certified": False,
            "all_order_SM_or_Higgs_composite_CAR_certified": False,
            "full_SM_higher_order_error_zero_assumed": False,
            "self_gravitating_Minkowski_solution_claimed": False,
            "n3_n4_entire_sources_equal_claimed": False,
            "new_axioms": 0,
        },
    })


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = calculate()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(data, f, ensure_ascii=False, sort_keys=True, indent=2)
            f.write("\n")
        receipt = "CREATED_EXCLUSIVELY"
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert saved == data, "saved result differs; no file was written"
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({
        "verdict": data["verdict"], "receipt": receipt,
        "probability_strict_lower": "961/18014398509481984",
        "relative_energy_strict_upper": 1450,
        "results_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()

