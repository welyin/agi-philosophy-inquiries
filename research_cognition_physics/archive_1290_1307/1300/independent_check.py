"""Independent finite-CAR and joint-budget diagnostics for round 1300.

No author code or result is imported. These checks do not prove the continuum
Hadamard, Ward, or apparatus-realization statements; the review audits those.
Default / --check is read-only. --write creates the result exclusively.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def unitary(h, time):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * time * vals)) @ vecs.conj().T


def unitary_derivative(h, dh, time):
    vals, vecs = np.linalg.eigh(h)
    fv = np.exp(-1j * time * vals)
    divided = np.empty((len(vals), len(vals)), dtype=complex)
    for i, x in enumerate(vals):
        for j, y in enumerate(vals):
            divided[i, j] = (-1j * time * fv[i] if abs(x-y) < 1e-12
                             else (fv[i]-fv[j])/(x-y))
    return vecs @ (divided * (vecs.conj().T @ dh @ vecs)) @ vecs.conj().T


def annihilator(mode, count=4):
    op = np.zeros((1 << count, 1 << count), complex)
    for ket in range(1 << count):
        if ket & (1 << mode):
            op[ket ^ (1 << mode), ket] = (-1) ** ((ket & ((1 << mode)-1)).bit_count())
    return op


def finite_car():
    eye2 = np.eye(2)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    gen = np.kron(sy, eye2)
    pminus = np.kron((eye2-sx)/2, eye2)
    fproj = np.diag([1., 1., 0., 0.])
    jcode = np.eye(4)[:, :2]
    aa = [annihilator(i) for i in range(4)]
    nn = [a.conj().T @ a for a in aa]

    def second_quantized(a):
        return sum(a[i, j] * aa[i].conj().T @ aa[j]
                   for i in range(4) for j in range(4))

    def objects(k):
        past = unitary(gen, k)
        p = past @ pminus @ past.conj().T
        dp = -1j * (gen @ p-p @ gen)
        h = np.kron(sx+k*sz, eye2)
        dh = np.kron(sz, eye2)
        u = unitary(h, time)
        du = unitary_derivative(h, dh, time)
        q = float(p[2, 2].real)
        complement = np.diag([(1-q)**2, q*(1-q), q*(1-q), q*q])
        c0 = (np.eye(4)-fproj) @ p @ (np.eye(4)-fproj)
        dc0 = (np.eye(4)-fproj) @ dp @ (np.eye(4)-fproj)
        return p, dp, h, u, du, complement, c0, dc0

    time = .7
    kappa = .2
    p, dp, h, u, du, complement, c0, dc0 = objects(kappa)
    residuals = [np.linalg.norm(p@p-p), np.linalg.norm(p-p.conj().T)]
    f_basis = np.eye(4)[:, [1, 2]]  # fixed odd-parity single-occupation code
    hf = second_quantized(h)
    uf = unitary(hf, time)
    probabilities = []
    for code in range(2):
        rf = np.outer(f_basis[:, code], f_basis[:, code])
        state = np.kron(complement, rf)
        out = uf @ state @ uf.conj().T
        cov = u @ (c0+np.outer(jcode[:, code], jcode[:, code])) @ u.conj().T
        actual_cov = np.array([[np.trace(out @ aa[j].conj().T @ aa[i])
                                for j in range(4)] for i in range(4)])
        residuals.append(np.linalg.norm(cov-actual_cov))
        probabilities.append(float(np.trace(out @ nn[2]).real))
    predicted_gap = np.sin(np.sqrt(1+kappa*kappa)*time)**2/(1+kappa*kappa)
    residuals.append(abs(probabilities[0]-probabilities[1]-predicted_gap))

    # Exact finite-code CP replacement with a Bell partner. Retain the full
    # field and partner; the complement is not silently reset after evolution.
    full = np.zeros((32, 32), complex)
    for i in range(2):
        for j in range(2):
            fij = np.outer(f_basis[:, i], f_basis[:, j])
            rij = np.outer(eye2[:, i], eye2[:, j])/2
            full += np.kron(np.kron(complement, fij), rij)
    initial_br = np.trace(full.reshape(4, 4, 2, 4, 4, 2), axis1=1, axis2=4).reshape(8, 8)
    residuals.append(np.linalg.norm(initial_br-np.kron(complement, eye2/2)))
    joint_u = np.kron(uf, eye2)
    evolved = joint_u @ full @ joint_u.conj().T
    pr = np.kron(nn[2], eye2)
    post = pr @ evolved @ pr+(np.eye(32)-pr) @ evolved @ (np.eye(32)-pr)
    r_after = np.trace(post.reshape(16, 2, 16, 2), axis1=0, axis2=2)
    residuals.append(np.linalg.norm(r_after-eye2/2))
    residuals.append(abs(np.trace(post)-1))

    # Differentiate a genuine one-body insertion with both source and state
    # dependence. This is a finite analogue, not renormalized stress data.
    rho = eye2/2
    cinit = c0+jcode @ rho @ jcode.conj().T
    c = u @ cinit @ u.conj().T
    dc = du @ cinit @ u.conj().T+u @ dc0 @ u.conj().T+u @ cinit @ du.conj().T
    insertion = np.kron(sz+(0.3+kappa)*sx, eye2)
    dinsertion = np.kron(sx, eye2)
    state_part = float(np.trace(insertion @ dc).real)
    contact_part = float(np.trace(dinsertion @ c).real)
    derivative = state_part+contact_part

    def source_value(k):
        _, _, _, uk, _, _, ck, _ = objects(k)
        evolved_c = uk @ (ck+jcode @ rho @ jcode.conj().T) @ uk.conj().T
        ak = np.kron(sz+(0.3+k)*sx, eye2)
        return float(np.trace(ak @ evolved_c).real)

    step = 2e-5
    fd = (source_value(kappa-2*step)-8*source_value(kappa-step)
          +8*source_value(kappa+step)-source_value(kappa+2*step))/(12*step)
    residuals.append(abs(fd-derivative))
    frozen_reference_dc = du @ cinit @ u.conj().T+u @ cinit @ du.conj().T
    frozen_reference_error = abs(float(np.trace(insertion @ (dc-frozen_reference_dc)).real))
    assert abs(contact_part) > 1e-3
    assert frozen_reference_error > 1e-3
    assert max(residuals) < 1e-9
    assert min(probabilities) > 0 and predicted_gap > .3
    return {
        "scope": "independent four-mode CAR diagnostic; not continuum field numerics",
        "reference_rank": int(round(np.trace(p).real)),
        "record_probabilities": [round(x, 13) for x in probabilities],
        "exact_transport_gap": round(float(predicted_gap), 13),
        "source_derivative": round(derivative, 13),
        "omitted_insertion_contact_error": round(abs(contact_part), 13),
        "frozen_reference_error": round(frozen_reference_error, 13),
        "initial_receiver_partner_factorization": True,
        "nonselective_partner_marginal_preserved": True,
        "full_field_partner_correlations_not_declared_unchanged": True,
        "maximum_residual_below_1e_minus_9": True,
    }


def fraction_string(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def budget_checks():
    # Independent rational instance of the joint inequality, not a numerical
    # assignment to a particular Dirac spectral density or apparatus.
    margins = [F(1, 64), F(1, 32), F(1, 128)]
    lipschitz = [F(3), F(5), F(7)]
    kappa = F(1, 65536)
    fixed_error = F(1, 1024)+F(1, 4096)
    residuals = [l*kappa+fixed_error for l in lipschitz]
    remaining = [m-e for m, e in zip(margins, residuals)]
    assert all(x > 0 for x in remaining)
    kmax = min((m-fixed_error)/(2*l) for m, l in zip(margins, lipschitz))
    assert 0 < kappa < kmax

    # Einstein residual κ² B + κ A becomes (κ B + A) after dividing by κ.
    # A is a fixed additional source; it cannot be hidden in a uniform O(κ²).
    a, b = F(1, 4096), F(3)
    ks = [F(1, 2**j) for j in (12, 16, 20)]
    balances = []
    for k in ks:
        einstein = k*k*b+k*a
        thermal = einstein/k
        assert thermal == k*b+a
        assert thermal-k*b == a
        balances.append({"kappa": fraction_string(k),
                         "einstein_residual": fraction_string(einstein),
                         "thermal_residual_without_common_2pi": fraction_string(thermal),
                         "omitted_source_error": fraction_string(a)})
    assert (ks[-1]*ks[-1]*b+ks[-1]*a)/(ks[-1]*ks[-1]) > 80

    # Fixed finite test menus permit nonzero input budgets and finite records.
    # A normalized Gaussian effective kernel has a finite union-bound failure
    # probability; this does not infer a quantum POVM from its noise law.
    channels = 3
    variance = F(1, 2**40)
    noise_radius = F(1, 2**12)
    failure_bound = channels*variance/(noise_radius*noise_radius)
    assert 0 < failure_bound < F(1, 10000)
    total_energy_initial = F(2)+F(2)+F(4)+F(1)+F(1)
    total_energy_final = F(2)+F(19, 10)+F(37, 10)+F(11, 10)+F(6, 5)+F(1, 10)
    assert total_energy_initial == total_energy_final == 10
    return {
        "scope": "rational nonempty template; not actual spectral/source calibration data",
        "kappa": fraction_string(kappa),
        "strict_kappa_upper_bound": fraction_string(kmax),
        "positive_remaining_margins": [fraction_string(x) for x in remaining],
        "fixed_source_deletion": balances,
        "gaussian_union_failure_bound": fraction_string(failure_bound),
        "finite_complete_ordinary_energy": fraction_string(total_energy_initial),
        "energy_table_not_used_as_stress_projection": True,
    }


def payload():
    sources = [HERE / "../_shared/round1300_admission.md",
               HERE / "../1298/proof.md",
               HERE / "../../archive_1219_1233/1232/proof.md"]
    return {
        "status": "PASS_INDEPENDENT_FINITE_DIAGNOSTICS",
        "car": finite_car(),
        "joint_budget": budget_checks(),
        "continuum_Hadamard_PDE_proved_by_code": False,
        "actual_instrument_or_UV_qualification_proved_by_code": False,
        "source_bindings": {str(p.relative_to(HERE)).replace("\\", "/"):
                            hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "independent_check_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = payload()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("CREATED", result["status"])
    else:
        expected = json.loads(OUT.read_text(encoding="utf-8"))
        if expected != result:
            raise SystemExit("MISMATCH: saved independent diagnostic differs")
        print("READ_ONLY_MATCH", result["status"])


if __name__ == "__main__":
    main()
