"""Independent 1284 algebra checks. No author imports; default never writes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "independent_results.json"


def real(z):
    assert abs(complex(z).imag) < 1e-9
    return float(complex(z).real)


def hermitian(rng, n):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (z + z.conj().T) / 2


def qcharge(pi, phi, r):
    z = pi @ r @ phi
    return real(1j * (z - z.conjugate()))


def phase_hamiltonian(pi, phi, alpha, cmats, smats, kmats):
    c0 = sum(a * c for a, c in zip(alpha, cmats))
    s0 = sum(a * s for a, s in zip(alpha, smats))
    k00 = sum(alpha[a] * alpha[b] * kmats[a][b]
              for a in range(len(alpha)) for b in range(len(alpha)))
    q0 = qcharge(pi, phi, c0)
    ks = 2 * float(np.real(pi @ s0 @ phi))
    return real(pi @ pi.conj() + q0 - ks +
                phi.conj() @ ((s0 + 1j*c0) @ (s0 - 1j*c0) - k00) @ phi)


def calculate():
    rng = np.random.default_rng(12840031)
    legendre_errors, boundary_errors, primary_errors = [], [], []
    for _ in range(24):
        n, m = 3, 2
        cm = [hermitian(rng, n) for _ in range(m)]
        sm = [hermitian(rng, n) for _ in range(m)]
        km = [[None for _ in range(m)] for _ in range(m)]
        for a in range(m):
            for b in range(a, m):
                km[a][b] = km[b][a] = hermitian(rng, n)
        alpha, adot, da = (rng.normal(size=m) for _ in range(3))
        phi, vel, dphi = [rng.normal(size=n) + 1j*rng.normal(size=n)
                           for _ in range(3)]
        c0 = sum(a*c for a, c in zip(alpha, cm))
        s0 = sum(a*s for a, s in zip(alpha, sm))
        k00 = sum(alpha[a]*alpha[b]*km[a][b]
                  for a in range(m) for b in range(m))
        sq = np.array([real(phi.conj() @ s @ phi) for s in sm])
        lc = real(vel.conj() @ vel + 1j*(phi.conj() @ c0 @ vel -
                  vel.conj() @ c0 @ phi) + phi.conj() @ k00 @ phi)
        lold = lc - float(adot @ sq)
        bdot = float(adot @ sq) + 2*float(np.real(phi.conj() @ s0 @ vel))
        lnew = lold + bdot
        piold = vel.conj() + 1j * phi.conj() @ c0
        pinew = piold + phi.conj() @ s0
        p0old = -sq
        hold = 2*float(np.real(piold @ vel)) + float(p0old @ adot) - lold
        hnew = 2*float(np.real(pinew @ vel)) - lnew
        hc = phase_hamiltonian(pinew, phi, alpha, cm, sm, km)
        legendre_errors += [abs(hold-hnew), abs(hnew-hc)]
        delta_b = float(da @ sq) + 2*float(np.real(phi.conj() @ s0 @ dphi))
        theta_diff = 2*float(np.real((pinew-piold) @ dphi)) - float(p0old @ da)
        boundary_errors.append(abs(theta_diff-delta_b))
        # Independent real finite difference in A0-dot verifies the old primary.
        for a in range(m):
            step = 1e-5
            avp, avm = adot.copy(), adot.copy()
            avp[a] += step
            avm[a] -= step
            deriv = ((lc-float(avp@sq))-(lc-float(avm@sq)))/(2*step)
            primary_errors.append(abs(deriv-p0old[a]))

    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.diag([1., -1.]).astype(complex)
    ident, zero = np.eye(2, dtype=complex), np.zeros((2, 2), complex)
    r, phi = z/2, np.array([1., 0.], complex)
    shift = 0.25
    pi_c = np.array([0., -1j], complex)
    pi_s = np.array([0., 1.], complex)
    defect_c = 2*float(np.real(pi_c @ (1j*shift*x) @ phi))
    defect_s = 2*float(np.real(pi_s @ (-shift*x) @ phi))
    defect_k = 2*real(phi.conj() @ (r@r-(r@r+shift*ident)) @ phi)
    assert defect_c == .5 and defect_s == -.5 and defect_k == -.5
    # An explicitly positive but noninvariant mass potential.
    phi_v = np.array([1., 1j], complex)/np.sqrt(2.)
    mass = 2*ident+x
    defect_v = real(1j * phi_v.conj() @ (mass@r-r@mass) @ phi_v)
    h = 1e-5
    up = np.diag(np.exp(1j*h*np.array([.5, -.5])))
    um = up.conj().T
    finite_v = (real((up@phi_v).conj() @ mass @ (up@phi_v)) -
                real((um@phi_v).conj() @ mass @ (um@phi_v)))/(2*h)
    assert np.min(np.linalg.eigvalsh(mass)) > 0

    # Correct two-generator temporal coefficients: independently evaluate at random A0.
    rs = [x/2, z/2]
    ks = [[(a@b+b@a)/2 for b in rs] for a in rs]
    correct_errors = []
    for _ in range(12):
        al = rng.normal(size=2)
        ph = rng.normal(size=2) + 1j*rng.normal(size=2)
        pi = rng.normal(size=2) + 1j*rng.normal(size=2)
        expected = real(pi @ pi.conj()) + sum(a*qcharge(pi, ph, rr) for a, rr in zip(al, rs))
        got = phase_hamiltonian(pi, ph, al, rs, [zero, zero], ks)
        correct_errors.append(abs(got-expected))

    # Local Gauss completion, D E = dE/dx - i[A,E], exact constant SU(2) charge.
    a, b, length = .7, .4, .6
    connection = a*z/2
    charge_matrix = b*x/2
    phg = np.array([1., 0.], complex)
    pig = np.array([0., -1j*b], complex)
    charge_components = [qcharge(pig, phg, t/2) for t in [x, y, z]]
    assert np.allclose(charge_components, [b, 0, 0], atol=1e-15)
    nodes, weights = np.polynomial.legendre.leggauss(48)
    grid = (nodes+1)*length/2
    errors, density = [], []
    for xx in grid:
        electric = b/(2*a)*(np.sin(a*xx)*x+(np.cos(a*xx)-1)*y)
        de = b/2*(np.cos(a*xx)*x-np.sin(a*xx)*y)
        de_cov = de-1j*(connection@electric-electric@connection)
        errors.append(float(np.linalg.norm(de_cov-charge_matrix, ord=2)))
        density.append(real(np.trace(electric@electric)))
    energy_quad = float(np.dot(weights, density)*length/2)
    energy_exact = b*b/(a*a)*(length-np.sin(a*length)/a)
    assert energy_exact > 0
    # This has a nonzero endpoint flux; no claim of fixed-zero/periodic boundaries.
    endpoint_flux_norm = float(np.linalg.norm(
        b/(2*a)*(np.sin(a*length)*x+(np.cos(a*length)-1)*y), ord=2))
    assert endpoint_flux_norm > 0

    assert max(legendre_errors) < 1e-10
    assert max(boundary_errors) < 1e-11
    assert max(primary_errors) < 1e-7
    assert abs(defect_v-1) < 1e-14 and abs(finite_v-defect_v) < 1e-9
    assert max(correct_errors) < 1e-10
    assert max(errors) < 1e-13
    assert abs(energy_quad-energy_exact) < 1e-14
    return {
        "status": "PASS",
        "independence": "No author imports or numerical payload; direct original-L Legendre and boundary checks.",
        "legendre": {"random_instances": 24,
                      "max_direct_legendre_error": max(legendre_errors),
                      "max_boundary_one_form_error": max(boundary_errors),
                      "max_old_primary_difference_error": max(primary_errors)},
        "finite_witnesses": {"C_minus_R": defect_c, "nonzero_S": defect_s,
                             "wrong_K": defect_k, "noninvariant_positive_V": defect_v,
                             "V_finite_difference_error": abs(finite_v-defect_v),
                             "correct_temporal_random_instances": 12,
                             "max_correct_temporal_error": max(correct_errors)},
        "local_gauss": {"A_z_coefficient": a, "Q_x": b, "box_length": length,
                        "charge_components": charge_components,
                        "max_covariant_divergence_residual": max(errors),
                        "electric_energy_exact": energy_exact,
                        "electric_energy_quadrature_error": abs(energy_quad-energy_exact),
                        "endpoint_flux_operator_norm": endpoint_flux_norm},
        "scope": {"representation_input": True,
                  "classical_local_initial_data_only": True,
                  "periodic_or_fixed_flux_extension_certified": False,
                  "unphysical_quantum_preparations_granted": False,
                  "complete_quantum_or_microscopic_implementation_certified": False,
                  "global_Ward_from_Gauss_state_identity_claimed": False}
    }


def compare(saved, fresh):
    if isinstance(fresh, dict):
        return saved.keys() == fresh.keys() and all(compare(saved[k], v) for k, v in fresh.items())
    if isinstance(fresh, list):
        return len(saved) == len(fresh) and all(compare(a, b) for a, b in zip(saved, fresh))
    if isinstance(fresh, float):
        return bool(np.isclose(saved, fresh, rtol=1e-10, atol=1e-12))
    return saved == fresh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = calculate()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("PASS: exclusively created independent_results.json")
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert compare(saved, data), "Independent payload mismatch"
        print("PASS: independent result matches; no files written")


if __name__ == "__main__":
    main()
