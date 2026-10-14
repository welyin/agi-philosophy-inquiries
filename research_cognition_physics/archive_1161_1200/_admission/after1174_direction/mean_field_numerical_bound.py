"""Finite-time Hermite residual enclosure for the fixed Hartree preview.

No files are written. Uses outward basic-operation intervals and a standard
gamma_n reduction bound; this is a transparent floating-point error model,
not a formally verified IEEE implementation. See the accompanying note.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np

UROUND = 2.0 ** -53


def down(x):
    return np.nextafter(np.asarray(x, float), -np.inf)


def up(x):
    return np.nextafter(np.asarray(x, float), np.inf)


def pt(x):
    a = np.asarray(x, float)
    return a, a


def input_interval(x):
    """Includes the exact displayed decimal parameter around its float value."""
    return down(x), up(x)


def add(a, b):
    return down(a[0] + b[0]), up(a[1] + b[1])


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    values = np.array([a[0] * b[0], a[0] * b[1],
                       a[1] * b[0], a[1] * b[1]])
    return down(values.min(axis=0)), up(values.max(axis=0))


def scale(a, b):
    return mul(a, pt(b))


def divpos(a, b):
    assert np.all(b[0] > 0)
    return mul(a, (down(1.0 / b[1]), up(1.0 / b[0])))


def take(a, key):
    return a[0][key], a[1][key]


def stack(a):
    return np.stack([x[0] for x in a]), np.stack([x[1] for x in a])


def sum0(a):
    """Arbitrary-association float sum enclosed by a gamma_(n-1) model."""
    n = a[0].shape[0]
    gamma = float(up((max(n - 1, 0) * UROUND)
                     / (1 - max(n - 1, 0) * UROUND)))
    lo, hi = np.sum(a[0], axis=0), np.sum(a[1], axis=0)
    ab = np.sum(np.maximum(np.abs(a[0]), np.abs(a[1])), axis=0)
    ab = up(ab / (1 - max(n - 1, 0) * UROUND))
    err = up(gamma * ab)
    return down(lo - err), up(hi + err)


def sqrt_upper(a):
    assert np.all(a >= 0)
    return up(np.sqrt(a))


def abs_upper(a):
    return np.maximum(np.abs(a[0]), np.abs(a[1]))


def norm2_upper(a, weights=None):
    v = abs_upper(a)
    squares = up(v * v)
    if weights is not None:
        squares = mul(pt(squares), weights)[1]
    return float(sqrt_upper(sum0(pt(squares.reshape(-1)))[1]))


def hermite(y0, y1, f0, f1, dt):
    """The exact Hermite curve of stored binary endpoint/tangent data."""
    zero, one, d0, d1 = pt(y0), pt(y1), scale(pt(f0), dt), scale(pt(f1), dt)
    diff = sub(one, zero)
    return stack([zero, d0, sub(scale(diff, 3), add(scale(d0, 2), d1)),
                  add(scale(diff, -2), add(d0, d1))])


def derivative(c, dt):
    rows = [divpos(scale(take(c, j + 1), j + 1), pt(dt)) for j in range(3)]
    rows += [pt(np.zeros_like(c[0][0])) for _ in range(4)]
    return stack(rows)


def product(a, b):
    rows = []
    for n in range(a[0].shape[0] + b[0].shape[0] - 1):
        terms = [mul(take(a, j), take(b, n - j))
                 for j in range(a[0].shape[0])
                 if 0 <= n - j < b[0].shape[0]]
        total = terms[0]
        for term in terms[1:]:
            total = add(total, term)
        rows.append(total)
    return stack(rows)


def bernstein(c):
    degree = c[0].shape[0] - 1
    rows = []
    for j in range(degree + 1):
        total = pt(np.zeros_like(c[0][0]))
        for k in range(j + 1):
            ratio = divpos(pt(math.comb(j, k)), pt(math.comb(degree, k)))
            total = add(total, mul(take(c, k), ratio))
        rows.append(total)
    return stack(rows)


def pad7(c):
    return stack([take(c, j) if j < c[0].shape[0]
                  else pt(np.zeros_like(c[0][0])) for j in range(7)])


def sin_interval(angle):
    # |angle| <= pi. Twenty-two terms, next Taylor term < 4e-34.
    term, total = angle, angle
    squared = mul(angle, angle)
    for j in range(1, 22):
        term = divpos(neg(mul(term, squared)), pt((2 * j) * (2 * j + 1)))
        total = add(total, term)
    remainder = 4e-34
    return sub(total, pt(remainder))[0], add(total, pt(remainder))[1]


def exact_source_intervals(grid, length):
    assert grid == 65 and length == 16.0, "Fixed candidate only."
    # These binary endpoints straddle mathematical pi.
    pi = (np.array(3.141592653589793), np.array(3.1415926535897936))
    idx = np.arange(33, dtype=float)
    angle = divpos(scale(mul(pi, pt(idx)), 2), pt(grid))
    sins = sin_interval(angle)
    dx = divpos(pt(length), pt(grid))
    sval, counts = [], []
    for trip in itertools.combinations_with_replacement(range(33), 3):
        sq = pt(0.0)
        for t in trip:
            v = take(sins, t)
            lo = max(0.0, float(v[0]))
            hi = max(abs(float(v[0])), abs(float(v[1])))
            sq = add(sq, (down(lo * lo), up(hi * hi)))
        root = (down(np.sqrt(max(0.0, float(sq[0])))),
                up(np.sqrt(max(0.0, float(sq[1])))))
        sval.append(divpos(root, dx))
        multiplicity = math.factorial(3)
        for t in set(trip):
            multiplicity //= math.factorial(trip.count(t))
        multiplicity *= 2 ** sum(t != 0 for t in trip)
        counts.append(multiplicity)
    assert sum(counts) == grid ** 3
    weights = divpos(pt(np.asarray(counts)), pt(grid ** 3))
    return stack(sval), weights, sins, dx


def bloch_residual(r0, r1, g0, g1, fr0, fr1, fg0, fg1,
                   dt, ss, weights, eps, nu, mass, omega):
    rr = hermite(r0, r1, fr0, fr1, dt)
    gg = hermite(g0, g1, fg0, fg1, dt)
    gz = take(gg, (slice(None), 2))
    zz = (gz[0][:, None], gz[1][:, None])
    eps_iv, nu_iv = input_interval(eps), input_interval(nu)
    dx = mul(mul(zz, nu_iv), (ss[0][None, :], ss[1][None, :]))
    dx = (dx[0].copy(), dx[1].copy())
    dx0 = add(take(dx, 0), ss)
    dx[0][0], dx[1][0] = dx0
    dz = mul(zz, scale(eps_iv, mass))
    dz = (np.broadcast_to(dz[0], dx[0].shape).copy(),
          np.broadcast_to(dz[1], dx[1].shape).copy())
    dz0 = add(take(dz, 0), pt(mass))
    dz[0][0], dz[1][0] = dz0
    rx = take(rr, (slice(None), slice(None), 0))
    ry = take(rr, (slice(None), slice(None), 1))
    rz = take(rr, (slice(None), slice(None), 2))
    rhs_r = stack([
        scale(product(dz, ry), -2),
        scale(sub(product(dz, rx), product(dx, rz)), 2),
        scale(product(dx, ry), 2),
    ])
    rhs_r = (np.moveaxis(rhs_r[0], 0, -1), np.moveaxis(rhs_r[1], 0, -1))
    bcoef = []
    for j in range(4):
        f_dot_r = add(mul(take(rz, j), scale(eps_iv, mass)),
                      mul(mul(ss, nu_iv), take(rx, j)))
        bcoef.append(scale(sum0(mul(weights, f_dot_r)), 2))
    bc = stack(bcoef)
    gx = take(gg, (slice(None), 0))
    gy = take(gg, (slice(None), 1))
    rhs_g = stack([
        scale(product(bc, gy), -2),
        scale(add(product(bc, gx), scale(pad7(gz), omega)), 2),
        scale(pad7(gy), -2 * omega),
    ])
    rhs_g = (rhs_g[0].T, rhs_g[1].T)
    residual_r = bernstein(sub(derivative(rr, dt), rhs_r))
    residual_g = bernstein(sub(derivative(gg, dt), rhs_g))
    residual = 0.0
    for j in range(7):
        rg = norm2_upper(take(residual_g, j))
        rv = abs_upper(take(residual_r, j))
        local_sq = sum0(pt(up(rv * rv).T))[1]
        weighted_sq = sum0(mul(weights, pt(local_sq)))[1]
        squared_bound = add(mul(pt(rg), pt(rg)), scale(pt(weighted_sq), 2))[1]
        residual = max(residual, float(sqrt_upper(squared_bound)))
    rb = abs_upper(bernstein(rr))
    rnorm = sqrt_upper(sum0(pt(up(rb * rb).transpose(2, 0, 1)))[1])
    return residual, float(rnorm.max()), rr, gg


def positive_factors(z):
    """Upper exp(z), phi1(z), phi2(z), using positive Taylor tails."""
    assert 0 <= z < 0.1
    values = []
    for offset in range(3):
        term = float(up(1 / math.factorial(offset)))
        total = term
        for n in range(1, 13):
            term = float(up(up(term * z) / (n + offset)))
            total = float(up(total + term))
        following = float(up(up(term * z) / (13 + offset)))
        tail = float(up(following / down(1 - up(z / (14 + offset)))))
        values.append(float(up(total + tail)))
    return values


def cmatmul(ar, ai, br, bi):
    def realmm(a, b):
        shape = np.broadcast_shapes(a[0].shape[:-2], b[0].shape[:-2]) + (4, 4)
        result = pt(np.zeros(shape))
        for j in range(4):
            av = (a[0][..., :, j, None], a[1][..., :, j, None])
            bv = (b[0][..., None, j, :], b[1][..., None, j, :])
            result = add(result, mul(av, bv))
        return result
    return sub(realmm(ar, br), realmm(ai, bi)), add(realmm(ar, bi), realmm(ai, br))


def probe_residual(u0, u1, fu0, fu1, gg, dt, kin, masses, beta, eps, nu):
    # All matrix data supplied as real intervals, except stored Hermite nodes.
    ur = hermite(u0.real, u1.real, fu0.real, fu1.real, dt)
    ui = hermite(u0.imag, u1.imag, fu0.imag, fu1.imag, dt)
    kr, ki = kin
    mass_iv = input_interval(masses[None, :, None, None])
    mass_r = mul(mass_iv, pt(beta[None, None, :, :].real))
    mass_i = mul(mass_iv, pt(beta[None, None, :, :].imag))
    eps_iv, nu_iv = input_interval(eps), input_interval(nu)
    h0r, h0i = add(kr, mass_r), add(ki, mass_i)
    hzr = add(mul(kr, nu_iv), mul(mass_r, eps_iv))
    hzi = add(mul(ki, nu_iv), mul(mass_i, eps_iv))
    hr, hi = [], []
    for j in range(4):
        z = take(gg, (j, 2))
        real, imag = mul(z, hzr), mul(z, hzi)
        if j == 0:
            real, imag = add(real, h0r), add(imag, h0i)
        hr.append(real)
        hi.append(imag)
    pr, pi = [], []
    for n in range(7):
        real = pt(np.zeros_like(u0.real))
        imag = pt(np.zeros_like(u0.real))
        for j in range(4):
            if 0 <= n - j < 4:
                re, im = cmatmul(hr[j], hi[j], take(ur, n - j), take(ui, n - j))
                real, imag = add(real, re), add(imag, im)
        # i H U = -Im(HU) + i Re(HU)
        pr.append(neg(imag))
        pi.append(real)
    residual_r = bernstein(add(derivative(ur, dt), stack(pr)))
    residual_i = bernstein(add(derivative(ui, dt), stack(pi)))
    ar, ai = abs_upper(residual_r), abs_upper(residual_i)
    squares = up(up(ar * ar) + up(ai * ai))
    # Frobenius majorizes each 4x4 block's operator norm.
    sums = sum0(pt(squares.reshape(7, 4, 2, 16).transpose(3, 0, 1, 2)))[1]
    return float(sqrt_upper(sums).max())


def run(steps=2048):
    assert steps > 0 and steps & (steps - 1) == 0, "Dyadic fixed time mesh."
    path = Path(__file__).parent / "mean_field_candidate" / "preview.py"
    spec = importlib.util.spec_from_file_location("fixed_mean_field_preview", path)
    preview = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(preview)
    ss, weights, sins, dx = exact_source_intervals(preview.GRID, preview.LENGTH)
    assert np.all(preview.SVAL >= ss[0]) and np.all(preview.SVAL <= ss[1])
    eps_iv, nu_iv = input_interval(preview.EPS), input_interval(preview.NU)
    mass_eps = scale(eps_iv, preview.MS)
    f2 = add(mul(mul(ss, ss), mul(nu_iv, nu_iv)), mul(mass_eps, mass_eps))
    F = float(sqrt_upper(sum0(mul(weights, f2))[1]))
    # Four modes are 0 and the three positive coordinate unit Fourier modes.
    kval = divpos(take(sins, 1), dx)
    k_lo = np.zeros_like(preview.KIN.real)
    k_hi = np.zeros_like(k_lo)
    ki_lo, ki_hi = k_lo.copy(), k_hi.copy()
    for j in range(1, 4):
        kmat = scale(kval, preview.ALPHA[j - 1].real)
        kimat = scale(kval, preview.ALPHA[j - 1].imag)
        k_lo[j], k_hi[j] = kmat
        ki_lo[j], ki_hi[j] = kimat
    kin = ((k_lo[:, None], k_hi[:, None]), (ki_lo[:, None], ki_hi[:, None]))
    r = np.zeros((len(ss[0]), 3))
    r[:, 2] = 1
    state = (r, np.array([1., 0., 0.]), np.zeros(2),
             np.tile(np.eye(4, dtype=complex), (4, 2, 1, 1)))
    dt = preview.TIME / steps
    error, integral_error = 0., 0.
    max_residual, max_rnorm, max_L = 0., 1., 0.
    probe_defect_integral = 0.
    gz_integral = pt(0.)
    for _ in range(steps):
        k1 = preview.derivative(state, preview.NU, False)
        k2 = preview.derivative(preview.plus(state, k1, dt / 2), preview.NU, False)
        k3 = preview.derivative(preview.plus(state, k2, dt / 2), preview.NU, False)
        k4 = preview.derivative(preview.plus(state, k3, dt), preview.NU, False)
        new = tuple(a + dt * (b + 2*c + 2*d + e) / 6
                    for a, b, c, d, e in zip(state, k1, k2, k3, k4))
        endderivative = preview.derivative(new, preview.NU, False)
        rr, rmax, _, gg = bloch_residual(
            state[0], new[0], state[1], new[1], k1[0], endderivative[0],
            k1[1], endderivative[1], dt, ss, weights,
            preview.EPS, preview.NU, preview.MS, preview.OMEGA)
        # sqrt2 and every positive accumulation are rounded upward.
        L = float(up(up(up(math.sqrt(2)) * F) * up(1 + rmax)))
        ex, ph1, ph2 = positive_factors(float(up(L * dt)))
        dt_phi1 = float(up(dt * ph1))
        dt2_phi2 = float(up(up(dt * dt) * ph2))
        integral_error = float(up(integral_error + up(
            up(error * dt_phi1) + up(rr * dt2_phi2))))
        error = float(up(up(error * ex) + up(rr * up(dt * ph1))))
        pr = probe_residual(state[3], new[3], k1[3], endderivative[3],
                            gg, dt, kin, preview.MP, preview.BETA, preview.EPS, preview.NU)
        probe_defect_integral = float(up(probe_defect_integral + up(dt * pr)))
        interval_int = pt(0.)
        for j in range(4):
            interval_int = add(interval_int, divpos(take(gg, (j, 2)), pt(j + 1)))
        gz_integral = add(gz_integral, scale(interval_int, dt))
        max_residual, max_rnorm, max_L = max(max_residual, rr), max(max_rnorm, rmax), max(max_L, L)
        state = new
    pmax = float(up(2 * np.pi / preview.LENGTH))
    massmax = float(input_interval(float(preview.MP.max()))[1])
    coefficient = float(add(mul(eps_iv, pt(massmax)), mul(nu_iv, pt(pmax)))[1])
    eta_bg = float(up(coefficient * integral_error))
    tau_interval = add(pt(preview.TIME), mul(gz_integral, eps_iv))
    tau_stored = float(state[2][0])
    tau_interpolation_difference = max(abs(tau_stored - float(tau_interval[0])),
                                       abs(tau_stored - float(tau_interval[1])))
    tau_error = float(add(mul(eps_iv, pt(integral_error)),
                          pt(tau_interpolation_difference))[1])
    xi_interval = add(pt(preview.TIME), mul(gz_integral, nu_iv))
    xi_stored = float(state[2][1])
    xi_interpolation_difference = max(abs(xi_stored - float(xi_interval[0])),
                                      abs(xi_stored - float(xi_interval[1])))
    xi_error = float(add(mul(nu_iv, pt(integral_error)),
                         pt(xi_interpolation_difference))[1])
    flat_h_bound = float(up(massmax + pmax))
    eta_flat_clock = float(up(flat_h_bound * tau_error))
    return {
        "status": "OUTWARD_RESIDUAL_ENCLOSURE_UNDER_DECLARED_IEEE_BASIC_OPERATION_MODEL",
        "steps": steps, "F_upper": F, "max_L": max_L,
        "max_Hermite_r_norm": max_rnorm, "max_Bloch_residual": max_residual,
        "Bloch_terminal_error_bound": error,
        "integrated_Bloch_error_bound": integral_error,
        "probe_background_coefficient_error": eta_bg,
        "probe_own_interpolation_residual_error": probe_defect_integral,
        "probe_total_numerical_error": float(up(eta_bg + probe_defect_integral)),
        "clock_tau_error": tau_error,
        "clock_xi_error": xi_error,
        "massless_active_plus_flat_numerical_operator_error":
            float(up(pmax * up(xi_error + tau_error))),
        "flat_clock_calibration_operator_error": eta_flat_clock,
        "active_plus_flat_numerical_operator_error":
            float(up(eta_bg + probe_defect_integral + eta_flat_clock)),
        "stored_tau": tau_stored,
        "stored_xi": float(state[2][1]),
        "not_claimed": "formal IEEE implementation verification or certified hardware/library behavior",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=2048)
    args = parser.parse_args()
    print(json.dumps(run(args.steps), ensure_ascii=False, indent=2))
