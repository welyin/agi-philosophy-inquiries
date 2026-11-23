"""Spinor-level Bhabha cross-check and single-box norm bound for a pure axial X."""

import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
EW = json.loads((HERE / CFG["predecessor"]).read_text(encoding="utf-8"), parse_float=F)["couplings"]
RESULT = HERE / "independent_results.json"
QED = EW["g2"]**2 * EW["gY"]**2 / (EW["g2"]**2 + EW["gY"]**2)
GX2 = EW["g2"]**2 + EW["gY"]**2
SW = EW["gY"]**2 / GX2
ZSHAPE = (SW-F(1, 2), SW)
ZV = np.sqrt(float(GX2)) * np.array(list(map(float, ZSHAPE)))


def sqrt_bound(x, upper):
    units = 10**12
    scaled = x * units**2
    root = math.isqrt(scaled.numerator // scaled.denominator)
    return F(root + int(upper and root**2 < scaled), units)


def norm_certificate():
    k0, k1 = map(F, CFG["certified_k_envelope"])
    emin, emax = CFG["energy_GeV"]
    mmin, mmax = CFG["mass_GeV"]
    gmin, gmax = map(F, CFG["gF"])
    qlo, qhi = gmin**2/2, gmax**2/2
    r0, r1 = F(mmin, emax)**2, F(mmax, emin)**2
    rz0, rz1 = (EW["mZ"]/emax)**2, (EW["mZ"]/emin)**2
    im_x = r1 * F(CFG["Gamma_over_mass"][1])
    p0 = qlo*(1-r1)/((1-r0)**2 + im_x**2)
    p1, t1 = qhi/(1-r1), qhi/r0
    z1 = 1/(1-rz1)
    ll, rr = [GX2 * a**2 for a in ZSHAPE]
    lr = abs(GX2 * ZSHAPE[0] * ZSHAPE[1])
    same = [QED + val*(z1+1/rz0) for val in (ll, rr)]
    bn = sqrt_bound(sum(b**2 for b in same)
                    + 2*(k1*(QED + lr*z1))**2 + 2*(lr/rz0)**2, True)
    dn = sqrt_bound(2*(p1+t1)**2 + 2*(k1*p1)**2 + 2*t1**2, True)
    slow = QED*sqrt_bound(2*(1-k1)**2+2, False) - k1*bn
    shigh = 2*QED + k1*bn
    lo = -4*QED*(p1+2*k1*t1) - 2*k1*bn*dn
    up = -4*QED*(1-k1)**2*p0 + k1*(2*bn*dn + dn**2)
    assert slow > 0 and lo < up < 0
    return [math.floor(10**6*k1*lo/slow**2), math.ceil(10**6*k0*up/shigh**2)]


# All four momenta lie in the xz plane.  Each column below has ubar or vbar
# completeness p-slash, with E=1/2 so s=1.  The gamma matrices use the Weyl basis.
I2 = np.eye(2)
ZERO = np.zeros((2, 2))
PAULI = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
G0 = np.block([[ZERO, I2], [I2, ZERO]])
GAMMA = [G0] + [np.block([[ZERO, sig], [-sig, ZERO]]) for sig in PAULI]
METRIC = np.array([1, -1, -1, -1])


def spinors(theta):
    c, s = math.cos(theta/2), math.sin(theta/2)
    return np.array([[-s, 0], [c, 0], [0, c], [0, s]])


def scattering(k, energy, mass, g, width, photon_only=False):
    angle = 2*math.asin(math.sqrt(k))
    states = [spinors(t) for t in (0, math.pi, angle, angle-math.pi)]
    a = math.sqrt(float(QED))
    bosons = [(a, a, 0., 0.)]
    if not photon_only:
        bosons.append((*ZV, float(EW["mZ"]/energy)**2,
                        float(EW["mZ"]*EW["GammaZ"]/energy**2)))
    if g:
        h = g/math.sqrt(2)
        r = (mass/energy)**2
        bosons.append((h, -h, r, r*width))
    amp = np.zeros((2, 2, 2, 2), dtype=complex)
    for left, right, r, iw in bosons:
        vertex = np.diag([left, left, right, right])
        ops = [G0 @ mu @ vertex for mu in GAMMA]

        def current(a_index, b_index):
            return np.array([states[a_index].T @ op @ states[b_index] for op in ops])

        annihilation = np.einsum("uba,ucd,u->abcd", current(1, 0), current(2, 3), METRIC)
        scattering_t = np.einsum("uca,ubd,u->abcd", current(2, 0), current(1, 3), METRIC)
        # Relative exchange sign for Wick contractions of identical electron fields.
        amp += annihilation / (1-r + 1j*iw) - scattering_t / (-k-r)
    return np.vdot(amp, amp).real / 4


def integrated_ratio(point, n):
    lo, hi = [math.sin(t/2000)**2 for t in CFG["narrow_arm_mrad"]]
    nodes, weights = np.polynomial.legendre.leggauss(n)
    ks = (lo + hi + nodes*(hi-lo))/2
    energy, mass, g, width = point
    full = [scattering(k, *point) for k in ks]
    sm = [scattering(k, energy, mass, 0, width) for k in ks]
    return np.dot(weights, full) / np.dot(weights, sm)


def run():
    interval = norm_certificate()
    assert interval[0] >= -2600 and interval[1] < 0
    samples = ((189, 56, .62, 0), (200, 60, .632, .05), (207, 64, .65, .1))
    values = []
    for e, m, g, w in samples:
        ratio = integrated_ratio((e,m,g,w), 24)
        assert abs(ratio-integrated_ratio((e,m,g,w), 12)) < 2e-13
        assert interval[0] <= (ratio-1)*1e6 <= interval[1]
        r = (m/e)**2
        slope = -(g*g/2) * (1-r) / ((1-r)**2+(r*w)**2) / float(QED)
        k = 1e-7
        numerical = (scattering(k,e,m,g,w)/scattering(k,e,m,0,w) - 1)/k
        assert abs(slope-numerical) < 2e-4
        values.append({"sqrt_s_m_g_width": [e,m,g,w], "cross_section_ratio":round(ratio,12),
                       "forward_derivative_in_k":round(slope,9)})
    for k in (.00037, .00060, .01, .2, .7):
        u = -(1-k)
        textbook = 2*float(QED**2)*((u*u+k*k) + (u*u+1)/k**2 - 2*u*u/k)
        assert abs(scattering(k, 200, 60, 0, 0, photon_only=True)/textbook - 1) < 1e-13
    return {"uniform_single_box_norm_bound_ppm":interval,"two_body_spinor_witnesses":values,
            "qed_crossing_and_spin_average_checked":True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(RESULT.read_text(encoding="utf-8"))
    else:
        with RESULT.open("x", encoding="utf-8") as out:
            out.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
