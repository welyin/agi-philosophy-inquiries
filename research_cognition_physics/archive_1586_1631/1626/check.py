"""Rational enclosure of the unpolarized OPAL two-body luminosity kernel."""
import argparse
import importlib.util
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
OLD = HERE / INPUT["predecessor"]
EW = json.loads(OLD.read_text(encoding="utf-8"), parse_float=F)["couplings"]


class Box:
    """Exact rational interval arithmetic, including sign changes."""
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        assert self.lo <= self.hi

    @classmethod
    def wrap(cls, value):
        return value if isinstance(value, cls) else cls(value)

    def __add__(self, other):
        other = self.wrap(other)
        return Box(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.wrap(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = self.wrap(other)
        values = [a * b for a in (self.lo, self.hi)
                  for b in (other.lo, other.hi)]
        return Box(min(values), max(values))

    __rmul__ = __mul__

    def square(self):
        upper = max(self.lo ** 2, self.hi ** 2)
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo ** 2, self.hi ** 2)
        return Box(lower, upper)

    def inv(self):
        assert self.lo * self.hi > 0
        return Box(1 / self.hi, 1 / self.lo)

    def divide(self, n):
        step = (self.hi - self.lo) / n
        return [Box(self.lo + i * step, self.lo + (i + 1) * step) for i in range(n)]


def prop(real, imaginary):
    norm_inv = (real.square() + imaginary.square()).inv()
    return real * norm_inv, -imaginary * norm_inv


def transfer_delta(k, energy, mass, g, width):
    x, y = EW["gY"] ** 2, EW["g2"] ** 2
    e2, sw = x * y / (x + y), x / (x + y)
    ll, rr = (x + y) * (sw - F(1, 2)) ** 2, (x + y) * sw ** 2
    lr = (x + y) * (sw - F(1, 2)) * sw
    inv_s = energy.square().inv()
    rz = EW["mZ"] ** 2 * inv_s
    wz = EW["mZ"] * EW["GammaZ"] * inv_s
    zre, zim = prop(1 - rz, wz)
    r = mass.square() * inv_s
    xre, xim = prop(1 - r, width * r)
    h = g.square() * F(1, 2)
    sr, si = h * k * xre, h * k * xim
    t = h * k * (k + r).inv()
    w = (1 - k).square()
    ktz = k * (k + rz).inv()
    pt = -lr * ktz
    ars, ais = k * (e2 + lr * zre), k * lr * zim
    den = 2 * (pt - e2).square() + 2 * k.square() * (ars.square() + ais.square())
    # Combine photon-t interference across helicities BEFORE interval bounding.
    delta = -4 * e2 * (w * sr + k * (2 - k) * t)
    delta += 2 * w * ((sr - t).square() + si.square()) + 4 * pt * t + 2 * t.square()
    delta += 2 * k.square() * (sr.square() + si.square() - 2 * (ars * sr + ais * si))
    for chiral in (ll, rr):
        pr, pi = k * (e2 + chiral * zre) - chiral * ktz, k * chiral * zim
        den += w * ((pr - e2).square() + pi.square())
        delta += 2 * w * (pr * (sr - t) + pi * si)
    return delta * den.inv()


def check():
    energy = Box(*INPUT["energy_GeV"])
    m, g = Box(*INPUT["mass_GeV"]), Box(*INPUT["gF"])
    width, k = Box(*INPUT["Gamma_over_mass"]), Box(*INPUT["certified_k_envelope"])
    theta = [F(a, 2000) for a in INPUT["narrow_arm_mrad"]]
    # For 0 < x < .03, x - x^3/6 <= sin(x) <= x.
    assert k.lo <= (theta[0] - theta[0] ** 3 / 6) ** 2
    assert theta[1] ** 2 <= k.hi
    values = [transfer_delta(kc, ec, mc, g, width)
              for kc in k.divide(INPUT["k_subdivisions"])
              for ec in energy.divide(2) for mc in m.divide(2)]
    scale = 1_000_000
    bound = [math.floor(scale * min(d.lo for d in values)),
             math.ceil(scale * max(d.hi for d in values))]
    assert -3000 < bound[0] < bound[1] < 0
    return bound


def old_integral(energy, mass, g, width, nodes):
    spec = importlib.util.spec_from_file_location("bhabha_1623", OLD.with_name("check.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    low, high = [math.sin(a / 2000) ** 2 for a in INPUT["narrow_arm_mrad"]]
    x, weights = np.polynomial.legendre.leggauss(nodes)
    k = (low + high + x * (high - low)) / 2
    terms = model.components(energy, model.I(1 - 2 * k),
                             model.I((mass / energy) ** 2), model.I(width))
    bsm = model.density(terms, model.I(g ** 2 / 2)).lo
    sm = model.density(terms, model.I(0)).lo
    return float(weights @ bsm / (weights @ sm))


def run():
    bound = check()
    samples = []
    for energy, mass, g, width in ((189,56,.62,0), (200,60,.632,.05),
                                    (207,64,.65,.1)):
        ratio = old_integral(energy, mass, g, width, 32)
        assert abs(ratio - old_integral(energy, mass, g, width, 64)) < 2e-12
        assert bound[0] < (ratio - 1) * 1e6 < bound[1]
        samples.append({"energy_GeV":energy, "mass_GeV":mass, "gF":g,
                        "Gamma_over_mass":width, "monitor_over_SM":ratio})
    assert old_integral(200,60,0,.05,32) == 1
    return {"pointwise_fractional_change_ppm_outward":bound,
            "uniform_common_acceptance_examples":samples,
            "method":"exact rational enclosure, 64 parameter cells; full 1623 helicity quadrature"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        assert result == json.loads((HERE / "results.json").read_text(encoding="utf-8"))
    print(json.dumps(result, indent=2))
