"""Independent once-subtracted spectral integral for a degenerate scalar threshold."""
import argparse
import json
from decimal import Decimal, localcontext
from math import hypot, pi, sqrt
from pathlib import Path

import numpy as np


def calculate():
    x, weights = np.polynomial.legendre.leggauss(64)
    speed, weights = (1+x)/2, weights/2
    def spectral_shape(z):
        if z == 4:
            return 2/9
        return z*np.sum(weights*speed**4 / (1-z*(1-speed**2)/4))/6

    def ratio_of_polarizations(mass):
        return spectral_shape((1500/mass)**2) / spectral_shape((500/mass)**2)

    target = ratio_of_polarizations(1000)
    def fit(r):
        low, high = 750.0, 2500.
        assert ratio_of_polarizations(low) > r > ratio_of_polarizations(high)
        for _ in range(64):
            test = (low+high)/2
            if ratio_of_polarizations(test) > r:
                low = test
            else:
                high = test
        return (low+high)/2
    interval = [fit(target*1.01/.99), fit(target*.99/1.01)]
    q_septet = np.arange(-3,4) + 2
    q_doublets = np.tile([0,1],56)
    assert sum(q_septet**2) == sum(q_doublets**2) == 56
    assert sum(q_septet**4) == 980 and sum(q_doublets**4) == 56
    factor = 56 / (4*pi*137.035999178)
    opal_sigma = hypot(58, 43, 30)*1e-5
    def spacelike(m):
        return factor*(spectral_shape(-1.81/m**2)-spectral_shape(-6.07/m**2))
    with localcontext() as ctx:
        ctx.prec = 72
        def series_shape(z):
            return z/30+z*z/420+z**3/3780
        mass = Decimal("2e13")
        ah, al = Decimal(1500)**2/mass**2, Decimal(500)**2/mass**2
        relative_excess = series_shape(ah)/series_shape(al)/9-1
        np.testing.assert_allclose(float(relative_excess),float((ah-al)/14),rtol=1e-15,atol=0.)
    return {"threshold_ratio_limit": spectral_shape(4)/spectral_shape(4/9),
            "timelike_ratio_1TeV":target, "one_percent_mass_bounds_GeV": interval,
            "two_traces_Q2": [int(sum(q_septet**2)),int(sum(q_doublets**2))],
            "two_traces_Q4": [int(sum(q_septet**4)),int(sum(q_doublets**4))],
            "opal_endpoint_one_TeV": spacelike(1e3),
            "opal_sigma_fractions_one_TeV_and_2e13": [spacelike(m)/opal_sigma for m in (1e3,2e13)],
            "relative_curvature_of_ratio_at_2e13": float(relative_excess)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    vals = calculate()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8") as out:
            json.dump(vals, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert saved.keys() == vals.keys()
        for key, val in vals.items():
            np.testing.assert_allclose(val, saved[key], rtol=2e-12, atol=0)
        print("1378: spectral density, inverse at finite precision, OPAL scale and trace degeneracy verified.")
