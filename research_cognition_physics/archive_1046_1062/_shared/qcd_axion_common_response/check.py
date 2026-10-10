"""Count0 identities and an explicitly ideal rectangular Ramsey filter only."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys
import numpy as np


def sinc(x):
    return 1.0 if x == 0.0 else math.sin(x)/x


def calculate():
    # Arbitrary dimensionless rational substitutions, not physical axion points.
    chi, rho, kappa = F(9), F(1, 8), F(2, 5)
    sqrt_chi, sqrt_2rho = F(3), F(1, 2)
    assert sqrt_chi**2 == chi and sqrt_2rho**2 == 2*rho
    cases = []
    for fa in (F(3), F(6)):
        mass = sqrt_chi/fa
        amplitude = sqrt_2rho/mass
        theta_amplitude = amplitude/fa
        edm_amplitude = kappa*theta_amplitude
        assert mass**2*fa**2 == chi
        assert mass**2*amplitude**2/2 == rho
        assert theta_amplitude == F(1, 6)
        assert edm_amplitude == F(1, 15)
        cases.append({"fa": str(fa), "mass": str(mass),
                      "field_amplitude": str(amplitude),
                      "theta_amplitude": str(theta_amplitude),
                      "edm_amplitude": str(edm_amplitude)})
    # Rectangular time average: changing mass changes the filter, not bare D.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    time0, duration, phase, D = 0.4, 1.7, 0.23, F(1, 15)
    rows, max_error = [], 0.0
    for omega in (0.0, 0.5, 1.0, 2*math.pi/duration):
        times = time0 + duration*(nodes+1)/2
        quadrature = float(D)*float(np.dot(weights, np.cos(omega*times+phase)))/2
        analytic = float(D)*sinc(omega*duration/2)*math.cos(omega*(time0+duration/2)+phase)
        error = abs(analytic-quadrature)
        max_error = max(max_error, error)
        assert error < 1e-13
        rows.append({"omega": omega, "rectangular_average": analytic,
                     "quadrature_error": error})
    assert abs(rows[-1]["rectangular_average"]) < 1e-15
    assert abs(rows[1]["rectangular_average"]-rows[2]["rectangular_average"]) > 1e-3
    # H=-d E sigma_z has an energy splitting of 2 d E.
    electric, d = F(7, 3), F(2, 11)
    levels = (-d*electric, d*electric)
    assert levels[1]-levels[0] == 2*d*electric
    return {"scope": "count0 algebra and ideal filter; no QCD or empirical recomputation",
            "dimensionless_inputs": {"chi": str(chi), "rho": str(rho), "kappa": str(kappa)},
            "rational_cases": cases,
            "ramsey_ideal_window": {"time0": time0, "duration": duration, "phase": phase},
            "filter_checks": rows,
            "max_filter_quadrature_error": max_error,
            "spin_energy_splitting_exact": str(levels[1]-levels[0]),
            "unchanged_bare_amplitude_does_not_imply_unchanged_finite_readout": True,
            "all_assertions_passed": True}


if __name__ == "__main__":
    result = calculate()
    target = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments for read-only checks, or --save-exclusive once.")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print("PASS: common theta/mass/amplitude and ideal finite-window identities; count0.")
