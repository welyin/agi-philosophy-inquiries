"""Opposite-transverse WZ, finite energy, beam folding and unmatched triple-vertex directions."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(number):
    return json.loads((BASE.parent/str(number)/"inputs.json").read_text(encoding="utf-8-sig"))


def two_body(energy, mw, mz):
    s = energy**2
    rw, rz = mw*mw/s, mz*mz/s
    alpha = 1-rw-rz
    beta = math.sqrt(alpha**2-4*rw*rz)
    return alpha, beta, beta/(1+rw-rz), beta/(1+rz-rw)


def x_from_rapidity(gap, vw, vz):
    if gap >= math.atanh(vw)+math.atanh(vz):
        return 1.
    lo, hi = 0., 1.
    for _ in range(60):
        x = (lo+hi)/2
        if math.atanh(vw*x)+math.atanh(vz*x) < gap:
            lo = x
        else:
            hi = x
    return (lo+hi)/2


def finite_point(energy, mw, mz, gu, gd):
    alpha, beta, vw, vz = two_body(energy, mw, mz)
    xzero = alpha*(gu+gd)/(beta*(gu-gd))
    def even_weight(x):
        return (1-x**4)/(alpha**2-beta**2*x*x)**2
    def raw(x):
        # Common normalisation is immaterial; compute separate helicities from t and u.
        factor = gu/(alpha+beta*x)+gd/(alpha-beta*x)
        return sum((1-x*x)*(lam-x)**2*factor**2 for lam in [-1,1])
    for x in np.linspace(-.95,.95,31):
        normalisation = 2*(beta*(gu-gd))**2
        assert np.isclose(raw(x), normalisation*even_weight(x)*(x-xzero)**2)
        for f in (.5,.76,1.):
            lhs = (f*raw(x)+(1-f)*raw(-x))/(normalisation*even_weight(x))
            rhs = (x-(2*f-1)*xzero)**2+4*f*(1-f)*xzero**2
            assert np.isclose(lhs,rhs)
    nodes, weights = np.polynomial.legendre.leggauss(48)
    def integral(lo,hi):
        x = lo+(hi-lo)*(nodes+1)/2
        return (hi-lo)/2*np.dot(weights, even_weight(x)*(x*x+xzero**2))
    a,b = [x_from_rapidity(k,vw,vz) for k in (.5,1.5)]
    return {"cos_scattering_zero":xzero, "unsigned_gap_at_tagged_zero":
            abs(math.atanh(vw*xzero)+math.atanh(vz*xzero)),
            "unsigned_gap_bin_edge_cos": [a,b],
            "Born_opposite_TT_depth_two_fixed_bins":1-2*integral(0,a)/integral(a,b),
            "reduced_mistag_floor_at_correct_fraction_0_76":4*.76*.24*xzero**2}


def audit():
    pole = read(1469)["Z_flavor_specific"]
    seff = (1-pole["gV"][0]/pole["gA"][0])/4
    gu, gd = .5-2*seff/3,-.5+seff/3
    mw = read(1404)["const_GeV"]["MW"]
    mz = read(1475)["Z_mass_GeV"][0]
    rz = (gu+gd)/(gu-gd)
    onset = math.sqrt(mw**2+mz**2+2*mw*mz/math.sqrt(1-rz**2))
    assert math.isclose(abs(two_body(onset,mw,mz)[0]/two_body(onset,mw,mz)[1]*rz),1)
    c3w = {"with_int_only_95_TeV_neg2":[-.38,.17],"with_square_also_95_TeV_neg2":[-.062,.052]}
    return {"central_sine_squared_from_electron_pole":seff, "mass_W_Z_GeV":[mw,mz],
            "WZ_threshold_GeV":mw+mz,"opposite_TT_physical_internal_zero_begins_above_GeV":onset,
            "u_dbar_W_angle_high_energy_zero_WZ":rz,
            "u_dbar_W_angle_zero_Wgamma":1/3,
            "fixed_partonic_energies_GeV":{str(e):finite_point(e,mw,mz,gu,gd) for e in (200.,350.,500.,1000.)},
            "CMS_138fb_Wgamma_C3W_pT_and_decay_angle":c3w,
            "conditional_95_interval_width_ratio":np.ptp(c3w["with_int_only_95_TeV_neg2"])/
              np.ptp(c3w["with_square_also_95_TeV_neg2"])}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
