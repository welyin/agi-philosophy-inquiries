"""Sharp triplet matching diagnostics at a gauge-fixed mass product."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def weak_squared(mass, mu, gw, first_mass=None):
    active = 0. if first_mass is None else max(0.,math.log(mass/first_mass))
    return 1/(1/gw**2+((19/6)*math.log(mass/mu)-(4/3)*active)/(8*math.pi**2))


def scalar_mass_ledger(masses, mu, gw):
    a,b = sorted(masses)
    return 6*(a*a*weak_squared(a,mu,gw)**2+b*b*weak_squared(b,mu,gw,a)**2)/LOOP**2


def degeneracy_value(mean, mu, gw):
    return 12*mean**2*weak_squared(mean,mu,gw)**2/LOOP**2


def spectral_trace(light_eV, heavy, angle, v):
    c,s = np.cos(angle),np.sin(angle)
    r = np.array([[c,s],[-s,c]],dtype=complex)
    y = 1j*math.sqrt(2)/v * np.diag(np.sqrt(np.array(light_eV)*1e-9)) @ r @ np.diag(np.sqrt(heavy))
    mass = -v*v/2*y@np.diag(1/np.array(heavy))@y.T
    assert np.allclose(mass*1e9,np.diag(light_eV),rtol=1e-9,atol=1e-13)
    return float(np.linalg.norm(y@np.diag(heavy))**2)


def evaluate():
    cfg = load(HERE/"inputs.json")
    gauge = load(HERE/cfg["gauge_product"])
    light = load(HERE/cfg["neutrino_products"])["neutrino"]
    old = load(HERE/cfg["reference"])
    mu, gw, v = old["matching_scale"],old["reference"]["g2"],cfg["vev_GeV"]
    mh = old["reference"]["mh"]
    mh_max = mh + 2*old["marginal_width"]["mh"]
    lower,upper = gauge["whole_triplet_geometric_mean_GeV"]
    mean = gauge["central_triplet_geometric_mean_GeV"]
    snap = {}
    for order,m in light.items():
        a,b = m["central_masses_eV"]
        ratio = (b/a)**(1/3)
        poles = [mean/math.sqrt(ratio),mean*math.sqrt(ratio)]
        trace = spectral_trace((a,b),poles,math.pi/2,v)
        sharp = 4*mean**3*math.sqrt(a*b)*1e-9/v**2
        assert math.isclose(trace,sharp,rel_tol=1e-12)
        for magn in (1.,ratio,10.,300.):
            trial = [mean/math.sqrt(magn),mean*math.sqrt(magn)]
            assert scalar_mass_ledger(trial,mu,gw)>=degeneracy_value(mean,mu,gw)*(1-1e-12)
            for z in (0.,math.pi/2,.3+.5j,1.2+1.3j):
                assert spectral_trace((a,b),trial,z,v)>=sharp*(1-1e-12)
        ys = m["full_product_range"]
        snap[order] = {"optimum_heavy_mass_ratio":ratio,"optimum_poles_GeV":poles,
                       "central_trace_floor_GeV2":trace,
                       "central_own_threshold_abs_muH2_floor_GeV2":6*trace/LOOP,
                       "full_box_abs_muH2_floor_envelope_GeV2":
                           [12*p*p*s/LOOP for p,s in [(lower,ys[0]),(upper,ys[1])]]}
    ledger_min = degeneracy_value(lower,mu,gw)
    left,right = mu,lower
    for _ in range(70):
        mid = (left+right)/2
        if degeneracy_value(mid,mu,gw)>mh_max**2:
            right = mid
        else:
            left = mid
    # With all colored inputs fixed, alpha1=alpha3 still fixes the meeting
    # scale. Changing the weak mass product changes its inverse coupling alone.
    threshold_shift = math.log(lower/right)/(3*math.pi**2)
    assert ledger_min/mh_max**2 > 9e4 and threshold_shift>0.15
    return {"triplet_yukawa":snap,"two_loop_gauge_abs_m2_central_GeV2":degeneracy_value(mean,mu,gw),
            "full_box_min_gauge_abs_m2_GeV2":ledger_min,
            "lowest_ratio_to_Mh2_budget":ledger_min/mh_max**2,
            "single_largest_threshold_absolute_m2_lower_GeV2":ledger_min/2,
            "maximum_geometric_mean_for_one_total_Mh2_GeV":right,
            "minimum_weak_matching_shift_inverse_g2":threshold_shift,
            "same_shift_inverse_alpha2":4*math.pi*threshold_shift,
            "mu_over_M_where_pure_gauge_finite_term_vanishes":math.exp(-1/12)}


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    opts = cli.parse_args()
    output = HERE/"results.json"
    data = evaluate()
    if opts.check:
        assert data==load(output)
    else:
        with output.open("x",encoding="utf-8") as file:
            file.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    print("1612: both sharp trace bounds, MS matching and compensating gauge threshold verified")
