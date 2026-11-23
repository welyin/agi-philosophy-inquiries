"""Historical W calibration and line-shape sensitivity; masses and widths in GeV."""
import argparse
import cmath
from itertools import combinations, product
import json
from math import hypot, isclose, log, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = {"CMS": (80.3602, .0099), "ATLAS": (80.3665, .0159), "CDF": (80.4335, .0094)}
MID = (125.11, 170.7, 91.1876, .11852, .0276)
ERR = (.11, 1.5, .0021, .00084, .0001)
ALPHA_LEPT = .0314977


def fit(h, t, z, strong, had):
    x = log(h/100.)
    y = (h/100.)**2
    u = (t/174.3)**2-1
    return (80.3779 + x*(-.05427 + x*(-.008931 + x*x*.0000882)) + .000161*(y-1)
            + u*(.5237 - .0679*u - .00179*x + .0000664*y)
            - 1.070*(had+ALPHA_LEPT-.05907)/.05907
            - .0795*(strong-.119)/.119 + 114.9*(z-91.1875)/91.1875)


def pole_names(m, width):
    s = m*m/complex(1, width/m)
    fixed = sqrt(s.real)
    energy = cmath.sqrt(s)
    assert isclose(-s.imag/fixed, width/hypot(1, width/m), rel_tol=1e-13)
    return fixed, energy.real, -2*energy.imag


def mass_intersection(data):
    scale = max(abs(a-b)/(sa+sb) for (a,sa), (b,sb) in combinations(data.values(), 2))
    low = max(a-scale*sa for a,sa in data.values())
    high = min(a+scale*sa for a,sa in data.values())
    assert abs(low-high) < 1e-12
    return {"factor": scale, "touch_mass_GeV": low}


def model_box(k):
    values = [fit(*(v+k*i*s for v,i,s in zip(MID, signs, ERR))) for signs in product((-1,1), repeat=5)]
    return [min(values), max(values)]


def run():
    bench = [100., 174.3, 91.1875, .119, .05907-ALPHA_LEPT]
    baseline = fit(*bench)
    # Compare published table 2 changes, which are calculated relative to the old baseline.
    shifts = []
    for i, delta in enumerate((100, 5.1, .0021, .0027, .00036)):
        a = bench.copy()
        a[i] += delta
        shifts.append(1000*(fit(*a)-baseline))
    tabulated = (-41.4,31.,2.6,-1.7,-6.5)
    assert all(abs(i-j)<.15 for i,j in zip(shifts,tabulated))
    assert abs(baseline-80.3779)<1e-10
    initial = fit(*MID)
    assert abs(initial-80.3438192137)<1e-9
    central_delta = fit(MID[0], MID[1]+1, *MID[2:]) - initial
    compensate = central_delta*.05907/1.070
    assert abs(fit(MID[0], MID[1]+1, *MID[2:4], MID[-1]+compensate) - initial) < 1e-12
    assert compensate > 2*ERR[-1]

    def line(m, g):
        fixed, energy, decay_width = pole_names(m,g)
        return {"run_m_GeV": m, "run_width_GeV": g, "sqrt_Re_s_GeV": fixed,
                "Re_sqrt_s_GeV": energy, "width_energy_GeV": decay_width,
                "run_minus_constant_width_MeV": (m-fixed)*1000,
                "run_minus_energy_mass_MeV": (m-energy)*1000}

    # ATLAS free (m,width) measurement. Transform the reported correlated inputs locally.
    m, g, sm, sg, rho = 80.3548, 2.198, .0161, .049, -.30
    x = (g/m)**2
    dm = (1+2*x)/(1+x)**1.5
    dg = -(g/m)/(1+x)**1.5
    pole_sigma = sqrt((dm*sm)**2+(dg*sg)**2+2*rho*dm*dg*sm*sg)
    joint = line(m, g)
    joint.update({"constant_width_local_error_MeV": 1000*pole_sigma,
                  "measured_running_mass_change_relative_to_nominal_MeV": 1000*(m-DATA["ATLAS"][0]),
                  "linear_regression_dM_over_dWidth": rho*sm/sg})
    free = DATA.copy()
    free["ATLAS"] = (m,sm)
    boxes = {}
    for k in (1,2):
        lo, hi = model_box(k)
        boxes[str(k)] = {"polynomial_GeV": [lo, hi], "window_GeV": [lo-.0045, hi+.0045]}
    return {"v3_origin_GeV": baseline,
            "v1_table2_narrative_minus_v3_origin_MeV": 1000*(80.3799-baseline),
            "table2_relative_shifts_MeV_H_t_Z_alphaS_DeltaAlpha": shifts,
            "calibrated_central_GeV": initial, "registered_rectangles": boxes,
            "increase_had_for_one_GeV_heavier_top_fixed_W": compensate,
            "W_CMS_running_convention": line(DATA["CMS"][0], 2.085),
            "Z_LEP_running_convention": line(91.1876, 2.4952),
            "ATLAS_free_width_from_same_events": joint,
            "sharp_common_width_factor": mass_intersection(DATA),
            "sharp_factor_replacing_ATLAS_with_free_width": mass_intersection(free)}


def verify(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            verify(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x,y in zip(a,b):
            verify(x,y)
    else:
        assert isclose(a,b,abs_tol=2e-10,rel_tol=2e-10), (a,b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    target = HERE/"independent_results.json"
    result = run()
    if args.write:
        with target.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(result, json.loads(target.read_text(encoding="utf-8")))
        print("1394 W calibration and pole conventions independently verified")
