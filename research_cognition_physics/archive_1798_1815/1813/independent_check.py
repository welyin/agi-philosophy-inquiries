"""Signed two-loop threshold projection in the non-supersymmetric five-and-adjoint model."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = runpy.run_path(str(HERE.parent/"1808"/"check.py"))
GAUGE = json.loads((HERE.parent/"1807"/"independent_results.json").read_bytes())
FITS = json.loads((HERE.parent/"1811"/"results.json").read_bytes())
SM_B = np.array([41/10, -19/6, -7])
SM_B2 = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12], [11/10, 9/2, -26]])
CA = np.array([0., 2., 3.])
W = np.array([5., -3., -2.])
LOOP = 16*math.pi**2
MPL = 2.435e18
CAP = 1.0  # Each indicated (squared) Yukawa Frobenius norm separately.


def coefficients(indices, casimirs, kind):
    t, c = np.array(indices), np.array(casimirs)
    if kind == "Dirac":
        return 4*t/3, 4*np.outer(t, c)+20*np.diag(t*CA)/3
    factor = {"complex": 1., "real": .5}[kind]
    return factor*t/3, factor*(4*np.outer(t, c)+2*np.diag(t*CA)/3)


FIELDS = {
    "D": coefficients([3/10, 1/2, 0], [3/20, 3/4, 0], "Dirac"),
    "T": coefficients([1/5, 0, 1/2], [1/15, 0, 4/3], "Dirac"),
    "S3": coefficients([0, 2, 0], [0, 2, 0], "real"),
    "S8": coefficients([0, 0, 3], [0, 0, 3], "real"),
    "HC": coefficients([1/5, 0, 1/2], [1/15, 0, 4/3], "complex"),
}


def log_integrals(a, b, dt):
    """With a=1/alpha, integrate g^2 dlog(mu) along each one-loop segment."""
    end = a - b*dt/(2*math.pi)
    assert min(end) > 0
    integrals = np.empty(3)
    for i, slope in enumerate(b):
        integrals[i] = (4*math.pi*dt/a[i] if abs(slope) < 1e-12 else
                        -8*math.pi**2*math.log1p(-slope*dt/(2*math.pi*a[i]))/slope)
    return end, integrals


def integral(data, events, endpoint):
    a = data["a"].copy()
    b, matrix = SM_B.copy(), SM_B2.copy()
    t0 = corr = 0.
    for t, field in sorted(events+[(endpoint, None)], key=lambda pair:pair[0]):
        if t > endpoint:
            break
        a, ig = log_integrals(a, b, t-t0)
        corr += W@matrix@ig/LOOP
        if field:
            db, db2 = FIELDS[field]
            b += db
            matrix += db2
        t0 = t
    return corr


def central(tag):
    data = BASE["inputs"](tag)
    old = GAUGE[tag]["two_split_complete_pairs"]
    light = data["MZ"]
    q = old["matching_scale_GeV"]
    xratio = FITS[tag]["MX_Q_lower_if_A2_1_same_g"]
    endpoint = math.log(q*xratio/light)
    events = [(math.log(m/light), name) for name in ["D", "T"] for m in old["M_"+name+"_GeV"]]
    rows = {}
    for f3, f8, fC in itertools.product([0,.5,1.], repeat=3):
        extra = [(f*endpoint, n) for f, n in zip((f3, f8, fC), ("S3", "S8", "HC"))]
        rows[tuple([f3, f8, fC])] = integral(data, events+extra, endpoint)
    return dict(physical_X_test_GeV=q*xratio, one_loop_reference_Q_GeV=q,
                exact_on_one_loop_path_gauge_first_correction_at_all_light=rows[(0.,0.,0.)],
                exact_on_one_loop_path_gauge_first_correction_at_all_X=rows[(1.,1.,1.)],
                sample_27_min_max=[min(rows.values()),max(rows.values())],
                smallest_first_order_deficit_before_gauge_Yukawa_and_finite=44*math.log(xratio)-2*endpoint)


def upper_integral(g0, slope, t):
    # d(1/g^2)/dt >= -slope/(8 pi^2) supplies a comparison solution.
    return -8*math.pi**2/slope*np.log1p(-slope*g0*t/(8*math.pi**2))


def uniform(tag):
    data = BASE["inputs"](tag)
    # All b additions and all delta B entries are positive; SM B33 remains signed.
    contents = ["D", "D", "T", "T", "S3", "S8", "HC"]
    bmax = SM_B+sum(FIELDS[n][0] for n in contents)
    matrix_max = SM_B2+sum(FIELDS[n][1] for n in contents)
    t_planck = math.log(MPL/data["MZ"])
    g0 = 4*math.pi/data["al"]
    box = np.array([.7, .44, 1.6])
    assert (g0 < box).all()
    beff = bmax + np.maximum(matrix_max, 0.)@box/LOOP
    assert max(beff[1:])<0
    end = g0/(1-g0*beff*t_planck/(8*math.pi**2))
    assert end[0]<box[0]
    assert np.isclose(W@SM_B,44)
    # Q = MZ exp{2pi w.a0/44} for the two full 5s; individual thresholds are w-blind.
    time_q_max = 2*math.pi*(W.clip(0)@data["ah"]+W.clip(max=0)@data["al"])/44
    ratio = min(r["MX_Q_lower_if_A2_1_same_g"] for r in FITS.values())
    t = time_q_max+math.log(ratio)

    def harmful_gauge(time):
        ig = upper_integral(g0, beff, time)
        return (81*ig[1]+24*ig[2])/LOOP

    def harmful_yukawa(time, cap):
        # 6*l + 12*ue + 6*ud; the top, new adjoint Yukawas and positive terms drop.
        return 24*cap*time/LOOP

    # Explicitly allow either sign of the MS vector matching constant, |c|=2.
    budget = 2.
    def remaining(logratio, cap):
        time = time_q_max+logratio
        return 44*logratio-2*time-harmful_gauge(time)-harmful_yukawa(time, cap)-budget
    def ratio_limit(cap):
        lo, hi = 0., math.log(ratio)
        assert remaining(hi, cap)>0
        for _ in range(80):
            mid = (lo+hi)/2
            if remaining(mid, cap)>0:
                hi=mid
            else:
                lo=mid
        return math.exp(hi)
    # At fixed X/Q, using maximal Q maximizes all harmful terms.
    need_before_yuk = 44*math.log(ratio)-2*t-harmful_gauge(t)-budget
    needed_average = need_before_yuk*LOOP/t
    max_colored_yuk_square_sum = 3*4*math.pi
    # Rank <=3 for both u^c e^c and u^c d^c. Color couplings turn on only above HC.
    max_color_mass = {}
    for key, ecap in [("trace_e_at_most_1",1.),("all_singular_values_at_most_sqrt4pi",max_colored_yuk_square_sum)]:
        minimum_log = (LOOP*need_before_yuk-6*ecap*t)/(18*max_colored_yuk_square_sum)
        max_color_mass[key] = data["MZ"]*math.exp(t-minimum_log)
    # Normalized trace insertion: delta L=-1/2 Tr[Z(Sigma) F^2],
    # with additional delta Z1 for hypercharge; delta alpha1^-1=4pi delta Z1.
    return dict(Q_high_two_width_GeV=data["MZ"]*math.exp(time_q_max),
                maximum_active_b_i=bmax.tolist(), two_loop_resummed_b_eff_upper=beff.tolist(),
                g_squared_tube=box.tolist(), upper_g_squared_at_reduced_Planck=end.tolist(),
                tested_X_over_Q=ratio, max_log_MX_over_MZ_at_ratio=t,
                harmful_two_loop_gauge_upper=harmful_gauge(t),
                full_range_gauge_upper_through_Planck=harmful_gauge(t_planck),
                favorable_finite_constant_abs=budget,
                required_negative_Yukawa_weighted_mean_after_gauge=needed_average,
                necessary_max_Frobenius_squared_of_three_sectors=needed_average/24,
                max_X_over_Q_with_all_three_Frobenius_squared_at_most_1=ratio_limit(CAP),
                max_X_over_Q_with_all_three_Frobenius_squared_at_most_4pi=ratio_limit(4*math.pi),
                max_HC_mass_GeV_if_colored_singular_values_at_most_sqrt4pi=max_color_mass,
                positive_deficit_units_at_proton_test_if_caps_1=remaining(math.log(ratio), CAP),
                min_hypercharge_kinetic_shift_if_caps_1=remaining(math.log(ratio), CAP)/(40*math.pi**2))


def calculate():
    fields = {"SM": (SM_B, SM_B2), **FIELDS}
    coeffs = {key:[float(W@db), (W@mat).tolist()] for key,(db,mat) in fields.items()}
    for key, expected in {"SM":[15,-13,60],"D":[0,-20,0],"T":[0,0,-20],
                          "S3":[0,-28,0],"S8":[0,0,-42],"HC":[0,0,-2]}.items():
        np.testing.assert_allclose(coeffs[key][1],expected,atol=3e-14)
    result = dict(projected_one_and_two_loop=coeffs)
    for tag in FITS:
        result[tag]=dict(central_first_correction=central(tag),uniform_bound=uniform(tag))
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    path=HERE/"independent_results.json"
    observed=json.loads(json.dumps(calculate(),allow_nan=False),parse_float=lambda x:float(f"{float(x):.10g}"))
    if p.parse_args().check:
        assert observed==json.loads(path.read_bytes())
    else:
        path.write_bytes((json.dumps(observed,indent=2)+"\n").encode())
    print("1813 independent",hashlib.sha256(path.read_bytes()).hexdigest()[:12])
