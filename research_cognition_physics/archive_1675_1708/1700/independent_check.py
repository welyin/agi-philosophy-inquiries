"""Trace identity for Higgs orientation and squark/slepton mass endpoints."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def beta_identity(families):
    # H_u,H_d,tr Q,tr U,tr L,tr E. Q,U are third-family or three-family traces.
    # Columns: four soft Yukawa sums, gauge terms g_i²|M_i|², g_1² S.
    k = families
    matrix = np.array([
        [6, 0, 0, 2, -F(6,5), -6, 0, F(3,5)],
        [0, 6, 2, 0, -F(6,5), -6, 0, -F(3,5)],
        [2, 2, 0, 0, -k*F(2,15), -6*k, -k*F(32,3), k*F(1,5)],
        [4, 0, 0, 0, -k*F(32,15), 0, -k*F(32,3), -k*F(4,5)],
        [0, 0, 2, 2, -F(18,5), -18, 0, -F(9,5)],
        [0, 0, 4, 0, -F(72,5), 0, 0, F(18,5)],
    ], dtype=object)
    return list(np.array([-1,1,-3,3,1,-1]) @ matrix)


def compute():
    assert beta_identity(1) == [0,0,0,0,F(24,5),0,0,-F(48,5)]
    assert beta_identity(3) == [0,0,0,0,-F(36,5),36,0,-F(78,5)]
    old = read("../1699/independent_results.json")
    low = read("../1699/results.json")["soft_interface"]
    ms = old["common_scale_GeV"]
    mh = read("../1372/inputs.json")["atlas"]["mass_GeV"][0]
    mz = read("../1369/inputs.json")["mu_ref_GeV"]
    vev = read("../1492/inputs.json")["v246_GeV"]
    dt = math.log(read("../1699/inputs.json")["comparison_scale_GeV"]/ms)
    gap = ms**2+mh**2/2
    sw2 = low["gY"]**2/(low["gY"]**2+low["g2"]**2)
    d_diff = abs(0.5-4*sw2/3)*mz**2
    # At all studied charged thresholds b1 cannot exceed MSSM plus two 5+5bar pairs.
    g1sq_max = 1/(1/low["g1_GUT"]**2 - F(43,5)*dt/(8*math.pi**2))
    loop = 16*math.pi**2
    alpha = float(F(24,5) * g1sq_max * dt / loop)
    delta = old["results_by_other_quartic_allowance"]["0"]["mechanisms"]["actual_neutral_soft_CW"][
        "max_abs_Hu_mass_squared_GeV2"]
    mt = low["yt"]*vev/math.sqrt(2)
    soft_result = {}
    for quartic in [0, .003, .01]:
        ratio = math.sqrt((low["lambda"]-low["stop_max"]-quartic)/low["D"])
        needed = gap*ratio - delta
        qmin = needed/3
        d = max(0.,qmin-d_diff)/ms**2
        r = (d+math.sqrt(4+d*d))/2
        # Directly diagonalize a concrete stop submatrix saturating the soft endpoint
        # relation with the neutral remainder accepted at its full allowance.
        # It only tests this algebraic measurement map.
        cos2 = -ratio
        diagdiff = qmin+(0.5-4*sw2/3)*mz**2*cos2
        off = math.sqrt(6)*mt*ms
        center = math.sqrt(ms**4+(diagdiff/2)**2+off**2)
        evals = np.linalg.eigvalsh([[center+diagdiff/2, off], [off, center-diagdiff/2]])
        mlo,mhi = np.sqrt(evals)
        np.testing.assert_allclose(mlo*mhi,ms**2,rtol=1e-14)
        np.testing.assert_allclose(evals[1]-evals[0],math.hypot(diagdiff,2*off),rtol=1e-13)
        assert mhi/mlo>=r
        # Direct exclusion for an interval of actual stop eigenmasses using the
        # nondegenerate pure-y_t^4 threshold with arbitrary real or complex mixing.
        cmax = (mt**2+mz**2)/ms**2
        def lambda_ceiling(trial):
            tmax = .5*math.log((trial+cmax)/(1/trial-cmax))
            logbase = math.log(math.cosh(math.log(trial))+cmax)
            shape = 3/(1-tmax**2/15)
            cosine = (3*((trial-1/trial)*ms**2+d_diff)+delta)/gap
            return low["D"]*min(1.,cosine)**2 + low["stop_max"]*(logbase+shape)/3 + quartic
        lo, hi = 1., 1.3
        assert lambda_ceiling(lo)<low["lambda"]<lambda_ceiling(hi)
        for _ in range(48):
            mid = (lo+hi)/2
            if lambda_ceiling(mid)<low["lambda"]:
                lo = mid
            else:
                hi = mid
        safe_lo = math.floor(lo*1e6)/1e6
        assert lambda_ceiling(safe_lo)<low["lambda"]
        soft_result[str(quartic)] = dict(minimum_abs_cos2beta=ratio, absolute_spectrum_budget_GeV2=needed,
                minimum_stop_L_minus_R_squared_if_equal_slepton_traces_GeV2=qmin,
                minimum_common_fractional_spread_in_soft_squared_masses=needed/(12*ms**2),
                minimum_stop_mass_ratio_at_fixed_geometric_mean_all_Xt=r,
                ratio_excluded_up_to_with_actual_nondegenerate_stop_threshold=safe_lo,
                quartic_ceiling_at_excluded_ratio=lambda_ceiling(safe_lo),
                stop_masses_GeV_if_max_mixing_and_budget_saturated=[mlo,mhi],
                small_tan_and_equal_IR_requires_some_abs_bino_over_Ms_at_least=math.sqrt(needed/ms**2/alpha))
    return dict(running_top_mass_proxy_GeV=mt, tree_stop_D_difference_bound_GeV2=d_diff,
                g1_squared_upper_from_largest_beta=float(g1sq_max),
                max_unsigned_J_over_Ms2_per_squared_bino_bound=alpha,
                inherited_total_quadratic_remainder_sensitivity_GeV2=delta,
                coefficients_third_quark_family=[str(x) for x in beta_identity(1)],
                coefficients_full_flavor_trace=[str(x) for x in beta_identity(3)],
                registered_quartic_sensitivities=soft_result)


def tidy(x):
    if isinstance(x,dict):
        return {key:tidy(value) for key,value in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,(float,np.floating)) else x


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result = tidy(compute())
    if p.parse_args().check:
        assert result == read("independent_results.json")
        print("1700 independent: exact soft traces, both signs, physical stop eigenvalues verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
