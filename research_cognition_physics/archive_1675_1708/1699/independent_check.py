"""Independent two-doublet eigenvalue and soft-boundary tests at the 30 TeV interface."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def gauge_at_soft(sh, common):
    scale = sh["soft_GeV"]
    ref = common["mu_ref_GeV"]
    g2 = np.square([common[key][0] for key in ("g1_GUT", "g2", "g3")])
    logs = np.array([41/10, -19/6, -7]) * math.log(scale / ref)
    for poles, coeff in [
        (sh["colored_fermions_GeV"], [4/15, 0, 2/3]),
        (sh["weak_fermions_GeV"], [2/5, 2/3, 0]),
    ]:
        logs += np.asarray(coeff) * sum(max(0, math.log(scale / mass)) for mass in poles)
    return 1 / np.sqrt(1/g2 - logs/(8*math.pi**2))


def envelopes(quartic, yt, dterm, ms, light_coef, cw_coef, allowance):
    stop_max = 9*yt**4/(16*math.pi**2)
    lower_cos = math.sqrt(max(0, (quartic-stop_max-allowance)/dterm))
    gap = ms**2-light_coef
    assert 0 < lower_cos < 1
    sine2_min = yt**2/(4*math.pi)
    required = gap*lower_cos
    diagonal_data = {}
    for key, projected in (("calibrated_fermion_only", light_coef), ("actual_neutral_soft_CW", cw_coef)):
        max_diag = abs(projected)/sine2_min
        full_drop = max_diag/gap
        missing = required - max_diag
        assert missing > 0
        diagonal_data[key] = dict(max_abs_Hu_mass_squared_GeV2=max_diag,
                                  min_additional_diagonal_difference_GeV2=missing,
                                  equal_bare_diagonals_quartic_shortfall_at_max_stop=
                                  quartic-stop_max-dterm*full_drop**2)
    # This point saturates the analytic bound. Directly diagonalize light/heavy sectors.
    ang = 0.5*math.acos(-lower_cos)
    v = np.array([math.cos(ang), math.sin(ang)])
    w = np.array([-v[1], v[0]])
    q = light_coef*np.outer(v,v) + ms**2*np.outer(w,w)
    eig, directions = np.linalg.eigh(q)
    np.testing.assert_allclose(eig, [light_coef, ms**2], rtol=2e-10)
    np.testing.assert_allclose(abs(v @ directions[:,0]), 1, atol=1e-13)
    np.testing.assert_allclose(abs(q[0,0]-q[1,1]), required, rtol=1e-13)
    np.testing.assert_allclose(abs(q[0,1]), gap*math.sin(2*ang)/2, rtol=1e-13)
    np.testing.assert_allclose(dterm*math.cos(2*ang)**2 + stop_max+allowance, quartic)
    return dict(stop_max=stop_max, minimum_abs_cos2beta=lower_cos, gap_GeV2=gap,
                required_diagonal_difference_GeV2=required,
                attainable_quadratic_matrix_GeV2=q.tolist(), mechanisms=diagonal_data)


def calculate():
    shared = read("../1698/results.json")["shared_witness"]
    soft_cw = read("../1698/independent_results.json")
    frozen = read("results.json")
    soft = frozen["soft_interface"]
    quartic, yt, gY, g2 = (soft[key] for key in ("lambda", "yt", "gY", "g2"))
    factors = gauge_at_soft(shared, read("../1369/inputs.json"))
    np.testing.assert_allclose([gY, g2, soft["g3"]],
                               [factors[0]*math.sqrt(.6), factors[1], factors[2]], rtol=3e-7)
    dterm = (gY*gY+g2*g2)/8
    light = soft_cw["target_potential_coefficient_GeV2"]
    cw = soft_cw["frozen_30TeV_universal_soft"]["soft_scale"]["neutral_CW_coefficient_GeV2"]
    ms = shared["soft_GeV"]
    result = {str(extra): envelopes(quartic, yt, dterm, ms, light, cw, extra)
              for extra in (0, .01)}
    # Check the stop-threshold maximum on the entire real mixing axis by completing a square.
    for x in (-5, -math.sqrt(6), -1, 0, 2, math.sqrt(6), 6):
        np.testing.assert_allclose(x*x - x**4/12, 3 - (x*x-6)**2/12, atol=1e-13)
    return dict(common_scale_GeV=ms, effective_lambda=quartic, effective_top=yt,
                independently_integrated_gauge=factors.tolist(), D_term_factor=dterm,
                results_by_other_quartic_allowance=result)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    options = cli.parse_args()
    result = calculate()
    # Discard platform-level differences in eigensolver and elementary functions.
    from decimal import Decimal
    record = json.loads(json.dumps(result), parse_float=lambda s: float(f"{Decimal(s):.10g}"))
    path = HERE / "independent_results.json"
    if options.check:
        assert record == json.loads(path.read_bytes())
        print("PASS: gauge thresholds, sharp two-doublet tilt, maximum stop mixing, two independent neutral budgets")
    else:
        with path.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(record, indent=2) + "\n")
