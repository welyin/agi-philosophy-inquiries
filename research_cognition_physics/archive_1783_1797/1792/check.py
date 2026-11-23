"""Tree-level charged-partner norm and finite-width certificates."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def charged_weights(mass, p):
    v = p["vev_GeV"]
    top_mass = p["top_kinematic_GeV"]
    top_yuk = math.sqrt(2) * p["top_matching_GeV"] / v
    weights = []
    for down_mass in p["down_kinematic_GeV"]:
        t, b = (top_mass / mass) ** 2, (down_mass / mass) ** 2
        a = 1 - t - b
        phase = a * a - 4 * t * b
        assert mass > top_mass + down_mass
        weights.append(3 * top_yuk**2 * phase**1.5 / (16 * math.pi * a))
    w, h = p["W_GeV"], p["higgs_GeV"]
    assert mass > w + h
    phase = (1 - ((w + h) / mass) ** 2) * (1 - ((w - h) / mass) ** 2)
    wh = mass**2 * phase**1.5 / (16 * math.pi * v**2)
    return np.array(weights), wh


def measure_delta(phase_degrees, c, magnitude):
    return abs(magnitude * np.exp(1j * math.radians(phase_degrees)) - c)


def width_min(delta, cap, t, w):
    if delta == 0:
        return 0.0
    u = min(cap, delta * math.sqrt(t / w))
    return t * delta**2 / u + w * u


def delta_max(budget, cap, t, w):
    u = min(cap, budget / (2 * w))
    return math.sqrt(u * (budget - w * u) / t)


def flavor_weights(p):
    ref = HERE / p["CKM_source"]
    point = json.loads(ref.read_text(encoding="utf-8"))["rows"]["archived_time_integrated_angle"]
    uv, cv = point["moduli"]["Vub"], point["moduli"]["Vcb"]
    bv2 = 1 - uv**2 - cv**2
    ratio = point["td_ts"] ** 2
    return np.array([(1 - bv2) * ratio / (1 + ratio), (1 - bv2) / (1 + ratio), bv2])


def phase_limit(delta, c, mag):
    if delta < abs(mag - c):
        return None
    return math.degrees(math.acos(max(-1, min(1, (mag**2 + c**2 - delta**2) / (2 * mag * c)))))


def calculate(p):
    flavor = flavor_weights(p)
    c, mag = p["c_trial"], p["kappa_magnitude"]
    cap = 1 - c**2
    def evaluate(mass):
        t, w = charged_weights(mass, p)
        inv = flavor / t
        tc = 1 / sum(inv)
        light_only = 1 / sum(inv[:2])
        widths = {str(deg): width_min(measure_delta(deg, c, mag), cap, tc, w)
                  for deg in p["phase_trials_deg"]}
        budgets = []
        for budget in p["pool_budget"]:
            for fraction in p["relative_born_width_floor"]:
                available = budget / fraction
                maximum = delta_max(available, cap, tc, w)
                no_b = delta_max(available, cap, light_only, w)
                b_min_at_real = max(0, abs(mag-c)-no_b) / math.sqrt(flavor[2])
                assert abs(width_min(maximum, cap, tc, w)-available) < 2e-13
                budgets.append({"Gamma_pool":budget,"born_floor":fraction,
                                "delta_max":maximum,"alpha_max_deg":phase_limit(maximum, c, mag),
                                "delta_max_if_no_tb":no_b,
                                "sum_abs_r_zb_min_at_positive_real":b_min_at_real})
        return dict(Q_GeV=mass, T_by_flavor=t.tolist(), T_effective=tc, U_wh=w,
                    min_Born_Gamma_over_M_sum_by_alpha_deg=widths, thresholds=budgets)
    def maximum_min_mass(phase, budget):
        delta = measure_delta(phase,c,mag)
        def difference(q):
            t, w = charged_weights(q,p)
            return width_min(delta, cap, 1 / sum(flavor / t), w) - budget
        low = max(p["higgs_GeV"] + p["W_GeV"], p["top_kinematic_GeV"] + max(p["down_kinematic_GeV"])) + 1e-4
        high = low * 2
        if difference(low) > 0:
            return {"below_open_Wh_threshold": True}
        while difference(high) < 0:
            high *= 2
        for _ in range(90):
            mid = (low + high) / 2
            if difference(mid) < 0:
                low = mid
            else:
                high = mid
        return {"upper_bound_GeV":high}
    direct = []
    for mass, upper in zip(p["atlas_charged_Wh"]["mass_range_GeV"], p["atlas_charged_Wh"]["limit_range_pb"]):
        row = evaluate(mass)
        budget = 0.10
        cut = next(b for b in row["thresholds"] if b["Gamma_pool"] == budget and b["born_floor"] == 1)
        for phase in [0, 15, 30]:
            d = measure_delta(phase, c, mag)
            if d > cut["delta_max"]:
                continue
            pb = max(0, d - cut["delta_max_if_no_tb"]) / math.sqrt(flavor[2])
            lower_factor = row["U_wh"] * pb**2 / budget
            direct.append(dict(mass_GeV=mass, phase_deg=phase, limit_pb=upper,
                               signal_factor_lower=lower_factor,
                               needed_unit_chiral_production_times_transfer_pb=upper/lower_factor))
    return dict(ckm_top_row_moduli_squared=flavor.tolist(),
                width_rows=[evaluate(m) for m in p["mass_tests_GeV"]],
                one_pole_charged_Wh_recast_if_narrow=direct,
                max_allowed_lightest_pole_if_all_open_at_pool_0p1={
                    str(phi):maximum_min_mass(phi,.1) for phi in p["phase_trials_deg"]})


def check_saturation(p):
    # Exactly sharp lower envelope for the two displayed partial widths at fixed Q.
    flavor = flavor_weights(p)
    t, w = charged_weights(500, p)
    c = p["c_trial"]
    d = measure_delta(15, c, p["kappa_magnitude"])
    factor = sum(flavor/t)
    u = min(1-c*c, d / math.sqrt(w * factor))
    z = d * np.sqrt(flavor) / (math.sqrt(u) * t * factor)
    assert abs(math.sqrt(u)*np.dot(np.sqrt(flavor),z) - d) < 1e-13
    assert abs(np.dot(t,z**2) + w*u - width_min(d,1-c*c,1/factor,w)) < 1e-14
    for j, down_mass in enumerate(p["down_kinematic_GeV"]):
        top = p["top_kinematic_GeV"]/500
        down = down_mass/500
        a = 1 - top**2 - down**2
        r = 2*top*down/a*z[j]
        k = math.sqrt(a*a-4*(top*down)**2)
        bare = 3/(16*math.pi) * 2*(p["top_matching_GeV"]/p["vev_GeV"])**2
        assert abs(bare*k*(a*(z[j]**2+r*r)-4*top*down*z[j]*r)-t[j]*z[j]**2) < 1e-13


def distinct_charged_masses(p):
    v, h, c = p["vev_GeV"], p["higgs_GeV"], p["c_trial"]
    r = math.sqrt(1-c*c)
    mix = np.array([[c,-r],[r,c]])
    neutral = mix @ np.diag([h**2, (2*h)**2]) @ mix.T
    quartic = neutral/v**2
    assert np.all(np.linalg.eigvalsh(quartic)>0)
    assert np.allclose(neutral @ mix[:,0],h*h*mix[:,0])
    top = (np.exp(1j*np.deg2rad(15))-c)/r
    assert abs(c+r*top-np.exp(1j*np.deg2rad(15)))<1e-13
    return dict(re_Yextra_tt_over_yt=float(top.real),im_Yextra_tt_over_yt=float(top.imag),
                K=quartic.tolist(),charged_quartic_for_250_and_500_GeV=[2*m*m/v**2 for m in [250,500]])


def readable(x):
    if isinstance(x, dict):
        return {k:readable(v) for k,v in x.items()}
    if isinstance(x, list):
        return [readable(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x, float) else x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    inputs = HERE / "inputs.json"
    p = json.loads(inputs.read_text(encoding="utf-8"))
    check_saturation(p)
    answer = readable(calculate(p))
    answer["independent_charged_mass_witness"] = readable(distinct_charged_masses(p))
    alternative = dict(p, top_matching_GeV=p["alternative_top_matching_GeV"])
    alt = calculate(alternative)
    answer["alternative_top_matching"] = readable({
        "mass_GeV":alternative["top_matching_GeV"],
        "phase_cap_deg_at_500GeV_pool_0p1":alt["width_rows"][2]["thresholds"][2]["alpha_max_deg"],
        "phase_cap_deg_at_500GeV_pool_0p1_85percent_floor":alt["width_rows"][2]["thresholds"][3]["alpha_max_deg"],
        "upper_lightest_mass_at_30deg_pool_0p1":alt["max_allowed_lightest_pole_if_all_open_at_pool_0p1"]["30"]})
    answer["inputs_sha256"] = hashlib.sha256(inputs.read_bytes()).hexdigest()
    raw = json.dumps(answer, ensure_ascii=False, indent=2) + "\n"
    results = HERE / "results.json"
    if args.check:
        assert json.loads(results.read_text(encoding="utf-8")) == answer
        print("1792 OK", hashlib.sha256(raw.encode()).hexdigest()[:12])
    else:
        results.write_text(raw,encoding="utf-8")


if __name__ == "__main__":
    main()
