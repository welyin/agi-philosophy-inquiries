"""Running-width SM W interpolation, inverse thresholds, and observable consistency."""
import argparse
import itertools
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def predict(point, model, coefficient_key="eq8"):
    reference = model["reference"]
    mass_h = point["m_H"] / 100
    log_h = math.log(mass_h)
    top = (point["m_t"] / reference["m_t"])**2 - 1
    z = point["m_Z"] / reference["m_Z"] - 1
    alpha = (point["delta_alpha_had"] + model["delta_alpha_lept"]) / reference["delta_alpha"] - 1
    strong = point["alpha_s"] / reference["alpha_s"] - 1
    basis = (1, -log_h, -log_h**2, log_h**4, mass_h**2 - 1,
             -alpha, top, -top**2, -log_h * top, mass_h**2 * top, -strong, z)
    return sum(weight * value for weight, value in zip(model[coefficient_key], basis, strict=True))


def delta_r(mass_w, mass_z, model):
    norm = math.sqrt(2) * model["Gmu_GeV_minus2"] * model["alpha_inverse"] / math.pi
    return norm * mass_w**2 * (1 - (mass_w / mass_z)**2) - 1


def corners(data, factor, skip_top=False):
    entries = {k: (v[0] - factor*v[1], v[0] + factor*v[1])
               for k, v in data.items() if not (skip_top and k == "m_t")}
    for values in itertools.product(*entries.values()):
        yield dict(zip(entries.keys(), values))


def inverse_top(wanted, point, model):
    # Check only inside the one-parameter marginal validation envelope of the paper.
    mid = model["reference"]["m_t"]
    span = 2 * model["historical_marginal_width"]["m_t"]
    low, high = mid - span, mid + span
    if not predict({**point, "m_t": low}, model) <= wanted <= predict({**point, "m_t": high}, model):
        return None
    for _ in range(56):
        probe = (low + high) / 2
        if predict({**point, "m_t": probe}, model) < wanted:
            low = probe
        else:
            high = probe
    return (low + high) / 2


def historical_radius(point, model):
    total = 0
    for k, scale in model["historical_marginal_width"].items():
        value = (point["delta_alpha_had"] + model["delta_alpha_lept"]
                 if k == "delta_alpha" else point[k])
        total += ((value - model["reference"][k]) / scale)**2
    return math.sqrt(total)


def calculate(inputs):
    model = inputs["formula"]
    constants = inputs["calibration"]
    center = {k:v[0] for k, v in constants.items()}
    limit_top = model["reference"]["m_t"] + 2 * model["historical_marginal_width"]["m_t"]
    budget = model["interpolation_budget_GeV"] + model["missing_orders_budget_GeV"]
    baseline = predict(center, model)
    alternate = predict({**center, "m_t": inputs["alternative_top"][0]}, model)
    increments = {k: 1000*(predict({**center, k:v[0]+v[1]}, model)-baseline)
                  for k, v in constants.items()}
    holdouts = {}
    windows = {}
    for label, (w, error) in inputs["measured"].items():
        holdouts[label] = {"observed_minus_SM_MeV":1000*(w - baseline),
                           "delta_r_eff":delta_r(w, center["m_Z"], model),
                           "delta_r_minus_SM":delta_r(w, center["m_Z"], model)-delta_r(baseline, center["m_Z"], model),
                           "nominal_required_top_GeV_or_null":inverse_top(w, center, model)}
    for width in (1,2):
        vertices = list(corners(constants, width))
        preds = [predict(p, model) for p in vertices]
        min_w, max_w = min(preds)-budget, max(preds)+budget
        key = str(width)
        windows[key] = {"polynomial_mass_GeV":[min(preds),max(preds)],
                        "mass_with_budget_GeV":[min_w,max_w],
                        "max_historical_scaled_radius_if_uncorrelated":max(historical_radius(p, model) for p in vertices)}
        for label, (w, error) in inputs["measured"].items():
            wlo, whi = w-width*error, w+width*error
            intervals = []
            for p in vertices:
                for prediction in [predict(p,model)-budget,predict(p,model)+budget]:
                    for measured in (wlo,whi):
                        intervals.append(delta_r(measured,p["m_Z"],model)-delta_r(prediction,p["m_Z"],model))
            others = list(corners(constants, width, skip_top=True))
            envelope = lambda t, side: (max if side == "high" else min)(
                predict({**p, "m_t":t},model) for p in others)
            # A larger prediction at a given t lowers the mass required for the lower edge.
            best = max(others, key=lambda p:predict({**p, "m_t":limit_top},model))
            worst = min(others, key=lambda p:predict({**p, "m_t":limit_top},model))
            for t in (164.1,174.3,184.5):
                assert abs(predict({**best, "m_t":t},model) - envelope(t,"high")) < 1e-12
                assert abs(predict({**worst,"m_t":t},model) - envelope(t,"low")) < 1e-12
            holdouts[label][key] = {
                "excess_MeV_range":[1000*(wlo-max_w),1000*(whi-min_w)],
                "effective_delta_r_range":[min(intervals),max(intervals)],
                "required_top_lower_and_upper_within_old_marginal_or_null":[
                    inverse_top(wlo-budget,best,model),inverse_top(whi+budget,worst,model)]}
    observations = list(inputs["measured"].items())
    pairwise = {}
    for (name1,(mass1,e1)),(name2,(mass2,e2)) in itertools.combinations(observations,2):
        pairwise[name1+"_"+name2] = abs(mass1-mass2)/(e1+e2)
    joint_factor = max(pairwise.values())
    joint_low = max(m - joint_factor*e for _,(m,e) in observations)
    joint_high = min(m + joint_factor*e for _,(m,e) in observations)
    assert abs(joint_low-joint_high) < 1e-10
    reference_test = {**center,"m_H":100,"m_t":174.3,"m_Z":91.1875,"alpha_s":.119,
                      "delta_alpha_had":.05907-model["delta_alpha_lept"]}
    assert abs(predict(reference_test,model)-80.3779) < 1e-12
    return {"central_W_GeV":baseline,"delta_r_effective_SM":delta_r(baseline,center["m_Z"],model),
            "alternate_CMS_top_W_GeV":alternate,
            "narrow_formula_minus_main_MeV":1000*(predict(center,model,"eq9")-baseline),
            "plus_one_width_shifts_MeV":increments,"boxes":windows,"holdouts":holdouts,
            "central_others_W_at_historical_upper_top_GeV":predict({**center,"m_t":limit_top},model),
            "pairwise_minimax_factors":pairwise,"all_minimax_factor":joint_factor,
            "all_minimax_mass_GeV":(joint_low+joint_high)/2,
            "running_to_pole_W_shift_example_MeV":1000*baseline*(1-1/math.sqrt(
                1+(inputs["width_examples"]["W"]/baseline)**2))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = calculate(json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig")))
    saved = ROOT/"results.json"
    def stable(x):
        if isinstance(x,dict):
            return {k:stable(v) for k,v in x.items()}
        if isinstance(x,list):
            return [stable(item) for item in x]
        return round(x,10) if isinstance(x,float) else x
    result = stable(result)
    if args.write_new:
        with saved.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(saved.read_text(encoding="utf-8-sig"))
        print("1394 W-mass interpolation, inverse boundaries and common residual verified")
