"""Round 1619. Author check for a conditional charged-lepton pole mass prediction."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = json.loads((HERE / "inputs.json").read_text(encoding="utf8"))


def roots(electron, muon):
    a, b = math.sqrt(electron), math.sqrt(muon)
    offset = math.sqrt(3 * (a * a + 4 * a * b + b * b))
    return [(2 * (a + b) - offset) ** 2, (2 * (a + b) + offset) ** 2]


def tau_at(exponent, electron, muon, log_offset=0.0):
    # log_offset is delta_heavy minus the equal electron and muon offsets.
    return math.exp((math.log(roots(electron**exponent, muon**exponent)[1])
                     - log_offset) / exponent)


def inverse_slope(tau, electron, muon, log_offset):
    low, high = .98, 1.02
    assert tau_at(low, electron, muon, log_offset) > tau
    assert tau_at(high, electron, muon, log_offset) < tau
    for _ in range(64):
        mid = (low + high) / 2
        if tau_at(mid, electron, muon, log_offset) > tau:
            low = mid
        else:
            high = mid
    return (low + high) / 2 - 1


def compare_sample(sample, em, mm, slope_qed):
    center = sample["center"]
    low = center - 2 * (sample["stat"][0] + sample["syst"][0])
    high = center + 2 * (sample["stat"][1] + sample["syst"][1])
    e_lo, e_hi = em[0] - 2 * em[1], em[0] + 2 * em[1]
    mu_lo, mu_hi = mm[0] - 2 * mm[1], mm[0] + 2 * mm[1]
    out = {}
    for drift in INPUT["threshold_log_oscillations"]:
        lower = inverse_slope(high, e_lo, mu_lo, drift)
        upper = inverse_slope(low, e_hi, mu_hi, -drift)
        out[str(drift)] = {"log_slope": [lower, upper],
                           "alphaF_over_4alpha": [1 - upper/slope_qed, 1 - lower/slope_qed]}
    # Minimum possible osc(delta) with alpha_F=0. The failing QED-only prediction is low.
    highest_qed = tau_at(1 + slope_qed, e_hi, mu_hi)
    assert highest_qed < low
    min_drift = (1 + slope_qed)*math.log(low/highest_qed)
    prediction = roots(em[0], mm[0])[1]
    side = 0 if center > prediction else 1
    return {"two_width_MeV": [low, high], "center_minus_predicted_MeV": center - prediction,
            "difference_over_quadrature_error": (center - prediction) /
            math.hypot(sample["stat"][side], sample["syst"][side]),
            "slope_windows": out, "min_log_osc_to_admit_QED_only": min_drift}


def certificate():
    em, mm = (INPUT["pole_MeV"][key] for key in ("electron", "muon"))
    slope_qed = 3/(2*math.pi*INPUT["alpha_inverse"])
    small, heavy = roots(em[0], mm[0])
    e_lo, e_hi = em[0] - 2 * em[1], em[0] + 2 * em[1]
    m_lo, m_hi = mm[0] - 2 * mm[1], mm[0] + 2 * mm[1]
    samples = {name: compare_sample(s, em, mm, slope_qed)
               for name, s in INPUT["tau_samples"].items()}
    values = [em[0], mm[0], heavy]
    assert abs(sum(values)/sum(math.sqrt(v) for v in values)**2-2/3) < 1e-14
    return {"roots_MeV": [small, heavy], "heavy_two_width_MeV":
            [roots(e_lo, m_lo)[1], roots(e_hi, m_hi)[1]], "pure_QED_slope": slope_qed,
            "QED_only_dictionary_prediction_MeV": tau_at(1+slope_qed, em[0], mm[0]),
            "family_vector_M12_over_M11": math.sqrt((em[0]+mm[0])/(2*em[0])),
            "family_vector_M23_over_M12": math.sqrt((mm[0]+heavy)/(em[0]+mm[0])),
            "measurements": samples}


def near(given, calculated):
    if isinstance(calculated, dict):
        return given.keys() == calculated.keys() and all(near(given[k], v) for k,v in calculated.items())
    if isinstance(calculated, list):
        return len(given) == len(calculated) and all(near(x,y) for x,y in zip(given,calculated))
    return math.isclose(given, calculated, rel_tol=5e-12, abs_tol=5e-14)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    options = parser.parse_args()
    assert not (options.write and options.check)
    results = certificate()
    path = HERE / "results.json"
    if options.write:
        with path.open("x", encoding="utf8") as stream:
            json.dump(results, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    else:
        assert near(json.loads(path.read_text(encoding="utf8")), results)
        print("PASS 1619")


if __name__ == "__main__":
    main()
