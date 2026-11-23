"""Infer one-loop threshold logarithms from two electroweak current correlators."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
D = Decimal
PI = D("3.1415926535897932384626433832795028841971693993751058209749445923")


def data(path):
    return json.loads((HERE / path).read_bytes(), parse_float=D)


def moments(doublets, triplets):
    readings = []
    for order in (1, 2):
        light = sum(m ** (-2 * order) for m in doublets)
        color = sum(m ** (-2 * order) for m in triplets)
        readings.append((light + color / 3, light / 2))
    return readings


def sum_log_relative(m1, m2, count, reference):
    # m1 and m2 are positive sums of the first and second powers of M^{-2}.
    assert m1 > 0
    if count == 1:
        product = m1
    else:
        assert count == 2 and m2 > 0
        assert m1**2 * (1 - D("1e-40")) <= 2 * m2
        product = (m1**2 - m2) / 2
        assert product > 0
    return -(product * reference ** (2 * count)).ln() / 2


def strong_prediction(readings, counts, reference, ew):
    doublets, triplets = counts
    weak_moments = [2 * w for _, w in readings]
    color_moments = [3 * (a - 2 * w) for a, w in readings[:triplets]]
    logs_d = sum_log_relative(*weak_moments, doublets, reference)
    logs_t = sum_log_relative(color_moments[0], color_moments[-1], triplets, reference)
    doublet_indices = [D(3)/5, D(1), D(0)]
    triplet_indices = [D(2)/5, D(0), D(1)]
    beta = [D(41)/10, -D(19)/6, -D(7)]
    j, k, h = [], [], []
    for u, v, b in zip(doublet_indices, triplet_indices, beta):
        weight = doublets * u + triplets * v
        j.append(6 + weight)
        k.append(b/(2*PI) + weight/(3*PI))
        h.append((logs_d*u + logs_t*v)/(3*PI))
    rhs = [a + shift for a, shift in zip(ew, h[:2])]
    det = j[0]*k[1] - j[1]*k[0]
    c = (rhs[0]*k[1] - rhs[1]*k[0])/det
    t = (j[0]*rhs[1] - j[1]*rhs[0])/det
    assert c > 0 and t > 0
    return 1 / (j[2]*c + k[2]*t - h[2]), reference*t.exp()


def result():
    inp = data("../1369/inputs.json")
    old = data("../1807/independent_results.json")["common_MZ"]
    ref = inp["mu_ref_GeV"]
    ew = [4*PI/inp[key][0]**2 for key in ("g1_GUT", "g2")]
    # Integrate the two Taylor coefficients independently over the Feynman parameter.
    import numpy as np
    nodes, weights = np.polynomial.legendre.leggauss(8)
    x = (nodes + 1)/2
    for n, normalization in ((1, 60), (2, 560)):
        val = weights @ ((x*(1-x))**(n+1))/(4*n)
        assert abs(val - 1/normalization) < 2e-16
    answer = {}
    for key in ("two_doublets_one_triplet", "two_split_complete_pairs"):
        sample = old[key]
        masses_d = sample["M_D_GeV"]
        masses_t = sample["M_T_GeV"]
        responses = moments(masses_d, masses_t)
        counts = (len(masses_d), len(masses_t))
        prediction, scale = strong_prediction(responses, counts, ref, ew)
        calibrated = inp["g3"][0]**2/(4*PI)
        assert abs(prediction / calibrated - 1) < D("5e-14")
        assert abs(scale / sample["matching_scale_GeV"] - 1) < D("5e-14")
        if counts[1] == 1:
            three_inputs = [responses[0], (None, responses[1][1])]
            assert strong_prediction(three_inputs, counts, ref, ew) == (prediction, scale)
        precision = [(a - 2*w)/(3*a) for a, w in responses]
        # If both normalized components carry independent errors <= epsilon*a,
        # a-2w can be moved all the way to zero at the epsilon value below.
        answer[key] = {
            "normalized_photon_and_weak_moments_GeV": [[float(a), float(w)] for a, w in responses],
            "relative_error_allowing_color_signal_to_vanish_per_moment": [float(p) for p in precision],
            "reconstructed_alpha3_from_EW_and_synthetic_moments": float(prediction),
            "reconstructed_matching_scale_GeV": float(scale),
        }
        if counts == (2, 1):
            orig = sum(m**(-2) for m in masses_d)
            frac = D("0.998")
            altered = [(f*orig)**(-D("0.5")) for f in (frac, 1-frac)]
            new = moments(altered, masses_t)
            assert max(abs(new[0][i]/responses[0][i]-1) for i in (0, 1)) < D("1e-70")
            new_prediction, new_scale = strong_prediction(new, counts, ref, ew)
            assert min(altered + masses_t) > ref and max(altered + masses_t) < new_scale
            assert abs(new_prediction-prediction) > D("0.001")
            answer["same_two_dimension_six_readings_competitor"] = {
                "doublet_masses_GeV": list(map(float, altered)),
                "matching_scale_GeV": float(new_scale),
                "new_alpha3": float(new_prediction),
                "alpha3_shift": float(new_prediction-prediction),
                "weak_second_moment_ratio": float(new[1][1]/responses[1][1])}
    return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    with localcontext() as ctx:
        ctx.prec = 110
        raw = json.dumps(result(), indent=2, allow_nan=False) + "\n"
    out = HERE / "independent_results.json"
    if parser.parse_args().check:
        assert json.loads(raw) == json.loads(out.read_bytes())
    else:
        out.write_text(raw, encoding="utf-8")
    print("1809 independent", hashlib.sha256(raw.encode()).hexdigest()[:12])
