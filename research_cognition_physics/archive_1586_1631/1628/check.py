"""A held-out atmospheric angle from one independently attached flavor column."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def limits(electron_ratio, muon_ratio, box):
    p, q = electron_ratio, muon_ratio
    norm = math.hypot(1, p, q)
    x = box["sin2_theta12"][0]
    y = box["sin2_theta13"][1]
    s12, c12 = math.sqrt(x), math.sqrt(1-x)
    s13, c13 = math.sqrt(y), math.sqrt(1-y)
    assert q > p*s13/c13 and q*c12*s13 < s12
    cot12 = c12/s12

    # Solve the triangle inequalities before comparing to atmospheric data.
    A, B = 1-q*cot12*s13, q+cot12*s13
    C = p*cot12*c13
    no_upper = math.atan2(B, A) + math.asin(C/math.hypot(A, B))
    io_lower = math.acos(p*s13/(c13*math.hypot(1,q))) - math.atan(q)

    z_low, z_high = box["sin2_theta23"]
    s_low, c_low = math.sqrt(z_low), math.sqrt(1-z_low)
    s_high, c_high = math.sqrt(z_high), math.sqrt(1-z_high)
    assert s_low > q*c_low
    gaps = {
        "normal": s12*(s_low-q*c_low) - c12*(s13*(c_low+q*s_low)+p*c13),
        "inverted": c13*(c_high-q*s_high) - p*s13
    }
    abs31_min = box["abs31_eV2"][0]
    eigengap_min = {"normal": abs31_min*box["solar_over_abs31"][0], "inverted": abs31_min}
    output = {}
    for ordering, numerator in gaps.items():
        overlap = numerator / norm
        assert 0 < overlap < 1
        # Positivity of D-d_min I also bounds its second moment by its first.
        strong_factor = 2*overlap*math.sqrt(1-overlap**2) if overlap < 1/math.sqrt(2) else 1.0
        output[ordering] = {
            "min_abs_overlap": overlap, "angle_to_nearest_orthogonal_ray_deg": math.degrees(math.asin(overlap)),
            "min_m_light_squared_over_source_if_background_PSD": overlap**2,
            "basic_background_spread_over_gap_min": overlap, "min_gap_eV2": eigengap_min[ordering],
            "variance_background_spread_over_gap_min": strong_factor,
            "min_background_spread_eV2": eigengap_min[ordering]*strong_factor
        }
    return {"normal_sin2_23_upper": math.sin(no_upper)**2,
            "inverted_sin2_23_lower": math.sin(io_lower)**2, "comparison": output}


def nonuniversal_witness(p, q, spec):
    """Fit contrasting two-column examples; the oscillation angles here are inputs."""
    import numpy as np
    x, y, z = spec["sin2_12_13_23"]
    s12, s13, s23 = np.sqrt([x,y,z])
    c12, c13, c23 = np.sqrt([1-x,1-y,1-z])
    U = np.array([[c12*c13, s12*c13, s13],
                  [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
                  [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13]])
    ratios = np.diag([p*p, q*q, 1.0])
    fits = {}
    for name, entry in spec["orderings"].items():
        spectrum = np.asarray(entry["relative_eigenvalues_eV2"])
        h = U@np.diag(spectrum)@U.T
        v = np.asarray(entry["relative_signs"]) * np.array([p,q,1])
        v /= np.linalg.norm(v)
        def overlap(slope):
            _, vectors = np.linalg.eigh(h-slope*ratios)
            n = vectors[:,0]
            return float(n@v)*np.sign(n[0])
        left, right = entry["bisection_eV2"]
        fl = overlap(left)
        assert fl * overlap(right) < 0
        for _ in range(54):
            mid = (left+right)/2
            if fl*overlap(mid)>0:
                left = mid
            else:
                right = mid
        slope = (left+right)/2
        values, vectors = np.linalg.eigh(h-slope*ratios)
        assert abs(overlap(slope)) < 1e-11 and values[1]>values[0]
        a = 1/np.sum((vectors[:,1:].T@v)**2/(values[1:]-values[0]))
        rem = h-slope*ratios-values[0]*np.eye(3)-a*np.outer(v,v)
        eig, cols = np.linalg.eigh(rem)
        assert a>0 and eig[2]>0 and np.max(abs(eig[:2]))<1e-13
        w = np.sqrt(eig[2])*cols[:,-1]
        if w[0]<0:
            w = -w
        # Choose the common shift so that the smallest diagonal entry is zero.
        diagonal = values[0]+slope*np.diag(ratios)
        shift = max(0.,-min(diagonal))
        diagonal += shift
        rebuilt = np.diag(diagonal) + a*np.outer(v,v) + np.outer(w,w)
        residual = float(np.linalg.norm(rebuilt-h-shift*np.eye(3),ord=2))
        assert min(diagonal)>=-1e-14 and residual<1e-14
        fits[name] = {"slope_eV2":float(slope), "fitted_source_weight_eV2":float(a),
                      "common_mass_squared_shift_eV2":float(shift),
                      "diagonal_eV2":diagonal.tolist(),"second_column_eV":w.tolist(),
                      "operator_residual_below_eV2":1e-14}
    return fits


def compute():
    cfg = read(HERE/"inputs.json")
    masses = read(HERE/cfg["charged_calibration"])["pole_MeV"]
    prediction = read(HERE/cfg["charged_prediction"])
    box = read(HERE/cfg["oscillation_windows"])["box"]
    tau = max(prediction["roots_MeV"])
    p, q = (math.sqrt(masses[i][0]/tau) for i in ("electron", "muon"))
    nominal = limits(p, q, box)
    # Outward rectangular enclosure of the frozen two-width charged calibration.
    tau_min = prediction["heavy_two_width_MeV"][0]
    p_hi, q_hi = (math.sqrt((masses[i][0]+2*masses[i][1])/tau_min) for i in ("electron","muon"))
    assert nominal["normal_sin2_23_upper"] < box["sin2_theta23"][0]
    assert nominal["inverted_sin2_23_lower"] > box["sin2_theta23"][1]
    return {"predicted_tau_MeV":tau,"source_ratios_e_mu_over_tau":[p,q],
            "two_width_upper_ratios":[p_hi,q_hi], "nominal":nominal,"outward":limits(p_hi,q_hi,box),
            "nonuniversal_witness":nonuniversal_witness(p,q,cfg["contrast_witness"])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    target = HERE/"results.json"
    result = compute()
    if args.check:
        if read(target) != result:
            raise SystemExit("1628 results differ")
        print("1628 PASS: angle holdout, two orders, charged calibration and spectral width")
    else:
        with target.open("x",encoding="utf-8") as file:
            file.write(json.dumps(result,indent=2)+"\n")
