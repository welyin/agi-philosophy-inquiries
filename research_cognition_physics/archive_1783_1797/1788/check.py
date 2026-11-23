"""Two independent projections of actual OPAL tables. Run with --check for replay."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

import numpy as np


HERE = Path(__file__).resolve().parent


def angular(data, luminosity_weighted=False):
    e = np.array(data["energy_GeV"])
    y = np.array(data["measured_pb_per_unit_abs_cos"])
    err = np.array(data["stat"])
    edges = np.linspace(0.0, 1.0, 11)
    x2 = (edges[1:]**2 + edges[:-1] * edges[1:] + edges[:-1]**2) / 3
    basis = np.column_stack((.75 * (1 + x2), 1.5 * (1 - x2)))
    proj = np.linalg.solve(basis.T @ basis, basis.T)
    np.testing.assert_allclose(proj @ basis, np.eye(2), atol=1e-13)
    exposure = np.array(data["table5_hadronic_luminosity_inverse_pb"]) if luminosity_weighted else np.ones(len(e))
    fractions = exposure / sum(exposure)
    scale = (e/e[0])**2 * fractions
    weights = proj[:, :, None] * scale
    values = np.sum(weights * y, axis=(1, 2))
    stat = np.sqrt(np.sum((weights * err)**2, axis=(1, 2)))
    shape = np.sum(np.abs(weights * y) * np.array(data["table11_sys_percent"])/100, axis=(1, 2))
    lumi = np.abs(proj @ y * scale) @ (np.array(data["table3_lumi_percent"])/100)
    sig = np.hypot(stat, shape + lumi)
    upper = values + 2*sig
    mu = data["born_mu_pb_at_first_energy"] * data["electron_interference_floor_round1785"]
    beta = np.sqrt(1 - 4*data["reference_mass_GeV"]**2/e**2)
    ordinary = np.sum(scale * data["table5_ZFITTER_reference_pb"] /
                      (1+2*data["ordinary_mass_cap_GeV"]**2/(data["minimum_sprime_over_s"]*e**2)))
    def mass_reach(w, target, scalar=False):
        def strength(mass):
            speed = np.sqrt(np.maximum(0.0, 1 - 4*mass**2/e**2))
            return mu*w*(fractions@(speed**3/4 if scalar else speed))
        if strength(0) <= target:
            return None
        lo, hi = 0., max(e)/2
        for _ in range(80):
            mid = (lo+hi)/2
            if strength(mid) > target:
                lo = mid
            else:
                hi = mid
        return lo
    scenarios = {}
    for w in (1/3, 1., 4., 9., 16.):
        f, s = mu*w*(fractions@beta), mu*w*(fractions@beta**3)/4
        scenarios[str(w)] = {"born_T_F_pb": f, "born_L_S_pb": s,
                             "unit_transfer_fermion_mass_reach_with_reference_GeV": mass_reach(w,upper[0]-ordinary),
                             "unit_transfer_scalar_mass_reach_no_reference_GeV": mass_reach(w,upper[1],True),
                             "eta_F_required_if_reference_exact_zero_leak": (upper[0]-ordinary)/f,
                             "eta_F_required_no_reference": upper[0]/f,
                             "eta_S_required_no_reference": upper[1]/s}
    return {"T_L_per_energy": (proj@y).T.tolist(), "projector": proj.tolist(),
            "scaled_T_L_pb": values.tolist(), "stat_pb": stat.tolist(),
            "arbitrary_shape_sys_pb": shape.tolist(), "arbitrary_lumi_pb": lumi.tolist(),
            "two_width_upper_pb": upper.tolist(), "ordinary_T_reference_pb": ordinary,
            "stat_absolute_correlation_upper_pb": np.sum(np.abs(weights)*err,axis=(1,2)).tolist(),
            "cases_by_charge_weight": scenarios}


def least_visible(response, stat):
    """Closest point to zero in whitened convex hull; witness also gives positive dual."""
    size = response.shape[1]
    wresponse = response/stat[:, None]
    gram = wresponse.T @ wresponse
    for count in range(1, size + 1):
        for support in combinations(range(size), count):
            ix = list(support)
            sub = np.linalg.solve(gram[np.ix_(ix, ix)], np.ones(count))
            mix = np.zeros(size)
            mix[ix] = sub/sub.sum()
            gram_mix = gram @ mix
            distance2 = float(mix @ gram_mix)
            if np.min(mix) < -1e-10 or min(gram_mix) < distance2-1e-9:
                continue
            w = (response @ mix)/(stat**2*distance2)
            assert np.min(w) >= -1e-10 and np.min(w @ response) >= 1-1e-8
            return mix, w / min(w @ response)
    raise RuntimeError("No active face satisfied the projection conditions")


def recoil(data):
    unit = data["bin_width_GeV"]*data["table_multiplier_per_GeV"]
    obs, stat, sys_err = (np.array(data[k])*unit for k in ("obs","stat","sys"))
    response = np.maximum(0, np.array(data["resolution"]) - data["rounding_absolute"])
    rows = []
    for size in range(1, 11):
        mix, w = least_visible(response[:, :size], stat)
        value = w @ obs
        stat_width = np.linalg.norm(w*stat)
        sys_width = w @ sys_err
        sigma = np.hypot(stat_width, sys_width)
        rows.append({"upper_edge_GeV": 10*size, "min_response": min(w @ response[:,:size]),
                     "worst_mix": mix.tolist(), "weights": w.tolist(), "value": value,
                     "stat": stat_width, "sys_sup": sys_width, "two_width_upper": value+2*sigma,
                     "all_stat_corr_two_width_upper": value+2*np.hypot(w@stat,sys_width)})
    return rows


def potentials(data):
    v = data["v_GeV"]
    rows = []
    for fraction in data["light_ZZ_shares"]:
        rot = np.array([[np.sqrt(fraction), np.sqrt(1-fraction)],
                        [np.sqrt(1-fraction), -np.sqrt(fraction)]])
        hessian = rot @ np.diag([data["light_GeV"]**2, data["heavy_GeV"]**2]) @ rot.T
        vals, vec = np.linalg.eigh(hessian)
        assert min(vals) > 0
        rows.append({"even_mass_GeV": np.sqrt(vals).tolist(), "light_ZZ_share": vec[0,0]**2,
                     "positive_quadratic_form_abc": (hessian/v**2).tolist(),
                     "d": 2*(data["charged_GeV"]/v)**2, "e": 2*(data["odd_GeV"]/v)**2,
                     "f": 1.0})
    return rows


def canonical(value):
    if isinstance(value, dict):
        return {k:canonical(v) for k,v in value.items()}
    if isinstance(value, list):
        return [canonical(v) for v in value]
    return round(float(value), 10) if isinstance(value, (float,np.floating)) else value


if __name__ == "__main__":
    inputs = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    output = canonical({"round":1788, "angular":angular(inputs["angular"]),
                        "angular_luminosity_weighted":angular(inputs["angular"], True),
                        "recoil":recoil(inputs["recoil"]), "potentials":potentials(inputs["potential_example"])})
    body = json.dumps(output,ensure_ascii=False,indent=2) + "\n"
    path = HERE/"results.json"
    if "--check" in sys.argv:
        assert json.loads(path.read_text(encoding="utf-8")) == output, "results changed"
    else:
        path.write_text(body, encoding="utf-8")
    print("1788", hashlib.sha256(body.encode()).hexdigest()[:12])
