"""Round 1543: independently construct generator masses, gauge match and dilepton trace."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PI = math.pi


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def run_g(start, initial, beta, scale):
    return (initial**-2 - beta * math.log(scale / start) / (8 * PI**2))**-.5


def mass_ratio(g_y, g_s, defect=0):
    return (g_s**2 / g_y**2 - 2 / 3 + defect)**-.5


def couplings(base, scale, signs=(0, 0)):
    g1, g3 = (base[k][0] + s * base[k][1]
              for k, s in zip(("g1_GUT", "g3"), signs))
    start = base["mu_ref_GeV"]
    return (run_g(start, math.sqrt(3/5) * g1, 41/6, scale),
            run_g(start, g3, -7, scale))


def action(vacuum, gen):
    response = gen @ vacuum
    return float(2 * np.vdot(response, response).real)


def generators(n, dim_two):
    flip = np.zeros((n, n))
    flip[0, -1] = flip[-1, 0] = .5
    charge = np.diag([1/6] * 3 + [-.5]) if n == 4 else np.diag([.5, -.5])
    vacuum = np.eye(n)[:, -1 if n == 4 else 0]
    if dim_two:
        def tensor_sum(matrix):
            return np.kron(matrix, np.eye(n)) + np.kron(np.eye(n), matrix)
        flip, charge = tensor_sum(flip), tensor_sum(charge)
        vacuum = np.kron(vacuum, vacuum)
    return flip, charge, vacuum


def generator_test():
    coefficients = {}
    for title, tensor in (("fundamental", False), ("symmetric", True)):
        u, x, left = generators(4, tensor)
        w, t3, right = generators(2, tensor)
        v = np.kron(left, right)
        light = np.kron(x, np.eye(len(right))) + np.kron(np.eye(len(left)), t3)
        assert np.linalg.norm(light @ v) < 1e-13
        coeff = [action(left, u), action(right, w), action(right, t3)]
        assert np.allclose(coeff, [1, 1, 2] if tensor else [.5, .5, .5])
        coefficients[title] = coeff
    return coefficients


def benchmark(inputs, base):
    b = inputs["benchmark"]
    f, rho, v = b["f_GeV"], b["rho"], b["v_EW_GeV"]
    scale = b["matching_scale_GeV"]
    gy, gs = couplings(base, scale)
    gl = run_g(base["mu_ref_GeV"], base["g2"][0], -19/6, scale)
    gr = gs * mass_ratio(gy, gs)
    gn = math.hypot(gr, math.sqrt(1.5) * gs)
    gz = math.hypot(gy, gl)
    w10 = f * math.sqrt(rho)
    u2, wr2, zr2 = (gs*f)**2, (gr*f)**2, (gn*f)**2*(1+rho)
    neutral = np.array([gz, -gr**2 / gn])
    poles = np.linalg.eigvalsh(np.diag([0., zr2]) + np.outer(neutral, neutral)*v**2/4)
    neutral_corrected = sum(poles) - (gz**2+gr**4/gn**2)*v**2/4
    wr_physical_sq = wr2 + gr**2*v**2/4
    reconstructed = neutral_corrected / gn**2 - (wr_physical_sq-gr**2*v**2/4)/gr**2
    assert math.isclose(reconstructed, w10**2, rel_tol=1e-11)
    strengths = np.array(b["y_singular"])
    omega = np.exp(2j*PI/3)
    q = omega ** np.outer(np.arange(3), np.arange(3)) / math.sqrt(3)
    y = q.T @ np.diag(strengths) @ q
    width_channels = [abs(y[i, j])**2/(16*PI*(1+(i==j)))
                      for i in range(3) for j in range(i,3)]
    trace = np.linalg.norm(y)**2
    assert math.isclose(32*PI*sum(width_channels), trace, rel_tol=1e-12)
    neutrino_trace = sum(np.linalg.svd(w10*y, compute_uv=False)**2)
    assert math.isclose(neutrino_trace, 32*PI*sum(width_channels)*reconstructed, rel_tol=1e-11)
    eps = inputs["pole_squared_rel_error"]
    # For true x in x_observed [1-eps,1+eps], infer w10² from N minus F.
    error_w10_sq = eps*(zr2/gn**2+wr2/gr**2)
    width_eps = inputs["independent_width_relative_error"]
    trace_lo = (1-width_eps)*trace*max(0, w10**2-error_w10_sq)
    trace_hi = (1+width_eps)*trace*(w10**2+error_w10_sq)
    extra_u2 = u2  # a free SU4 adjoint raises physical U mass squared alone.
    wrong_class_ratio = zr2/(1.5*(u2+extra_u2)+wr2)
    assert wrong_class_ratio < 1
    return {"w4_GeV":f*math.sqrt(2*(1-rho)), "w10_GeV":w10,
            "PS_masses_GeV_U_WR_ZR":list(np.sqrt([u2, wr2, zr2])),
            "EW_poles_GeV_WR_Zlight_Zheavy":list(np.sqrt([wr_physical_sq, *poles])),
            "EW_unsubtracted_ratio":poles[1]/(1.5*u2+wr_physical_sq),
            "corrected_ratio":zr2/(1.5*u2+wr2), "normalized_total_dilepton_width":sum(width_channels),
            "neutrino_singular_GeV":list(w10*strengths),
            "sum_neutrino_mass_squared_GeV2":neutrino_trace,
            "mass_squared_sum_working_interval_GeV2":[trace_lo,trace_hi],
            "independent_SU4_adjoint_ratio":wrong_class_ratio}


def compute():
    inputs = read(HERE/"inputs.json")
    base = read(HERE/inputs["baseline"])
    low, high = inputs["matching_window_GeV"]
    delta = inputs["matching_defect_abs"]
    min_args = couplings(base, low, (-2, 2))
    max_args = couplings(base, high, (2, -2))
    lo, hi = mass_ratio(*min_args), mass_ratio(*max_args)
    table = {}
    for name, row in inputs["table2_comparison"].items():
        gy = math.sqrt(3/5) * row["g1"][0]
        table[name] = {"SM_table_central":mass_ratio(gy, row["g3"][0]),
                       "one_loop_central":mass_ratio(*couplings(base, float(name)))}
    eps = inputs["pole_squared_rel_error"]
    return {"coefficients_per_unit_vev2_U_WR_neutral":generator_test(),
            "WR_U_two_width_box":[lo,hi], "WR_U_with_signed_threshold":[
                    mass_ratio(*min_args, delta),mass_ratio(*max_args, -delta)],
            "ZR_U_two_width_box":[math.sqrt(1.5+lo*lo),math.sqrt(2*(1.5+hi*hi))],
            "ZR_WR_two_width_box":[math.sqrt(1+1.5/hi**2),math.sqrt(2*(1+1.5/lo**2))],
            "minimum_abs_inverse_ratio_shift_for_equal_WR_U":1/hi**2-1,
            "observed_ratio_min_for_positive_tensor_at_pole_precision":(1+eps)/(1-eps),
            "higher_order_same_data":table, "synthetic":benchmark(inputs,base)}


def rounded(value):
    if isinstance(value, dict):
        return {k:rounded(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [rounded(i) for i in value]
    return round(float(value),10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = rounded(compute())
    path = HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
    elif args.check:
        assert read(path)==result, "Stored output differs from independent computation"
        print("1543 tensors, matching, bounded windows, pole subtraction and width trace verified.")
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2))
