"""Verify one-loop charged-lepton texture, spectra and distance certificates."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def window(readout, widths):
    centre, minus, plus = readout
    return centre - widths * minus, centre + widths * plus


def edge_ratio(masses):
    squared = np.square(masses)
    return (squared[1] - squared[0]) / (squared[2] - squared[0])


def texture(r):
    return np.array([[0, r, -1], [r, 0, 1-r], [-1, 1-r, 0]])


def spectral_test(r, amplitude_upper, no_gap_upper):
    matrix = texture(r)
    masses, vec = np.linalg.eigh(matrix)
    order = np.argsort(abs(masses))
    light, middle, high = abs(masses[order])
    p0 = float(vec[0, order[0]] ** 2)
    assert np.allclose(np.sort(abs(masses)), np.linalg.svd(matrix, compute_uv=False)[::-1])
    pmax = (1-math.sqrt(1-amplitude_upper))/2
    angles = [abs(math.sqrt(p0*(1-p))-math.sqrt(p*(1-p0))) for p in (pmax, 1-pmax)]
    angle_sine = min(angles)
    io_lower = (middle*angle_sine-light)/(2*high)
    low, mid = light/high, middle/high
    no_ratio_at_io = ((mid-io_lower)**2-(low+io_lower)**2)/(1+io_lower)**2
    assert mid > io_lower and no_ratio_at_io > no_gap_upper
    return {
        "edge_ratio": float(r), "singular_values": [float(t) for t in (light,middle,high)],
        "light_flavour_weights": (vec[:, order[0]]**2).tolist(),
        "reactor_amplitude": 4*p0*(1-p0),
        "small_to_total_squared_span": float((high**2-middle**2)/(high**2-light**2)),
        "angle_sine_lower": angle_sine, "all_orderings_extra_relative_norm_lower": float(io_lower),
        "NO_solar_to_total_min_at_that_radius": float(no_ratio_at_io),
        "box_row_factor": float(np.linalg.norm(matrix, ord=np.inf)/high)
    }


def calculate(data):
    widths = data["widths"]
    upper_a = window(data["Daya_Bay_A13"], widths)[1]
    solar = window(data["KamLAND_only_dm21_eV2"], widths)
    no_gap = window(data["T2K_mu_disapp_dm32_NO_eV2"], widths)
    io_gap = window(data["T2K_mu_disapp_dm13_IO_eV2"], widths)
    solar_to_span = [solar[0]/(solar[0]+no_gap[1]), solar[1]/(solar[1]+no_gap[0])]
    io_solar_to_span = [solar[0]/(solar[0]+io_gap[1]), solar[1]/(solar[1]+io_gap[0])]
    charged = np.asarray(data["charged_yukawas_mz"])
    rref = edge_ratio(charged[:,0])
    nominal = spectral_test(rref, upper_a, solar_to_span[1])
    pole = spectral_test(edge_ratio(data["pole_comparison_MeV"]), upper_a, solar_to_span[1])
    limits = [edge_ratio(charged[:,0]+widths*charged[:,1]*np.asarray(signs))
              for signs in product((-1,1), repeat=3)]
    dr = max(abs(r-rref) for r in limits)
    # Changing r changes just two edges. Their matrix has norm sqrt(2) |dr|.
    charged_error = math.sqrt(2)*dr/nominal["singular_values"][2]
    bound = nominal["all_orderings_extra_relative_norm_lower"]
    envelope = (bound-charged_error)/(1+charged_error)
    # The loop mass difference is an average over ln(M_charged^2), bounded
    # by the Feynman-parameter integral at the lowest M and highest m_lepton.
    z = (data["charged_lepton_mass_upper_GeV"]/
         data["charged_scalar_mass_lower_for_finite_loop_GeV"])**2
    finite_single = 1-(1-z+z*math.log(z))/(1-z)**2
    edges = [(0,1),(0,2),(1,2)]
    edge_factors = [
        (1+w)/(1-w) for i,j in edges for w in
        [((charged[i,0]+widths*charged[i,1])/(charged[j,0]-widths*charged[j,1]))**2]
    ]
    edge_loop = finite_single*max(edge_factors)
    width_delta = 6*data["charged_tag_absolute_fraction_budget"]
    x_low, x_high = math.sqrt(1-width_delta), math.sqrt(1+width_delta)
    edge_shape = max(1-x_low, x_high-1)
    transfer = data["extra_relative_mass_edge_transfer"]
    # Triangle inequality allows complex matching errors after edge rephasing.
    shape_and_loop = edge_shape + x_high*edge_loop
    total_edge_error = shape_and_loop + x_high*(1+edge_loop)*transfer
    # The exact common-scale mass matrix at each r differs by charged_error;
    # its row norm is at most the central row norm plus 2|dr|.
    row_factor = nominal["box_row_factor"]+2*dr/nominal["singular_values"][2]
    total = charged_error+row_factor*total_edge_error
    assert 0 < total < bound
    # Unequal f can cancel every charged-lepton squared difference analytically.
    inverse_edges = [1/rref, 1/(1-rref), 1.]
    alternative = texture(rref)*np.array([[0,inverse_edges[0],1],
                                         [inverse_edges[0],0,inverse_edges[1]],
                                         [1,inverse_edges[1],0]])
    assert np.allclose(np.linalg.svd(alternative, compute_uv=False), [2,1,1])
    return {"common_scale":nominal,"pole_shape_comparison":pole,
            "NO_solar_to_total_data":solar_to_span,"IO_solar_to_total_data":io_solar_to_span,
            "charged_mz_r_interval":[min(limits),max(limits)],"relative_charged_input_budget":charged_error,
            "two_width_charged_box_extra_norm_over_actual_nominal":envelope,
            "finite_charged_scalar_single_line_upper":finite_single,
            "finite_charged_scalar_edge_upper":edge_loop,
            "charged_width_relative_amplitude_uncertainty":edge_shape,
            "relative_total_budget_to_central":total,
            "finite_contract_extra_norm_over_central_nominal":bound-total,
            "finite_contract_extra_norm_over_actual_loop":(bound-total)/(1+total),
            "free_edge_equal_matrix_singular_values":[2,1,1]}


def close(saved, fresh):
    if isinstance(fresh, dict):
        return saved.keys() == fresh.keys() and all(close(saved[k], fresh[k]) for k in fresh)
    return np.allclose(saved, fresh, rtol=2e-9, atol=1e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate(json.loads((ROOT/"inputs.json").read_text(encoding="utf-8")))
    out = ROOT/"results.json"
    if args.write:
        with out.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
    elif not close(json.loads(out.read_text(encoding="utf-8")), result):
        raise SystemExit("Mismatch with saved results")
    else:
        print("1537: charged-flavour reads, both gap orderings and finite corrections verified")
