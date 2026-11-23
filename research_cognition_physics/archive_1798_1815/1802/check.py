"""Two common normalized parent Yukawas, doublet fractions and Majorana scales."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
CAP = math.sqrt(4 * math.pi)


def spec(matrix):
    return np.linalg.norm(matrix, ord=2)


def complete_row(row):
    # Columns of V are orthonormal. The first row of V* spans the light projection.
    _, _, vh = np.linalg.svd(np.array([row]), full_matrices=True)
    phase = row @ vh[0].conj()
    result = vh.copy()
    result[0] *= phase / abs(phase)
    np.testing.assert_allclose(row, result[0], atol=1e-13)
    np.testing.assert_allclose(result @ result.conj().T, np.eye(3), atol=1e-13)
    return result


def point(coords, api, data, vu):
    coords = np.asarray(coords)
    up, down, _, f_eff, _, _, _, tagged = api["matrices"](coords, data)
    _, _, rh, rf, _, type_i = api["symmetric_sources"](coords, data)
    h_eff = down - f_eff
    h, f = spec(h_eff), spec(f_eff)
    ph, pf = abs(rh)**2, abs(rf)**2
    assert ph > pf > 1
    required = math.sqrt(max(h*h + f*f, ph*h*h + pf*f*f))
    min_x, min_y = (h/CAP)**2, (f/CAP)**2
    max_y = min(1-min_x, (1-ph*min_x)/pf)
    assert min_y < max_y
    # ph > pf: at the maximal total down fraction, x has its smallest value.
    spectator_min = 1-min_x-max_y
    # Check the minimax parent norm is reached with normalized up and down rows.
    x, y = (h/required)**2, (f/required)**2
    np.testing.assert_allclose(ph*x+pf*y, 1.0, atol=1e-12)
    neutral = api["mixed"](coords, data)
    light = neutral["low"]["common_light_scale_eV"]*1e-9
    zeta = complex(coords[7], coords[8])
    right_scale = vu**2*np.linalg.norm(type_i)/(2*light*abs(zeta))
    left_scale_eV = light*1e9/np.linalg.norm(f_eff)
    # Work in the down Takagi frame. The ordinary down coupling is diagonal.
    fd = tagged.conj().T @ f_eff @ tagged.conj()
    hd = tagged.conj().T @ h_eff @ tagged.conj()
    ud = tagged.conj().T @ up @ tagged.conj()
    np.testing.assert_allclose(hd[1,2], -fd[1,2], atol=1e-14)
    np.testing.assert_allclose(ud[1,2], (rf-rh)*fd[1,2], atol=1e-14)
    # Saturate the Cauchy bound even at the registered finite parent norm cap.
    x = 1 / (abs(rh)*(abs(rh)+abs(rf)))
    y = 1 / (abs(rf)*(abs(rh)+abs(rf)))
    assert x >= min_x and y >= min_y and x+y < 1
    down_row = np.sqrt([x, y, 1-x-y])
    up_row = np.array([rh*math.sqrt(x), rf*math.sqrt(y), 0])
    np.testing.assert_allclose(np.linalg.norm(up_row), 1, atol=3e-14)
    complete_row(up_row)
    partners = complete_row(down_row.astype(complex)) @ np.array(
        [hd[1,2]/math.sqrt(x), fd[1,2]/math.sqrt(y), 0])
    charged_partner_norm = (abs(rh)+abs(rf))*abs(fd[1,2])
    np.testing.assert_allclose(abs(partners[0]), 0, atol=1e-14)
    np.testing.assert_allclose(np.linalg.norm(partners[1:]), charged_partner_norm, atol=1e-13)
    return {
        "rH_abs":abs(rh), "rF_abs":abs(rf), "required_minimax_parent_norm":required,
        "down_spectator_fraction_min_at_cap":spectator_min,
        "down_spectator_fraction_infimum_without_cap":1-1/pf,
        "allowed_down_F_modulus":[math.sqrt(min_y),math.sqrt(max_y)],
        "vR_effective_GeV":[right_scale*math.sqrt(min_y),right_scale*math.sqrt(max_y)],
        "vL_effective_eV":[left_scale_eV*math.sqrt(min_y),left_scale_eV*math.sqrt(max_y)],
        "fixed_abs_vL_eV_per_vR_GeV":left_scale_eV/right_scale,
        "min_inclusive_bs_parent_Yukawa":charged_partner_norm,
        "min_inclusive_bs_squared":charged_partner_norm**2,
        "observed_up_in_down_frame_bs":abs(ud[1,2]),
        "norm_cap_at_inclusive_min":[h/math.sqrt(x),f/math.sqrt(y)],
        "max_mR_singular_GeV":right_scale*f,
    }


def bs_clock():
    earlier = HERE.parent
    matrix_inputs = json.loads((earlier/"1795/inputs.json").read_bytes())
    had = matrix_inputs["mesons"]["Bs"]["matrix_elements"]
    width = matrix_inputs["matrix_element_relative_half_window"]
    anomalous = np.array([[6., -40.], [1/6, 16+34/3]])
    gap = np.ptp(np.linalg.eigvals(anomalous))
    scalar_max = (1+width)*max(abs(had[2]), abs(had[3])/(6*gap))
    kernel_floor = (1-width)*had[1]-scalar_max
    assert kernel_floor > 0 and had[2]*had[3] > 0
    measured = json.loads((earlier/"1795/results.json").read_bytes())["observed_norm_working_caps_GeV"]["Bs"]
    top_budget = json.loads((earlier/"1799/independent_results.json").read_bytes())["physical_mesons"]["Bs"]["other_abs_budget_from_1581_inputs_GeV"]
    return {"positive_mixing_floor_GeV3":kernel_floor,"total_permitted_scalar_amplitude_GeV":measured+top_budget,
            "GeV_per_unit_flavour_norm":math.sqrt(kernel_floor/(measured+top_budget))}


def calculate():
    api = runpy.run_path(str(HERE / CFG["matching_code"]))
    points = json.loads((HERE / CFG["source_points"]).read_bytes())
    tan_beta = runpy.run_path(str(HERE / CFG["tan_beta_code"]))["TAN"]
    vu = CFG["vev_GeV"]*tan_beta/math.hypot(1, tan_beta)
    data = api["backend"]()
    # Reverse triangle bound over all unmeasured Takagi phases, for the same CKM and masses.
    terms = data[0]*abs(data[3][:,1]*data[3][:,2])
    lo = max(0., 2*max(terms)-sum(terms))
    cases = {k:point(v, api, data, vu) for group in ("mixed_type_I_II","fixed_delta_minus70_mass_readout")
             for k,v in points[group].items()}
    assert all(row["min_inclusive_bs_parent_Yukawa"] >= lo for row in cases.values())
    clock = bs_clock()
    scale = clock["GeV_per_unit_flavour_norm"]
    for row in cases.values():
        mass_floor = scale*row["min_inclusive_bs_parent_Yukawa"]
        row["Bs_lower_GeV_on_harmonic_partner_mass_over_transfer"] = mass_floor
        row["fraction_bs_norm_squared_at_or_below_1TeV_upper_if_unit_transfer"] = (1000/mass_floor)**2
    return {"tan_beta":tan_beta, "vu_GeV":vu, "assigned_parent_norm_cap":CAP,
            "all_family_rephases_interval_on_abs_up_bs":[float(lo),float(sum(terms))],
            "Bs_clock":clock,"all_embeddings_min_harmonic_mass_over_transfer_GeV":lo*scale,
            "cases":cases}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    results = calculate()
    blob = json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    path = HERE / "results.json"
    if parser.parse_args().check:
        assert json.loads(path.read_bytes()) == results
    else:
        path.write_text(blob, encoding="utf-8")
    print("1802", hashlib.sha256(blob.encode()).hexdigest()[:12])
