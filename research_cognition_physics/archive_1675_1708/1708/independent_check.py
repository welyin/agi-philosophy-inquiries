"""Hold the registered 9:4 terminal ratio fixed and calibrate one neutrino frequency."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def load_model(name):
    return runpy.run_path(str(HERE / name))


def interval(measurement):
    centre, lower, upper = measurement
    return np.array([centre - lower, centre + upper])


def fixed_shape(sim, bdata, reactor, golden, clock, vertices, dz, k, role):
    d, rep, trace, _, _, _, common = sim["graph"](k, dz, role, .20)
    basis = sim["zer"]["unpack"](rep["charged_left_MZ"])
    # Sum the probes at the two physical poles before evolving. This uses no fitted H ratio.
    pair = sim["pair_sources"](d).sum(axis=1)
    matrix = basis.T @ sim["lower_tensors"](pair, trace) @ basis
    _, masses, vh = np.linalg.svd(matrix)
    U = vh.conj().T[:, ::-1]
    light = masses[::-1]
    assert light[0] < 1e-12 * light[2]
    probabilities = abs(U)**2
    s13 = probabilities[0, 2]
    s12, s23 = probabilities[0, 1]/(1-s13), probabilities[1, 2]/(1-s13)
    angles = np.array([4*s13*(1-s13), s12, s23])
    tau = math.sqrt(s12*(1-s12)*s23*(1-s23)*s13)
    sin_delta = (U[0, 0]*U[1, 1]*(U[0, 1]*U[1, 0]).conj()).imag/((1-s13)*tau)
    cos_delta = (probabilities[1, 0]-s12*(1-s23)-(1-s12)*s23*s13)/(2*tau)
    np.testing.assert_allclose(sin_delta**2 + cos_delta**2, 1., atol=1e-12)
    delta = math.atan2(sin_delta, cos_delta)
    solar = reactor["juno_NO_dm21"][0]
    # Positive common gain scales mass linearly. s is its square times the original solar gap.
    mass_ratio_sq = (light[2]/light[1])**2
    windows = np.array([interval(reactor["juno_NO_dm21"]),
                        interval(reactor["juno_NO_dm31"])/mass_ratio_sq,
                        interval(reactor["accelerator_NO_dm32"])/(mass_ratio_sq-1)])
    left, right = max(windows[:, 0]), min(windows[:, 1])
    g = rep["gamma_MZ_deg"]
    fl = sim["fl"]
    V = fl["low_ckm"](g, fl["angles"]("direct"))
    beta = np.angle(-V[1, 0]*V[1, 2].conj()/(V[2, 0]*V[2, 2].conj()))
    phi_s = -2*np.angle(-V[2, 1]*V[2, 2].conj()/(V[1, 1]*V[1, 2].conj()))
    ratio = abs(V[2, 0]/V[2, 1])
    S = math.sin(2*beta)
    visible = bdata["from_tree"](abs(V[0, 1]), abs(V[0, 2]), abs(V[1, 2]), g)
    np.testing.assert_allclose([S, ratio], visible, rtol=1e-13)
    poly = read("../1705/independent_results.json")["official_coordinates_rad"]
    f68 = vertices(poly["conf68_0"], s23, delta)
    f90 = any(vertices(poly[name], s23, delta) for name in ("conf90_0", "conf90_1"))
    checks = {"gamma_time_integrated": 60.1 <= g <= 64.8, "golden_Bd": golden[0] <= S <= golden[1],
              "clock_ratio_4sea": clock[0] <= ratio <= clock[1],
              "three_angles": bool(((sim["ANG"][:, 0] <= angles) & (angles <= sim["ANG"][:, 1])).all()),
              "FC_joint_90": f90, "mass_frequencies_at_solar_centre": left <= solar <= right,
              "phi_s_LHCb_Run2": abs(phi_s+.039) <= math.hypot(.022, .006)}
    return dict(gamma_deg=g, golden_Bd=S, ratio_B=ratio, phi_s_rad=phi_s,
                cms_combined_phi_s_marginal_widths=abs(phi_s+.075)/.023,
                solar_input_eV2=solar, solar_scale_common_interval_eV2=[left,right],
                juno_atmospheric_eV2=solar*mass_ratio_sq, beam_atmospheric_eV2=solar*(mass_ratio_sq-1),
                no_tuning_two_source_ratio=2.25, common_terminal_multiplier=math.sqrt(solar)/light[1],
                angles=angles.tolist(), delta_deg=math.degrees(delta), FC68=f68,
                T2K2023_1D90=(-3.01 <= delta <= -.52), checks=checks, joint=all(checks.values()),
                matched_Higgs_moment_GeV2=common["pole_Yukawa_Higgs_moment_GeV2"])


def compute():
    sim = load_model("../1706/check.py")
    bdata = load_model("../1707/independent_check.py")
    vertices = load_model("../1705/independent_check.py")["interior"]
    reactor = read("../1668/inputs.json")
    rad = bdata["DATA"]["golden_S"]
    golden = [rad[0]+sgn*math.hypot(*rad[1:]) for sgn in (-1,1)]
    clock = [bdata["ratio_at"](t)[0] for t in (-1,1)]
    expected = read("../1689/inputs.json")["bare_H0_H6_masses_in_sterile_units"]
    H = sim["H0"]
    np.testing.assert_allclose(H[1,1]/H[0,0], expected[1]/expected[0])
    assert expected == [4, 9]
    vub = sim["fl"]["cfg"]["direct_semileptonic_moduli"]["Vub"]
    vub[0] -= np.linalg.norm(vub[1:])
    cases = {f"d{dz}_k{k}_{role}": fixed_shape(sim, bdata, reactor, golden, clock, vertices, dz, k, role)
             for dz in (0,-1) for k in (1,2) for role in ("unit_third","shared_third")}
    passing = [key for key, result in cases.items() if result["joint"]]
    assert passing == ["d-1_k1_unit_third"], passing
    return {"fixed_Vub":vub[0],"Bd_window":golden,"Bd_Bs_clock_interval":clock,
            "passes_declared_joint_windows":passing,"cases":cases}


def tidy(obj):
    if isinstance(obj, dict):
        return {key: tidy(item) for key,item in obj.items()}
    if isinstance(obj, list):
        return [tidy(item) for item in obj]
    if isinstance(obj, (bool,np.bool_)):
        return bool(obj)
    if isinstance(obj, (float,np.floating)):
        return float(f"{obj:.11g}")
    return obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    value = tidy(compute())
    if parser.parse_args().check:
        assert value == read("independent_results.json")
        print("Verified fixed two-source shape, joint quark and lepton observables, one solar calibration.")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as out:
            out.write(json.dumps(value,indent=2)+"\n")
