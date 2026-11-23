"""Test charged, basis-invariant sources for a rank-two neutrino texture."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


cfg = read("inputs.json")
sim = runpy.run_path(str(HERE / cfg["runner"]))
state = sim["graph"].__globals__
fl = sim["fl"]
charged = read(cfg["charged_data"])["direct_semileptonic_moduli"]
spectra = read(cfg["same_spectrum"])["shared_witness"]
vacuum = read(cfg["flavor_vacuum"])
sterile = read(cfg["sterile_couplings"])
alignment = read(cfg["alignment_vectors"])
a, b = (np.array(alignment[key], float) for key in ("source_a", "source_b"))
m, g, q = (sterile[key] for key in ("m_N0H0", "coupling_to_s2", "coupling_to_x"))
x, rho = vacuum["real_X_over_cutoff"], vacuum["scalar_magnitude_over_cutoff"]
s = rho * np.exp(1j * np.pi * vacuum["scalar_angle_over_pi"])
phases = (s/rho)**np.array(alignment["su5_generation_charges_5bar"])
portal = cfg["primary_charged_endpoint"]
quarks = runpy.run_path(str(HERE / "../1708/check.py"))["quarks"]
measure = read(cfg["gaps"])
tables = read(cfg["tables"])
clocks = read(cfg["clocks"])


def strength(rule, us):
    if rule == "two_port_Cabibbo":
        return (a @ a) * us**2
    if rule == "one_port_Cabibbo":
        return us**2
    if rule == "normalized_both":
        return (a @ a) / (b @ b)
    if rule == "equal_raw":
        return 1.
    raise ValueError(rule)


def new_dirac(r):
    one = math.sqrt(r) * b
    two = s**2 * (x / rho**2 * a + g / m * one)
    shape = phases[:, None] * np.column_stack((one, two))
    unit_c = np.array([[m, 0], [g*s**2, q*x]])
    expected = phases[:, None] * np.column_stack((one/m, (s/rho)**2 * a/q))
    np.testing.assert_allclose(np.linalg.solve(unit_c, shape.T).T, expected, atol=5e-15)
    # Start the outer solver in GeV using only a measured Higgs mass moment
    # and the fixed physical neutral poles, without the old light-mass ratio.
    yn = (shape @ state["U"].conj()) / sim["vev"]
    mu = spectra["neutral_only_Higgs_mass_squared_GeV2"]
    gain = math.sqrt(8*np.pi**2*mu / np.sum(abs(yn * state["m"])**2))
    return gain * shape


def upper_cabibbo(rep, step):
    """Run up and down Yukawas to exactly the same upper neutral mass."""
    y = np.zeros((4, 3, 4), complex)
    y[0, :, :3] = np.diag(fl["central"]("up"))
    theta = fl["angles"]("direct")
    down = fl["central"]("down")
    down[0] += portal["yd_widths"] * fl["low"]["yukawa"]["down"][0][1]
    y[1, :, :3] = fl["low_ckm"](rep["gamma_MZ_deg"], theta) @ np.diag(down)
    basis = sim["zer"]["unpack"](rep["charged_left_MZ"])
    y[2, :, :3] = basis @ np.diag(fl["central"]("charged_leptons"))
    poles = spectra["neutral_poles_GeV"]
    for start, susy, dt, gauge in fl["segments"](step):
        if start >= poles[2]:
            up, _ = fl["basis"](y[0])
            dw, _ = fl["basis"](y[1])
            return abs(np.vdot(up[:, 0], dw[:, 1]))
        assert not susy
        for jj, pole in enumerate(poles):
            if start == pole:
                y[3, :, jj] = fl["physics"]["Y0"][:, jj]
        for i in range(0, len(gauge)-1, 2):
            beta = fl["beta"]
            k1 = beta(y, gauge[i], False)
            k2 = beta(y+dt*k1/2, gauge[i+1], False)
            k3 = beta(y+dt*k2/2, gauge[i+1], False)
            k4 = beta(y+dt*k3, gauge[i+2], False)
            y += dt*(k1+2*k2+2*k3+k4)/6
    raise ValueError("upper neutral threshold missing")


def case(rule, shift, equal_terminal=False, step=.2):
    us_mz = charged["Vus"][0] + shift * np.linalg.norm(charged["Vus"][1:])
    values = fl["cfg"]["direct_semileptonic_moduli"]
    values["Vus"][0] = us_mz
    values["Vub"][0] = charged["Vub"][0] + portal["Vub_widths"]*np.linalg.norm(charged["Vub"][1:])
    us = us_mz
    for _ in range(6):
        r = strength(rule, us)
        state["D0"] = new_dirac(r)
        # The absolute seed (1 keV) drops out of the common solar calibration.
        weights = np.ones(2) if equal_terminal else np.array([m*m, q*q])
        state["H0"] = 1e-6*np.diag(weights)
        graph = sim["graph"](portal["phase_k"], portal["yd_widths"], portal["third_port"], step)
        d, rep, trace, _, err, _, moment = graph
        updated = upper_cabibbo(rep, step)
        if abs(strength(rule, updated)-r) < 2e-11:
            break
        us = updated
    else:
        raise RuntimeError("charged source matching did not converge")
    light = sim["lower_tensors"](sim["pair_sources"](d).sum(axis=1), trace)
    basis = sim["zer"]["unpack"](rep["charged_left_MZ"])
    mass = np.linalg.svd(light, compute_uv=False)[::-1]
    solar = measure["juno_NO_dm21"][0]
    gain = math.sqrt(solar) / mass[1]
    obs = sim["observables"](light*gain, basis)
    fc = sim["accelerator"](obs)
    mr2 = (mass[2]/mass[1])**2
    qobs = quarks(rep["gamma_MZ_deg"], fl["angles"]("direct"))

    def interval(row):
        return np.array([row[0]-row[1], row[0]+row[2]])

    def inside(v, w):
        return bool(w[0] <= v <= w[1])

    allowed = np.stack([interval(measure["juno_NO_dm21"]), interval(measure["juno_NO_dm31"])/mr2,
                        interval(measure["accelerator_NO_dm32"])/(mr2-1)])
    lo, hi = max(allowed[:, 0]), min(allowed[:, 1])
    tests = {"gamma": inside(qobs["gamma_deg"], tables["time_integrated_gamma_deg"]["68"]),
             "golden": inside(qobs["golden_S"], clocks["one_width_golden_S_range"]),
             "clocks": inside(qobs["Vtd_over_Vts"], clocks["one_width_clock_ratio_range"]),
             "lepton_angles": all(obs["windows_pass"]),
             "joint_FC90": fc["conf90_0"] or fc["conf90_1"],
             "all_frequencies_fixed_solar": bool(lo <= solar <= hi),
             "all_frequencies_for_some_allowed_solar": bool(lo <= hi)}
    r_eff = r * (weights[0]/weights[1]) * (q/m)**2
    relation = (4+121*r_eff**2-16*r_eff)/(6*r_eff)
    tree_ratio = (relation+math.sqrt(relation**2-4))/2
    return {"Vus_MZ": us_mz, "Vus_at_upper_neutral_cluster": updated, "Dirac_ratio_squared": r,
            "terminal_H6_over_H0": weights[1]/weights[0], "tree_rank_one_ratio": r_eff,
            "mass_ratio_tree": tree_ratio, "mass_ratio_after_transport": math.sqrt(mr2),
            "angles": obs["reactor_solar_atm"], "J": obs["J"], "delta_deg": np.degrees(fc["delta_rad"]),
            "accel": fc, "quarks": qobs, "delta_m31_squared_eV2": mr2*solar,
            "delta_m32_squared_eV2": (mr2-1)*solar, "allowed_solar_eV2": [lo, hi],
            "common_H0_H6_keV": (weights*gain).tolist(),
            "weak_weight": moment["heavy_weak_weight_leading"],
            "Higgs_moment_GeV2": moment["pole_Yukawa_Higgs_moment_GeV2"],
            "common_boundary_residual": err, "tests": tests}


def terminal_exchange(us):
    # H and Z' have identical charges. Z is their conjugate. M couples Z to
    # a normalized heavy combination. A single LNV source couples both Z' legs.
    mixing = .001
    heavy = 1e7
    B = 1.
    weights = np.array([m, q])
    p = mixing*weights
    c = np.sqrt(1-p*p)
    result = []
    for proj, orth in zip(p, c):
        # Original order is Z,H,Z'. Columns are Z, retained H, heavy Y.
        rot = np.array([[1., 0., 0.], [0., orth, proj], [0., -proj, orth]])
        old = np.array([[0., heavy*proj, heavy*orth],
                        [heavy*proj, 0., 0.], [heavy*orth, 0., B]])
        mass = rot.T @ old @ rot
        indices = [0, 2]
        reduced = mass[1, 1] - mass[1, indices] @ np.linalg.solve(
            mass[np.ix_(indices, indices)], mass[indices, 1])
        np.testing.assert_allclose(reduced, B*proj**2, rtol=1e-11)
        result.append(reduced)
    u = np.array([1., 0., 0.])
    d = np.array([us, math.sqrt(1-us*us), 0.])
    P, T = np.outer(u, u), np.outer(d, d)
    for ell in range(1, 5):
        np.testing.assert_allclose(np.trace(np.linalg.matrix_power(P@T, ell)), us**(2*ell))
    return {"common_heavy_neutral_GeV": heavy, "source_B_GeV": B, "transfer_amplitudes": p.tolist(),
            "induced_H_GeV": result, "H6_over_H0": result[1]/result[0],
            "closed_trace_word_lengths": [2,4,6,8]}


def main():
    results = {rule: case(rule, 0) for rule in cfg["shape_candidates"]}
    cases = {str(n): case("two_port_Cabibbo", n) for n in cfg["Vus_marginal_widths"] if n}
    other = case("two_port_Cabibbo", 0, equal_terminal=True)
    fine = case("two_port_Cabibbo", 0, step=.1)
    base = results["two_port_Cabibbo"]
    obs = ["Dirac_ratio_squared", "delta_m31_squared_eV2", "J"]
    convergence = {key:abs(fine[key]-base[key]) for key in obs}
    # One common attenuation and explicit missing-port Gram complete the two
    # success columns to an isometry. Total weights on these ports need not be equal.
    us = base["Vus_at_upper_neutral_cluster"]
    z = np.column_stack((a/np.linalg.norm(a), us*b))
    scale = 1/np.linalg.norm(z, 2)
    eig, basis = np.linalg.eigh(np.eye(2)-(scale*z).T @ (scale*z))
    assert min(eig) > -1e-14
    fail = np.diag(np.sqrt(np.maximum(eig, 0))) @ basis.T
    ok = np.vstack((scale*z, fail))
    np.testing.assert_allclose(ok.T@ok, np.eye(2), atol=1e-14)
    assert max(convergence.values()) < 2e-8
    return {"same_1708_point": results, "alternative_terminal_at_same_point": other,
            "additional_Vus_marginal_endpoints": cases, "half_step_differences": convergence,
            "terminal_and_charged_port_checks": terminal_exchange(us),
            "successful_ports_common_attenuation": scale, "isometry_error": np.linalg.norm(ok.T@ok-np.eye(2)),
            "cluster_scales_GeV": state["m"].tolist()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    output = sim["tidy"](main())
    file = HERE / "results.json"
    if parser.parse_args().check:
        assert json.loads(file.read_bytes()) == output
        print("1709 verified")
    else:
        file.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
