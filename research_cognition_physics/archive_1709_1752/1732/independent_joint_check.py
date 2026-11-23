"""Rebuild the lower-Vus witness and check it against the frozen joint result."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
POINT = "registered_with_lower_us"
# Explicit witness at the precision published by R1732. None is optimized here.
OFFSETS = {"us": -.95, "ub": -.93, "cb": .93,
           "strange": -.66128807485, "bottom": .72091507347}


def read(path):
    return json.loads((HERE / path).read_bytes())


def inside(value, bounds):
    return bool(bounds[0] <= value <= bounds[1])


def winding(point, ring):
    coords = np.asarray(ring)
    vertices = coords[:, 0] + 1j*coords[:, 1] - point
    # Distance to the nearest edge in native (sin²θ23, rad) coordinates.
    edges = np.roll(vertices, -1)-vertices
    nonzero = abs(edges) > 0
    line = edges[nonzero]
    origin = vertices[nonzero]
    blend = np.clip(-(origin*line.conjugate()).real/abs(line)**2, 0, 1)
    assert min(abs(origin+blend*line)) > 1e-6
    turns = np.sum(np.angle(np.roll(vertices, -1)/vertices))/(2*math.pi)
    assert abs(turns-round(turns)) < 1e-10
    return bool(round(turns))


def check():
    config = read("inputs.json")
    saved = read("results.json")
    charged = saved["exploratory_k3"][POINT]
    neutral = saved["same_point_neutral_fixed_gain_solar_only"][POINT]
    for key in ("us", "ub", "cb"):
        assert OFFSETS[key] == config["exploratory_weak_offsets"][POINT][key]

    engine = runpy.run_path(str(HERE / config["neutral_transport"]))
    point = engine["endpoint"]("explicit_open", config["log_step"]/2, explicit_offsets=OFFSETS)
    row = point["row"]
    spectrum = np.array(read("../1386/inputs.json")["yukawa"]["down"])
    down = np.array(row["MZ_down_yukawas"])
    pulls = (down-spectrum[:,0])/spectrum[:,1]
    np.testing.assert_allclose(down, charged["down_strange_bottom_MZ"], rtol=3e-7, atol=0)
    for reference in (charged["down_strange_bottom_marginal_widths"],neutral["charged_mass_widths"]):
        np.testing.assert_allclose(pulls, reference, atol=1e-5, rtol=0)
    assert max(abs(pulls)) < 1 and abs(row["return_relative"]) < 1e-9
    assert charged["k"] == 3
    assert abs(row["x"] - row["Vud_10PeV"]**3) < 1e-8
    np.testing.assert_allclose([row["x"],row["Vud_10PeV"]],
                               [charged["positive_x"],charged["Vud_10PeV"]],atol=1e-8,rtol=0)
    high_ratio = math.prod(row["high_e_mu"])*row["tau_10PeV"]/math.prod(row["high_down_yukawas"])
    np.testing.assert_allclose(high_ratio,charged["high_mass_product_ratio"],atol=1e-8,rtol=0)
    gamma = row["gamma_MZ_deg"]
    assert abs(gamma-charged["gamma_MZ_deg"]) < 1e-5
    gamma_windows = read("../1707/inputs.json")["time_integrated_gamma_deg"]
    assert inside(gamma,gamma_windows["68"]) == charged["gamma_in_time_integrated_68"] == neutral["charged_gamma_in_68"] == True
    assert inside(gamma,gamma_windows["95"]) == charged["gamma_in_time_integrated_95"] == True

    roots = read("../1718/results.json")["original_basis_verification"]["direct_forward_both_roots"]
    gain = min(r["raw_real_b_coefficient"] for r in roots)
    assert gain == neutral["frozen_raw_gain"]
    assert engine["CFG"]["up_monomials"] == neutral["up_monomials"] == [1,2,0]
    mass = point["C"] + gain*np.outer(point["b"],point["b"])
    np.testing.assert_allclose(mass,mass.T,atol=1e-12)
    squared, u = np.linalg.eigh(mass.conj().T @ mass)
    assert abs(squared[0]) < squared[2]*1e-13
    inputs = read("../1668/inputs.json")
    solar_frequency = inputs["juno_NO_dm21"][0]
    gain_eV = math.sqrt(solar_frequency/(squared[1]-squared[0]))
    masses = np.sqrt(np.maximum(0.,squared))*gain_eV
    masses[0] = 0.   # common, independently verified exact zero line
    phases = np.angle(np.diag(u.T @ mass @ u))
    u *= np.exp(-.5j*phases)
    np.testing.assert_allclose(u.T @ mass @ u*gain_eV,np.diag(masses),atol=1e-13)
    prob = abs(u)**2
    z = prob[0,2]
    solar = prob[0,1]/(1-z)
    atm = prob[1,2]/(1-z)
    angles = [4*z*(1-z),solar,atm]
    np.testing.assert_allclose(angles,neutral["neutral_angles"],atol=2e-7,rtol=0)
    windows = read("../1659/results.json")["joint_windows"]
    assert [inside(x,box) for x,box in zip(angles,windows)] == neutral["angles_in_box"] == [True]*3
    sc, se = read("../1485/inputs.json")["solar_sin2"]
    solar_pull = (solar-sc)/se
    assert 1 < solar_pull < 2
    splittings = (squared[2]-squared[:2])*gain_eV**2
    np.testing.assert_allclose(splittings,neutral["dm31_and_32_eV2"],atol=2e-9,rtol=0)
    frequency_status = []
    for x,key in zip(splittings,("juno_NO_dm31","accelerator_NO_dm32")):
        center,lo,hi = inputs[key]
        frequency_status.append(inside(x,(center-lo,center+hi)))
    assert frequency_status == neutral["dm31_and_32_in_box"] == [True,True]
    j = (u[0,0]*u[1,1]*(u[0,1]*u[1,0]).conjugate()).imag
    delta = math.atan2(j/(1-z),.5*(prob[1,0]-solar*(1-atm)-(1-solar)*atm*z))
    assert abs(math.degrees(delta)-neutral["delta_deg"]) < 1e-4
    curves = read("../1705/independent_results.json")["official_coordinates_rad"]
    tags = [any(winding(atm+1j*delta,ring) for key,ring in curves.items()
                if key.startswith(f"conf{level}_")) for level in (68,90)]
    assert tags == neutral["SK_T2K_FC68_90"] == [False,True]
    mbb = 1000*abs(np.sum(masses*u[0]**2))
    mbeta = 1000*math.sqrt(prob[0]@(masses**2))
    assert abs(mbb-1000*abs(mass[0,0])*gain_eV) < 1e-10
    np.testing.assert_allclose([mbb,mbeta],[neutral["mbb_meV"],neutral["mbeta_meV"]],
                               atol=2e-5,rtol=0)
    print(f"1732 {POINT} verified: d,s,b pulls {pulls.round(6).tolist()}; gamma {gamma:.6f} deg")
    print(f"dm31, dm32 = {splittings[0]:.10g}, {splittings[1]:.10g} eV^2; "
          f"three angles {np.round(angles,8).tolist()}; solar pull {solar_pull:.5f}")
    print(f"delta {math.degrees(delta):.6f} deg, FC 68/90 {tags}; "
          f"m_bb {mbb:.6f} meV, m_beta {mbeta:.6f} meV")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", required=True)
    parser.parse_args()
    check()
