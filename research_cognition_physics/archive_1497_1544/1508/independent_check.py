"""Neutral-atom Abelian norm, MICROSCOPE to Washington transfer, and a rotated kinetic basis."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def eta(source, f_a, f_b):
    return -source*(f_a-f_b)/(1-source*(f_a+f_b)/2)


def inner(a, b, kinetic):
    return float(a @ np.linalg.solve(kinetic, b))


def calculate():
    # Neutron charge per atomic mass unit. Orbit uses the actual alloys, not pure Ti/Pt.
    orbital = [.54043, .59668]
    ground = [5/9.01218, .54146]
    sat_positive = (-1.5+2*(2.3+1.5))*1e-15
    lower_lab = (.3-2*1.8)*1e-13
    upper_lab = (.3+2*1.8)*1e-13
    delta_orbit = orbital[1]-orbital[0]
    mean_orbit = sum(orbital)/2
    sat_limit = sat_positive/(delta_orbit+sat_positive*mean_orbit)
    y_ground = eta(sat_limit, *ground)
    # A separate illustrative transfer envelope covers composition and terrestrial source changes.
    net_limit = sat_limit*1.1
    conservative = net_limit*.015/(1-net_limit*.56)
    # Invert a signed Be-minus-Ti result, keeping the denominator of the Eotvos ratio.
    lab_only = -lower_lab/(ground[0]-ground[1]-lower_lab*sum(ground)/2)
    assert lower_lab < -conservative < y_ground < 0 < upper_lab
    assert math.isclose(eta(sat_limit, *orbital), sat_positive, rel_tol=1e-12)
    assert math.isclose(eta(lab_only, *ground), lower_lab, rel_tol=1e-12)

    # q_i are covectors of two currents; the kinetic matrix is positive definite.
    k = np.array([[1.9, -.47],[-.47,1.2]])
    u = np.array([.88,-.11])
    w = np.array([-.3,.47])
    rotation = np.array([[1.2, .8],[-.17,1.4]])
    def project(kinetic, em, neutron):
        n2 = inner(neutron, neutron, kinetic)
        p2 = n2 - inner(em, neutron, kinetic)**2 / inner(em, em, kinetic)
        assert 0 <= p2 <= n2
        return np.array([n2, p2])
    initial = project(k, u, w)
    transformed = project(rotation.T@k@rotation, rotation.T@u, rotation.T@w)
    assert np.allclose(initial, transformed, rtol=1e-14)
    for z, neutrons in [(4,5),(22,26),(78,117)]:
        up, down, electron = 2*u/3+w/3, -u/3+w/3, -u-w
        atom = z*(2*up+down+electron)+neutrons*(up+2*down)
        assert np.allclose(atom, neutrons*w, rtol=1e-14)
    gn_const = 4*math.pi*6.7088e-39*.9314941**2
    earth_min = .48
    return {"N_over_mu_orbit_TiPt": orbital, "N_over_mu_ground_BeTi": ground,
            "ground_per_orbital_leading_ratio": -(ground[0]-ground[1])/delta_orbit,
            "orbital_sign_allowed_eta_upper_two_width": sat_positive,
            "orbital_source_factor_max": sat_limit,
            "ground_eta_if_orbital_saturates": y_ground,
            "ground_abs_eta_max_broad_transfer": conservative,
            "ground_measurement_box": [lower_lab, upper_lab],
            "ground_alone_source_factor_max": lab_only,
            "massless_n_gauge_norm_max_at_fEarth_min_0_48": math.sqrt(sat_limit*gn_const/earth_min),
            "covariant_metric_example_neutron_norm2_perp_norm2": initial.tolist(),
            "inputs_scope": "Same unscreened massless repulsive current and a terrestrial monopole; extra ground source budget 10%, material difference <= .015, mean <= .56. Natural Be and Ti for central transfer; errors are working boxes.",
            "sources": ["https://arxiv.org/abs/2209.15487","https://doi.org/10.1103/PhysRevD.97.055039",
                        "https://arxiv.org/abs/0712.0607"]}


def agrees(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(agrees(v,b[k]) for k,v in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(agrees(x,y) for x,y in zip(a,b))
    if isinstance(a, float):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-42)
    return a == b


if __name__ == "__main__":
    obj = calculate()
    if sys.argv[1:] == ["--write-new"]:
        with OUT.open("x", encoding="utf-8") as handle:
            json.dump(obj, handle, indent=2)
            handle.write("\n")
    elif sys.argv[1:] in ([],["--check"]):
        assert agrees(obj,json.loads(OUT.read_text(encoding="utf-8")))
        print("1508 independent: charge basis, signed terrestrial transfer and terrestrial-only bound verified.")
    else:
        raise SystemExit("Use --write-new or --check.")
