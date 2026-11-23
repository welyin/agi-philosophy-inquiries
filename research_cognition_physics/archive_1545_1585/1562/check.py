"""Oriented flavor cofactor, complete mass-label check and held-out up spectrum."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def cofactors(spec):
    a, b, c = spec
    return np.array([b*c, a*c, a*b])


def classify(spec, target, widths):
    # Every output is bounded by the largest allowed mass. Identify the original
    # heavy direction, then use the second mass to certify the middle direction.
    d, s, b = spec
    low, high = target[:, 0]-widths*target[:, 1], target[:, 0]+widths*target[:, 1]
    cap_heavy = (high[2]+d)/(s*b)
    heavy_floor = b-cap_heavy*d*s
    assert heavy_floor > high[1]
    cap_middle = (high[1]+d)/(s*b)
    middle_floor = s-cap_middle*d*b
    assert middle_floor > high[0]
    coefficient = [max(0.,low[0]-d,d-high[0])/(s*b),(high[0]+d)/(s*b)]
    return dict(heavy_floor=float(heavy_floor), middle_floor=float(middle_floor),
                coefficient_range=coefficient)


def jacobian_checks():
    def adjoint_polynomial(m):
        return np.column_stack([np.cross(m[:,(i+1)%3], m[:,(i+2)%3]) for i in range(3)]).T
    rng = np.random.default_rng(1562)
    def unitary():
        q, r = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        return q @ np.diag(np.diag(r)/abs(np.diag(r)))
    l, r, y = unitary(), unitary(), np.diag([0., .13, .61])
    moved = l @ y @ r.conj().T
    adj = adjoint_polynomial(y.conj().T)
    weight = np.linalg.det(r)/np.linalg.det(l)
    assert np.allclose(adjoint_polynomial(moved.conj().T),weight*l@adj@r.conj().T)
    assert np.linalg.norm(adj[:,0]) > 0 and np.linalg.matrix_rank(adj)==1
    assert np.allclose(np.linalg.det(moved)*adjoint_polynomial(moved.conj().T),
                       moved@adjoint_polynomial(moved.conj().T@moved))


def calculate():
    config = read(HERE/"inputs.json")
    targets = read(HERE/config["charged"])["high_scale"]
    points = read(HERE/config["tensors"])
    n = config["widths"]
    edm = read(HERE/config["neutron"])["edm"]
    slope_floor = edm["abs_slope_e_fm"][0]-sum(edm["abs_slope_e_fm"][1:])
    theta_budget = edm["limit_90_ecm"]/(edm["fm_in_cm"]*slope_floor)
    outputs = {}
    jacobian_checks()
    for tag in config["branches"]:
        v = dict(zip(points["labels"],points["candidates"][tag]["values"]))
        up = np.array([v[x] for x in ("yu","yc","yt")])
        down = np.array([v[x] for x in ("yd","ys","yb")])
        d,s,b = down
        low, high = classify(down,np.array(targets["down"]),n), classify(up,np.array(targets["up"]),n)
        tar = targets["down"][0][0]
        weights = cofactors(down)
        # z parametrizes all complex couplings consistent with the *one* central d.
        positive = (tar-d)/weights[0]
        ranges = {}
        for i, name in enumerate(("s","b"), start=1):
            center, radius = down[i]-weights[i]*d/weights[0], weights[i]*tar/weights[0]
            ranges[name] = [center-radius,center+radius]
            assert center > radius
        out = down + positive*weights
        shift = out**2-down**2
        gap_ratio = (1/s**2-1/b**2)/(1/d**2-1/s**2)
        assert math.isclose((shift[1]-shift[2])/(shift[0]-shift[1]),gap_ratio,rel_tol=2e-6)
        cpslope = sum(d*tar/(h*h-d*d-d*tar) for h in (s,b))
        same = low["coefficient_range"]
        # Independent phases and permutations. Reverse triangle in each up entry.
        up_floor = np.maximum(up-cofactors(up)*same[1], cofactors(up)*same[0]-up)
        assert min(up_floor) > targets["up"][0][0]+n*targets["up"][0][1]
        assert np.all(abs(out-np.array(targets["down"])[:,0])<=n*np.array(targets["down"])[:,1])
        outputs[tag] = dict(
            down_permutation_certificate=low, up_permutation_certificate=high,
            cofactor_coefficient_real_positive=positive, calibrated_positive_down=out.tolist(),
            other_down_masses_any_phase=ranges, squared_shift_gap_ratio=gap_ratio,
            min_up_mass_if_equal_strength=float(min(up_floor)),
            up_entries_floors_if_equal_strength=up_floor.tolist(),
            maximum_up_to_down_coefficient=high["coefficient_range"][1]/same[0],
            radial_kappa_d_if_cubic_positive=float(3-2*d/tar),
            aligned_strong_angle_slope_budget=cpslope,
            maximum_light_phase_rad_on_small_theta_branch=theta_budget/(1-cpslope))
    return dict(scale_GeV=config["scale_GeV"],neutron_small_angle_working_budget=theta_budget,branches=outputs)


def compare(new, old):
    if isinstance(new, dict):
        assert new.keys()==old.keys()
        for key, value in new.items():
            compare(value,old[key])
    else:
        np.testing.assert_allclose(new,old,rtol=2e-10,atol=1e-15)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    dest=HERE/"results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write("\n")
    else:
        compare(result,read(dest))
        print("1562 cofactor, all-phase quark spectra and physical alternatives verified")
