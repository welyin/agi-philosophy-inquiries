"""Exact weak-current minors in a five-state neutral-lepton sector."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def mass_and_angular_windows(cfg, box):
    xlo, xhi = box["sin2_theta12"]
    ylo, yhi = box["sin2_theta13"]
    rlo, rhi = box["solar_over_abs31"]
    dlo, dhi = box["abs31_eV2"]
    windows = {
        "normal": ((1-xhi)*(1-yhi), (1-xlo)*(1-ylo),
                   dlo*math.sqrt(rlo), dhi*math.sqrt(rhi)),
        "inverted": (ylo, yhi, dlo*math.sqrt(1+rlo), dhi*math.sqrt(1+rhi)),
    }
    floor_eV = cfg["heavy_mass_floor_GeV_benchmark"]*1e9
    return {
        order: dict(mutau_minor_share=[lo, hi], any_one_heavy_e_fraction_upper=1-lo,
                    massive_light_product_eV2=[mlo, mhi],
                    sum_weak_current_minors_upper_at_mass_floor=mhi/floor_eV**2)
        for order, (lo, hi, mlo, mhi) in windows.items()
    }


def direct_mass_matrices():
    """Diagonalize symmetric mass matrices, independently of any mixing ansatz."""
    rng = np.random.default_rng(1629)
    largest_error = 0.0
    determinant_factors = []
    for _ in range(40):
        active_singlet = rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2))
        sterile_seed = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
        sterile = 6*np.eye(2) + sterile_seed + sterile_seed.T
        mass = np.block([[np.zeros((3,3)), active_singlet], [active_singlet.T, sterile]])
        left, values, _ = np.linalg.svd(mass)
        # M = V* m V^dagger, so conjugates of the left singular
        # vectors give the weak-current columns up to column phases.
        ids = np.argsort(values)
        values = values[ids]
        active_rows = left[:3, ids].conj()
        n, heavy = active_rows[:, 0], active_rows[:, 3:]
        defect = heavy @ heavy.conj().T
        minors = []
        for i, j in ((1, 2), (0, 2), (0, 1)):
            minors.append(float(np.linalg.det(defect[np.ix_([i,j],[i,j])]).real))
        total = sum(minors)
        norm = float(np.linalg.det(np.eye(3)-defect).real)
        err = max(float(np.max(abs(np.array(minors)/total-abs(n)**2))),
                  abs(values[3]*values[4]*total/(values[1]*values[2]*norm)-1),
                  abs(np.linalg.norm(n)-1), float(values[0]/values[-1]))
        largest_error = max(largest_error, err)
        determinant_factors.append(norm)
    assert largest_error < 1e-9
    assert min(determinant_factors) < 0.8
    return {"trials":40, "all_residuals_below":1e-9, "includes_det_NN_dagger_below":0.8}


def heavy_split_witness(cfg):
    spec = cfg["majorana_phase_witness"]
    h = np.asarray(spec["h"])
    root = np.diag(np.sqrt((1-h)/h))
    angle = spec["active_plane_rotation_rad"]
    q = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
    pair = spec["massive_light_masses_eV"]
    outputs = []
    for phase in spec["phases_rad"]:
        a = q @ np.diag([pair[0], pair[1]*np.exp(1j*phase)]) @ q.T
        heavy = np.sort(np.linalg.svd(root @ a @ root, compute_uv=False))
        projected = math.prod(pair)*np.prod(1-h)/np.prod(h)
        assert abs(np.prod(heavy)/projected-1) < 1e-12
        outputs.append([round(float(v)/1e6,6) for v in heavy])
    assert outputs[0][0] != outputs[1][0]
    return {"phase_pairs_heavy_masses_MeV":outputs, "common_product_MeV2":round(float(projected)/1e12,6)}


def extra_pair(cfg):
    s = cfg["extra_Dirac_pair_witness"]
    h = [s["d"]**2/(s["d"]**2+s["M"]**2)]
    for d, majorana in s["known_two_seesaw_blocks"]:
        h.append((1-majorana/math.sqrt(majorana**2+4*d*d))/2)
    minors = [h[1]*h[2], h[0]*h[2], h[0]*h[1]]
    return {"heavy_current_diagonal":h, "mutau_minor_share":minors[0]/sum(minors),
            "massless_polar_column_e_fraction":1.0}


def run():
    cfg = json.loads((HERE/"inputs.json").read_bytes())
    box = json.loads((HERE/cfg["oscillation_box"]).read_bytes())["box"]
    return {"windows":mass_and_angular_windows(cfg, box), "direct_Takagi_tests":direct_mass_matrices(),
            "undetermined_heavy_splitting":heavy_split_witness(cfg),
            "four_singlet_counterexample":extra_pair(cfg)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE/"results.json"
    if args.check:
        assert result == json.loads(path.read_bytes())
    else:
        path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("1629 PASS: exact weak-current minors, heavy product, split witnesses and four-singlet counterexample")
