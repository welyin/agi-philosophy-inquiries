"""Independent Jensen moment bound and neutrino-mass cost outside the muon port."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def jensen(b, c, min_mass, cutoff):
    # At fixed total weight the single geometric mass saturates both Jensen bounds.
    assert 0 < min_mass < cutoff and b > 0 and c > 0
    mass = max(min_mass, min((b / c)**0.25, cutoff / math.sqrt(math.e)))
    weight = min(b / mass**2, c * mass**2)
    return mass, weight, weight * math.log(cutoff / mass)


def central_block(point):
    h = np.diag(point["H_diag"])
    f = np.zeros((3, 3), complex)
    upper = np.triu_indices(3)
    f[upper] = np.array(point["F_real_upper"]) + 1j*np.array(point["F_imag_upper"])
    f[upper[::-1]] = f[upper]
    w = h - 3 * complex(point["s_real"], point["s_imag"]) * f
    k = point["t_eV"] * w @ np.linalg.solve(f, w.T)
    _, _, vh = np.linalg.svd(h - 3*f)
    # Right singular vectors are the physical left-doublet axes in this convention.
    charged = vh.conj().T[:, ::-1]
    block = charged[:, [0, 2]].T @ k @ charged[:, [0, 2]]
    return float(np.sum(np.linalg.svd(block, compute_uv=False)))


def two_widths(row):
    return row[0] - 2*row[1], row[0] + 2*row[-1]


def universal_block(data, ordering):
    cfg = data[ordering]
    sa, sb = two_widths(cfg["solar_sine_squared"])
    qlo, qhi = two_widths(cfg["atmospheric_sine_squared"])
    amp_lo, amp_hi = two_widths(data["DayaBay_amplitude"])
    assert 0 < amp_lo < amp_hi < 1
    s13_hi = (1-math.sqrt(1-amp_hi)) / 2  # physical low branch
    # Bound the muon fraction of the lightest mass state for any Dirac phase.
    if ordering == "NO":
        project = max(0, math.sqrt(sa*(1-qhi))-math.sqrt((1-sa)*qhi*s13_hi))**2
    else:
        project = qlo*(1-s13_hi)
    solar_min, solar_max = two_widths(cfg["solar_gap"])
    atm_min, atm_max = two_widths(cfg["atmospheric_gap"])
    # All physical lightest masses greater than cut already give norm >= cut.
    if ordering == "NO":
        min_pair = math.sqrt(solar_min*(solar_min + atm_min))
    else:
        min_pair = math.sqrt(atm_min*(atm_min + solar_min))
    def determinant_floor(cut):
        largest = math.sqrt(cut**2 + solar_max + atm_max)
        return min_pair * project - cut * largest * (1-project)
    lo, hi = 0., .1
    for _ in range(90):
        cut = (lo+hi)/2
        if determinant_floor(cut) >= (cut/2)**2:
            lo = cut
        else:
            hi = cut
    return project, lo


def compute():
    cfg = read(HERE/"inputs.json")
    u, cutoff = cfg["v_over_sqrt2_GeV"], cfg["cutoff_GeV"]
    mh, err = cfg["higgs_mass_GeV"]
    budget = 8*math.pi**2 * (mh+2*err)**2 / cfg["higgs_transfer_lower"]
    output = {"moments": {}, "oscillation_projection": {}, "actual_tensors": {}}
    for key in ("mu_mu", "trace"):
        c = 2*cfg["eta_"+key+"_upper"] / u**2
        m, weight, integral = jensen(budget, c, cfg["weak_scale_GeV"], cutoff)
        a = math.log(cutoff/m)
        # Fraction of total optimal weight outside [M*/2,2M*] at 99% saturation.
        def gap(x):
            return a*math.cosh(2*x)-.5*math.sinh(2*x) - (a-x)
        tail = .01*a / min(gap(math.log(2)),gap(-math.log(2)))
        assert abs(weight*m**2-budget) < 1e-8*budget
        output["moments"][key] = dict(mass_GeV=m, weight=weight, integral=integral,
                                       tail_in_optimal_weight_at_99pct=tail)
    osc = read(HERE/"../1386/inputs.json")
    for ordering in ("NO", "IO"):
        p, bound = universal_block(osc, ordering)
        output["oscillation_projection"][ordering] = dict(
            mu_lightest_probability_lower=p, transverse_mass_norm_lower_eV=bound,
            all_heavy_mass_floor_must_be_below_GeV=(budget*u**2/(bound*1e-9))**(1/3))
    candidates = read(HERE/cfg["quark_lepton_points"])["candidates"]
    for key in ("release_down_A", "release_down_B"):
        transverse = central_block(candidates[key]["tensors"])
        output["actual_tensors"][key] = dict(
            transverse_mass_norm_eV=transverse,
            heavy_mass_floor_ceiling_GeV=(budget*u**2/(transverse*1e-9))**(1/3))
    primary = read(HERE/"results.json")
    for key in output["moments"]:
        assert abs(primary["bounds"][key]["upper"]["log_moment_upper"] -
                   output["moments"][key]["integral"]) < 1e-12
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = compute()
    target = HERE/"independent_results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    else:
        previous = read(target)
        for group, entries in result.items():
            assert previous[group].keys() == entries.keys()
            for key, row in entries.items():
                for name, value in row.items():
                    assert math.isclose(previous[group][key][name], value, rel_tol=1e-10, abs_tol=1e-13)
        print("1559 independent moments and transverse mass costs verified")
