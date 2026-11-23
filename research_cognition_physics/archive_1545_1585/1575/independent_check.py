"""Reconstruct matching directly, then test the Higgs quadratic budget."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def run():
    cfg = read("inputs.json")
    raw = read("../1572/inputs.json")
    g = read(cfg["gauge"])
    osc = read(cfg["oscillation"])
    transport = read("../1492/inputs.json")["positive_Higgs_threshold_transport"][0]
    hmass, hwidth = read("../1372/inputs.json")["atlas"]["mass_GeV"]
    ceiling = (hmass + 2*hwidth)**2
    v = cfg["vEW_GeV"]
    mz = g["mu_ref_GeV"]
    bt, bo, bq = (np.array(raw["adjoint_beta"][k]) for k in ("weak3", "color8", "dirac32"))
    # 2pi*A = const + (b+bt+bo)*L - bo*o - bt*t + bq*vq,
    # t=log(MT/MZ), o=log(MO/MZ), vq=log(MX/MQ).
    matrix = np.column_stack([np.ones(3), raw["SM_b"] + bt + bo, -bo])
    slope_t = np.linalg.solve(matrix, bt)[2] - 1
    slope_vq = np.linalg.solve(matrix, -bq)[2]
    assert min(slope_t, slope_vq) > 0
    table = {}
    for offset in cfg["inverse_alpha_offsets"]:
        samples = []
        for signs in itertools.product((-1,1), repeat=6):
            alphainv = []
            for i, label in enumerate(("g1_GUT","g2","g3")):
                value, error = g[label]
                alphainv.append(4*math.pi/(value+signs[i]*2*error)**2 + signs[i+3]*offset)
            _, L, o = np.linalg.solve(matrix, 2*math.pi*np.array(alphainv))
            samples.append((o,L))
        o,L = min(samples)
        assert 0 < o < L
        ratio = (2*math.exp(o)-3)/5
        ref = read("results.json")["global_matching_boxes"][str(offset)]["mass"]["S_over_T_lower"]
        assert math.isclose(ratio, ref, rel_tol=1e-9)
        rows = {}
        for name, floor in [("NO",math.sqrt(osc["NO"]["solar"][1])),
                            ("IO",math.sqrt(osc["IO"]["atm_abs"][1]-osc["IO"]["solar"][2]))]:
            # In the finite matching convention of 1558 and the canonical singlet portal.
            ms = mz*ratio
            min_shift = transport * (floor*1e-9)*ms**3/(4*math.pi**2*v**2)
            mt_cap = mz*(ceiling/min_shift)**(1/3)
            assert mt_cap < mz and 0.5*min_shift > ceiling
            rows[name] = {"light_mass_floor_eV": floor, "minimum_MS_GeV": ms,
                          "singlet_Higgs_mass_parameter_sq_floor_GeV2": min_shift,
                          "multiple_of_ATLAS_pole_budget": min_shift/ceiling,
                          "MT_ceiling_if_sole_positive_budget_GeV": mt_cap}
        table[str(offset)] = rows
    return {"rge_log_ratio_slopes_for_t_vq": [float(slope_t), float(slope_vq)],
            "ATLAS_comparison_mass_parameter_ceiling_GeV2": ceiling,
            "minimum_positive_transport": transport, "windows": table}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new", action="store_true")
    flags.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    dest = HERE / "independent_results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    else:
        assert json.loads(dest.read_text()) == result
    print(json.dumps(result))
