"""Determine the physical hypercharge index from a color/weak crossing at one loop."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def matching(raw, exposure, mz):
    g1, g2, g3 = raw
    ay, a2, a3 = 4*np.pi/np.array([g1**2*3/5, g2**2, g3**2])
    log = 6*(a2-a3)/(23+exposure)
    ax = a3+7*log
    ay_high = ay - (41+exposure)/6*log
    k = ay_high/ax
    a2_high = a2 - (-19+exposure)/6*log
    assert np.isclose(a2_high, ax) and np.isclose(k, 1+(ay-a2-10*log)/ax)
    # Reconstruct the three weak-scale couplings by running in two actual segments.
    rows = np.array([[41/6, -19/6, -7], [42/6, -18/6, -7]])
    restored = np.array([k*ax, ax, ax]) + log*np.array([1-exposure, exposure]) @ rows
    assert np.allclose(restored, [ay, a2, a3])
    return float(k), float(np.log10(mz)+2*np.pi*log/np.log(10))


def compute():
    info = json.loads((HERE.parent/"1369"/"inputs.json").read_text(encoding="utf-8"))
    g, err = np.array([info[key] for key in ["g1_GUT", "g2", "g3"]]).T
    endpoints = [matching(g, e, info["mu_ref_GeV"]) for e in (0, 1)]
    candidates = [matching(g+2*err*np.array(s), e, info["mu_ref_GeV"])
                  for s in product([-1, 1], repeat=3) for e in (0, 1)]
    k, scales = np.array(candidates).T
    upper = float(max(k))
    canonical = 2*(3/9+2/4)
    assert upper < canonical and endpoints[0][0] < endpoints[1][0]
    return {"round": 1417, "central_k_log10_scale_no_new_doublet_exposure": endpoints[0],
            "central_k_log10_scale_full_new_doublet_exposure": endpoints[1],
            "two_width_k_min_max": [float(min(k)), upper],
            "two_width_log10_MX_GeV_min_max": [float(min(scales)), float(max(scales))],
            "one_weak_doublet_embedded_five_block_trace_floor": canonical,
            "minimum_gap_in_k": canonical-upper}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.check:
        saved = json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        assert result.keys() == saved.keys() and all(np.allclose(result[k],saved[k],rtol=1e-10) for k in saved)
        print("1417 inverse hypercharge index and segmented running verified")
    else:
        print(json.dumps(result, indent=2))
