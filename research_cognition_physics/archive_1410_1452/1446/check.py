"""A common two-channel seed versus two neutrino oscillation frequencies."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def box(values, width):
    return values[0]-width*values[1], values[0]+width*values[-1]


def evaluate():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf8"))
    old = json.loads((HERE/cfg["prior"]).read_text(encoding="utf8"))
    juno = cfg["juno_207_days"]
    t = math.sqrt(2)-1
    r = t*t
    # q=m_middle/m_heavy; normal ordering with m1=0
    n = cfg["widths"]
    solar = box(juno["solar_21"], n)
    normal = box(juno["normal_atmospheric_31"], n)
    inverted = box(juno["inverted_abs_31"], n)
    q_lo, q_hi = math.sqrt(solar[0]/normal[1]), math.sqrt(solar[1]/normal[0])
    io_min = math.sqrt(inverted[0]/(inverted[0]+solar[1]))
    atm32 = old["NO"]["atmospheric_gap"]
    factor = r**2 / (1-r**2)
    predicted_solar_t2k = [factor*x for x in box(atm32, n)]
    assert q_lo < r < q_hi < t < io_min
    shape_lo = q_lo + 1/q_lo
    shape_seed = r + 1/r
    max_overlap = math.sqrt((shape_lo-shape_seed)/(shape_lo-2))
    # Complex source vectors with adjustable Hermitian overlap and Majorana phase.
    v = np.array([1.,0.,0.])
    tests = 0
    for h in [0.,0.04,max_overlap,0.4,0.9]:
        for angle in [0.,0.51,1.7,math.pi]:
            w = np.array([h,math.sqrt(1-h*h)*np.exp(0.37j),0j])
            kernel = np.outer(v,v) + r*np.exp(1j*angle)*np.outer(w,w)
            heavy,middle,zero = np.linalg.svd(kernel,compute_uv=False)
            q = middle/heavy
            trace_ratio = (shape_seed + 2*h*h*math.cos(angle))/(1-h*h)
            assert abs(q+1/q-trace_ratio) < 1e-10 and q <= r + 1e-13 and zero < 1e-14
            if h == max_overlap and angle == math.pi:
                assert abs(q-q_lo) < 1e-12
            tests += 1
    return {"seed_amplitude_ratio":t, "two_vertex_mass_ratio":r,
            "predicted_normal_gap21_over_gap31":r*r,
            "JUNO_central_gap_ratio":juno["solar_21"][0]/juno["normal_atmospheric_31"][0],
            "normal_mass_ratio_box":[q_lo,q_hi], "inverted_mass_ratio_floor":io_min,
            "T2K_predicted_gap21_central":factor*atm32[0],
            "T2K_predicted_gap21_box":predicted_solar_t2k, "JUNO_gap21_box":list(solar),
            "JUNO_atm_implied_solar_center":r*r*juno["normal_atmospheric_31"][0],
            "orthogonal_relative_transfer_allowed":[q_lo/r,q_hi/r],
            "sharp_source_overlap_upper_if_equal_transfer":max_overlap,
            "departure_from_orthogonality_angle_deg_max":math.degrees(math.asin(max_overlap)),
            "mass_sum_eV_using_JUNO_atm31":math.sqrt(juno["normal_atmospheric_31"][0])*(1+r),
            "complex_matrix_checks":tests}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=evaluate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        assert out==json.loads(path.read_text(encoding="utf8"))
    print("1446 OK",out["predicted_normal_gap21_over_gap31"])
