"""Cross-check the mass moment and its contribution to weak and scalar running."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOOP = 16 * math.pi**2


def boson_part(lam, g2, gy):
    return (24*lam**2 - (9*g2**2 + 3*gy**2)*lam
            + 3*(2*g2**4 + (g2**2+gy**2)**2)/8)


def bounds(masses, v, floor, lam, g2, gy):
    lowest = np.column_stack((masses, [floor, floor]))
    y2 = 2*(lowest/v)**2
    negative = 6 * np.sum(np.maximum(y2-lam, 0.0)**2)
    boson = boson_part(lam, g2, gy)
    demand = (negative-boson)/(4*lam**2)
    n_min = max(0, math.ceil(demand))
    b_min = -43/6 + n_min/3
    a_low = 4*math.pi/g2**2
    finite = []
    for scale_ratio, max_growth in [(2., 1.10), (10., 1.10)]:
        max_b = 2*math.pi*a_low*(1-1/max_growth)/math.log(scale_ratio)
        n_max = math.floor(3*(max_b + 43/6))
        finite.append(dict(scale_ratio=scale_ratio, max_doublets=n_max,
                           beta_lambda_upper=(boson+4*lam**2*n_max-negative)/LOOP))
    return dict(assumed_lam=lam, assumed_g2=g2, assumed_gY=gy, mass_floor_GeV=floor,
                min_weighted_doublets_for_nonnegative_beta=n_min, min_b2=b_min,
                max_scale_ratio_before_ten_percent_alpha_growth=
                    math.exp(2*math.pi*a_low*(1-1/1.1)/b_min) if b_min>0 else None,
                two_point_comparisons=finite)


def random_check(seed):
    rng = np.random.default_rng(seed)
    lam = rng.uniform(.05, .3)
    # Quark and charged-lepton matrices include both direct and mirror blocks.
    ranks = []
    t2 = t4 = slots = 0
    for colour, n_plus, n_minus in [(3,4,2), (1,4,3)]:
        slots += colour*(n_plus+n_minus)
        for sector in range(2):
            if colour == 1 and sector == 1:
                shapes = [(n_plus+n_minus, 19)]
            else:
                shapes = [(n_plus, 11), (n_minus, 7)]
            for shape in shapes:
                y = (rng.normal(size=shape) + 1j*rng.normal(size=shape)) * rng.uniform(.02,.6)
                s = np.linalg.svd(y, compute_uv=False)
                gram = y @ y.conj().T
                t2 += colour*np.trace(gram).real
                t4 += colour*np.trace(gram @ gram).real
                ranks.extend([i for i in s**2 for _ in range(colour)])
    assert len(ranks) <= 2*slots
    direct = 4*lam*t2 - 2*t4
    completed = 2*lam**2*len(ranks)-2*np.sum((np.array(ranks)-lam)**2)
    np.testing.assert_allclose(completed, direct, rtol=1e-11, atol=1e-10)
    assert direct <= 4*lam**2*slots + 1e-10


def run():
    main = json.loads((HERE/"results.json").read_bytes())
    old = json.loads((HERE/"../1767/results.json").read_bytes())["matching"]
    v = old["vev_GeV"]
    masses = np.array([np.array(old["high_yu"])*v/math.sqrt(2),
                       [old["gate_d"]["light_GeV"],old["daughter_yd"][1]*v/math.sqrt(2),
                        old["gate_b"]["light_GeV"]]])
    summary = []
    for key, row in main["moments"]["floor_rows"].items():
        m = np.column_stack((masses, [float(key)]*2))
        for power in [2,4]:
            np.testing.assert_allclose(row[f"T{power}_min"], 3*np.sum((math.sqrt(2)*m/v)**power),
                                       rtol=2e-9)
    for seed in range(72):
        random_check(1778+seed)
    # Test the bound against an explicit anomaly-free, purely chiral four-family
    # assignment and any number of extra vectorlike leptons with active singlets.
    lam,g2,gy = .1,.63,.36
    for floor in [200.,300.,400.]:
        summary.append(bounds(masses,v,floor,lam,g2,gy))
    y2_quark = 2*(np.column_stack((masses, [300.,300.]))/v)**2
    b0 = boson_part(lam,g2,gy)+3*np.sum(4*lam*y2_quark-2*y2_quark**2)
    # Four standard lepton doublets, including Dirac singlets, here have zero
    # Yukawas; each added vectorlike lepton doublet has four couplings sqrt(lam).
    n_pairs = math.ceil(-b0/(8*lam**2))
    assert b0 + n_pairs*8*lam**2 >= 0 > b0 + (n_pairs-1)*8*lam**2
    n_doublets = 16+2*n_pairs
    assert n_doublets >= summary[1]["min_weighted_doublets_for_nonnegative_beta"]
    b_2 = -43/6+n_doublets/3
    a = 4*math.pi/g2**2
    return dict(benchmark_inputs_are_conditional=True, six_masses_at_registered_scale_GeV=masses.tolist(),
                coefficient_checks=72, independent_bounds=summary,
                compensation_example=dict(vectorlike_lepton_doublet_pairs=n_pairs,
                    total_colour_weighted_Weyl_doublets=n_doublets,
                    one_loop_b2=b_2, instantaneous_beta_lambda=(b0+n_pairs*8*lam**2)/LOOP,
                    scale_factor_for_ten_percent_weak_growth=math.exp(2*math.pi*(a-a/1.1)/b_2)))


def rounded(x):
    if isinstance(x, dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,float) else x


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(run())
    payload=json.dumps(result,indent=2)+"\n"
    out=HERE/"independent_results.json"
    if args.check:
        assert json.loads(out.read_bytes())==result
    else:
        out.write_text(payload,encoding="utf-8")
    print("1778 weak and scalar coefficients",hashlib.sha256(payload.encode()).hexdigest()[:12])
