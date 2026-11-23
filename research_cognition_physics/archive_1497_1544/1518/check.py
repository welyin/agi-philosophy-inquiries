"""One-loop matching and singular-value split budget for the 1517 matter table."""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
PI = math.pi
SM = np.array([4.1, -19/6, -7])
MSSM = np.array([33/5, 1, -3])
COLOR = np.array([.4, 0, 1])
WEAK = np.array([.6, 1, 0])
SCHEME = np.array([0, 2, 3]) / (12*PI)

def read(path):
    return json.loads((BASE / path).read_text(encoding="utf-8-sig"))

def match(couplings, mz, soft, ratio):
    inverse = 4*PI / np.array(couplings)**2 - SCHEME
    log_soft = math.log(soft/mz)
    h = np.dot([-1, 3, -2], inverse)
    p = np.dot([5, -3, -2], inverse)
    mx = mz * math.exp((2*PI*p - 8*log_soft)/36)
    eff = ratio*mz * math.exp((2*PI*h + 2*log_soft)*5/12)
    return dict(MX_GeV=mx, even_effective_mass_GeV=eff, effective_over_MX=eff/mx)

def min_split(target, floor):
    low, high = 0., target
    for _ in range(100):
        mid = (low+high)/2
        if mid*(1+mid/floor)**2 >= target:
            high = mid
        else:
            low = mid
    return high

def odd_masses(cfg, mx, ratio):
    spec = importlib.util.spec_from_file_location("round1516", BASE / cfg["flavor_program"])
    flavor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(flavor)
    core = flavor.parent_core(read(cfg["flavor_input"]))
    phi = cfg["phase_degrees"]
    _, roots, _ = flavor.mixing_core(core, phi)
    phase = np.exp(1j*math.radians(phi))
    real = core[3][1:] / abs(1 + 2*phase*roots)
    unit = np.linalg.norm(real*roots)
    odd_color = mx / (5*unit)
    np.testing.assert_allclose(np.linalg.norm(np.r_[core[3][0], real*(1-3*phase*roots)]), ratio)
    return unit, odd_color, ratio*odd_color

def rotation(i, j, angle):
    v = np.eye(3)
    v[i,i] = v[j,j] = math.cos(angle)
    v[i,j], v[j,i] = -math.sin(angle), math.sin(angle)
    return v

def example(cfg, couplings, mz, soft, ratio, central):
    mx, eff = central["MX_GeV"], central["even_effective_mass_GeV"]
    k = cfg["split_over_X"]*mx
    d = k/(math.sqrt(eff/k)-1)
    ell, odd_t, odd_d = odd_masses(cfg, mx, ratio)
    a = 4*PI / np.array(couplings)**2 - SCHEME
    inverse_from_each = a - (SM*math.log(soft/mz) + MSSM*math.log(mx/soft)
        + COLOR*(math.log(mx/k)+2*math.log(mx/(d+k))+math.log(mx/odd_t))
        + WEAK*(2*math.log(mx/d)+math.log(mx/odd_d)))/(2*PI)
    assert np.ptp(inverse_from_each) < 1e-12
    g = math.sqrt(4*PI/inverse_from_each[0])
    vev = mx/(5*g)
    cutoff = cfg["cutoff_GeV"]
    left, right = rotation(0,1,.49), rotation(0,2,-.43)
    # At the physical threshold the light higgsino pair has mass soft.
    doublets = left @ np.diag([soft,d,d]) @ right.T
    triplets = left @ np.diag([k,d+k,d+k]) @ right.T
    c = (doublets-triplets)/5
    mu = doublets-9*c
    np.testing.assert_allclose(mu+4*c, triplets, atol=.1)
    assert np.linalg.norm(triplets-doublets,2) <= k*(1+1e-13)
    assert np.linalg.det(triplets)>0
    np.testing.assert_allclose(np.linalg.det(triplets)/d**2, eff, rtol=1e-12)
    small_tolerance = cfg["light_mass_perturbation_GeV"]/(9*(k-soft)/5)
    return dict(soft_GeV=soft, gauge_inverse_from_each=inverse_from_each.tolist(),
        gauge_unified=g, Sigma_vev_GeV=vev, cutoff_over_vev=cutoff/vev,
        norm_odd_flavor_source_in_color_mass_units=ell,
        odd_color_GeV=odd_t, odd_weak_GeV=odd_d, odd_over_X=[odd_t/mx,odd_d/mx],
        split_GeV=k, even_doublets_GeV=[soft,d,d], even_triplets_GeV=[k,d+k,d+k],
        norm_dimensionless_quadratic=cutoff*np.linalg.norm(c,2)/vev**2,
        light_left_right_overlap=float(left[:,0]@right[:,0]),
        fractional_quadratic_only_shift_for_one_extra_TeV=small_tolerance)

def calculate():
    cfg = read("inputs.json")
    table, odd = read(cfg["gauge_table"]), read(cfg["odd_ratio"])
    gauge = np.array([table[key] for key in ["g1_GUT", "g2", "g3"]])
    mz, soft = table["mu_ref_GeV"], cfg["susy_threshold_GeV"]
    ratio = odd["R_center"]
    central = match(gauge[:,0], mz, soft, ratio)
    data = []
    for signs in itertools.product((-1,1),repeat=3):
        couplings = gauge[:,0]+cfg["gauge_half_widths"]*np.array(signs)*gauge[:,1]
        for r in odd["R_two_width"]:
            data.append(match(couplings, mz, soft, r))
    bounds = {key:[min(row[key] for row in data),max(row[key] for row in data)] for key in central}
    floor = cfg["heavy_doublet_floor_over_X"]
    target = bounds["effective_over_MX"][0]
    k0 = cfg["split_over_X"]
    assert target > k0
    return dict(central=central, two_width_rectangle=bounds,
        necessary_split_over_X_for_both_doublets_above_floor=min_split(target,floor),
        smallest_doublet_ceiling_over_X_if_split_at_most_k0=k0/(math.sqrt(target/k0)-1),
        alternatives_soft_central={str(v):match(gauge[:,0],mz,v,ratio)
            for v in cfg["alternative_susy_threshold_GeV"]},
        witness=example(cfg,gauge[:,0],mz,soft,ratio,central))

def main():
    p = argparse.ArgumentParser()
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args = p.parse_args()
    result = calculate()
    if args.write:
        (BASE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    else:
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a:
                    compare(a[key],b[key])
            else:
                np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-14)
        compare(result,read("results.json"))
        print("1518 OK: two independent combinations and all three matched couplings")

if __name__ == "__main__":
    main()
